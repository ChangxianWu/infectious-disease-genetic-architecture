# -*- coding: utf-8 -*-
""" 12_matrix_runner.py — 主矩阵：munge 全部 FinnGen 终点 → 两两 bivariate LDSC


用法：python 12_matrix_runner.py [--munge-only]
输出：work/output/ldsc_rg/ldsc_matrix_<DATE>.csv（300 行，含 status + BH/Bonferroni）
"""
import csv
import io
import sys
import time
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, norm

sys.path.insert(0, "/path/to/project/work/code")
sys.path.insert(0, "/path/to/project/work/code/ldsc_py3")
from constants import DATA, OUT, DATE, FINNGEN_ENDPOINTS, GWAS_DATA  # noqa: E402
import regressions as reg  # noqa: E402  (移植版；seterr 语义与官方一致)

MUNGED = DATA / "munged"
RGDIR = OUT / "ldsc_rg"
RGDIR.mkdir(parents=True, exist_ok=True)
MANIFEST = DATA / "finngen_R12_pheno_n.tsv"
LDSC_REF = GWAS_DATA / "ldsc_ref"
MTB_BLOCKED = {"AB1_TUBERCULOSIS", "AB1_RESP_TUBERCU_CONF", "AB1_RESP_TUBERCU_NONCONF",
               "AB1TUBERCU_MILIARY", "AB1_TUBERCU_OTH_ORG", "AB1_SEQULAE_TUBERCU", "TBC_RESP"}
N_BLOCKS = 200


def load_manifest_counts():
    rows = {}
    with open(MANIFEST, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[r["phenocode"]] = (int(r["num_cases"]), int(r["num_controls"]))
    return rows


def munge_missing(counts):
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "munge", "/path/to/project/work/code/10_munge_finngen.py")
    mg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mg)
    todo = []
    for label, ph in FINNGEN_ENDPOINTS.items():
        if ph in MTB_BLOCKED or (MUNGED / ("%s.sumstats" % ph)).exists():
            continue
        if ph not in counts:
            raise KeyError("官方 manifest 无 %s" % ph)
        nc, nk = counts[ph]
        if nc < 500:
            print("[munge] %s (%s) 病例 %d <500，功效门剔除" % (label, ph, nc))
            continue
        todo.append((label, ph, nc, nk))
    if not todo:
        print("[munge] 无需补跑")
        return
    print("[munge] 加载 frq/panel（一次）...")
    frq = mg.load_1000g_frq()
    panel_ids = set()
    for i in range(1, 23):
        d = pd.read_csv(mg.LDSC_REF / "LDscore" / ("LDscore.%d.l2.ldscore.gz" % i),
                        sep="\t", usecols=["SNP"])
        panel_ids.update(d["SNP"])
    for label, ph, nc, nk in todo:
        t0 = time.time()
        c = mg.munge_one(ph, nc, nk, frq, panel_ids)
        buf = io.StringIO()
        w = csv.DictWriter(buf, fieldnames=list(c.keys()), delimiter="\t")
        cnt_path = MUNGED / "munge_counts.csv"
        if not cnt_path.exists():
            w.writeheader()
        w.writerow(c)
        mode = "a" if cnt_path.exists() else "w"
        with open(cnt_path, mode, encoding="utf-8-sig", newline="") as fh:
            fh.write(buf.getvalue())
        print("[munge] %s=%s 完成 %s SNPs mean_chi2=%.2f (%.0fs)"
              % (label, ph, c["final_snps"], c["mean_chisq"], time.time() - t0))


def load_master(label_to_ph):
    """一次载入共用 SNP 集 + 每性状 Z/A1/A2/N（v2 提速核心）。"""
    l2s, wls, m550 = [], [], []
    for i in range(1, 23):
        l2s.append(pd.read_csv(LDSC_REF / "LDscore" / ("LDscore.%d.l2.ldscore.gz" % i),
                               sep="\t", usecols=["SNP", "L2"]))
        wls.append(pd.read_csv(LDSC_REF / "1000G_Phase3_weights_hm3_no_MHC"
                               / ("weights.hm3_noMHC.%d.l2.ldscore.gz" % i),
                               sep="\t", usecols=["SNP", "L2"]))
        m550.append(float(open(LDSC_REF / "LDscore" / ("LDscore.%d.l2.M_5_50" % i)).read().strip()))
    l2 = pd.concat(l2s, ignore_index=True).drop_duplicates("SNP", keep="first")
    wl = pd.concat(wls, ignore_index=True).drop_duplicates("SNP", keep="first").rename(
        columns={"L2": "W_L2"})
    base = l2.merge(wl, on="SNP", how="inner")
    meta = {"SNP": base["SNP"].values, "L2": base["L2"].values.astype(np.float64),
            "W_L2": base["W_L2"].values.astype(np.float64)}
    idx = pd.Index(meta["SNP"])
    traits = {}
    for label, ph in label_to_ph.items():
        d = pd.read_csv(MUNGED / ("%s.sumstats" % ph), sep="\t",
                        usecols=["SNP", "A1", "A2", "N", "Z"])
        d = d.drop_duplicates("SNP", keep="first").set_index("SNP").reindex(idx)
        traits[label] = d
    return meta, traits, float(np.sum(m550))


def align(t1, t2):
    """两性状等位对齐：A1 同向保留；A1 互换（回文）标记剔除。"""
    same = (t1["A1"] == t2["A1"]) & (t1["A2"] == t2["A2"])
    swap = (t1["A1"] == t2["A2"]) & (t1["A2"] == t2["A1"])
    ok = same | swap
    return ok, same


def run_pair(lab1, lab2, t1, t2, meta, M):
    ok, same = align(t1, t2)
    z1 = t1["Z"].values[ok.values]
    z2raw = t2["Z"].values[ok.values]
    z2 = np.where(same.values[ok.values], z2raw, -z2raw)
    N1 = t1["N"].values[ok.values].astype(np.float64)
    N2 = t2["N"].values[ok.values].astype(np.float64)
    x = meta["L2"][ok.values].reshape((-1, 1))
    w = meta["W_L2"][ok.values].reshape((-1, 1))
    n = len(z1)
    row = {"axis1": lab1, "axis2": lab2, "n_snps": n, "M_5_50": M,
           "n_allele_swap": int((~same.values[ok.values]).sum()),
           "mean_chisq_1": float(np.mean(z1 ** 2)), "mean_chisq_2": float(np.mean(z2 ** 2))}
    t0 = time.time()
    try:
        rg = reg.RG(z1.reshape((n, 1)), z2.reshape((n, 1)), x, w,
                    N1.reshape((n, 1)), N2.reshape((n, 1)), np.array([[M]]), n_blocks=N_BLOCKS)
    except FloatingPointError as exc:
        # 官方 ldsc.py 对每对 try/except 后继续（sumstats.estimate_rg）；同规记录不中止
        row.update({"status": "fp_error:%s" % str(exc)[:40], "rg": np.nan, "rg_se": np.nan,
                    "rg_p": np.nan, "rg_z": np.nan, "runtime_s": round(time.time() - t0, 1)})
        return row
    row["runtime_s"] = round(time.time() - t0, 1)
    if rg._negative_hsq is not None:
        row.update({"status": "negative_hsq", "rg": np.nan, "rg_se": np.nan,
                    "rg_p": np.nan, "rg_z": np.nan,
                    "h2_1": float(rg.hsq1.tot), "h2_1_se": float(rg.hsq1.tot_se),
                    "h2_2": float(rg.hsq2.tot), "h2_2_se": float(rg.hsq2.tot_se)})
    else:
        row["rg_raw"] = rg.rg_ratio          # 原始估计（审计用；越界值在此保留）
        row.update({"status": "ok", "rg": rg.rg_ratio, "rg_se": rg.rg_se,
                    "rg_p": float(rg.p), "rg_z": float(rg.z),
                    "h2_1": float(rg.hsq1.tot), "h2_1_se": float(rg.hsq1.tot_se),
                    "h2_2": float(rg.hsq2.tot), "h2_2_se": float(rg.hsq2.tot_se)})
    # 截距 QC 列（2026-10-25 增补，解读防御：样本重叠/分层混杂的标准检查）。
    # gencov 截距≈1 表示 cross-trait 回归无超出 LD 的系统膨胀。确定量（jackknife 估计，
    # 非墙钟），不受 NONDET_COLS 指纹口径影响；fp_error 分支提前 return，无此三列。
    for att, key in ((getattr(rg, "hsq1", None), "int_hsq_1"),
                     (getattr(rg, "hsq2", None), "int_hsq_2"),
                     (getattr(rg, "gencov", None), "int_gencov")):
        try:
            row[key] = float(att.intercept)
        except (TypeError, ValueError, AttributeError):
            row[key] = np.nan
    # 方向：rg_ratio 越界（|rg|>1.2）官方判 nan（"rg out of bounds"）。
    # 口径纪律：**rg 列只在可解释时为数值，其余一律 nan**（与官方输出语义一致）；
    # 原始越界值保留在 rg_raw，避免"把 nan 写成数值"或"把数值写成 nan 丢证据"两个错误。
    if row.get("status") == "ok" and abs(row["rg"]) > 1.2:
        row["status"] = "rg_out_of_bounds"
        row["rg"] = np.nan
        row["rg_se"] = np.nan
        row["rg_p"] = np.nan
        row["rg_z"] = np.nan
    return row


def main():
    t_start = time.time()
    counts = load_manifest_counts()
    munge_missing(counts)
    if "--munge-only" in sys.argv:
        print("[matrix] munge-only，退出")
        return
    label_to_ph = {lab: ph for lab, ph in FINNGEN_ENDPOINTS.items()
                   if ph not in MTB_BLOCKED and (MUNGED / ("%s.sumstats" % ph)).exists()}
    print("[matrix] 载入 %d 性状主表 ..." % len(label_to_ph))
    meta, traits, M = load_master(label_to_ph)
    print("[matrix] 共用 SNP 集 %d | M_5_50=%.0f" % (len(meta["SNP"]), M))
    rows = []
    pairs = list(combinations(sorted(traits.keys()), 2))
    for k, (lab1, lab2) in enumerate(pairs, 1):
        row = run_pair(lab1, lab2, traits[lab1], traits[lab2], meta, M)
        rows.append(row)
        if k % 25 == 0 or k == len(pairs):
            print("[matrix %d/%d] 最新 %s x %s: rg=%s status=%s"
                  % (k, len(pairs), lab1, lab2,
                     ("%.3f" % row["rg"]) if pd.notna(row.get("rg")) else "NA", row["status"]))
    df = pd.DataFrame(rows)
    # 功效门（规则 89 扩充）：两性状 h² 单侧 p<0.05 且 mean χ²≥1.02
    for i in (1, 2):
        se = df["h2_%d_se" % i]
        df["h2_%d_p" % i] = norm.sf(df["h2_%d" % i] / se)
        df["h2_%d_p" % i] = df["h2_%d_p" % i].where(df["h2_%d" % i].notna(), np.nan)
    df["power_ok"] = (df["status"] == "ok") & (df["h2_1_p"] < 0.05) & (df["h2_2_p"] < 0.05) \
        & (df["mean_chisq_1"] >= 1.02) & (df["mean_chisq_2"] >= 1.02)
    ok = df["power_ok"] == True  # noqa: E712
    df["rg_q_bh"] = np.nan
    df["rg_p_bonf"] = np.nan
    n_ok = int(ok.sum())
    if n_ok:
        df.loc[ok, "rg_q_bh"] = false_discovery_control(df.loc[ok, "rg_p"].values, method="bh")
        df.loc[ok, "rg_p_bonf"] = np.minimum(1.0, df.loc[ok, "rg_p"].values * n_ok)
    out = RGDIR / ("ldsc_matrix_%s.csv" % DATE)
    df.to_csv(out, index=False)
    print()
    print("[matrix] 完成 %d 对（%.1f min）" % (len(df), (time.time() - t_start) / 60))
    print("[matrix] status 分布：")
    print(df["status"].str.split(":").str[0].value_counts().to_string())
    print("[matrix] power_ok=%d | BH q<0.05: %d | Bonferroni p<0.05: %d"
          % (n_ok, int((df["rg_q_bh"] < 0.05).sum()), int((df["rg_p_bonf"] < 0.05).sum())))
    print("[matrix] -> %s" % out.name)


if __name__ == "__main__":
    main()
