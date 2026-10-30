# -*- coding: utf-8 -*-
""" 45_autoimmune_mhc_sens.py — E3 附加：自免×感染共享的 MHC 包含敏感性
（2026-10-25）

问题：主 E3 层（43 号）与主矩阵同口径用 `weights.hm3_noMHC` → "非 MHC 架构"。
      RA×细菌预演信号可能由 MHC 驱动（RA 是典型 MHC-II 病）。须给出 MHC **包含**
      时的对照，才能区分"自免×感染共享=泛 MHC 效应"与"非 MHC 真共享"。

口径（完全继承 23_mhc_sensitivity.py 的"含 MHC 近似"）：
  - 参考 LD 与权重都用同一份官方 `LDscore.*.l2.ldscore.gz`（含 MHC）自加权；
  - 只跑主 E3 层 **过功效门** 的对（其余对无解读意义）；
  - 与 43 号主口径（noMHC 权重）逐对并列：rg 差 Δ、符号翻转、p 值变化；
  - 出界策略与主层一致：|rg|>1.2 判 nan（rg_raw 保留证据，不静默）；
  - 诚实性声明同 23 号：非官方严格对照，解读上限为"方向/显著性是否改变"，
    不得给出 MHC 定量贡献百分比。

前置：43_autoimmune_runner.py 已跑完（读其主层 CSV 过门对）。
用法：python 45_autoimmune_mhc_sens.py
输出：work/output/autoimmune_mhc_sens.csv + autoimmune_mhc_sens_summary.csv
"""
import importlib.util
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, "/path/to/project/work/code")
sys.path.insert(0, "/path/to/project/work/code/ldsc_py3")
from constants import DATE, GWAS_DATA, OUT  # noqa: E402
import regressions as reg  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "mr", "/path/to/project/work/code/12_matrix_runner.py")
mr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mr)

LDSC_REF = GWAS_DATA / "ldsc_ref"
PRIMARY = OUT / ("ldsc_rg_autoimmune_%s.csv" % DATE)
OUTCSV = OUT / "autoimmune_mhc_sens.csv"
SUMCSV = OUT / "autoimmune_mhc_sens_summary.csv"


def load_ref_full():
    l2s, m550 = [], []
    for i in range(1, 23):
        l2s.append(pd.read_csv(LDSC_REF / "LDscore" / ("LDscore.%d.l2.ldscore.gz" % i),
                               sep="\t", usecols=["SNP", "L2"]))
        m550.append(float(open(LDSC_REF / "LDscore" / ("LDscore.%d.l2.M_5_50" % i)).read().strip()))
    l2 = pd.concat(l2s, ignore_index=True).drop_duplicates("SNP", keep="first")
    return l2, float(np.sum(m550))


def main():
    prim = pd.read_csv(PRIMARY)
    pw = prim[prim["power_ok"] == True]  # noqa: E712
    if not len(pw):
        raise SystemExit("[sens] 主 E3 层无过门对——先跑 43 号")
    pairs = list({(r["axis1"], r["axis2"]) for _, r in pw.iterrows()})
    labels = sorted({l for p in pairs for l in p})
    print("[sens] 过门对 %d，涉及性状 %d" % (len(pairs), len(labels)))

    l2full, M = load_ref_full()
    idx_full = pd.Index(l2full["SNP"])
    z_l2 = l2full["L2"].values.astype(np.float64)

    full = {}
    # label→phenocode：自免 label=phenocode；感染 label 走 FINNGEN_ENDPOINTS 映射
    # （munge 产物以 phenocode 命名——43 号 load_master 同规）
    label_to_ph = dict(mr.FINNGEN_ENDPOINTS)
    for lab in labels:
        ph = label_to_ph.get(lab, lab)
        d = pd.read_csv(mr.MUNGED / ("%s.sumstats" % ph), sep="\t",
                        usecols=["SNP", "A1", "A2", "N", "Z"]).drop_duplicates("SNP")
        full[lab] = d.set_index("SNP").reindex(idx_full)

    rows = []
    for k, (l1, l2lab) in enumerate(sorted(pairs), 1):
        a, b = full[l1], full[l2lab]
        same = (a["A1"] == b["A1"]) & (a["A2"] == b["A2"])
        swap = (a["A1"] == b["A2"]) & (a["A2"] == b["A1"])
        ok = (same | swap).values
        ok = ok & np.isfinite(a["Z"].values) & np.isfinite(b["Z"].values)
        z1 = a["Z"].values[ok]
        z2raw = b["Z"].values[ok]
        z2 = np.where(same.values[ok], z2raw, -z2raw)
        N1 = a["N"].values[ok].astype(float).reshape((-1, 1))
        N2 = b["N"].values[ok].astype(float).reshape((-1, 1))
        x = z_l2[ok].reshape((-1, 1))
        w = x
        n = len(z1)
        row = {"trait1": l1, "trait2": l2lab, "n_snp_incl_mhc": n}
        try:
            rg = reg.RG(z1.reshape((n, 1)), z2.reshape((n, 1)), x, w, N1, N2,
                        np.array([[M]]), n_blocks=200)
            status = "ok" if rg._negative_hsq is None else "negative_hsq"
        except FloatingPointError as exc:
            status = "fp_error:%s" % str(exc)[:30]
            rg = None
        row["status_incl_mhc"] = status
        if rg is not None and status == "ok":
            r = float(rg.rg_ratio)
            row.update({"rg_raw_incl_mhc": r,
                        "rg_incl_mhc": r if abs(r) <= 1.2 else np.nan,
                        "rg_se_incl_mhc": float(rg.rg_se) if abs(r) <= 1.2 else np.nan,
                        "p_incl_mhc": max(float(rg.p), 1e-300) if abs(r) <= 1.2 else np.nan})
        else:
            row.update({"rg_raw_incl_mhc": np.nan, "rg_incl_mhc": np.nan,
                        "rg_se_incl_mhc": np.nan, "p_incl_mhc": np.nan})
        rows.append(row)
        if k % 25 == 0 or k == len(pairs):
            print("[sens %d/%d] %s x %s" % (k, len(pairs), l1, l2lab))
    sens = pd.DataFrame(rows)

    # 与主口径并列（43 号：axis1/axis2 同序）
    j = pw.merge(sens, left_on=["axis1", "axis2"], right_on=["trait1", "trait2"],
                 how="left", suffixes=("", "_sens"))
    j["delta_rg"] = j["rg_incl_mhc"] - j["rg"]
    j["sign_flip"] = (np.sign(j["rg_incl_mhc"]) != np.sign(j["rg"])) & \
        j["rg_incl_mhc"].notna() & j["rg"].notna()
    j[["axis1", "axis2", "rg", "rg_se", "rg_p", "rg_q_bh_layer",
       "rg_incl_mhc", "rg_se_incl_mhc", "p_incl_mhc", "delta_rg", "sign_flip",
       "status_incl_mhc", "n_snp_incl_mhc"]].to_csv(OUTCSV, index=False)

    cmp_ok = j["rg_incl_mhc"].notna() & j["rg"].notna()
    auto_labels = sorted({r["axis1"] for _, r in j.iterrows()})
    by_auto = {a: j[(j["axis1"] == a) & cmp_ok] for a in auto_labels}
    summ_rows = [
        {"key": "n_pairs_compared", "value": str(int(cmp_ok.sum())),
         "note": "主层过门且含 MHC 口径 ok 的对（口径=23 号自加权近似）"},
        {"key": "n_sign_flip", "value": str(int(j["sign_flip"].sum())),
         "note": "含 MHC 后符号翻转数"},
        {"key": "median_delta_rg", "value": "%.3f" % j.loc[cmp_ok, "delta_rg"].median(),
         "note": "Δrg 中位（正=含 MHC 后共享变大）"},
        {"key": "n_delta_gt_0p1", "value": str(int((j.loc[cmp_ok, "delta_rg"].abs() > 0.1).sum())),
         "note": "|Δrg|>0.1 的对数"},
        {"key": "honesty_note", "value": "approx",
         "note": "非官方严格对照（23 号口径）；解读上限=方向/显著性是否改变"},
    ]
    for a in auto_labels:
        sub = by_auto[a]
        if len(sub):
            summ_rows.append({"key": "median_delta_%s" % a,
                              "value": "%.3f" % sub["delta_rg"].median(),
                              "note": "%s 各对 Δrg 中位（%d 对）" % (a, len(sub))})
    pd.DataFrame(summ_rows).to_csv(SUMCSV, index=False)
    print("[sens] 比较 %d 对；符号翻转 %d；Δrg 中位 %.3f"
          % (cmp_ok.sum(), j["sign_flip"].sum(), j.loc[cmp_ok, "delta_rg"].median()))


if __name__ == "__main__":
    main()
