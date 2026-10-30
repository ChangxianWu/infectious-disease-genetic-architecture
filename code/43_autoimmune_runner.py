# -*- coding: utf-8 -*-
""" 43_autoimmune_runner.py — E3 自身免疫锚层：munge 自免终点 → 自免×感染 rg 批跑
（2026-10-25）

设计：
  10 个自免终点（T1D/RA/MS/乳糜泻/UC-strict/Crohn-strict/银屑病/强直/干燥/Graves）
  × 25 个感染终点（与主矩阵同面板、同 munge、同 SNP 集）两两 bivariate LDSC。
  目的：检验"感染×自免"的 rg 分布与 MHC 占比是否同构于"感染×感染"——
  共享架构=抗原呈递臂主张的第三个证据腿（三种结局都有信息量，见方案文档 §三 E3）。

口径（全部继承主矩阵，不另立口径）：
  1. munge 复用 10_munge_finngen.munge_one（同 1000G frq 对齐 / 同 MAF / 同守卫）
  2. rg 复用 12_matrix_runner v2 的 load_master/align/run_pair（同共用 SNP 集、
     同 n_blocks=200、同回文剔除、同 try/except 逐对 status 不静默、同 |rg|>1.2 出界判 nan）
  3. 功效门同主矩阵：两性状 h² 单侧 p<0.05 且 mean χ²≥1.02
  4. BH 校正：本层 **power_ok** 的 rg_p 内部 BH（与主矩阵同域；2026-10-28 统一）（记 rg_q_bh_layer 列，不与主矩阵 300 对混池）
  5. 病例数/对照数取自官方 manifest（finngen_R12_pheno_n.tsv），禁硬编码

前置：38_download_autoimmune.sh 已完成（gz 落 gwas_data/finngen_r12/）。
用法：python 43_autoimmune_runner.py [--munge-only]
输出：work/output/ldsc_rg_autoimmune_<DATE>.csv + autoimmune_anchor_summary.csv
"""
import csv
import io
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, norm

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATA, OUT, DATE, FINNGEN_ENDPOINTS, GWAS_DATA  # noqa: E402,F401

AUTOIMMUNE = ["T1D", "M13_RHEUMA", "G6_MS", "K11_COELIAC", "K11_UC_STRICT2",
              "K11_CD_STRICT2", "L12_PSORIASIS", "M13_ANKYLOSPON", "M13_SJOGREN",
              "E4_GRAVES_STRICT"]
OUTCSV = OUT / ("ldsc_rg_autoimmune_%s.csv" % DATE)
SUMCSV = OUT / "autoimmune_anchor_summary.csv"


def load_spec(spec_path, name):
    spec = importlib.util.spec_from_file_location(name, spec_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mr = load_spec("/path/to/project/work/code/12_matrix_runner.py", "matrix_runner_20b2")


def munge_autoimmune():
    mg = load_spec("/path/to/project/work/code/10_munge_finngen.py", "munge")
    counts = mr.load_manifest_counts()
    todo = []
    for ph in AUTOIMMUNE:
        if (mr.MUNGED / ("%s.sumstats" % ph)).exists():
            print("[munge] %s 已有 sumstats，跳过" % ph)
            continue
        if ph not in counts:
            raise KeyError("官方 manifest 无 %s" % ph)
        nc, nk = counts[ph]
        todo.append((ph, nc, nk))
    if not todo:
        return
    print("[munge] 加载 frq/panel（一次）...")
    frq = mg.load_1000g_frq()
    panel_ids = set()
    for i in range(1, 23):
        d = pd.read_csv(mg.LDSC_REF / "LDscore" / ("LDscore.%d.l2.ldscore.gz" % i),
                        sep="\t", usecols=["SNP"])
        panel_ids.update(d["SNP"])
    for ph, nc, nk in todo:
        t0 = time.time()
        c = mg.munge_one(ph, nc, nk, frq, panel_ids)
        cnt_path = mr.MUNGED / "munge_counts.csv"
        buf = io.StringIO()
        w = csv.DictWriter(buf, fieldnames=list(c.keys()), delimiter="\t")
        if not cnt_path.exists():
            w.writeheader()
        w.writerow(c)
        with open(cnt_path, "a", encoding="utf-8-sig", newline="") as fh:
            fh.write(buf.getvalue())
        print("[munge] %s 完成 %s SNPs mean_chi2=%.2f (%.0fs)"
              % (ph, c["final_snps"], c["mean_chisq"], time.time() - t0))


def main():
    munge_autoimmune()
    if "--munge-only" in sys.argv:
        print("[anchor] munge-only，退出")
        return

    inf_labels = {lab: ph for lab, ph in FINNGEN_ENDPOINTS.items()
                  if ph not in mr.MTB_BLOCKED and (mr.MUNGED / ("%s.sumstats" % ph)).exists()}
    auto_labels = {ph: ph for ph in AUTOIMMUNE
                   if (mr.MUNGED / ("%s.sumstats" % ph)).exists()}
    missing = [ph for ph in AUTOIMMUNE if ph not in auto_labels]
    if missing:
        print("[anchor] 缺 munge 产物（下载未完成？）：%s" % ", ".join(missing))
    print("[anchor] 感染 %d 性状 × 自免 %d 性状" % (len(inf_labels), len(auto_labels)))

    all_labels = dict(inf_labels)
    all_labels.update(auto_labels)
    # 分批载入（每批 5 自免 + 全部感染 ≈ 主矩阵 25 性状的内存峰值 +20%；16 GB 机器实测安全）
    CHUNK = 5
    auto_sorted = sorted(auto_labels)
    rows = []
    n_total_pairs = len(auto_sorted) * len(inf_labels)
    done = 0
    for c0 in range(0, len(auto_sorted), CHUNK):
        chunk = auto_sorted[c0:c0 + CHUNK]
        batch_labels = dict(inf_labels)
        batch_labels.update({a: auto_labels[a] for a in chunk})
        print("[anchor] 批次 %s（自免 %d + 感染 %d）载入 ..."
              % (",".join(chunk), len(chunk), len(inf_labels)))
        meta, traits, M = mr.load_master(batch_labels)
        print("[anchor] 共用 SNP 集 %d | M_5_50=%.0f" % (len(meta["SNP"]), M))
        for a in chunk:
            for lab2 in sorted(inf_labels):
                row = mr.run_pair(a, lab2, traits[a], traits[lab2], meta, M)
                rows.append(row)
                done += 1
                if done % 25 == 0 or done == n_total_pairs:
                    print("[anchor %d/%d] %s x %s: rg=%s status=%s"
                          % (done, n_total_pairs, a, lab2,
                             ("%.3f" % row["rg"]) if pd.notna(row.get("rg")) else "NA",
                             row["status"]))
        del traits, meta
    df = pd.DataFrame(rows)
    for i in (1, 2):
        se = df["h2_%d_se" % i] if "h2_%d_se" % i in df.columns else np.nan
        df["h2_%d_p" % i] = norm.sf(df["h2_%d" % i] / se) if "h2_%d" % i in df.columns else np.nan
        df["h2_%d_p" % i] = df["h2_%d_p" % i].where(df["h2_%d" % i].notna(), np.nan)
    df["power_ok"] = (df["status"] == "ok") & (df["h2_1_p"] < 0.05) & (df["h2_2_p"] < 0.05) \
        & (df["mean_chisq_1"] >= 1.02) & (df["mean_chisq_2"] >= 1.02)
    # 2026-10-28 修正（五评委核查）：推断域统一为**过功效门**的对，与主矩阵 12 号同规
    # （非过门对按正文声明"excluded from inference"，不应计入 FDR 复现集）。
    okm = df["power_ok"] == True  # noqa: E712
    df["rg_q_bh_layer"] = np.nan
    if okm.any():
        df.loc[okm, "rg_q_bh_layer"] = false_discovery_control(
            df.loc[okm, "rg_p"].values, method="bh")
    df.to_csv(OUTCSV, index=False)

    # ---- 汇总（对照读数：主矩阵感染×感染的同口径数字） ----
    n_ok = int(okm.sum())
    n_power = int(df["power_ok"].sum())
    pw = df[df["power_ok"]]
    n_pos = int((pw["rg"] > 0).sum()) if n_power else 0
    med = float(pw["rg"].median()) if n_power else float("nan")
    n_bh = int(((df["rg_q_bh_layer"] < 0.05).sum()))
    main_matrix = pd.read_csv(OUT / "ldsc_rg" / "ldsc_matrix_2026-10-11.csv")
    mm_pw = main_matrix[main_matrix["power_ok"] == True]  # noqa: E712
    auto_pass = sorted(pw["axis1"].unique()) if n_power else []
    summ = pd.DataFrame([
        {"key": "auto_endpoints_munged", "value": str(len(auto_labels)),
         "note": "10 计划；缺 munge：%s" % (", ".join(missing) or "无")},
        {"key": "n_pairs_attempted", "value": str(len(df)), "note": "自免×感染全对"},
        {"key": "n_status_ok", "value": str(n_ok), "note": "status==ok（可解释 rg）"},
        {"key": "n_power_ok", "value": str(n_power),
         "note": "功效门（两性状 h² 单侧 p<0.05 且 mean χ²≥1.02，同主矩阵口径）"},
        {"key": "n_rg_positive", "value": "%d/%d" % (n_pos, n_power), "note": "过门对中 rg>0"},
        {"key": "median_rg_power_ok", "value": ("%.3f" % med) if n_power else "NA",
         "note": "过门对 rg 中位（对照：主矩阵感染×感染 BH 集中位 0.531）"},
        {"key": "n_bh_layer", "value": str(n_bh), "note": "本层 BH<0.05（rg_q_bh_layer 列）"},
        {"key": "main_matrix_power_ok_ref", "value": str(int(mm_pw.shape[0])),
         "note": "主矩阵感染×感染过门对数（对照）"},
        {"key": "auto_endpoints_passing_gate", "value": ", ".join(auto_pass) or "无",
         "note": "过功效门的自免终点（其余如实报告为不可估计）"},
    ])
    summ.to_csv(SUMCSV, index=False)
    print("[anchor] 全对 %d / ok %d / 过门 %d / rg>0 %d / BH %d / 中位 rg %.3f"
          % (len(df), n_ok, n_power, n_pos, n_bh, med))
    print("[anchor] 落盘 %s / %s" % (OUTCSV.name, SUMCSV.name))


if __name__ == "__main__":
    main()
