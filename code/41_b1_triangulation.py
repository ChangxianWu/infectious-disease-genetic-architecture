# -*- coding: utf-8 -*-
""" 41_b1_triangulation.py — E4 与 20B1 蛋白层（pQTL-MR）三角互证（2026-10-25）

逻辑：两稿共享同一批 FinnGen R12 感染终点。20B1 在蛋白组层（cis-pQTL 工具变量）
对若干终点找到 MR 因果信号；project 在汇总统计层给出终点间 rg 连通度。
若"蛋白层有信号"的终点恰是"rg 网络中连通度高"的终点，则两层证据互证同一
感染易感生物学核心（非同一分析的两个名字）。

口径：
  1. 20B1 侧（只读其产物，不重算其多重校正）：
     - sens_multiple_testing_by_endpoint.csv：endpoint × n_tests × n_bh05
     - mr_all_endpoints.csv：端点级 min(pval_out)（全部蛋白-工具检验的原始最小 p）
     - s7_canonical_numbers.tsv：20B1 canonical 六个 BH 存活对（hit_*_pbh 键），
       逐条列出 protein × endpoint（来源=20B1 冻结数字，不复制方法论）
  2. project 侧：ldsc_matrix_2026-10-11.csv 中 power_ok=True 的行 → 每端点 mean rg
     （端点连通度）+ 该端点参与过的 BH 存活对数
  3. 统计：两表共有端点上 Spearman(20B1 -log10(min p), project mean rg)
     + "20B1 BH 端点在 project 连通度排序中的位置"
  4. 端点码映射：20B1 PNEUMONIA ↔ project PNEUMONIA_MIXED（同一 FinnGen phenocode
     J10_PNEUMONIA；以 20B1 mr_main 文件名+project table1 双侧核对后映射，非记忆）

用法：python 41_b1_triangulation.py
输出：work/output/b1_triangulation.csv（端点级）+ b1_triangulation_summary.csv
"""
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, "/path/to/project/work/code")
from constants import OUT  # noqa: E402

B1 = "/path/to/project/20B1/work/output"
CODE_MAP = {"PNEUMONIA": "PNEUMONIA_MIXED"}  # J10_PNEUMONIA 两稿码名差异（见 docstring 4）


def main():
    # ---- 20B1 侧 ----
    mt = pd.read_csv(B1 + "/sens_multiple_testing_by_endpoint.csv")
    allmr = pd.read_csv(B1 + "/mr_all_endpoints.csv", usecols=["endpoint", "pval_out"],
                        low_memory=False)
    minp = allmr.groupby("endpoint")["pval_out"].min().rename("b1_min_p")
    b1 = mt.merge(minp, left_on="endpoint", right_index=True, how="left")
    canon = pd.read_csv(B1 + "/s7_canonical_numbers.tsv", sep="\t",
                        names=["key", "value", "note"], skiprows=1)
    # 20B1 canonical 对每个 BH 对登记 5 个键（hit_<P>_<E>_{OR,CI,p,pbh,ncase}）——
    # 只取 _pbh 键才是"对数"（首版取全部 hit_ 键 → 30 个键被当成 30 对，2026-10-25 修）
    bh_hits = canon[canon["key"].str.startswith("hit_")
                    & canon["key"].str.endswith("_pbh")].copy()
    bh_hits["pair"] = bh_hits["key"].str.replace("hit_", "", regex=False).str.replace(
        "_pbh$", "", regex=True)
    bh_hits[["protein", "endpoint_b1"]] = bh_hits["pair"].str.rsplit("_", n=1, expand=True)
    bh_pairs = [(r["protein"], r["endpoint_b1"]) for _, r in bh_hits.iterrows()]
    b1_endpoints_bh = sorted({CODE_MAP.get(e, e) for _, e in bh_pairs})
    b1["endpoint_b2"] = b1["endpoint"].map(CODE_MAP).fillna(b1["endpoint"])

    # ---- project 侧 ----
    m = pd.read_csv(OUT / "ldsc_rg" / "ldsc_matrix_2026-10-11.csv")
    ok = m[m["power_ok"] == True]  # noqa: E712
    conn = {}
    bh_count = {}
    for _, r in ok.iterrows():
        for a, b in ((r["axis1"], r["axis2"]), (r["axis2"], r["axis1"])):
            conn.setdefault(a, []).append(r["rg"])
            if pd.notna(r.get("rg_q_bh")) and r["rg_q_bh"] < 0.05:
                bh_count[a] = bh_count.get(a, 0) + 1
    b2 = pd.DataFrame({"endpoint_b2": sorted(conn),
                       "b2_mean_rg_gate_ok": [np.mean(conn[e]) for e in sorted(conn)],
                       "b2_n_bh_pairs": [bh_count.get(e, 0) for e in sorted(conn)]})

    # ---- 连接与统计 ----
    j = b1.merge(b2, on="endpoint_b2", how="inner")
    # B1 端点级强度：min p 大量下溢为 0（<1e-300，超双精度）→ 用端点分层 BH 存活数
    # 作强度统计（scope=20B1 sens_multiple_testing_by_endpoint.csv 的端点内校正，
    # 与 20B1 canonical 的全局 BH 6 对不同口径，两者并报不混用）
    rho, p_sp = spearmanr(j["n_bh05"], j["b2_mean_rg_gate_ok"])
    j["b1_has_bh_hit"] = j["endpoint"].isin(b1_endpoints_bh) | j["endpoint_b2"].isin(b1_endpoints_bh)
    j["b2_connectivity_rank"] = j["b2_mean_rg_gate_ok"].rank(ascending=False).astype(int)
    j.to_csv(OUT / "b1_triangulation.csv", index=False)

    bh_rows = j[j["b1_has_bh_hit"]]
    summ = pd.DataFrame([
        {"key": "b1_endpoints_total", "value": str(len(b1)),
         "note": "20B1 受检端点数（sens_multiple_testing_by_endpoint.csv）"},
        {"key": "b1_bh_survivors_canonical", "value": str(len(bh_pairs)),
         "note": "全局 BH 口径（20B1 canonical）；" + "；".join("%s×%s" % (p, e) for p, e in bh_pairs)},
        {"key": "b1_strength_stat", "value": "n_bh05",
         "note": "端点分层 BH 存活数（20B1 sens 文件 scope）；min-p 下溢为 0 不可用"},
        {"key": "joined_endpoints", "value": str(len(j)),
         "note": "两稿共有且 project 功效门通过的端点数（PNEUMONIA→PNEUMONIA_MIXED 映射后）；"
                 "20B1 的 ENDOCARD/MENINGBACT/COVID19 因未过 project 功效门（无可估 rg）不入连接，如实记录"},
        {"key": "unmatched_b1_endpoints", "value": ", ".join(sorted(set(b1["endpoint"]) - set(j["endpoint"]))),
         "note": "20B1 端点不在 project 连通度表内（未过功效门或面板外）"},
        {"key": "spearman_nbh05_vs_meanrg_rho", "value": ("%.3f" % rho) if pd.notna(rho) else "NA",
         "note": "端点级 20B1 n_bh05 × project mean-rg 连通度；Spearman"},
        {"key": "spearman_nbh05_vs_meanrg_p", "value": ("%.3g" % p_sp) if pd.notna(p_sp) else "NA",
         "note": "同上"},
        {"key": "bh_endpoints_all_in_b2_core", "value": str(bool(len(bh_rows)) and
         (bh_rows["b2_n_bh_pairs"] >= 1).all()),
         "note": "20B1 BH 端点是否全部落在 project BH 共享网络内（定性互证）"},
        {"key": "bh_endpoints_mean_rg", "value": ("%.3f" % bh_rows["b2_mean_rg_gate_ok"].mean())
         if len(bh_rows) else "", "note": "20B1 BH 端点的 project 平均 rg"},
    ])
    summ.to_csv(OUT / "b1_triangulation_summary.csv", index=False)
    print("[tri] 共有端点 %d；Spearman ρ=%s (p=%s)" %
          (len(j), ("%.3f" % rho) if pd.notna(rho) else "NA",
           ("%.3g" % p_sp) if pd.notna(p_sp) else "NA"))
    if len(bh_rows):
        print("[tri] 20B1 BH 端点全部在 project 共享网络内：%s；mean rg %.3f" %
              (bool((bh_rows["b2_n_bh_pairs"] >= 1).all()), bh_rows["b2_mean_rg_gate_ok"].mean()))


if __name__ == "__main__":
    main()
