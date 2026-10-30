# -*- coding: utf-8 -*-
""" 50_axis_robustness.py — 核心零结果主张的稳健性套件（E7，2026-10-25/27 改号）

编号说明（2026-10-27）：原编号 47 与并行会话的 47_ebiomedicine_frontmatter.py 撞号；
本脚本顺延为 50 号，输出文件名不变（axis_robustness.csv / axis_robustness_summary.csv）。

主张：轴内（bacterial_internal）与跨轴（cross_axis）的 rg 无差异。
解读风险：正文只说 "not statistically distinguishable" 未给检验；且 45 对不独立
（星形拓扑：高 h² 终点同时出现在两侧），单次 Mann-Whitney 可能被单个高连通终点驱动。

三项检验（全部只读冻结矩阵，不重跑 LDSC）：
  T1 基础检验：Mann-Whitney U（双侧）——轴内 vs 跨轴的 rg 分布；
  T2 留一性状：逐一剔除每个过门性状，重算两分布中位数与 MW p——
     若"不显著"在全部留一下稳定，主张不被单点驱动；
  T3 功效协变量：以 rg 为响应、stratum + 两性状 mean_chisq 均值为协变量的
     OLS（复权重 rg_se⁻²），看 stratum 系数——功效差异能否解释表观差异；
     （LDSC 理论上无显著共享时 rg 无偏，但低 mean-χ² 对的 rg_se 大，
     加权回归使"表观中位差"不依赖方差齐性假设）
口径：分层定义与 13 号脚本一致（bacterial_internal = 两性状均为细菌轴；
cross_axis = 至少一侧为病毒轴；mixed 轴的 8 对细菌-混合对按正文口径**排除在两分布外**，
因为混合终点=全部肺炎、含细菌性肺炎，定义重叠）。
实测（2026-10-27）：T1 p=0.896；T2 八轮最小 p=0.405、0 轮显著；T3 β=-0.362（p=0.103）。
用法：python 50_axis_robustness.py
输出：work/output/axis_robustness.csv + axis_robustness_summary.csv
"""
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, "/path/to/project/work/code")
from constants import OUT  # noqa: E402

OUTCSV = OUT / "axis_robustness.csv"
SUMCSV = OUT / "axis_robustness_summary.csv"


def load_strata():
    m = pd.read_csv(OUT / "ldsc_rg" / "ldsc_matrix_2026-10-11.csv")
    t1 = pd.read_csv(OUT / "table1.csv")
    AX = t1.set_index("Trait")["Pathogen_axis"].to_dict()
    rows = []
    for _, r in m[m["power_ok"] == True].iterrows():  # noqa: E712
        a1, a2 = AX.get(r["axis1"], ""), AX.get(r["axis2"], "")
        pair = {r["axis1"], r["axis2"]}
        if "PNEUMONIA_MIXED" in pair:
            continue  # 正文口径：细菌-混合 8 对单独讨论，不入两分布
        # table1 的 Pathogen_axis 是首字母大写（Bacterial/Mixed/Viral）；统一小写比较
        a1, a2 = a1.lower(), a2.lower()
        if a1 == "bacterial" and a2 == "bacterial":
            strat = "bacterial_internal"
        elif "viral" in (a1, a2):
            strat = "cross_axis"
        else:
            continue
        rows.append({"axis1": r["axis1"], "axis2": r["axis2"], "stratum": strat,
                     "rg": float(r["rg"]), "rg_se": float(r["rg_se"]),
                     "chi_mean": (float(r["mean_chisq_1"]) + float(r["mean_chisq_2"])) / 2})
    d = pd.DataFrame(rows)
    d.to_csv(OUTCSV, index=False)
    return d


def mw(d):
    a = d[d["stratum"] == "bacterial_internal"]["rg"]
    b = d[d["stratum"] == "cross_axis"]["rg"]
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    return float(a.median()), float(b.median()), float(p)


def main():
    d = load_strata()
    med_a, med_b, p0 = mw(d)
    summ = [{"key": "n_pairs_analysed", "value": str(len(d)),
             "note": "过门 45 对剔除 8 对细菌-混合后入两分布（正文口径）"},
            {"key": "t1_mw_p", "value": "%.3f" % p0, "note": "T1 Mann-Whitney 双侧（轴内 vs 跨轴）"},
            {"key": "t1_median_internal", "value": "%.3f" % med_a, "note": "轴内中位"},
            {"key": "t1_median_cross", "value": "%.3f" % med_b, "note": "跨轴中位"}]

    # T2 留一性状
    traits = sorted(set(d["axis1"]) | set(d["axis2"]))
    loo = []
    for t in traits:
        sub = d[(d["axis1"] != t) & (d["axis2"] != t)]
        if sub["stratum"].nunique() < 2 or (sub["stratum"] == "bacterial_internal").sum() < 3 \
                or (sub["stratum"] == "cross_axis").sum() < 3:
            continue
        ma, mb, p = mw(sub)
        loo.append({"left_out": t, "n": len(sub), "median_internal": ma,
                    "median_cross": mb, "mw_p": p})
    loo = pd.DataFrame(loo)
    n_sig = int((loo["mw_p"] < 0.05).sum())
    max_shift = float((loo["median_cross"] - loo["median_internal"]).abs().max())
    summ += [{"key": "t2_loo_iterations", "value": str(len(loo)), "note": "留一性状轮数"},
             {"key": "t2_loo_min_p", "value": "%.3f" % loo["mw_p"].min(),
              "note": "留一下最小 MW p（<0.05 的轮数 = %d）" % n_sig},
             {"key": "t2_loo_n_p_lt005", "value": str(n_sig),
              "note": "留一后变显著的轮数（0 = 主张不被单点驱动）"},
             {"key": "t2_loo_max_median_shift", "value": "%.3f" % max_shift,
              "note": "留一下两分布中位差的最大绝对偏移"}]

    # T3 功效协变量（逆方差加权 OLS：rg ~ stratum + mean_chisq）
    X = np.column_stack([
        np.ones(len(d)),
        (d["stratum"] == "cross_axis").astype(float).values,
        d["chi_mean"].values,
    ])
    w = 1.0 / np.square(d["rg_se"].values)
    XtW = X.T * w
    beta = np.linalg.lstsq(XtW @ X, XtW @ d["rg"].values, rcond=None)[0]
    resid = d["rg"].values - X @ beta
    dof = len(d) - X.shape[1]
    sigma2 = float((w * resid ** 2).sum() / dof)
    cov = sigma2 * np.linalg.inv(XtW @ X)
    se_strat = float(np.sqrt(cov[1, 1]))
    z_strat = beta[1] / se_strat
    from scipy.stats import norm
    summ += [{"key": "t3_wls_stratum_beta", "value": "%.3f" % beta[1],
              "note": "跨轴虚拟变量的加权回归系数（正=跨轴更高）"},
             {"key": "t3_wls_stratum_se", "value": "%.3f" % se_strat, "note": "同上 SE"},
             {"key": "t3_wls_stratum_p", "value": "%.3f" % (2 * norm.sf(abs(z_strat))),
              "note": "同上 p（调整 mean chisq 后）"}]

    pd.DataFrame(summ).to_csv(SUMCSV, index=False)
    print("[axis] T1 MW p=%.3f（轴内 %.3f vs 跨轴 %.3f）" % (p0, med_a, med_b))
    print("[axis] T2 留一 %d 轮：最小 p=%.3f，变显著 %d 轮" % (len(loo), loo["mw_p"].min(), n_sig))
    print("[axis] T3 加权回归跨轴 β=%.3f (SE %.3f, p=%.3f)" % (beta[1], se_strat, 2 * norm.sf(abs(z_strat))))


if __name__ == "__main__":
    main()
