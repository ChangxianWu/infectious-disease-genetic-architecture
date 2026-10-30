# -*- coding: utf-8 -*-
""" 42_prs_transfer.py — E5 PRS 传递矩阵（派生层，2026-10-25）

把 rg 架构翻译成临床可读量：用 A 终点 GWAS 建立的（无限小模型收缩）PRS 预测
B 终点的期望 R² = rg²(A,B) × h²(B)（观测尺度；Dudbridge 2013 PLoS Genet
多基因评分预测准确率公式在 rg≠1 时的推广）。纯派生：不引入新数据、不训练模型。

口径（诚实限定，写入输出 note）：
  1. 观测尺度（表 1 h² 与矩阵 rg 同尺度，未做患病率换算 → 保守）
  2. 只报 BH 存活（rg_q_bh<0.05）且 power_ok=True 的对，双向各算一行
  3. 无个体数据 → 这是架构上限（无限小模型、或acles 权重），不是实测预测性能；
     用途 = 给 rg 一个"临床量级"的解读锚，防读者把 rg=0.5 误读为强可预测
用法：python 42_prs_transfer.py
输出：work/output/prs_transfer.csv + prs_transfer_summary.csv
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/path/to/project/work/code")
from constants import OUT  # noqa: E402


def main():
    m = pd.read_csv(OUT / "ldsc_rg" / "ldsc_matrix_2026-10-11.csv")
    t1 = pd.read_csv(OUT / "table1.csv")
    h2 = t1.set_index("Trait")["h2_observed"].to_dict()

    sig = m[(m["power_ok"] == True) & (m["rg_q_bh"].notna()) & (m["rg_q_bh"] < 0.05)]  # noqa: E712
    recs = []
    for _, r in sig.iterrows():
        for src, dst in ((r["axis1"], r["axis2"]), (r["axis2"], r["axis1"])):
            if dst not in h2:
                continue
            r2 = (r["rg"] ** 2) * h2[dst]
            recs.append({"prs_from": src, "predicted_endpoint": dst, "rg": r["rg"],
                         "rg_q_bh": r["rg_q_bh"], "h2_observed_target": h2[dst],
                         "expected_R2": r2, "same_axis": r["axis1_axis_same"] if
                         "axis1_axis_same" in r.index else ""})
    out = pd.DataFrame(recs).sort_values("expected_R2", ascending=False)
    out.to_csv(OUT / "prs_transfer.csv", index=False)

    within = out[out["prs_from"].map(lambda a: _axis(a)) == out["predicted_endpoint"].map(lambda a: _axis(a))]
    cross = out[out["prs_from"].map(lambda a: _axis(a)) != out["predicted_endpoint"].map(lambda a: _axis(a))]
    top = out.head(5)
    summ = pd.DataFrame([
        {"key": "n_bh_pairs", "value": str(len(sig)), "note": "BH 存活对（双向展开 %d 行）" % len(out)},
        {"key": "max_expected_R2", "value": "%.5f" % out["expected_R2"].max(),
         "note": "最大传递：%s→%s" % (out.iloc[0]["prs_from"], out.iloc[0]["predicted_endpoint"])
         if len(out) else ""},
        {"key": "median_R2_within_axis", "value": "%.5f" % within["expected_R2"].median()
         if len(within) else "", "note": "同轴传递 R² 中位（%d 行）" % len(within)},
        {"key": "median_R2_cross_axis", "value": "%.5f" % cross["expected_R2"].median()
         if len(cross) else "", "note": "跨轴传递 R² 中位（%d 行）" % len(cross)},
        {"key": "top5", "value": "；".join("%s→%s R²=%.4f%%" % (r["prs_from"], r["predicted_endpoint"],
                                                               100 * r["expected_R2"])
                                           for _, r in top.iterrows()),
         "note": "最大五个传递对（百分数口径）"},
        {"key": "scale_note", "value": "observed", "note":
         "观测尺度（未做患病率换算，保守）；无限小模型架构上限，非实测预测性能"},
    ])
    summ.to_csv(OUT / "prs_transfer_summary.csv", index=False)
    print("[prs] BH 对 %d（%d 行）；最大 R² %.5f（%s→%s）" %
          (len(sig), len(out), out["expected_R2"].max(),
           out.iloc[0]["prs_from"], out.iloc[0]["predicted_endpoint"]))
    print("[prs] 同轴中位 %.5f / 跨轴中位 %.5f" %
          (within["expected_R2"].median() if len(within) else float("nan"),
           cross["expected_R2"].median() if len(cross) else float("nan")))


def _axis(trait):
    if not hasattr(_axis, "map"):
        t1 = pd.read_csv(OUT / "table1.csv")
        _axis.map = t1.set_index("Trait")["Pathogen_axis"].to_dict()
    return _axis.map.get(trait, "")


if __name__ == "__main__":
    main()
