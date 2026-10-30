# -*- coding: utf-8 -*-
""" 48_intercept_qc.py — LDSC 截距 QC 汇总（E8，2026-10-25）

问题：样本重叠（FinnGen 共享对照）与群体分层的混杂此前只有文字声明；
读者标准问句是"截距如何？"。gencov 截距 = cross-trait 回归 z1·z2 的截距，
**零假设为 0**（区别于 h² 截距的零假设 1）；显著为正即提示重叠/分层膨胀。

口径：
  - 三个矩阵（主矩阵 / 自免锚层 / 对照层）的 status=ok 行的 int_gencov /
    int_hsq_1 / int_hsq_2 分布（中位、IQR、最大）；
  - 显著（BH<0.05）子集单独给 max（读者真正关心的是"显著对是否被膨胀抬高"）；
  - 跨源对（HGI×FinnGen，完全独立样本）单独给中位——应最接近 0（内部正对照）。
用法：python 48_intercept_qc.py
输出：work/output/intercept_qc_summary.csv
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATE, OUT  # noqa: E402

MATRICES = [
    ("main", OUT / "ldsc_rg" / ("ldsc_matrix_%s.csv" % DATE)),
    ("autoimmune", OUT / ("ldsc_rg_autoimmune_%s.csv" % DATE)),
    ("comparator", OUT / ("ldsc_rg_comparator_%s.csv" % DATE)),
]
OUTCSV = OUT / "intercept_qc_summary.csv"


def qstats(v):
    if not len(v):
        return {"n": 0, "median": np.nan, "q25": np.nan, "q75": np.nan, "max": np.nan}
    return {"n": len(v), "median": float(np.median(v)),
            "q25": float(np.percentile(v, 25)), "q75": float(np.percentile(v, 75)),
            "max": float(np.max(v))}


def main():
    rows = []
    for tag, path in MATRICES:
        if not path.exists():
            print("[intqc][WARN] 缺矩阵：%s" % path.name)
            continue
        d = pd.read_csv(path)
        ok = d[d["status"] == "ok"]
        if "int_gencov" not in ok.columns:
            raise SystemExit("[intqc] %s 无 int_gencov 列（run_pair 未含截距增补？）" % path.name)
        g = ok["int_gencov"].dropna().values
        bh = ok[ok.get("rg_q_bh_layer", pd.Series(dtype=float)).notna() \
                & (ok["rg_q_bh_layer"] < 0.05)] if "rg_q_bh_layer" in ok.columns \
            else ok[ok.get("rg_q_bh", pd.Series(dtype=float)).notna() & (ok["rg_q_bh"] < 0.05)]
        for name, s in (("all_ok", qstats(g)),
                        ("bh05", qstats(bh["int_gencov"].dropna().values))):
            rows.append({"matrix": tag, "subset": name, "stat_gencov_intercept": "",
                         **{("int_gencov_" + k): v for k, v in s.items()}})
        # 跨源子集：HGI 外部层对由 13 号自有循环产出、暂无截距列（50 行均空），
        # 该子集恒为空，不再产出 NaN 行（2026-10-27 修：空子集产出 NaN 行属规则 220 家族形态）。
        for side in ("int_hsq_1", "int_hsq_2"):
            s = qstats(ok[side].dropna().values)
            rows.append({"matrix": tag, "subset": "all_ok", "stat_gencov_intercept": side,
                         **{(side + "_" + k): v for k, v in s.items()}})
    df = pd.DataFrame(rows)
    df.to_csv(OUTCSV, index=False)
    for _, r in df[df["subset"] != "all_ok"].iterrows():
        print("[intqc] %-11s %-13s gencov int 中位 %.4f 最大 %.4f (n=%s)" %
              (r["matrix"], r["subset"], r["int_gencov_median"], r["int_gencov_max"], r["int_gencov_n"]))


if __name__ == "__main__":
    main()
