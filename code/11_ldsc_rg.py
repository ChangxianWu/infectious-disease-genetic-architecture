# -*- coding: utf-8 -*-
""" 11_ldsc_rg.py — bivariate LDSC 驱动（py3 移植 regressions.RG）
（Stage 1，2026-10-11）

数据流（与官方 ldsc.py --rg 一致）：
  x = 官方 l2 LD scores（1000G EUR，1.19M SNPs）∩ weights.hm3_noMHC（w_ld）
  y = z1*z2（等位对齐后）
  M = 官方 .l2.M_5_50 总和（22 染色体；rg 为比值，M 选择约掉）
  回归 SNP 集 = 两 sumstats ∩ l2 ∩ w_ld（MHC 由 noMHC weights 自动排除）
  RG 内部：Hsq(z1²) + Hsq(z2²) + Gencov(z1,z2)（固定 h1/h2/截距）→ rg = gencov/√(h1h2)
  jackknife n_blocks=200（官方默认）

用法：python 11_ldsc_rg.py TRAIT1 TRAIT2 [OUT_SUFFIX]
输出：work/output/ldsc_rg/<T1>_<T2>[suffix].csv（rg/se/p/z/h²/截距/mean χ²/n_snp）
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/path/to/project/work/code")
sys.path.insert(0, "/path/to/project/work/code/ldsc_py3")
from constants import DATA, OUT, GWAS_DATA  # noqa: E402
import regressions as reg  # noqa: E402  (移植版)

LDSC_REF = GWAS_DATA / "ldsc_ref"
RGDIR = OUT / "ldsc_rg"
RGDIR.mkdir(parents=True, exist_ok=True)
MUNGED = DATA / "munged"


def load_ref():
    l2s, wls, m550 = [], [], []
    for i in range(1, 23):
        l2s.append(pd.read_csv(LDSC_REF / "LDscore" / ("LDscore.%d.l2.ldscore.gz" % i),
                               sep="\t", usecols=["SNP", "L2"]))
        wls.append(pd.read_csv(LDSC_REF / "1000G_Phase3_weights_hm3_no_MHC"
                               / ("weights.hm3_noMHC.%d.l2.ldscore.gz" % i),
                               sep="\t", usecols=["SNP", "L2"]))
        m550.append(float(open(LDSC_REF / "LDscore" / ("LDscore.%d.l2.M_5_50" % i)).read().strip()))
    l2 = pd.concat(l2s, ignore_index=True).drop_duplicates("SNP", keep="first")
    wl = pd.concat(wls, ignore_index=True).drop_duplicates("SNP", keep="first")
    return l2, wl, float(np.sum(m550))


def run_rg(t1, t2, l2, wl, M, n_blocks=200, shuffle2=False, rng=None):
    s1 = pd.read_csv(MUNGED / ("%s.sumstats" % t1), sep="\t",
                     usecols=["SNP", "A1", "A2", "N", "Z"])
    s2 = pd.read_csv(MUNGED / ("%s.sumstats" % t2), sep="\t",
                     usecols=["SNP", "A1", "A2", "N", "Z"])
    m = s1.merge(s2, on="SNP", suffixes=("_1", "_2")).merge(l2, on="SNP") \
          .merge(wl, on="SNP", suffixes=("_l2", "_w"))
    # 等位对齐（A1_2==A1_1 同向；A1_2==A2_1 反向翻转 z2）
    same = m["A1_2"] == m["A1_1"]
    flip = m["A1_2"] == m["A2_1"]
    m = m[same | flip].copy()
    m.loc[m.index[same[m.index].values], "z2"] = m.loc[m.index[same[m.index].values], "Z_2"]
    m.loc[m.index[flip[m.index].values], "z2"] = -m.loc[m.index[flip[m.index].values], "Z_2"]
    m = m.dropna(subset=["Z_1", "z2"])
    n_snp = len(m)
    if shuffle2:
        m["z2"] = rng.permutation(m["z2"].values)
    z1 = m["Z_1"].values.reshape((n_snp, 1))
    z2 = m["z2"].values.reshape((n_snp, 1))
    x = m["L2_l2"].values.reshape((n_snp, 1))
    w = m["L2_w"].values.reshape((n_snp, 1))
    N1 = m["N_1"].values.reshape((n_snp, 1)).astype(float)
    N2 = m["N_2"].values.reshape((n_snp, 1)).astype(float)
    t0 = time.time()
    rghat = reg.RG(z1, z2, x, w, N1, N2, np.array([[M]]), n_blocks=n_blocks)
    runtime = time.time() - t0
    row = {
        "trait1": t1, "trait2": t2, "n_snps": n_snp, "M_5_50": M,
        "rg": float(rghat.rg_ratio) if rghat._negative_hsq is None else np.nan,
        "rg_se": float(rghat.rg_se) if rghat._negative_hsq is None else np.nan,
        "rg_p": float(rghat.p) if rghat._negative_hsq is None else np.nan,
        "rg_z": float(rghat.z) if rghat._negative_hsq is None else np.nan,
        "h2_1": float(rghat.hsq1.tot), "h2_1_se": float(rghat.hsq1.tot_se),
        "h2_2": float(rghat.hsq2.tot), "h2_2_se": float(rghat.hsq2.tot_se),
        "intercept_1": float(rghat.hsq1.intercept), "intercept_2": float(rghat.hsq2.intercept),
        "gencov_intercept": float(rghat.gencov.intercept),
        "mean_chisq_1": float(rghat.hsq1.mean_chisq), "mean_chisq_2": float(rghat.hsq2.mean_chisq),
        "runtime_s": round(runtime, 1), "shuffled": bool(shuffle2),
    }
    return row, rghat


def main():
    t1, t2 = sys.argv[1], sys.argv[2]
    suffix = sys.argv[3] if len(sys.argv) > 3 else ""
    print("[rg] 加载官方参考 ...")
    l2, wl, M = load_ref()
    print("[rg] l2=%d w_ld=%d M_5_50=%.0f" % (len(l2), len(wl), M))
    row, rghat = run_rg(t1, t2, l2, wl, M)
    print("[rg] %s x %s: rg=%.4f (%.4f) p=%.3e | h2=%.4f/%.4f | n_snps=%d | %.0fs"
          % (t1, t2, row["rg"], row["rg_se"], row["rg_p"], row["h2_1"], row["h2_2"],
             row["n_snps"], row["runtime_s"]))
    out = RGDIR / ("%s_%s%s.csv" % (t1, t2, suffix))
    pd.DataFrame([row]).to_csv(out, index=False)
    print("[rg] -> %s" % out.name)
    return row


if __name__ == "__main__":
    main()
