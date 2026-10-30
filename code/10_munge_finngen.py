# -*- coding: utf-8 -*-
""" 10_munge_finngen.py — FinnGen R12 汇总统计 → LDSC 官方 .sumstats 格式
（Stage 1，2026-10-11）

口径（全部计数落盘，规则 51/71）：
  输入列：#chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls
  1. rsID 对齐官方 l2 面板（1000G EUR，1,190,321 SNPs；panel 外丢弃，计数）
  2. 等位对齐 1000G frq（A1/A2 同向保留；反向翻转 beta 符号；双等位不一致丢弃，计数）
  3. MAF ∈ [0.01, 0.99]（1000G frq 口径，与官方 munge --maf 0.01 一致；计数）
  4. z = beta/sebeta（对齐后符号）；CHISQ = z²；se<=0/非有限丢弃（计数）
  5. p-Z 一致性守卫：chi2_z=(z)² vs chi2_p=chisq.isf(pval)；相对差 >30% 丢弃（计数）
  6. 重复 rsID：保留首条（计数）
  7. N = num_case + num_control（LDSC 官方 FAQ 口径：case-control 用总样本量）
  输出：work/data/munged/<phenocode>.sumstats（SNP A1 A2 N Z CHISQ）+ munge_counts.csv

用法：python 10_munge_finngen.py PHENOCODE CASE_N CONTROL_N
"""
import csv
import gzip
import io
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATA, FINNGEN, GWAS_DATA  # noqa: E402

LDSC_REF = GWAS_DATA / "ldsc_ref"
MUNGED = DATA / "munged"
MUNGED.mkdir(parents=True, exist_ok=True)

USECOLS = ["rsids", "ref", "alt", "pval", "beta", "sebeta", "af_alt"]


def load_1000g_frq():
    """1000G EUR frq 合并（22 染色体）：rsid → (a1, a2, maf)。"""
    frames = []
    for i in range(1, 23):
        f = LDSC_REF / ("1000G_Phase3_frq/1000G.EUR.QC.%d.frq" % i)
        frames.append(pd.read_csv(f, sep=r"\s+", usecols=["SNP", "A1", "A2", "MAF"]))
    frq = pd.concat(frames, ignore_index=True)
    frq = frq.drop_duplicates(subset="SNP", keep="first")
    return frq


def munge_one(phenocode, case_n, control_n, frq, panel_ids):
    t0 = time.time()
    gz = FINNGEN / ("finngen_R12_%s.gz" % phenocode)
    if not gz.exists():
        raise FileNotFoundError(gz)
    counts = {"phenocode": phenocode, "N": int(case_n) + int(control_n),
              "num_case": case_n, "num_control": control_n}
    chunks = []
    total = 0
    for ch in pd.read_csv(gz, sep="\t", usecols=USECOLS, chunksize=4_000_000,
                          dtype={"rsids": str, "ref": str, "alt": str,
                                 "beta": np.float64, "sebeta": np.float64,
                                 "pval": np.float64, "af_alt": np.float64},
                          compression="gzip"):
        total += len(ch)
        ch = ch[ch["rsids"].isin(panel_ids)]
        chunks.append(ch)
    counts["total_variants"] = total
    counts["in_l2_panel"] = sum(len(c) for c in chunks)
    dat = pd.concat(chunks, ignore_index=True)
    del chunks
    # 多 rsID 串（FinnGen 个别行 'rs1;rs2'）：取第一个
    multi = dat["rsids"].str.contains(";", na=False)
    counts["multi_rsids"] = int(multi.sum())
    dat["rsids"] = dat["rsids"].str.split(";").str[0]
    n_before = len(dat)
    dat = dat.drop_duplicates(subset="rsids", keep="first")
    counts["dup_rsids_dropped"] = n_before - len(dat)

    m = dat.merge(frq, left_on="rsids", right_on="SNP", how="inner")
    counts["merged_1000g"] = len(m)

    ref, alt = m["ref"].str.upper(), m["alt"].str.upper()
    a1, a2 = m["A1"].str.upper(), m["A2"].str.upper()
    same = (alt == a1) & (ref == a2)
    flip = (alt == a2) & (ref == a1)
    keep_idx = same | flip
    m = m[keep_idx].copy()
    counts["allele_mismatch_dropped"] = counts["merged_1000g"] - len(m)
    sign = np.where(same[m.index], 1.0, np.where(flip[m.index], -1.0, np.nan))
    m["z"] = sign * (m["beta"] / m["sebeta"])
    m["A1"] = a1[m.index]
    m["A2"] = a2[m.index]

    ok = np.isfinite(m["z"]) & (m["sebeta"] > 0)
    counts["se_or_nan_dropped"] = int((~ok).sum())
    m = m[ok].copy()

    maf_ok = (m["MAF"] >= 0.01) & (m["MAF"] <= 0.99)
    counts["maf_out_dropped"] = int((~maf_ok).sum())
    m = m[maf_ok].copy()

    m["chisq_z"] = m["z"] ** 2
    m["chisq_p"] = chi2_dist.isf(m["pval"].clip(1e-300, 1.0), 1)
    rel = (m["chisq_z"] - m["chisq_p"]).abs() / np.maximum(m[["chisq_z", "chisq_p"]].max(axis=1), 1.0)
    pz_ok = rel <= 0.30
    counts["pz_mismatch_dropped"] = int((~pz_ok).sum())
    m = m[pz_ok].copy()

    out = pd.DataFrame({
        "SNP": m["rsids"], "A1": m["A1"], "A2": m["A2"],
        "N": counts["N"], "Z": m["z"], "CHISQ": m["chisq_z"],
    })
    out_path = MUNGED / ("%s.sumstats" % phenocode)
    out.to_csv(out_path, sep="\t", index=False)
    counts["final_snps"] = len(out)
    counts["mean_chisq"] = float(np.mean(out["CHISQ"]))
    counts["runtime_s"] = round(time.time() - t0, 1)
    return counts


def main():
    phenocode, case_n, control_n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    print("[munge] 加载 1000G frq 与 l2 面板 ...")
    frq = load_1000g_frq()
    panel = pd.read_csv(LDSC_REF / "LDscore" / "LDscore.1.l2.ldscore.gz", sep="\t", nrows=5)
    panel_ids = set()
    for i in range(1, 23):
        d = pd.read_csv(LDSC_REF / "LDscore" / ("LDscore.%d.l2.ldscore.gz" % i), sep="\t",
                        usecols=["SNP"])
        panel_ids.update(d["SNP"])
    print("[munge] panel SNPs: %d | frq SNPs: %d" % (len(panel_ids), len(frq)))
    c = munge_one(phenocode, case_n, control_n, frq, panel_ids)
    cnt_path = MUNGED / "munge_counts.csv"
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(c.keys()), delimiter="\t")
    if not cnt_path.exists():
        w.writeheader()
    w.writerow(c)
    mode = "a" if cnt_path.exists() else "w"
    with open(cnt_path, mode, encoding="utf-8-sig", newline="") as fh:
        fh.write(buf.getvalue())
    print("[munge] %s 完成：%s" % (phenocode, {k: v for k, v in c.items() if k != "phenocode"}))


if __name__ == "__main__":
    main()
