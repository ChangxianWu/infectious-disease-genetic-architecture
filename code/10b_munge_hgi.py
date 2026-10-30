# -*- coding: utf-8 -*-
""" 10b_munge_hgi.py — COVID-19 HGI freeze7（EUR 分层）→ LDSC .sumstats
（Stage 2 外部层，2026-10-11）

列名**实查落盘**（禁猜，规则 49）：work/data/covidhgi 两个文件的 header 为
  #CHR POS REF ALT SNP all_meta_N all_inv_var_meta_beta all_inv_var_meta_sebeta
  all_inv_var_meta_p all_inv_var_meta_cases all_inv_var_meta_controls
  all_inv_var_meta_effective all_inv_var_het_p all_meta_AF rsid
口径：
  - SNP = rsid 列（缺失/'.' 丢弃）；A1 = ALT（beta 的效应等位），A2 = REF；
  - Z = beta/sebeta；N = cases + controls（逐 SNP，与 FinnGen 端 N=ncase+ncon 同规）；
  - 与 FinnGen 端完全相同的 1000G EUR 对齐 / MAF / p-Z 一致性过滤；
  - all_meta_N（首个文件首行=25）是**队列数**而非样本量，不得当 N 用（此处显式记录该坑）。
用法：python 10b_munge_hgi.py <C2|A2>
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATA, GWAS_DATA  # noqa: E402

HGI = DATA / "covidhgi"
MUNGED = DATA / "munged"
LDSC_REF = GWAS_DATA / "ldsc_ref"
FILES = {
    "C2": HGI / "COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz",
    "A2": HGI / "COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz",
}
USECOLS = ["REF", "ALT", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta",
           "all_inv_var_meta_p", "all_inv_var_meta_cases", "all_inv_var_meta_controls",
           "rsid"]


def load_1000g_frq():
    frames = []
    for i in range(1, 23):
        frames.append(pd.read_csv(LDSC_REF / ("1000G_Phase3_frq/1000G.EUR.QC.%d.frq" % i),
                                  sep=r"\s+", usecols=["SNP", "A1", "A2", "MAF"]))
    return pd.concat(frames, ignore_index=True).drop_duplicates(subset="SNP", keep="first")


def main():
    key = sys.argv[1]
    t0 = time.time()
    frq = load_1000g_frq()
    panel = set()
    for i in range(1, 23):
        d = pd.read_csv(LDSC_REF / "LDscore" / ("LDscore.%d.l2.ldscore.gz" % i),
                        sep="\t", usecols=["SNP"])
        panel.update(d["SNP"])
    counts = {"phenocode": "HGI_%s" % key}
    chunks, total = [], 0
    for ch in pd.read_csv(FILES[key], sep="\t", usecols=USECOLS, chunksize=4_000_000,
                          dtype={"REF": str, "ALT": str, "rsid": str},
                          compression="gzip"):
        total += len(ch)
        ch = ch[ch["rsid"].notna() & (ch["rsid"] != ".") & ch["rsid"].isin(panel)]
        chunks.append(ch)
    counts["total_variants"] = total
    counts["in_l2_panel"] = sum(len(c) for c in chunks)
    dat = pd.concat(chunks, ignore_index=True)
    del chunks
    n0 = len(dat)
    dat = dat.drop_duplicates(subset="rsid", keep="first")
    counts["dup_rsids_dropped"] = n0 - len(dat)
    m = dat.merge(frq, left_on="rsid", right_on="SNP", how="inner")
    counts["merged_1000g"] = len(m)
    ref, alt = m["REF"].str.upper(), m["ALT"].str.upper()
    a1, a2 = m["A1"].str.upper(), m["A2"].str.upper()
    same = (alt == a1) & (ref == a2)
    flip = (alt == a2) & (ref == a1)
    keep = same | flip
    m = m[keep].copy()
    counts["allele_mismatch_dropped"] = counts["merged_1000g"] - len(m)
    sign = np.where(same[m.index], 1.0, -1.0)
    m["z"] = sign * (m["all_inv_var_meta_beta"] / m["all_inv_var_meta_sebeta"])
    m["A1"] = a1[m.index]
    m["A2"] = a2[m.index]
    ok = np.isfinite(m["z"]) & (m["all_inv_var_meta_sebeta"] > 0)
    counts["se_or_nan_dropped"] = int((~ok).sum())
    m = m[ok].copy()
    maf_ok = (m["MAF"] >= 0.01) & (m["MAF"] <= 0.99)
    counts["maf_out_dropped"] = int((~maf_ok).sum())
    m = m[maf_ok].copy()
    m["chisq_z"] = m["z"] ** 2
    m["chisq_p"] = chi2_dist.isf(m["all_inv_var_meta_p"].clip(1e-300, 1.0), 1)
    rel = (m["chisq_z"] - m["chisq_p"]).abs() / np.maximum(
        m[["chisq_z", "chisq_p"]].max(axis=1), 1.0)
    pz_ok = rel <= 0.30
    counts["pz_mismatch_dropped"] = int((~pz_ok).sum())
    m = m[pz_ok].copy()
    m["N"] = m["all_inv_var_meta_cases"] + m["all_inv_var_meta_controls"]
    counts["N_median"] = int(np.median(m["N"]))
    counts["N_min"] = int(m["N"].min())
    counts["N_max"] = int(m["N"].max())
    out = pd.DataFrame({"SNP": m["rsid"], "A1": m["A1"], "A2": m["A2"],
                        "N": m["N"], "Z": m["z"], "CHISQ": m["chisq_z"]})
    out.to_csv(MUNGED / ("HGI_%s.sumstats" % key), sep="\t", index=False)
    counts["final_snps"] = len(out)
    counts["mean_chisq"] = float(np.mean(out["CHISQ"]))
    counts["runtime_s"] = round(time.time() - t0, 1)
    print("[hgi-munge] %s: %s" % (key, {k: v for k, v in counts.items() if k != "phenocode"}))
    cnt = MUNGED / "munge_counts.csv"
    hdr = open(cnt, encoding="utf-8-sig").readline().rstrip("\n").split("\t")
    row = "\t".join(str(counts.get(h, "")) for h in hdr)
    with open(cnt, "a", encoding="utf-8-sig", newline="") as fh:
        fh.write(row + "\n")


if __name__ == "__main__":
    main()
