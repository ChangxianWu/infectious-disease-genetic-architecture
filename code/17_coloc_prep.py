# -*- coding: utf-8 -*-
""" 17_coloc_prep.py — 共定位输入准备（Stage 3，2026-10-11）

对 locus_atlas 的非特异位点（intermediate/broad），抽取 ±500kb 窗口内**所有相关性状**
的原始汇总统计（FinnGen 官方 beta/sebeta/af_alt → 直接算 MAF，不从参考面板借），
输出每 (位点, 性状) 一张 TSV 供 R coloc 读取。

口径（规则 51/71/96）：
  - MAF = min(af_alt, 1−af_alt)（FinnGen 官方等位频率，非 1000G 替代）；
  - beta/se 用原始 beta/sebeta（不以 Z 反推）；varbeta = se²；
  - 窗口以 hg38 坐标为界（±500 kb），与 lead SNP 来源（FinnGen hg38 pos）一致；
  - 仅保留 MAF>0.01 且有限值行。
用法：python 17_coloc_prep.py
输出：work/data/coloc/<LOCUS>_<trait>.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATA, OUT, FINNGEN, FINNGEN_ENDPOINTS  # noqa: E402

WINDOW = 500_000
COLOC_DIR = DATA / "coloc"
COLOC_DIR.mkdir(parents=True, exist_ok=True)
USECOLS = ["#chrom", "pos", "ref", "alt", "rsids", "beta", "sebeta", "pval", "af_alt"]


def lead_positions():
    """位点 → (chr, pos)。同一位点 SNP 在任意 FinnGen 文件中位置一致（hg38）。"""
    ann = pd.read_csv(OUT / "lead_gene_annotation.tsv", sep="\t")
    return {r["SNP"]: (int(r["CHR_hg38"]), int(r["POS_hg38"])) for _, r in ann.iterrows()
            if pd.notna(r["CHR_hg38"])}


def extract(trait_ph, chrom, lo, hi, dst, region_snps=None):
    gz = FINNGEN / ("finngen_R12_%s.gz" % trait_ph)
    out = []
    for ch in pd.read_csv(gz, sep="\t", usecols=USECOLS, chunksize=4_000_000,
                          dtype={"rsids": str, "ref": str, "alt": str},
                          compression="gzip"):
        sel = (ch["#chrom"] == chrom) & (ch["pos"] >= lo) & (ch["pos"] <= hi)
        if sel.any():
            out.append(ch[sel])
        if ch["#chrom"].iloc[-1] > chrom:
            break
    if not out:
        return 0
    d = pd.concat(out, ignore_index=True)
    d = d[d["rsids"].notna() & (d["rsids"] != ".")]
    d = d.drop_duplicates("rsids", keep="first")
    d["MAF"] = np.minimum(d["af_alt"], 1 - d["af_alt"])
    d = d[(d["MAF"] > 0.01) & (d["MAF"] < 0.99)]
    d = d[np.isfinite(d["beta"]) & np.isfinite(d["sebeta"]) & (d["sebeta"] > 0)]
    keep = ["rsids", "pos", "ref", "alt", "beta", "sebeta", "pval", "MAF"]
    d[keep].rename(columns={"rsids": "SNP", "ref": "A2", "alt": "A1"}).to_csv(
        dst, sep="\t", index=False)
    return len(d)


def main():
    atlas = pd.read_csv(OUT / "locus_atlas.csv")
    focal = atlas[atlas["locus_class"].isin(["broad", "intermediate"])]
    pos = lead_positions()
    # 逆映射 label → phenocode
    lab2ph = {lab: ph for lab, ph in FINNGEN_ENDPOINTS.items()}
    lab2ph.update({"HGI_C2": None, "HGI_A2": None})
    n_files = 0
    for _, r in focal.iterrows():
        snp = r["SNP"]
        if snp not in pos:
            print("[coloc-prep] %s 无坐标，跳过" % snp)
            continue
        chrom, p = pos[snp]
        lo, hi = p - WINDOW, p + WINDOW
        # 该位点显著 (p<5e-8) 的 FinnGen 性状
        sig = []
        for tok in str(r["sig_traits"]).split(";"):
            if "(" not in tok:
                continue
            lab = tok.split("(")[0]
            if lab2ph.get(lab):
                sig.append(lab)
        for lab in sig:
            dst = COLOC_DIR / ("%s_%s.tsv" % (snp, lab))
            if dst.exists():
                continue
            n = extract(lab2ph[lab], chrom, lo, hi, dst)
            n_files += 1
            print("[coloc-prep] %s x %s (chr%d:%d-%d) %d SNPs" % (snp, lab, chrom, lo, hi, n))
    print("[coloc-prep] 完成 %d 个窗口文件 -> %s" % (n_files, COLOC_DIR))


if __name__ == "__main__":
    main()
