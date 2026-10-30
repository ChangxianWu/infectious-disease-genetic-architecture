# -*- coding: utf-8 -*-
""" 16_annotate_leads.py — lead SNP 基因注释（FinnGen 官方 nearest_genes，hg38）
（Stage 2b 辅助，2026-10-11）

口径：注释直接取 FinnGen R12 汇总统计自带的 nearest_genes 列（官方注释，非自建），
      避免自行按坐标猜基因；行匹配用 rsids 精确等值（本课题 munge 已证 FinnGen
      rsids 无分号多值串：0 条）。
用法：python 16_annotate_leads.py
输出：work/output/lead_gene_annotation.tsv
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATA, OUT, FINNGEN  # noqa: E402

CLUMP = OUT / "clump"
SCAN_SOURCES = ["AB1_INTESTINAL_INFECTIONS", "N14_PYELONEPHR", "AB1_OTHER_SEPSIS",
                "AB1_ERYSIPELAS", "J10_PNEUMONIA"]
USECOLS = ["#chrom", "pos", "rsids", "nearest_genes"]


def main():
    snps = set()
    for f in CLUMP.glob("*.leads_annot.tsv"):
        snps.update(pd.read_csv(f, sep="\t")["SNP"].astype(str))
    snps = sorted(snps)
    print("[annot] 待注释 lead SNP: %d" % len(snps))
    seen = {}
    for ph in SCAN_SOURCES:
        gz = FINNGEN / ("finngen_R12_%s.gz" % ph)
        if not gz.exists():
            print("[annot] 缺文件，跳过 %s" % ph)
            continue
        for ch in pd.read_csv(gz, sep="\t", usecols=USECOLS, chunksize=4_000_000,
                              dtype={"rsids": str}, compression="gzip"):
            hit = ch[ch["rsids"].isin(snps)]
            for _, r in hit.iterrows():
                seen.setdefault(r["rsids"], (r["#chrom"], r["pos"], r["nearest_genes"]))
            if len(seen) >= len(snps):
                break
        print("[annot] %s 后累计 %d/%d" % (ph, len(seen), len(snps)))
        if len(seen) >= len(snps):
            break
    out = pd.DataFrame(
        [{"SNP": s, "CHR_hg38": seen.get(s, (None, None, None))[0],
          "POS_hg38": seen.get(s, (None, None, None))[1],
          "nearest_genes": seen.get(s, (None, None, None))[2] or "NA"} for s in snps])
    out.to_csv(OUT / "lead_gene_annotation.tsv", sep="\t", index=False)
    miss = int((out["nearest_genes"] == "NA").sum())
    print("[annot] 落盘 %d 行（未命中 %d）-> %s" % (len(out), miss,
          OUT / "lead_gene_annotation.tsv"))


if __name__ == "__main__":
    main()
