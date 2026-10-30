# -*- coding: utf-8 -*-
""" 22_hgi_coloc_prep.py — COVID HGI 端共定位窗口抽取（跨源共定位，Stage 3c）

口径：
  - 坐标体系已实测对齐：同一 rsID 在 FinnGen(hg38) 与 HGI 默认 build 的 (CHR,POS)
    **逐位一致**（rs3116651/rs876537/rs492602/rs516246/rs495828 五点复核 True）→
    可直接沿用 FinnGen 侧 hg38 窗口边界，不需 liftover。
  - 复现与 FinnGen 侧同口径：MAF = min(af,1−af)（HGI 列 `all_meta_AF`）、
    beta/sebeta 原值、pval = `all_inv_var_meta_p`、去掉 MAF 越界与 se<=0。
  - 输出文件名与 FinnGen 侧一致（<LOCUS>_HGI_C2.tsv / <LOCUS>_HGI_A2.tsv），
    使 19_coloc_susie.R 的分组逻辑无需改动即可把 HGI 纳入配对。
用法：python 22_hgi_coloc_prep.py
输出：work/data/coloc/<LOCUS>_HGI_{C2,A2}.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATA, OUT  # noqa: E402

WINDOW = 500_000
HGI_DIR = DATA / "covidhgi"
COLOC = DATA / "coloc"
FILES = {"HGI_C2": "COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz",
         "HGI_A2": "COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz"}
USECOLS = ["#CHR", "POS", "REF", "ALT", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta",
           "all_inv_var_meta_p", "all_meta_AF", "rsid"]


def main():
    ann = pd.read_csv(OUT / "lead_gene_annotation.tsv", sep="\t")
    loci = ann[ann["SNP"].isin(["rs3116651", "rs492602", "rs495828"])]
    regions = [(r["SNP"], int(r["CHR_hg38"]), int(r["POS_hg38"])) for _, r in loci.iterrows()]
    print("[hgi-coloc] 窗口: %s" % regions)
    for lab, fn in FILES.items():
        f = HGI_DIR / fn
        if not f.exists():
            print("[hgi-coloc] 缺 %s，跳过" % fn)
            continue
        got = {s: [] for s, _, _ in regions}
        for ch in pd.read_csv(f, sep="\t", usecols=USECOLS, chunksize=4_000_000,
                              dtype={"rsid": str, "REF": str, "ALT": str},
                              compression="gzip"):
            for snp, c, p in regions:
                sel = (ch["#CHR"] == c) & (ch["POS"].between(p - WINDOW, p + WINDOW))
                if sel.any():
                    got[snp].append(ch[sel])
            if all(len(v) for v in got.values()):
                pass
        for snp, c, p in regions:
            if not got[snp]:
                print("[hgi-coloc] %s %s 无窗口行" % (lab, snp))
                continue
            d = pd.concat(got[snp], ignore_index=True)
            d = d[d["rsid"].notna() & (d["rsid"] != ".")].drop_duplicates("rsid", keep="first")
            d["MAF"] = np.minimum(d["all_meta_AF"], 1 - d["all_meta_AF"])
            d = d[(d["MAF"] > 0.01) & (d["MAF"] < 0.99)]
            d = d[np.isfinite(d["all_inv_var_meta_beta"])
                  & np.isfinite(d["all_inv_var_meta_sebeta"])
                  & (d["all_inv_var_meta_sebeta"] > 0)]
            out = pd.DataFrame({
                "SNP": d["rsid"], "pos": d["POS"], "A2": d["REF"], "A1": d["ALT"],
                "beta": d["all_inv_var_meta_beta"], "sebeta": d["all_inv_var_meta_sebeta"],
                "pval": d["all_inv_var_meta_p"], "MAF": d["MAF"]})
            dst = COLOC / ("%s_%s.tsv" % (snp, lab))
            out.to_csv(dst, sep="\t", index=False)
            print("[hgi-coloc] %s %s: %d SNP -> %s" % (lab, snp, len(out), dst.name))


if __name__ == "__main__":
    main()
