# -*- coding: utf-8 -*-
""" 17b_coloc_N.py — 写共定位所需的有效样本量映射（官方 manifest）
coloc.abf 的 N 用 FinnGen 官方病例+对照（有效样本量口径与主矩阵 N 一致）。
用法：python 17b_coloc_N.py
输出：work/data/coloc/trait_N.tsv（trait_label → N）
"""
import sys

import pandas as pd

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATA, FINNGEN_ENDPOINTS  # noqa: E402

man = pd.read_csv(DATA / "finngen_R12_pheno_n.tsv", sep="\t")
cnt = {r["phenocode"]: int(r["num_cases"]) + int(r["num_controls"]) for _, r in man.iterrows()}
rows = [{"trait": lab, "phenocode": ph, "N": cnt.get(ph)}
        for lab, ph in FINNGEN_ENDPOINTS.items() if cnt.get(ph)]
cnt2 = pd.read_csv(DATA / "munged/munge_counts.csv", sep="\t")
for _, r in cnt2.iterrows():
    ph = r["phenocode"]
    if ph.startswith("HGI_") and pd.notna(r.get("N_median")):
        rows.append({"trait": ph, "phenocode": ph, "N": int(r["N_median"])})
out = pd.DataFrame(rows)
out.to_csv(DATA / "coloc" / "trait_N.tsv", sep="\t", index=False)
print(out.to_string(index=False))
print("[coloc-N] %d 性状 -> trait_N.tsv" % len(out))
