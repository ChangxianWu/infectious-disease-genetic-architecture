# -*- coding: utf-8 -*-
""" 19_gene_sets.py — 基因集跨病原一致性检验（Stage 4，2026-10-11）

口径（设计文档 B2 第 5、6 步）：
  - 先验基因集**取自可核验来源**，不手写记忆清单：
      * 干扰素轴 / 补体 / Toll 样受体 / 白细胞介素 —— HGNC 官方 gene_group 字段
        
  - 检验对象 = 位点图谱 lead SNP 的最近基因并集（FinnGen 官方 nearest_genes）
  - 背景 = HGNC 全部 protein-coding 基因
  - 统计 = Fisher 精确检验（单侧富集）+ BH（多重检验）
  - IEI 对照（设计文档第 6 步）：**未做**——IUIS 2024 IEI 清单无稳定可取回的机器可读
    端点（genenames.org REST 无 IEI gene group；IUIS 表格为 PDF），按规则 96 不手抄
    记忆清单充数，列为本课题明确待办（人工取回 IUIS 表 1–8 后单列一层）。
用法：python 19_gene_sets.py
输出：work/output/gene_set_enrichment.csv
"""
import sys

import pandas as pd
from scipy.stats import false_discovery_control, fisher_exact

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATA, OUT  # noqa: E402

HGNC = DATA / "hgnc_complete_set.txt"
SET_PATTERNS = {
    "interferon": ["interferon"],
    "complement": ["complement"],
    "toll_like_receptor": ["toll like"],
    "interleukin": ["interleukin"],
    "chemokine": ["chemokine"],
    "HLA": ["histocompatibility complex"],
}


def main():
    d = pd.read_csv(HGNC, sep="\t", low_memory=False,
                    usecols=["symbol", "gene_group", "locus_group"])
    coding = d[d["locus_group"] == "protein-coding gene"].dropna(subset=["symbol"])
    background = set(coding["symbol"]) | {s.split("-AS1")[0] for s in coding["symbol"]}
    print("[geneset] 背景（protein-coding）: %d" % len(background))

    atlas = pd.read_csv(OUT / "locus_atlas.csv")
    lead_genes = set()
    for x in atlas["nearest_genes"].dropna():
        for gg in str(x).split(","):
            g = gg.strip()
            if g and g != "NA":
                lead_genes.add(g)
                lead_genes.add(g.split("-AS1")[0])
    lead_genes &= background
    print("[geneset] lead 基因（限 protein-coding）: %d -> %s"
          % (len(lead_genes), sorted(lead_genes)))

    rows = []
    for name, pats in SET_PATTERNS.items():
        sel = coding[coding["gene_group"].fillna("").str.lower().apply(
            lambda s: any(p in s for p in pats))]
        members = set(sel["symbol"])
        if not members:
            continue
        a = len(lead_genes & members)
        b = len(lead_genes - members)
        c = len(members - lead_genes)
        d_ = len(background - members - lead_genes)
        orr, p = fisher_exact([[a, b], [c, d_]], alternative="greater")
        rows.append({"gene_set": name, "n_members": len(members),
                     "n_lead_in_set": a, "odds_ratio": orr, "p_fisher": p,
                     "lead_in_set": ",".join(sorted(lead_genes & members))})
    out = pd.DataFrame(rows)
    if len(out):
        out["q_bh"] = false_discovery_control(out["p_fisher"].values, method="bh")
        out = out.sort_values("p_fisher")
        out.to_csv(OUT / "gene_set_enrichment.csv", index=False)
        print()
        print(out[["gene_set", "n_members", "n_lead_in_set", "odds_ratio", "p_fisher",
                   "q_bh", "lead_in_set"]].to_string(index=False))
    print()
    print("[geneset] IEI 对照层未做（IUIS 2024 无机器可读端点；严禁手抄记忆清单）")


if __name__ == "__main__":
    main()
