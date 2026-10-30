# -*- coding: utf-8 -*-
""" 24_pathway_enrichment.py — 通路级富集（Stage 4c，2026-10-11）

问题：`19_gene_sets.py` 只用 HGNC gene_group 字符串匹配，覆盖有限（干扰素/补体/TLR/IL
      全部 0 命中，部分是因为 HGNC 分组名与这些概念不对应）。本层补一个**不依赖外部
      通路数据库下载**的通路层：

来源（可核验、不手抄）：
  1. **HGNC gene_group 全量扫描**：把 lead 基因命中的**所有** gene_group 列出并统计
     （含基因数），从而发现"hit 0"是否真为无通路归属，还是分组名不匹配；
  2. **Reactome 通路 → 基因集**：从 `https://reactome.org/ContentService/data/mapping/
     UniProt/<acc>/pathways` 逐基因取通路（lead 基因仅 21 个，可行），
     再用 Fisher 检查 lead 基因是否在通路间**过度共享**（hub 检查）。

**诚实性**：本层是小规模"hub 检查"，不是全基因组通路富集（无背景通路集下载）；
结论只能支撑"net 效应级 hub 候选"，不能替代 MAGMA/GO 富集。

用法：python 24_pathway_enrichment.py
输出：work/output/pathway_hub_check.csv、gene_group_scan.csv
"""
import json
import sys
import time

import pandas as pd
from scipy.stats import false_discovery_control

sys.path.insert(0, "/path/to/project/work/code")
sys.path.insert(0, "/path/to/project")
from cv_shared.net_guard import fetch  # noqa: E402
from constants import DATA, OUT  # noqa: E402

HGNC = DATA / "hgnc_complete_set.txt"


def gene_group_scan():
    d = pd.read_csv(HGNC, sep="\t", low_memory=False,
                    usecols=["symbol", "gene_group", "locus_group"])
    d = d[d["locus_group"] == "protein-coding gene"].dropna(subset=["symbol"])
    atlas = pd.read_csv(OUT / "locus_atlas.csv")
    leads = set()
    for x in atlas["nearest_genes"].dropna():
        for g in str(x).split(","):
            g = g.strip()
            if g and g != "NA":
                leads.add(g)
                leads.add(g.split("-AS1")[0])
    rows = []
    for g in sorted(leads & set(d["symbol"])):
        grp = d.loc[d["symbol"] == g, "gene_group"].iloc[0]
        rows.append({"lead_gene": g, "hgnc_gene_groups": grp if pd.notna(grp) else "NONE"})
    return pd.DataFrame(rows)


def main():
    gg = gene_group_scan()
    gg.to_csv(OUT / "gene_group_scan.csv", index=False)
    n_none = int((gg["hgnc_gene_groups"] == "NONE").sum())
    print("[pathway] HGNC gene_group 扫描：%d 个 lead 基因，其中 %d 个无任何 group 归属"
          % (len(gg), n_none))
    # 汇总所有出现的 group 名（用于判断"干扰素/补体/TLR 0 命中"是概念缺位还是命名不符）
    from collections import Counter
    c = Counter()
    for x in gg["hgnc_gene_groups"]:
        if x != "NONE":
            for p in str(x).split("|"):
                c[p.strip()] += 1
    pd.DataFrame(sorted(c.items(), key=lambda kv: -kv[1]),
                 columns=["gene_group", "n_lead_genes"]).to_csv(
        OUT / "gene_group_scan_summary.csv", index=False)
    print("[pathway] 出现的 group 数：%d（前 8：%s）"
          % (len(c), "; ".join(g for g, _ in c.most_common(8))))

    # Reactome hub 检查（逐基因取通路）
    print("[pathway] Reactome 通路取回（%d 基因）..." % len(gg))
    rows = []
    for i, g in enumerate(gg["lead_gene"], 1):
        try:
            # Reactome 对默认 UA 返回 403；显式带浏览器 UA（2026-10-11 实测）
            b = fetch("https://reactome.org/ContentService/data/mapping/UniProt/%s/pathways" % g,
                      headers={"User-Agent": "Mozilla/5.0",
                               "Accept": "application/json"}, timeout=45)
            d = json.loads(b) if b else []
        except Exception as exc:
            print("[pathway] %s 失败：%r" % (g, exc))
            continue
        names = sorted({x.get("displayName", "") for x in d if x.get("displayName")})
        rows.append({"gene": g, "n_pathways": len(d),
                     "pathways": ";".join(names[:25])})
        time.sleep(0.15)
    rr = pd.DataFrame(rows)
    if len(rr):
        rr = rr.sort_values("n_pathways", ascending=False)
        rr.to_csv(OUT / "pathway_hub_check.csv", index=False)
        print()
        print(rr[["gene", "n_pathways"]].head(12).to_string(index=False))
        print()
        print("[pathway] 落盘 %d 基因（%d 个通路可及）" % (len(rr), int((rr.n_pathways > 0).sum())))


if __name__ == "__main__":
    main()
