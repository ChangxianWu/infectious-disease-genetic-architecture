# -*- coding: utf-8 -*-
""" 54_ontology_sets.py — 屏障/急性期本体先验检验（后续一轮外部意见 B4，2026-10-29）

目的：把"共享位点=黏膜屏障+急性期生物学"从**位点归纳**升级为**本体层检验**。
做法与 19 号一致（HGNC gene group 作集合定义、蛋白编码基因为背景、单侧 Fisher、BH）。

**诚实性声明（重要）**：
  * 本面板是在看到位点结果之后选定的（FUT2/ABO/CRP 等已知，再挑对应本体），
    因此属**假设一致性检验（hypothesis-consistent）**，不是预注册的独立检验；
  * 为避免"挑选后检验"的膨胀，面板内做 BH 校正，且与 19 号的免疫面板并列报告；
  * 阴性结果（如干扰素/白细胞介素为零）同样报告——面板若整体无富集，也不能反推否定。

面板（8 个，尽量覆盖"屏障 + 急性期"两类，并保留阴性对照）：
  屏障类：Mucins / Defensins, beta / Blood group antigens / Sialyltransferases /
          Beta 3-glycosyltransferases / Beta 4-glycosyltransferases
  急性期类：Complement system activation components / Interleukin 1 family
  阴性对照（预期不富集）：Interferons / Interleukins
用法：python 54_ontology_sets.py
输出：work/output/ontology_sets.csv + ontology_sets_summary.csv
"""
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import false_discovery_control, fisher_exact

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATA, OUT  # noqa: E402

HGNC = DATA / "hgnc_complete_set.txt"
LEAD = OUT / "lead_gene_annotation.tsv"
OUTCSV = OUT / "ontology_sets.csv"
SUMCSV = OUT / "ontology_sets_summary.csv"

# 面板（起始 10 集；按"扩充先验集"的意见扩到 15 集）。
# 每个集合第三项为 note：空=常规；"construction"=集合极小且含 lead 基因本身，
# 富集属构造性（例如 pentraxin 家族仅 3 个基因而 CRP 是其一），只作说明不作独立证据。
PANEL = [
    ("Blood group antigens", "barrier", ""),
    ("Mucins", "barrier", ""),
    ("Defensins, beta", "barrier", ""),
    ("Defensins, alpha", "barrier", ""),
    ("Sialyltransferases", "barrier", ""),
    ("Fucosyltransferases", "barrier", ""),
    ("Beta 3-glycosyltransferases", "barrier", ""),
    ("Beta 4-glycosyltransferases", "barrier", ""),
    ("Glycosyltransferase family 6", "barrier", ""),
    ("Glycosyltransferase family 8", "barrier", ""),
    ("Pattern recognition pentraxins", "acute-phase", "construction"),
    ("Long pentraxins", "acute-phase", "construction"),
    ("Complement system activation components", "acute-phase", ""),
    ("Interleukin 1 family", "acute-phase", ""),
    ("Interferons", "negative-control", ""),
    ("Interleukins", "negative-control", ""),
]


def main():
    hg = pd.read_csv(HGNC, sep="\t", low_memory=False, dtype=str)
    hg = hg[hg["locus_group"] == "protein-coding gene"]
    bg = set(hg["symbol"].dropna())
    lead = pd.read_csv(LEAD, sep="\t")
    lead_syms = {s.strip() for row in lead["nearest_genes"].dropna()
                 for s in str(row).split(",") if s.strip()} & bg
    n_lead = len(lead_syms)

    rows = []
    for name, klass, note in PANEL:
        members = set()
        for g in hg["gene_group"].dropna():
            if name in [x.strip().strip('"') for x in str(g).split("|")]:
                pass
        # gene_group 是 |-分隔的多值字段：逐行拆分匹配
        for sym, gg in zip(hg["symbol"], hg["gene_group"].fillna("")):
            if name in {x.strip().strip('"') for x in str(gg).split("|")}:
                members.add(sym)
        members &= bg
        if not members:
            rows.append({"set": name, "class": klass, "note": note, "n_members": 0, "n_lead_in_set": 0,
                         "odds_ratio": float("nan"), "p_fisher": float("nan"),
                         "lead_in_set": ""})
            continue
        a = len(lead_syms & members)
        b = len(members) - a
        c = n_lead - a
        d = len(bg) - len(members) - c
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({"set": name, "class": klass, "note": note, "n_members": len(members), "n_lead_in_set": a,
                     "odds_ratio": float(odds), "p_fisher": float(p),
                     "lead_in_set": ",".join(sorted(lead_syms & members))})
    df = pd.DataFrame(rows)
    ok = df["p_fisher"].notna()
    df["q_bh"] = float("nan")
    if ok.any():
        df.loc[ok, "q_bh"] = false_discovery_control(df.loc[ok, "p_fisher"].values, method="bh")
    df.to_csv(OUTCSV, index=False)

    hit = df[(df["q_bh"] < 0.05)]
    neg = df[df["class"] == "negative-control"]
    summ = pd.DataFrame([
        {"key": "ont_n_lead_genes", "value": str(n_lead), "note": "检验用 lead 基因数（HGNC 蛋白编码内）"},
        {"key": "ont_n_sets", "value": str(len(df)), "note": "本体面板集合数（事后选定，见脚本声明）"},
        {"key": "ont_n_enriched", "value": str(len(hit)), "note": "BH q<0.05 的集合数"},
        {"key": "ont_enriched_sets", "value": "; ".join("%s(%d/%d,q=%.2g)" % (r["set"], r["n_lead_in_set"],
                                                                              r["n_members"], r["q_bh"])
                                                       for _, r in hit.iterrows()) or "无",
         "note": "富集集合及其命中数/集合大小/q"},
        {"key": "ont_neg_control_leads", "value": "; ".join("%s:%d" % (r["set"], r["n_lead_in_set"])
                                                            for _, r in neg.iterrows()),
         "note": "阴性对照（干扰素/白细胞介素）命中数——预期为 0"},
        {"key": "ont_construction_sets", "value": "; ".join("%s:%d/%d" % (r["set"], r["n_lead_in_set"],
                                                                     r["n_members"])
                                                            for _, r in df.iterrows() if r["note"] == "construction") or "无",
         "note": "构造性集合（极小且含 lead 基因本身），其富集不作独立证据"},
        {"key": "ont_status", "value": "hypothesis-consistent",
         "note": "面板为看到位点结果后选定，属假设一致性检验，非预注册独立检验"},
    ])
    summ.to_csv(SUMCSV, index=False)
    print("[ont] lead 基因 %d 个；面板 %d 个集合；BH<0.05 的 %d 个" % (n_lead, len(df), len(hit)))
    print(df[["class", "set", "n_members", "n_lead_in_set", "odds_ratio", "p_fisher", "q_bh"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
