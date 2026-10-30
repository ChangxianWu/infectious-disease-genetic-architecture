# -*- coding: utf-8 -*-
""" 39_iei_layer.py — E1 IEI 对照层（设计文档第 6 步补交付）
（2026-10-25）

端点：PanelApp REST API Panel 398 "Primary immunodeficiency or monogenic
inflammatory bowel disease"（Immunology 组）。IUIS 2024 清单为 PDF 无机器可读端点
（19_gene_sets.py 当初跳过的原因），PanelApp 为版本化 REST API → 规则 96 合规：
取回即落盘原始 JSON + 记录 panel version/version_created/取回日期，禁手抄。

口径：
  1. 只取 confidence_level=3（green，经过同行评审的确认证据基因；PanelApp 官方分级）
  2. 重叠对象：(a) 27 lead 的最近基因（locus_atlas.csv，逗号拆分）；
     (b) 4 个跨轴区域的最近基因（table2.csv）；(c) coloc 共享（PP.H4>0.8）区域内基因
  3. 背景：HGNC complete set 全部 protein-coding 基因（与 19_gene_sets.py 背景一致）
  4. 统计：Fisher 精确检验（双边）——lead 基因集中 IEI 基因占比 vs 背景占比
  5. 全部计数落盘；重叠基因逐条列出（含其 PanelApp 实体名）

用法：python 39_iei_layer.py
输出：work/output/iei_overlap.csv（逐基因）+ iei_summary.csv（统计 + 版本溯源）
"""
import json
import sys
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATA, OUT  # noqa: E402

PANEL_ID = 398
API = "https://panelapp.genomicsengland.co.uk/api/v1/panels/%d/" % PANEL_ID
RAW_JSON = DATA / "panelapp_398_raw.json"
LEAD = OUT / "lead_gene_annotation.tsv"
T2 = OUT / "table2.csv"
COLOC = OUT / "coloc_results.csv"
HGNC = DATA / "hgnc_complete_set.txt"
OUT_LEAD = OUT / "iei_overlap.csv"
OUT_SUM = OUT / "iei_summary.csv"


def fetch_panel():
    """分页取回 Panel 398 全量 genes，原始 JSON 落盘（版本化溯源）。"""
    rows, meta = [], {}
    url = API
    while url:
        req = urllib.request.Request(url, headers={"User-Agent": "project-iei-layer/1.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            payload = json.loads(r.read().decode("utf-8"))
        if not meta:
            for k in ("id", "name", "version", "version_created", "disease_group"):
                meta[k] = payload.get(k)
        rows.extend(payload.get("results", [payload] if "genes" not in payload else payload["genes"]))
        url = payload.get("next")
    return meta, rows


def load_green_genes(rows):
    """confidence_level==3 → gene symbol 集合（逐条计数）。"""
    green, total = set(), 0
    for g in rows:
        total += 1
        gd = g.get("gene_data", {})
        sym = gd.get("gene_symbol")
        if sym and g.get("confidence_level") == "3":
            green.add(sym)
    return green, total


def main():
    # 网络失败时回退到早先实现原始 JSON 快照（版本化溯源不变：panel 版本与取回日期
    # 始终记录在快照与输出里；构建链不因一次性网络故障而整链失败）
    try:
        meta, rows = fetch_panel()
        RAW_JSON.write_text(json.dumps({"meta": meta, "fetched": date.today().isoformat(),
                                        "n_entries": len(rows), "rows": rows},
                                       ensure_ascii=False))
    except Exception as exc:
        if not RAW_JSON.exists():
            raise
        snap = json.loads(RAW_JSON.read_text(encoding="utf-8"))
        meta, rows = snap["meta"], snap["rows"]
        print("[iei][WARN] PanelApp 取回失败（%s），回退本地快照 %s（fetched=%s）"
              % (str(exc)[:70], RAW_JSON.name, snap.get("fetched")))
    green, n_entries = load_green_genes(rows)
    print("[iei] Panel %s v%s (%s)：条目 %d，green 基因 %d" %
          (meta.get("id"), meta.get("version"), meta.get("version_created"), n_entries, len(green)))

    hgnc = pd.read_csv(HGNC, sep="\t", low_memory=False, dtype=str)
    hgnc = hgnc[hgnc["locus_group"] == "protein-coding gene"]
    bg = set(hgnc["hgnc_id"])
    bg_syms = set(hgnc["symbol"].dropna())
    # 背景 = protein-coding 符号（PanelApp symbol 与 HGNC symbol 同源，对不上的如实计数丢弃）
    green_in_bg = green & bg_syms
    print("[iei] green 基因在 HGNC protein-coding 背景内：%d / %d（对不上 %d 个，丢弃计数）" %
          (len(green_in_bg), len(green), len(green - bg_syms)))

    lead = pd.read_csv(LEAD, sep="\t")
    lead_syms = sorted({s.strip() for row in lead["nearest_genes"].dropna()
                        for s in str(row).split(",") if s.strip()})
    t2 = pd.read_csv(T2)
    t2_syms = sorted({s.strip() for row in t2["Nearest_genes"].dropna()
                      for s in str(row).split(",") if s.strip()})
    coloc = pd.read_csv(COLOC)
    # 2026-10-28 修正（第二轮评审同型缺陷，规则 224）：列名是 PP_H4（下划线），且该表
    # **没有** nearest_genes 列——须经 locus（lead rsID）到 locus_atlas 的最近基因映射。
    # 旧写法用 "PP.H4"/"nearest_genes" 的 in 判断，静默产出 0（空集假读）。
    assert "PP_H4" in coloc.columns, "coloc_results 缺 PP_H4 列"
    lead_map = {str(r["SNP"]): str(r["nearest_genes"]) for _, r in lead.iterrows()}
    shared = coloc[coloc["PP_H4"] > 0.8]
    syms = set()
    for loc in shared["locus"].astype(str):
        for s in lead_map.get(loc, "").split(","):
            if s.strip():
                syms.add(s.strip())
    coloc_syms = sorted(syms & bg_syms)

    recs = []
    for layer, syms in (("lead", lead_syms), ("cross_axis_regions", t2_syms),
                        ("coloc_shared", coloc_syms)):
        for s in syms:
            hit = s in green_in_bg
            recs.append({"layer": layer, "symbol": s, "iei_green": hit,
                         "panel_version": meta.get("version"),
                         "panel_fetched": date.today().isoformat()})
    pd.DataFrame(recs).to_csv(OUT_LEAD, index=False)

    # Fisher：lead 基因（限 HGNC 背景内）vs 背景
    lead_in_bg = [s for s in lead_syms if s in bg_syms]
    k = sum(1 for s in lead_in_bg if s in green_in_bg)
    K, N = len(green_in_bg), len(bg_syms)
    n = len(lead_in_bg)
    table = [[k, n - k], [K - k, N - K - n + k]]
    odds, p = fisher_exact(table, alternative="two-sided")
    # 跨轴区域基因同一检验（区域基因是架构主张的载体）
    t2_in_bg = [s for s in t2_syms if s in bg_syms]
    k2 = sum(1 for s in t2_in_bg if s in green_in_bg)
    odds2, p2 = fisher_exact([[k2, len(t2_in_bg) - k2], [K - k2, N - K - len(t2_in_bg) + k2]],
                             alternative="two-sided")

    hits = [(r["layer"], r["symbol"]) for r in recs if r["iei_green"]]
    summ = pd.DataFrame([
        {"key": "panel_id", "value": str(meta.get("id")), "note": "PanelApp REST API（版本化端点）"},
        {"key": "panel_version", "value": str(meta.get("version")),
         "note": "version_created=%s，取回=%s" % (meta.get("version_created"), date.today().isoformat())},
        {"key": "panel_entries", "value": str(n_entries), "note": "含全部 confidence 级"},
        {"key": "iei_green_genes", "value": str(K), "note": "confidence_level==3（green）且在 HGNC 背景内"},
        {"key": "background_genes", "value": str(N), "note": "HGNC protein-coding（与 19 号脚本背景一致）"},
        {"key": "lead_genes_n", "value": str(len(lead_in_bg)), "note": "27 lead 最近基因（限背景内）"},
        {"key": "lead_iei_overlap", "value": str(k),
         "note": "重叠：%s" % (", ".join(s for l, s in hits if l == "lead") or "无")},
        {"key": "lead_iei_fisher_p", "value": "%.3g" % p, "note": "Fisher 双边；OR=%.2f" % odds},
        {"key": "crossaxis_genes_n", "value": str(len(t2_in_bg)), "note": "4 个跨轴区域最近基因"},
        {"key": "crossaxis_iei_overlap", "value": str(k2),
         "note": "重叠：%s" % (", ".join(s for l, s in hits if l == "cross_axis_regions") or "无")},
        {"key": "crossaxis_iei_fisher_p", "value": "%.3g" % p2, "note": "Fisher 双边；OR=%.2f" % odds2},
        {"key": "coloc_shared_genes_n", "value": str(len(coloc_syms)),
         "note": "coloc 共享 PP.H4>0.8 对涉及基因（背景内）"},
        {"key": "coloc_shared_iei_overlap", "value": str(sum(1 for l, s in hits if l == "coloc_shared")),
         "note": "重叠：%s" % (", ".join(s for l, s in hits if l == "coloc_shared") or "无")},
    ])
    summ.to_csv(OUT_SUM, index=False)
    print("[iei] lead 重叠 %d/%d（p=%.3g）；跨轴区域重叠 %d/%d（p=%.3g）" %
          (k, len(lead_in_bg), p, k2, len(t2_in_bg), p2))
    print("[iei] 落盘 %s / %s" % (OUT_LEAD.name, OUT_SUM.name))


if __name__ == "__main__":
    main()
