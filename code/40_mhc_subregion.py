# -*- coding: utf-8 -*-
""" 40_mhc_subregion.py — E2 MHC 亚区位置分解（2026-10-25）

主张升级：把"共享架构走 MHC-II 抗原呈递臂而非 I 型干扰素"从基因集富集证据
升级为位置层定量——MHC 邻近 lead / coloc 共享变异的最近基因按 HLA 类 I / 类 II /
类 III 区归类（基因家族归属法，非 Mb 边界切分，避免坐标近似争议）。

口径：
  1. 对象：locus_atlas.csv 中 in_mhc=True 的 lead + coloc_results.csv 中
     PP.H4>0.8 且 chr6 的共享变异
  2. 归类（最近基因名 → 家族）：
     - HLA class II：HLA-DR/HLA-DQ/HLA-DO/HLA-DM/HLA-DP 前缀
     - HLA class I：HLA-A/B/C/E/F/G/H 精确匹配
     - class III / extended：TNXB、C2、C4A/B、BF、CYP21A2、NFKBIL1、ATP6V1G2、
       AGER、RDBP、STK19、SKIV2L、CSNK2B、GABBR1、PBX2（仅列本项目实际出现的；
       全部为 6 号染色体 MHC class III/extended 区基因，GRCh38 ~31.3–32.3 Mb 邻域）
     - 其他 MHC：其余（如 MCC/YWHAZ 等非 HLA 基因）
  3. 每个 lead 的 axes 数、discovered trait 一并带出（表格可直接进补充材料）
  4. 计数落盘：n_mhc_leads、n_class1/2/3/other、coloc 同口径

用法：python 40_mhc_subregion.py
输出：work/output/mhc_subregion_leads.csv + mhc_subregion_summary.csv
"""
import sys
import re

import pandas as pd

sys.path.insert(0, "/path/to/project/work/code")
from constants import OUT  # noqa: E402

CLASS2_RE = re.compile(r"^HLA-(DR|DQ|DO|DM|DP)")
CLASS1_SET = {"HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-H"}
CLASS3_SET = {"TNXB", "C2", "C4A", "C4B", "BF", "CYP21A2", "NFKBIL1", "ATP6V1G2",
              "AGER", "RDBP", "STK19", "SKIV2L", "CSNK2B", "GABBR1", "PBX2"}


def classify_gene(sym):
    if CLASS2_RE.match(sym):
        return "HLA_class_II"
    if sym in CLASS1_SET:
        return "HLA_class_I"
    if sym in CLASS3_SET:
        return "class_III_extended"
    return "MHC_other"


def classify_first(nearest):
    """最近基因列（逗号分隔）：取第一个可归类基因的类别；全不可归类 → MHC_other。"""
    for s in str(nearest).split(","):
        c = classify_gene(s.strip())
        if c != "MHC_other":
            return c, s.strip()
    return "MHC_other", str(nearest).split(",")[0].strip()


def _nearest_of(snp, atlas):
    """coloc 行的 top_snp 是 'rs1;rs2;rs3' 复合串——经 atlas 反查该 lead 的最近基因。"""
    first = str(snp).split(";")[0].strip()
    hit = atlas[atlas["SNP"] == first]
    return str(hit.iloc[0]["nearest_genes"]) if len(hit) else first


def main():
    atlas = pd.read_csv(OUT / "locus_atlas.csv")
    mhc = atlas[atlas["in_mhc"] == True]  # noqa: E712
    recs = []
    for _, r in mhc.iterrows():
        cls, gene = classify_first(r["nearest_genes"])
        recs.append({"source": "lead", "snp": r["SNP"], "pos": r["BP"],
                     "nearest_genes": r["nearest_genes"], "classified_gene": gene,
                     "subregion": cls, "n_axes_sig": r["n_axes_sig"], "axes": r["axes"],
                     "locus_class": r["locus_class"]})
    coloc = pd.read_csv(OUT / "coloc_results.csv")
    # 列名是 PP_H4（下划线）；locus 列=lead rsID → 经 locus_atlas 的 in_mhc lead 集合判定
    # （首版误用 'PP.H4' 使分支静默跳过、0 为空集假读数——规则 221 家族陷阱，2026-10-25 修）
    if "PP_H4" in coloc.columns:
        col = coloc[coloc["PP_H4"] > 0.8]
        mhc_lead_ids = set(mhc["SNP"])
        col_mhc = col[col["locus"].isin(mhc_lead_ids)]
        for _, r in col_mhc.iterrows():
            # 归类按 **lead** 的最近基因（atlas 反查）；top_snp 是 credible 集变异 rsID，
            # 只随行记录不作分类依据（rsID→基因归属属记忆推断，禁用）
            lead_atlas = atlas[atlas["SNP"] == r["locus"]]
            ng_lead = str(lead_atlas.iloc[0]["nearest_genes"]) if len(lead_atlas) else ""
            cls, gene = classify_first(ng_lead)
            recs.append({"source": "coloc_shared", "snp": r["locus"],
                         "pos": "", "nearest_genes": ng_lead,
                         "classified_gene": gene, "subregion": cls,
                         "n_axes_sig": "", "axes": "%s x %s" % (r["trait1"], r["trait2"]),
                         "locus_class": "coloc_PP_H4>0.8 top_snp=%s" % r["top_snp"]})
    out = pd.DataFrame(recs)
    out.to_csv(OUT / "mhc_subregion_leads.csv", index=False)

    lead_out = out[out["source"] == "lead"]
    n = len(lead_out)
    cnt = lead_out["subregion"].value_counts().to_dict()
    summ = pd.DataFrame([
        {"key": "n_mhc_leads", "value": str(n), "note": "in_mhc=True 的 lead 变异数"},
        {"key": "n_mhc_class_II", "value": str(cnt.get("HLA_class_II", 0)),
         "note": "最近基因为 HLA-DR/DQ/DO/DM/DP（类 II 抗原呈递臂）"},
        {"key": "n_mhc_class_I", "value": str(cnt.get("HLA_class_I", 0)), "note": "HLA-A/B/C/E/F/G"},
        {"key": "n_mhc_class_III", "value": str(cnt.get("class_III_extended", 0)),
         "note": "class III/extended 区基因（TNXB/NFKBIL1/ATP6V1G2 等）"},
        {"key": "n_mhc_other", "value": str(cnt.get("MHC_other", 0)), "note": "其余"},
        {"key": "share_class_II", "value": ("%.3f" % (cnt.get("HLA_class_II", 0) / n)) if n else "",
         "note": "类 II 占 MHC lead 比例（位置层证据）"},
        {"key": "n_mhc_coloc_shared", "value": str(len(out[out["source"] == "coloc_shared"])),
         "note": "coloc 共享（PP_H4>0.8）中落在 MHC lead 的行数（variant 级；经 locus∈in_mhc lead 集判定）"},
    ])
    summ.to_csv(OUT / "mhc_subregion_summary.csv", index=False)
    print("[mhc-sub] MHC lead %d：类II %d / 类I %d / 类III %d / 其他 %d" %
          (n, cnt.get("HLA_class_II", 0), cnt.get("HLA_class_I", 0),
           cnt.get("class_III_extended", 0), cnt.get("MHC_other", 0)))
    print("[mhc-sub] 落盘 mhc_subregion_leads.csv / mhc_subregion_summary.csv")


if __name__ == "__main__":
    main()
