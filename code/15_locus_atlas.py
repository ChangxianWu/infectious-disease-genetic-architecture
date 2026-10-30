# -*- coding: utf-8 -*-
""" 15_locus_atlas.py — 广谱 vs 病原特异位点分类（Stage 2b 核心，2026-10-11）

设计文档 B2 第 4 步（"广谱 vs 特异位点分类"）的实现：
  1. 候选位点 = 功效门通过性状的 lead SNP 并集（14_clump_leads.py，plink r²0.01/500kb/5e-8）
  2. 每个 lead SNP 在**全部 25 个 FinnGen 终点 + HGI C2/A2** 中取 z（等位对齐，
     以各性状 A1 为基准）；显著 = p < 5e-8
  3. 位点分类（轴级）：
       broad      ≥3 个不同病原轴出现同向显著
       intermediate  2 个轴
       pathogen_specific  仅 1 个轴
       opposite  出现方向相反的显著（多效性方向冲突，单列）
  4. 输出位点×性状 z 矩阵（供稿件补充表）与分类计数。

用法：python 15_locus_atlas.py
输出：work/output/locus_atlas.csv、locus_by_trait_z.csv、locus_class_summary.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, "/path/to/project/work/code")
sys.path.insert(0, "/path/to/project/work/code/ldsc_py3")
import importlib.util  # noqa: E402
from constants import DATA, OUT, FINNGEN_ENDPOINTS  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "mr", "/path/to/project/work/code/12_matrix_runner.py")
mr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mr)
A2_SPEC = importlib.util.spec_from_file_location(
    "s2", "/path/to/project/work/code/13_stage2_analysis.py")
s2 = importlib.util.module_from_spec(A2_SPEC)
A2_SPEC.loader.exec_module(s2)
ATLAS_AXIS = dict(s2.AXIS)
ATLAS_AXIS["HGI_C2"] = "viral"     # SARS-CoV-2 infection
ATLAS_AXIS["HGI_A2"] = "viral"     # severe COVID-19（外部层，轴归属同病毒轴）
ATLAS_AXIS["other"] = "other"
LEADS_DIR = OUT / "clump"
P_SIG = 5e-8


def load_master_all():
    label_to_ph = {lab: ph for lab, ph in FINNGEN_ENDPOINTS.items()
                   if ph not in mr.MTB_BLOCKED and (mr.MUNGED / ("%s.sumstats" % ph)).exists()}
    label_to_ph["HGI_C2"] = "HGI_C2"
    label_to_ph["HGI_A2"] = "HGI_A2"
    meta, traits, M = mr.load_master(label_to_ph)
    return meta, traits, M


def load_leads_by_snp(all_leads, labels):
    """候选位点 × 性状 z 矩阵：直读各性状 .sumstats 取候选位点。

    不使用主矩阵的共用索引——该索引经 noMHC weights 口径滤波后不含 MHC
    （chr6:25–35 Mb），而 HLA 恰是感染宿主遗传最经典的共享位点，位点级图谱
    必须包含，否则系统性漏掉本题最强的共享位点证据。
    效应等位统一到 base 性状的 A1（a1 相同取 z；a2 相同取 −z；否则 NA）。
    """
    per_trait = {}
    print("[atlas] 直读 %d 个性状 sumstats 取候选位点（含 MHC 区）..." % len(labels))
    for lab in labels:
        ph = FINNGEN_ENDPOINTS.get(lab, lab)
        f = mr.MUNGED / ("%s.sumstats" % ph)
        d = pd.read_csv(f, sep="\t", usecols=["SNP", "A1", "A2", "Z"])
        per_trait[lab] = d[d["SNP"].isin(all_leads)].drop_duplicates("SNP").set_index("SNP")
    base_lab = "INTESTINAL" if "INTESTINAL" in labels else labels[0]
    zcols = {}
    for s in all_leads:
        base_allele = None
        for lab in [base_lab] + labels:
            t = per_trait.get(lab)
            if t is not None and s in t.index:
                v = t.loc[s, "A1"]
                if isinstance(v, str) and v in ("A", "C", "G", "T"):
                    base_allele = v
                    break
        if base_allele is None:
            continue
        row = {}
        for lab in labels:
            t = per_trait[lab]
            if s not in t.index:
                row[lab] = np.nan
                continue
            r = t.loc[s]
            a1v, a2v = r["A1"], r["A2"]
            zv = float(r["Z"]) if np.isfinite(r["Z"]) else np.nan
            if not np.isfinite(zv) or not isinstance(a1v, str):
                row[lab] = np.nan
            elif a1v == base_allele:
                row[lab] = zv
            elif a2v == base_allele:
                row[lab] = -zv
            else:
                row[lab] = np.nan
        if any(np.isfinite(v) for v in row.values()):
            zcols[s] = row
    print("[atlas] 可用位点 %d / %d（基准等位取 %s 的 A1）"
          % (len(zcols), len(all_leads), base_lab))
    return per_trait, zcols


def main():
    leads = {}
    for f in sorted(LEADS_DIR.glob("*.leads_annot.tsv")):
        lab = f.name.replace(".leads_annot.tsv", "")
        d = pd.read_csv(f, sep="\t")
        leads[lab] = d
    if not leads:
        print("[atlas] 无 lead 文件，先跑 14_clump_leads.py")
        return
    all_leads = sorted(set(s for d in leads.values() for s in d["SNP"]))
    print("[atlas] 候选位点（lead SNP 并集）: %d 个，来自 %d 个性状" % (len(all_leads), len(leads)))

    meta, traits, M = load_master_all()
    labels = list(traits.keys())
    per_trait, zcols = load_leads_by_snp(all_leads, labels)
    zdf = pd.DataFrame(zcols, index=labels).T if zcols else pd.DataFrame()
    pdf = pd.DataFrame(2 * norm.sf(np.abs(zdf.values)), index=zdf.index, columns=zdf.columns)

    # MHC 标注：位置取自 plink clump 注释文件（1000G EUR QC 参考，hg19）
    # hg19 MHC 常用定义 chr6:25–35 Mb（与 GRCh38 的 28.5–33.4 Mb 区间一致覆盖 HLA）
    posmap = {}
    for lab, d in leads.items():
        for _, r in d.iterrows():
            if pd.notna(r.get("CHR")) and pd.notna(r.get("BP")):
                posmap.setdefault(r["SNP"], (int(r["CHR"]), int(r["BP"])))
    mhc_flag = {s: (posmap.get(s, (0, 0))[0] == 6 and 25_000_000 <= posmap.get(s, (0, 0))[1] <= 35_000_000)
                for s in zdf.index}
    print("[atlas] MHC 区位点: %d / %d" % (sum(mhc_flag.values()), len(mhc_flag)))

    # 位点分类
    rows = []
    for s in zdf.index:
        row = {"SNP": s}
        sig_axes, sig_labels, opp = set(), [], False
        for lab in labels:
            z = zdf.loc[s, lab]
            p = pdf.loc[s, lab]
            if not np.isfinite(z) or not np.isfinite(p):
                continue
            if p < P_SIG:
                ax = ATLAS_AXIS.get(lab, "other")
                sig_axes.add(ax)
                sig_labels.append("%s(%.1f)" % (lab, z))
        # 方向冲突：跨性状 z 符号不一致且都显著
        zs = [zdf.loc[s, l] for l in labels
              if np.isfinite(pdf.loc[s, l]) and pdf.loc[s, l] < P_SIG]
        if zs and min(zs) < 0 < max(zs):
            opp = True
        row["n_axes_sig"] = len(sig_axes)
        row["axes"] = ",".join(sorted(sig_axes))
        row["sig_traits"] = ";".join(sig_labels)
        row["opposite_direction"] = bool(opp)
        row["in_mhc"] = bool(mhc_flag.get(s, False))
        if opp:
            row["locus_class"] = "opposite"
        elif len(sig_axes) >= 3:
            row["locus_class"] = "broad"
        elif len(sig_axes) == 2:
            row["locus_class"] = "intermediate"
        else:
            row["locus_class"] = "pathogen_specific"
        rows.append(row)
    atlas = pd.DataFrame(rows)
    # 归属性状（来自 clump）
    own = {}
    for lab, d in leads.items():
        for s in d["SNP"]:
            own.setdefault(s, []).append(lab)
    atlas["discovered_in"] = atlas["SNP"].map(lambda s: ",".join(own.get(s, [])))
    # 位点归并：同染色体且物理距离 ≤500 kb 的 lead SNP 视为同一区域（plink 已按
    # r²0.01/500kb 去冗余，但不同性状独立 clump 可能各自选出该区域内不同 index SNP，
    # 不归并会把同一区域重复计数为多个"独立共享位点"）。
    chr_bp = {s: posmap.get(s, (None, None)) for s in atlas["SNP"]}
    atlas["CHR"] = atlas["SNP"].map(lambda s: chr_bp[s][0])
    atlas["BP"] = atlas["SNP"].map(lambda s: chr_bp[s][1])
    region_ids = {}
    rid = 0
    for chrom in sorted(set(c for c in atlas["CHR"] if pd.notna(c))):
        sub = atlas[(atlas["CHR"] == chrom) & atlas["BP"].notna()].sort_values("BP")
        last_bp = None
        for _, r in sub.iterrows():
            if last_bp is None or (r["BP"] - last_bp) > 500_000:
                rid += 1
                last_bp = r["BP"]
            region_ids[r["SNP"]] = "R%02d_chr%d_%dMb" % (rid, int(chrom), int(r["BP"] // 1e6))
    atlas["region_id"] = atlas["SNP"].map(region_ids)
    # 基因注释（FinnGen nearest_genes，hg38；由 work/code/16_annotate_leads.py 落盘）
    ann_path = OUT / "lead_gene_annotation.tsv"
    if ann_path.exists():
        ann = pd.read_csv(ann_path, sep="\t")
        atlas = atlas.merge(ann[["SNP", "nearest_genes"]], on="SNP", how="left")
    atlas = atlas.sort_values(["locus_class", "n_axes_sig"], ascending=[True, False])
    atlas.to_csv(OUT / "locus_atlas.csv", index=False)
    zdf.to_csv(OUT / "locus_by_trait_z.csv")
    pdf.to_csv(OUT / "locus_by_trait_p.csv")
    cls = atlas["locus_class"].value_counts().rename_axis("locus_class").reset_index(name="n")
    cls.to_csv(OUT / "locus_class_summary.csv", index=False)
    print()
    print("[atlas] 位点分类：")
    print(cls.to_string(index=False))
    print()
    cols = ["SNP", "locus_class", "n_axes_sig", "axes", "in_mhc", "nearest_genes", "discovered_in"]
    cols = [c for c in cols if c in atlas.columns]
    print(atlas[cols].to_string(index=False))
    print()
    nonmhc = atlas[(atlas["locus_class"] != "pathogen_specific") & (~atlas["in_mhc"])]
    print("[atlas] 非 MHC 跨轴位点（跨病原共享的独立证据）: %d 个" % len(nonmhc))
    if len(nonmhc):
        print(nonmhc[cols].to_string(index=False))


if __name__ == "__main__":
    main()
