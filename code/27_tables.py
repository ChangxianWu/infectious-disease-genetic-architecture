# -*- coding: utf-8 -*-
""" 27_tables.py — 稿件表格构建（全部由产物 CSV 直读，规则 11/33：禁手抄）

表设计：
  Table 1  终点面板与功效（25 FinnGen 终点 + 2 HGI 层；轴、病例数、h²、mean χ²、是否过功效门）
  Table 2  跨病原共享区域（4 个区域；lead SNP、最近基因、覆盖轴、区内 LD、区内共定位证据）
  Table 3  共定位汇总（同源/跨源分层；share/suggest/not_support 计数 + 显著对明细）
  Table 4  MHC 敏感性（45 对中变化最大的 10 对；含/不含 MHC 的 rg 与 q 并列）
输出：work/output/table1.csv ... table4.csv
用法：python 27_tables.py
"""
import csv
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = Path(r"/path/to/project")
OUT = ROOT / "work" / "output"
DATA = ROOT / "work" / "data"

sys.path.insert(0, str(ROOT / "work" / "code"))
from constants import FINNGEN_ENDPOINTS  # noqa: E402

AXIS_OF = {
    "INTESTINAL": "Bacterial", "GASTRO_NOS": "Bacterial", "PYELONEPHR": "Bacterial",
    "ERYSIPELAS": "Bacterial", "BACT_NOS": "Bacterial", "PNEUMOBACT": "Bacterial",
    "SEPSIS": "Bacterial", "OSTEOMYEL": "Bacterial", "MENINGBACT": "Bacterial",
    "ENDOCARD": "Bacterial", "PNEUMONIA_MIXED": "Mixed",
    "COVID19": "Viral", "INFLUENZA": "Viral", "ZOSTER": "Viral", "VIRALPNEUMO": "Viral",
    "HERPES_SIMPLEX": "Viral", "VIRAL_WARTS": "Viral", "EBV": "Viral",
    "VIRAL_CNS": "Viral", "VARICELLA": "Viral",
    "MYCOSES": "Fungal", "CANDIDIASIS": "Fungal", "DERMATOPHYTOSIS": "Fungal",
    "INFESTATIONS": "Parasitic", "PROTOZOAL": "Parasitic",
    "HGI_C2": "Viral (external)", "HGI_A2": "Viral (external)",
}


def rd(p, sep=","):
    if not Path(p).exists():
        print("[tables] 缺文件 %s" % Path(p).name)
        return []
    with open(p, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=sep))


def wr(name, header, rows):
    p = OUT / name
    # 表 CSV 一律**无 BOM**：BOM 会被下游解析器（R read.csv 等）当作列名首字符，
    # 造成"列不存在"的假失败（本课题 data_audit 第 2 层实际踩到）。
    with open(p, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print("[tables] %s: %d 行" % (name, len(rows)))


# ---------- Table 1：终点面板与功效 ----------
man = {r["phenocode"]: r for r in rd(DATA / "finngen_R12_pheno_n.tsv", "\t")}
MAT = rd(OUT / "ldsc_rg" / "ldsc_matrix_2026-10-11.csv")
HGI = rd(OUT / "hgi_external_pairs.csv")

t1 = []
for lab, ph in FINNGEN_ENDPOINTS.items():
    m = man.get(ph)
    # 取该性状在任一对中出现时的 h2 / mean chisq（主矩阵内 h2 跨对恒定，见台账 §三）
    h2 = mc = h2p = ""
    for r in MAT:
        for a, b in (("axis1", "h2_1"), ("axis2", "h2_2")):
            if r.get(a) == lab and r.get(b):
                h2 = "%.5f" % float(r[b])
                mc = "%.3f" % float(r["mean_chisq_1" if b == "h2_1" else "mean_chisq_2"])
                h2p = "%.3g" % float(r["h2_1_p" if b == "h2_1" else "h2_2_p"])
                break
        if h2:
            break
    t1.append([lab, AXIS_OF.get(lab, ""), m["num_cases"] if m else "", m["num_controls"] if m else "",
               h2, mc, h2p,
               "yes" if h2 and float(h2p) < 0.05 and float(mc) >= 1.02 else "no"])
for lab in ("HGI_C2", "HGI_A2"):
    h2 = mc = ""
    for r in HGI:
        if r["axis2"] == lab and r.get("h2_2"):
            h2 = "%.5f" % float(r["h2_2"])
            mc = "%.3f" % float(r["mean_chisq_2"])
            break
    t1.append([lab, AXIS_OF[lab], "", "", h2, mc, "",
               "yes" if len([r for r in HGI if r["power_ok"] == "True" and r["axis2"] == lab]) else "no"])
wr("table1.csv", ["Trait", "Pathogen_axis", "Cases", "Controls", "h2_observed",
                  "Mean_chisq", "h2_p_one_sided", "Passed_power_gate"], t1)

# ---------- Table 2：跨病原共享区域 ----------
L = rd(OUT / "locus_atlas.csv")
C = rd(OUT / "coloc_results.csv")
REGIONS = {}
for r in L:
    if r["locus_class"] == "pathogen_specific":
        continue
    REGIONS.setdefault(r["region_id"], []).append(r)
t2 = []
for rid in sorted(REGIONS):
    rs = REGIONS[rid]
    snps = sorted(x["SNP"] for x in rs)
    genes = sorted(set(g.strip() for x in rs for g in str(x["nearest_genes"]).split(",") if g.strip()))
    axes = sorted(set(a for x in rs for a in str(x["axes"]).split(",") if a))
    mhc = any(x["in_mhc"] == "True" for x in rs)
    cls = "/".join(sorted(set(x["locus_class"] for x in rs)))
    traits = sorted(set(t.split("(")[0] for x in rs for t in str(x["sig_traits"]).split(";") if "(" in t))
    sub = [c for c in C if c["locus"] in snps]
    # **按 trait 对去重**：同一区域内多个 lead SNP 会对同一 trait 对重复测试（实测 CRP 区
    # 两个 lead 对 BACT_NOS×PNEUMOBACT 的 PP.H4 逐位相同 = 同一次分析）；不去重会把 4 对
    # 写成 8 条、把 1 对写成 4 条（MHC）。取每组最大 PP.H4（保守），存 (PP, t1, t2) 以便排序。
    best = {}
    for c in sub:
        k = tuple(sorted([c["trait1"], c["trait2"]]))
        v = float(c["PP_H4"])
        if k not in best or v > best[k][0]:
            best[k] = (v, c["trait1"], c["trait2"])
    share = [v for v in best.values() if v[0] > 0.8]
    share.sort(key=lambda x: -x[0])
    share_txt = "; ".join("%s x %s (%.3f)" % (t1, t2, v) for v, t1, t2 in share[:8])
    t2.append([rid, ",".join(snps), len(snps), ",".join(genes), ",".join(axes),
               "yes" if mhc else "no", cls, len(traits), len(best), len(sub),
               sum(1 for v in best.values() if v[0] > 0.8), share_txt,
               "%.3f" % max([float(c["PP_H4"]) for c in sub] or [0])])
wr("table2.csv", ["Region", "Lead_SNPs", "N_lead_SNPs", "Nearest_genes", "Pathogen_axes", "In_MHC",
                  "Classification", "N_traits_significant", "N_coloc_pairs_region_level",
                  "N_coloc_tests_variant_level", "N_coloc_share_region_level",
                  "Colocalised_pairs_examples", "Max_PP_H4"], t2)

# ---------- Table 3：共定位汇总 ----------
REGION_OF = {r["SNP"]: r["region_id"] for r in L}
REGNAME = {"R02_chr1_159Mb": "CRP", "R21_chr19_49Mb": "FUT2",
           "R14_chr9_136Mb": "ABO", "R07_chr6_32Mb": "MHC"}
# 跨源 = 恰一侧为外部性状（XOR）；与 25 号 canonical、30 号复算同口径（2026-10-29 统一）。
# 旧口径（任一 HGI）把 HGI_A2×HGI_C2 算作跨源，而 A2 是 C2 的严重度子集、属表型重叠。
XS = lambda c: c["trait1"].startswith("HGI") != c["trait2"].startswith("HGI")


def dedup_region(rows):
    """同一 (区域, trait 对) 只保留 PP.H4 最大的一行（同区多 lead = 同一次分析的重复）。"""
    best = {}
    for c in rows:
        k = (REGION_OF[c["locus"]], tuple(sorted([c["trait1"], c["trait2"]])))
        if k not in best or float(c["PP_H4"]) > float(best[k]["PP_H4"]):
            best[k] = c
    return list(best.values())


def band_counts(rows):
    return (len(rows),
            sum(1 for c in rows if float(c["PP_H4"]) > 0.8),
            sum(1 for c in rows if 0.5 < float(c["PP_H4"]) <= 0.8),
            sum(1 for c in rows if float(c["PP_H4"]) <= 0.5))


CR = dedup_region(C)
t3 = []
t3.append(["--- variant level (each lead variant counted separately) ---", "", "", "", ""])
for name, sub in (("All tests (coloc.abf)", C),
                  ("FinnGen-FinnGen (same source)", [c for c in C if not XS(c)]),
                  ("FinnGen-HGI (cross source)", [c for c in C if XS(c)])):
    n, sh, sg, ns = band_counts(sub)
    t3.append([name, n, sh, sg, ns])
t3.append(["--- region level (one test per region and trait pair; primary) ---", "", "", "", ""])
for name, sub in (("All tests (coloc.abf)", CR),
                  ("FinnGen-FinnGen (same source)", [c for c in CR if not XS(c)]),
                  ("FinnGen-HGI (cross source)", [c for c in CR if XS(c)])):
    n, sh, sg, ns = band_counts(sub)
    t3.append([name, n, sh, sg, ns])
t3.append(["--- region level, per region ---", "", "", "", ""])
for rid in sorted(REGNAME):
    for tag, sel in (("same source", lambda c: not XS(c)), ("cross source", XS)):
        sub = [c for c in CR if REGION_OF[c["locus"]] == rid and sel(c)]
        if not sub:
            continue
        n, sh, sg, ns = band_counts(sub)
        t3.append(["%s (%s)" % (REGNAME[rid], tag), n, sh, sg, ns])
wr("table3.csv", ["Stratum", "N_pairs", "PP_H4_gt_0.8", "PP_H4_0.5_to_0.8", "PP_H4_le_0.5"], t3)

# ---------- Table 4：MHC 敏感性（变化最大的 10 对） ----------
H = [r for r in rd(OUT / "mhc_sensitivity.csv") if r["status"] == "ok"]
H.sort(key=lambda r: -abs(float(r["mhc_effect"])))
t4 = [[r["pair"], "%.3f" % float(r["rg_excl_mhc"]), "%.3g" % float(r["q_excl_mhc"]),
       "%.3f" % float(r["rg_incl_mhc"]), "%.3g" % float(r["q_incl_mhc"]),
       "%+.3f" % float(r["mhc_effect"]), "yes" if r["flip_bh05"] == "True" else "no"]
      for r in H[:10]]
wr("table4.csv", ["Trait_pair", "rg_noMHC", "q_noMHC", "rg_withMHC", "q_withMHC",
                  "Delta_rg", "BH_conclusion_flipped"], t4)
