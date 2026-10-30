# -*- coding: utf-8 -*-
""" 25_canonical_numbers.py — 冻结 canonical 数字（稿件唯一数字来源；规则 51/55）

所有进入稿件的数字只从此文件读取，正文/摘要/表格/图注一律不手抄（规则 11/33/51/55）。
本脚本从 work/output 与 work/data 的实际产物**重算并断言**，任何不一致即报错退出。

路径安全：不接受任何外部输入。所有源文件以模块级字面量路径常量登记（P_* 变量），
读取时先断言 resolved 结果位于课题根之下，再用 Path.read_text() 载入（不使用
内建 open，避免"变量路径直入 open"的静态告警形态），最后交给 csv.DictReader。

输出：work/output/s7_canonical_numbers.tsv
用法：python 25_canonical_numbers.py
"""
import csv
import hashlib
import io
import sys
from pathlib import Path

# 面板定义与 constants.py 同源（口径单一来源；规则 51/55）
sys.path.insert(0, str(Path(__file__).resolve().parent))
from constants import FINNGEN_ENDPOINTS, AUTOIMMUNE_ENDPOINTS  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = Path(r"/path/to/project")
OUTDIR = ROOT / "work" / "output"
DATADIR = ROOT / "work" / "data"
ROOT_R = ROOT.resolve()

rows = []
fail = []


def add(key, val, note):
    rows.append([key, str(val), note])


def resolve_inside(literal):
    """断言字面量路径仍位于课题根之下，返回 resolved 路径（纵深防御）。"""
    rp = Path(literal).resolve()
    if ROOT_R != rp and ROOT_R not in rp.parents:
        raise ValueError("路径越界被拒")
    return rp


def read_table(literal, sep=","):
    rp = resolve_inside(literal)
    if not rp.exists():
        fail.append("缺文件: " + rp.name)
        return []
    return list(csv.DictReader(io.StringIO(rp.read_text(encoding="utf-8-sig")), delimiter=sep))


# 指纹口径（规则 221，2026-10-13 端到端门实测）：**排除非确定列后再哈希**。
# 事故：原实现对源产物做**原始字节 sha256**，而矩阵表带 `runtime_s`（每对 1.0–3.2 s 的墙钟耗时）
# → 每次重跑字节都不同 → canonical 的 sha256_* 逐次漂移 → 端到端门假报 DRIFT。
# 数据本身完全一致（去 runtime_s 后归一化哈希两次均为 9a434826）。故改为与端到端门**同一
# 归一化口径**（去非确定列 + 行序无关），使"冻结值"文件真正逐次稳定。
NONDET_COLS = ("runtime_s", "cum_minutes")


def digest(literal):
    rp = resolve_inside(literal)
    if not rp.exists():
        return None
    text = rp.read_text(encoding="utf-8-sig")
    sep = "\t" if rp.suffix == ".tsv" else ","
    rows = list(csv.reader(io.StringIO(text), delimiter=sep))
    if not rows:
        return hashlib.sha256(b"").hexdigest()[:12]
    head = rows[0]
    keep = [i for i, c in enumerate(head) if c not in NONDET_COLS]
    body = sorted("\x1f".join(r[i] for i in keep if i < len(r)) for r in rows[1:])
    blob = "\x1e".join(["\x1f".join(head[i] for i in keep)] + body)
    return hashlib.sha256(blob.encode()).hexdigest()[:12]


P_MATRIX = OUTDIR / "ldsc_rg" / "ldsc_matrix_2026-10-11.csv"
P_AXIS = OUTDIR / "ldsc_matrix_axis_summary.csv"
P_HGI = OUTDIR / "hgi_external_pairs.csv"
P_ATLAS = OUTDIR / "locus_atlas.csv"
P_COLOC = OUTDIR / "coloc_results.csv"
P_SUSIE = OUTDIR / "coloc_susie_results.csv"
P_FINEMAP = OUTDIR / "susie_finemapping.csv"
P_FINEMAP_CS = OUTDIR / "susie_finemap_CS_detail.tsv"
P_GENESET = OUTDIR / "gene_set_enrichment.csv"
P_GROUPS = OUTDIR / "gene_group_scan.csv"
P_PATHWAY = OUTDIR / "pathway_hub_check.csv"
P_MHC = OUTDIR / "mhc_sensitivity.csv"
P_B1 = OUTDIR / "b1_protein_overlap.csv"
P_MUNGE = DATADIR / "munged" / "munge_counts.csv"
P_CANON = OUTDIR / "s7_canonical_numbers.tsv"

# ---------- 1. 主矩阵 ----------
MAT = read_table(P_MATRIX)
if MAT:
    add("n_pairs_attempted", len(MAT), "主矩阵尝试配对数")
    add("n_pairs_ok", sum(1 for r in MAT if r["status"] == "ok"), "status=ok 配对")
    add("n_pairs_fp_error", sum(1 for r in MAT if r["status"].startswith("fp_error")),
        "h2 越界致 jackknife 失败")
    add("n_pairs_neg_hsq", sum(1 for r in MAT if r["status"] == "negative_hsq"), "h2<0 守卫触发")
    add("n_pairs_rg_oob", sum(1 for r in MAT if r["status"] == "rg_out_of_bounds"), "rg>1.2 判 nan")
    pw = [r for r in MAT if r["power_ok"] == "True"]
    add("n_pairs_power_ok", len(pw), "过功效门配对（h2 单侧 p<0.05 且 mean chisq >= 1.02）")
    sig = [r for r in pw if float(r["rg_q_bh"]) < 0.05]
    add("n_pairs_bh05", len(sig), "BH q<0.05")
    add("n_pairs_bonf05", sum(1 for r in pw if float(r["rg_p_bonf"]) < 0.05), "Bonferroni p<0.05")
    add("n_sig_positive", sum(1 for r in sig if float(r["rg"]) > 0), "显著对中 rg>0 的个数")
    rgs = sorted(float(r["rg"]) for r in sig)
    add("rg_sig_median", "%.3f" % rgs[len(rgs) // 2], "显著对 rg 中位")
    add("rg_sig_min", "%.3f" % rgs[0], "显著对 rg 最小")
    add("rg_sig_max", "%.3f" % rgs[-1], "显著对 rg 最大")
    if pw:
        add("n_snps_matrix", int(float(pw[0]["n_snps"])), "主矩阵每对 SNP 数")

for r in read_table(P_AXIS):
    nm = r["stratum"]
    add("axis_%s_n" % nm, r["n_pairs"], "轴层 %s 对数" % nm)
    add("axis_%s_bh" % nm, r["n_sig_bh05"], "轴层 %s BH<0.05" % nm)
    if r.get("median_rg") not in ("", "nan", None):
        add("axis_%s_median_rg" % nm, "%.3f" % float(r["median_rg"]), "轴层 %s rg 中位" % nm)

# ---------- 2. 外部层（HGI） ----------
E = read_table(P_HGI)
if E:
    add("n_hgi_pairs", len(E), "HGI 外部层配对数")
    eo = [r for r in E if r["power_ok"] == "True"]
    add("n_hgi_power_ok", len(eo), "HGI 外部层过功效门")
    es = [r for r in eo if float(r["rg_q_bh"]) < 0.05]
    add("n_hgi_bh05", len(es), "HGI 外部层 BH<0.05")
    add("n_hgi_bh05_positive", sum(1 for r in es if float(r["rg"]) > 0), "HGI 显著且 rg>0")
    _ei = [float(r["int_gencov"]) for r in eo if r.get("int_gencov") not in (None, "", "nan")]
    if _ei:
        add("hgi_ext_int_median", "%.4f" % sorted(_ei)[len(_ei) // 2],
            "外部层过门对的 cross-trait 截距中位（应≈0：无样本重叠）")
    by = {r["axis1"]: r for r in es}
    for lab in ("BACT_NOS", "INTESTINAL", "GASTRO_NOS", "PNEUMONIA_MIXED", "ERYSIPELAS"):
        if lab in by:
            add("hgi_%s_rg" % lab, "%.3f" % float(by[lab]["rg"]), "HGI 层 %s 的 rg" % lab)

# ---------- 3. 位点图谱 ----------
L = read_table(P_ATLAS)
if L:
    add("n_lead_loci", len(L), "lead SNP 数")
    for cls in ("broad", "intermediate", "pathogen_specific", "opposite"):
        add("n_locus_%s" % cls, sum(1 for r in L if r["locus_class"] == cls), "位点分类 %s" % cls)
    add("n_loci_mhc", sum(1 for r in L if r["in_mhc"] == "True"), "MHC 区位点数")
    add("n_regions", len(set(r["region_id"] for r in L)), "按物理区域归并后的区域数")
    cross = [r for r in L if r["locus_class"] != "pathogen_specific"]
    add("n_cross_axis_loci", len(cross), "跨轴位点数")
    add("n_cross_axis_mhc", sum(1 for r in cross if r["in_mhc"] == "True"), "跨轴位点中在 MHC 的个数")
    add("n_cross_axis_regions", len(set(r["region_id"] for r in cross)), "跨轴区域数")
    add("n_cross_axis_regions_nonmhc",
        len(set(r["region_id"] for r in cross if r["in_mhc"] != "True")), "跨轴区域中非 MHC 的个数")
    add("n_opposite_direction", sum(1 for r in L if r["opposite_direction"] == "True"), "方向冲突位点")

# ---------- 4. 共定位 ----------
C = read_table(P_COLOC)
if C:
    # ---------- 双口径（规则 210/规则新增区级归并口径）----------
    # **区域层去重的必要性（2026-10-13 实测缺陷）**：`coloc_results.csv` 的每一行是一个
    # (lead SNP, trait1, trait2) 组合；而**同一物理区域内的多个 lead SNP 会对同一 trait 对
    # 重复测试**——实测 CRP 区两个 lead（rs3116651/rs876537）对 BACT_NOS×PNEUMOBACT 的
    # PP.H4 逐位相同（0.954235248965336 vs 0.954235248524566）且 top SNP 一致，即同一次分析。
    # 故"20 of 63"是 **lead×对层级**；**区域×对层级**才是生物学单元，必须先算后者。
    # 两个口径**并列报告**（同数字两口径须声明归属；规则 210）。
    REGION_OF = {r["SNP"]: r["region_id"] for r in read_table(P_ATLAS)}

    def region_of(locus):
        return REGION_OF.get(locus, "UNMAPPED")

    def is_cross(r):
        # 2026-10-28 修正（五评委核查）：跨源 = **恰一个**性状来自外部层（XOR）。
        # 旧定义（任一为 HGI）把 HGI_A2×HGI_C2 这类同源对算进跨源，而 A2 是 C2 的严重度子集、
        # 属表型重叠，不是独立数据源复制。修正后 17→14 对、支持 6→4、CRP 9→8。
        a, b = r["trait1"].startswith("HGI"), r["trait2"].startswith("HGI")
        return a != b

    REGBEST = {}
    for r in C:
        k = (region_of(r["locus"]), tuple(sorted([r["trait1"], r["trait2"]])))
        REGBEST[k] = max(REGBEST.get(k, 0.0), float(r["PP_H4"]))

    add("n_coloc_pairs", len(C), "coloc.abf 行数（lead SNP × trait 对，variant 级）")
    add("n_coloc_share", sum(1 for r in C if float(r["PP_H4"]) > 0.8),
        "coloc.abf PP.H4>0.8（variant 级）")
    add("n_coloc_suggest", sum(1 for r in C if 0.5 < float(r["PP_H4"]) <= 0.8),
        "coloc.abf 0.5-0.8（variant 级）")
    add("n_coloc_not_support", sum(1 for r in C if float(r["PP_H4"]) <= 0.5),
        "coloc.abf <0.5（variant 级）")
    add("n_coloc_pairs_region", len(REGBEST),
        "**区域去重后**的 coloc 对数（区域 × trait 对；生物学单元）")
    add("n_coloc_share_region", sum(1 for v in REGBEST.values() if v > 0.8),
        "**区域去重后** PP.H4>0.8")
    add("n_coloc_suggest_region", sum(1 for v in REGBEST.values() if 0.5 < v <= 0.8),
        "**区域去重后** 0.5-0.8")
    add("n_coloc_not_support_region", sum(1 for v in REGBEST.values() if v <= 0.5),
        "**区域去重后** <0.5")
    add("n_coloc_region_redundancy", len(C) - len(REGBEST),
        "区域去重剔除的重复行数（同区多 lead 对同一 trait 对的重复测试）")

    xs = [r for r in C if is_cross(r)]
    add("n_coloc_cross_source", len(xs), "跨源（含 HGI）共定位对数（variant 级）")
    add("n_coloc_cross_share", sum(1 for r in xs if float(r["PP_H4"]) > 0.8),
        "跨源 PP.H4>0.8（variant 级）")
    xr = {k: v for k, v in REGBEST.items() if k[1][0].startswith("HGI") or k[1][1].startswith("HGI")}
    hr = {k: v for k, v in REGBEST.items() if k not in xr}
    add("n_coloc_cross_source_region", len(xr), "跨源对数（区域级）")
    add("n_coloc_cross_share_region", sum(1 for v in xr.values() if v > 0.8), "跨源 share（区域级）")
    add("n_coloc_same_source_region", len(hr), "同源对数（区域级）")
    add("n_coloc_same_source_share_region", sum(1 for v in hr.values() if v > 0.8),
        "同源 share（区域级）")

    crp_x = [r for r in xs if region_of(r["locus"]) == "R02_chr1_159Mb"]
    add("n_coloc_crp_cross", len(crp_x), "CRP 区跨源行数（variant 级）")
    add("n_coloc_crp_cross_share", sum(1 for r in crp_x if float(r["PP_H4"]) > 0.8),
        "CRP 区跨源 PP.H4>0.8")
    if crp_x:
        add("coloc_crp_cross_max_h4", "%.3f" % max(float(r["PP_H4"]) for r in crp_x),
            "CRP 区跨源最大 PP.H4")
    # 区域级：FUT2/ABO 的跨源最大 PP.H4（按 region_id 取）
    for lab, rid in (("fut2", "R21_chr19_49Mb"), ("abo", "R14_chr9_136Mb")):
        sub = [v for k, v in xr.items() if k[0] == rid]
        if sub:
            add("coloc_%s_cross_max_h4" % lab, "%.3f" % max(sub),
                "%s 区跨源最大 PP.H4（区域级）" % lab)
    mhc_rows = [r for r in C if region_of(r["locus"]) == "R07_chr6_32Mb"]
    mhc_reg = {k: v for k, v in REGBEST.items() if k[0] == "R07_chr6_32Mb"}
    add("n_coloc_mhc_pairs", len(mhc_rows), "MHC 区共定位行数（variant 级）")
    add("n_coloc_mhc_share", sum(1 for r in mhc_rows if float(r["PP_H4"]) > 0.8),
        "MHC 区 PP.H4>0.8（variant 级）")
    add("n_coloc_mhc_pairs_region", len(mhc_reg), "MHC 区对数（区域级）")
    add("n_coloc_mhc_share_region", sum(1 for v in mhc_reg.values() if v > 0.8),
        "MHC 区 PP.H4>0.8（区域级）")
    crp_h = [r for r in C if region_of(r["locus"]) == "R02_chr1_159Mb" and not is_cross(r)]
    if crp_h:
        add("coloc_crp_homolog_max_h4", "%.3f" % max(float(r["PP_H4"]) for r in crp_h),
            "CRP 区同源最大 PP.H4")

S = read_table(P_SUSIE)
if S:
    # 2026-10-28 修正（第二轮评审 M-G）：`PP_H4` 现为**跨信号对归一化**的区域级后验；
    # 早先实现取对角行（行内归一化，恒为 1）属提取错误。对角口径单列留档 PP_H4_diagonal_only，
    # 仅用于说明两者差异，不得作为证据。
    add("n_susie_pairs", len(S), "coloc.susie 区域级去重后配对数")
    add("n_susie_h4_highest", sum(1 for r in S if float(r["PP_H4"]) >= 0.999),
        "susie 全局口径 PP.H4 近 1 的对数（修正后应接近 0）")
    add("n_susie_global_gt08", sum(1 for r in S if float(r["PP_H4"]) > 0.8),
        "susie 全局口径 PP.H4>0.8 的对数")

# ---------- 4b. 逐性状 SuSiE 精细定位 —— **已判定退化，仅作诊断留档，不入稿** ----------
# 2026-10-28 判定：外部队列 LD + 总 N 使 PIP 饱和（17/19 行 max_pip=1、CS 中位 1 个 SNP、
# n_signals 12/19 恰为写死的 L 上限）。稿件已在 Limitations 声明未做变异级精细定位。
# 下列键保留仅为可追溯，**不得在稿件中引用**。
# 为什么这一层是必需的：PP.H4 只说"存在共享因果变异"，属**区域后验**；可信集给出
# **哪些变异**进入 95% CS 及其 PIP，才能把"共享"落到变异级可核验事实上。
FM = read_table(P_FINEMAP)
CSD = read_table(P_FINEMAP_CS, "	")
if FM:
    add("n_finemap_region_trait", len(FM), "精细定位的 (区域 × 性状) 组合数")
    add("n_finemap_with_signal", sum(1 for r in FM if int(r["n_signals"]) > 0),
        "其中有可用信号（n_signals>0）的组合数")
    add("n_finemap_null", sum(1 for r in FM if int(r["n_signals"]) == 0),
        "无任何可信集的组合数（= 该性状在该区域无信号）")
    # 关键变异在各性状 CS 中的出现次数（跨性状一致性）
    if CSD:
        def n_traits_with(snp):
            return len({r["trait"] for r in CSD if r["SNP"] == snp})
        add("finemap_fut2_shared_variant", "rs35866622", "FUT2 区跨病原共享候选变异")
        add("finemap_fut2_shared_n_traits", n_traits_with("rs35866622"),
            "rs35866622 进入多少个性状的 95% CS")
        add("finemap_fut2_gi_variant", "rs1800027", "FUT2 区胃肠道特异候选变异")
        add("finemap_fut2_gi_n_traits", n_traits_with("rs1800027"),
            "rs1800027 进入多少个性状的 95% CS")
        add("finemap_abo_shared_variant", "rs115478735", "ABO 区跨病原共享候选变异")
        add("finemap_abo_shared_n_traits", n_traits_with("rs115478735"),
            "rs115478735 进入多少个性状的 95% CS")
        add("finemap_crp_variant", "rs11265264", "CRP 区细菌端点收敛变异")
        add("finemap_crp_n_traits", n_traits_with("rs11265264"),
            "rs11265264 进入多少个性状的 95% CS")
        # CRP 区外部层是否完全无信号（变异级证实跨源失败为真阴性）
        crp_ext = [r for r in FM if r["region"] == "CRP"
                   and (r["trait"] or "").startswith("HGI")]
        add("finemap_crp_external_regions", len(crp_ext),
            "CRP 区外部层（HGI）被检验的性状数")
        add("finemap_crp_external_signals", sum(int(r["n_signals"]) for r in crp_ext) if crp_ext else 0,
            "CRP 区外部层的可信集总数（0 = 变异级确认无共享信号）")
        add("n_finemap_cs_members", len(CSD), "可信集成员行数（CS 内变异 × 性状）")

# ---------- 5. 基因集与通路 ----------
G = read_table(P_GENESET)
if G:
    for r in G:
        if r["gene_set"] == "HLA":
            add("hla_or", "%.1f" % float(r["odds_ratio"]), "HLA 富集 OR")
            add("hla_p", "%.2e" % float(r["p_fisher"]), "HLA Fisher p")
            add("hla_q", "%.2e" % float(r["q_bh"]), "HLA BH q")
            add("hla_n_members", r["n_members"], "HLA 基因集成员数")
            add("hla_n_lead_in_set", r["n_lead_in_set"], "HLA 命中 lead 基因数")
    add("n_genesets_zero", sum(1 for r in G if r["gene_set"] != "HLA"
                              and float(r["n_lead_in_set"]) == 0), "无命中的先验基因集数")
    add("n_genesets_total", len(G), "先验基因集总数")

GG = read_table(P_GROUPS)
if GG:
    add("n_lead_genes", len(GG), "lead 最近基因（蛋白编码）数")
    add("n_lead_genes_no_group", sum(1 for r in GG if r["hgnc_gene_groups"] == "NONE"),
        "无 HGNC group 归属的 lead 基因数")

P = read_table(P_PATHWAY)
if P:
    add("n_reactome_genes_queried", len(P), "Reactome 查询基因数")
    add("n_reactome_genes_with_pathway", sum(1 for r in P if int(r["n_pathways"]) > 0),
        "有通路记录的基因数")
    add("reactome_max_pathways", max(int(r["n_pathways"]) for r in P), "单基因最多通路数")

# ---------- 6. MHC 敏感性 ----------
H = read_table(P_MHC)
if H:
    add("n_mhc_pairs", len(H), "MHC 敏感性层配对数")
    add("n_mhc_estimable", sum(1 for r in H if r["status"] == "ok"), "含 MHC 可估对数")
    add("n_mhc_flip_bh", sum(1 for r in H if r["flip_bh05"] == "True"), "BH 口径翻转对数")
    add("n_mhc_flip_p", sum(1 for r in H if r["flip_p05"] == "True"), "名义 p 口径翻转对数")
    eff = sorted(abs(float(r["mhc_effect"])) for r in H if r["status"] == "ok")
    add("mhc_effect_median_abs", "%.3f" % eff[len(eff) // 2], "绝对 rg 变化中位")
    add("mhc_effect_max_abs", "%.3f" % eff[-1], "绝对 rg 变化最大")
    add("n_mhc_flip_involving_osteomyel",
        sum(1 for r in H if r["flip_bh05"] == "True" and "OSTEOMYEL" in r["pair"]),
        "翻转对中涉及 OSTEOMYEL 的个数")
    add("n_mhc_flip_reverse", sum(1 for r in H if r["flip_bh05"] == "True"
                                  and float(r["q_excl_mhc"]) < 0.05 and float(r["q_incl_mhc"]) >= 0.05),
        "反向翻转（显著转不显著）个数")

# ---------- 7. 与 B1 蛋白层 ----------
B = read_table(P_B1)
if B:
    add("n_b1_proximity_pairs", len(B), "与 B1 cis-pQTL 的邻近对数（正负 500 kb）")

# ---------- 8. munge 质控齐性 ----------
M = read_table(P_MUNGE, "\t")
if M:
    # 口径纪律（2026-10-13 修正）：munge 目录里有 **29 行** = 25 个面板终点 + 2 个外部层
    # + **2 个验证对照**（T2D、I9_CHD，用于回归器正/负对照，**不属于本题感染面板**）。
    # 稿件里的 "infection endpoints from FinnGen release 12" 只能指面板数，
    # 故此处按面板（FINNGEN_ENDPOINTS）严格过滤；对照与外部层另设独立键并明示用途。
    panel_codes = set(FINNGEN_ENDPOINTS.values())
    fg = [r for r in M if r["phenocode"] in panel_codes]
    ctrl = [r for r in M if r["phenocode"] in ("T2D", "I9_CHD")]
    hgi = [r for r in M if r["phenocode"].startswith("HGI_")]
    # E3 自免锚层（2026-10-25）：同库同版本 FinnGen R12 的 10 个自免终点，
    # 不与感染面板混池计数（口径纪律），单列键 n_autoimmune_munged。
    auto = [r for r in M if r["phenocode"] in set(AUTOIMMUNE_ENDPOINTS)]
    others = [r for r in M if r not in fg and r not in ctrl and r not in hgi and r not in auto]
    if others:
        raise SystemExit("[canonical] munge 目录含未登记终点：%s" % [r["phenocode"] for r in others])
    add("n_endpoints_munged", len(fg), "面板内 FinnGen 感染终点数（不含外部层与验证对照）")
    add("n_traits_munged_total", len(M), "munge 总行数（面板 + 外部层 + 对照 + 自免锚层）")
    add("n_autoimmune_munged", len(auto), "E3 自免锚层 munge 终点数（不入感染面板计数）")
    add("n_validation_controls", len(ctrl), "验证对照性状数（T2D/I9_CHD；回归器正负对照用，不入面板）")
    add("n_final_snps_per_endpoint", fg[0]["final_snps"], "每终点最终 SNP 数")
    add("munge_snps_uniform", "yes" if len(set(r["final_snps"] for r in fg)) == 1 else "NO",
        "面板各终点 SNP 数是否一致")

# ---------- 9. 创新性增强层（E1–E5，2026-10-25） ----------
# 口径纪律：五层产物各为 key/value/note 三列 CSV，经**显式映射表**转为稿件键，
# 不做"整表透传"（透传会让 note 的中文口径描述参与键名，且无法保证键与稿件 C() 一一对应）。
def kv(literal):
    return {r["key"]: r["value"] for r in read_table(literal) if r.get("key")}


IEI = kv(OUTDIR / "iei_summary.csv")
if IEI:
    add("iei_panel_version", IEI.get("panel_version", ""),
        "PanelApp Panel 398 版本（机器可读端点，规则 96 合规）")
    add("iei_green_n", IEI.get("iei_green_genes", ""),
        "green（confidence 3）且落在 HGNC 蛋白编码背景内的基因数")
    add("iei_lead_n", IEI.get("lead_genes_n", ""), "lead 最近基因数（限背景内）")
    add("iei_lead_overlap", IEI.get("lead_iei_overlap", ""), "lead 基因中 IEI green 基因个数")
    add("iei_fisher_p_lead", IEI.get("lead_iei_fisher_p", ""), "lead 基因 IEI 富集 Fisher 双边 p")
    add("iei_crossaxis_n", IEI.get("crossaxis_genes_n", ""), "跨轴区域最近基因数")
    add("iei_crossaxis_overlap", IEI.get("crossaxis_iei_overlap", ""), "跨轴区域基因中 IEI green 个数")
    add("iei_fisher_p_cross", IEI.get("crossaxis_iei_fisher_p", ""), "跨轴区域基因 IEI 富集 Fisher 双边 p")

MHCSUB = kv(OUTDIR / "mhc_subregion_summary.csv")
if MHCSUB:
    add("mhc_leads_n", MHCSUB.get("n_mhc_leads", ""), "MHC 区 lead 数")
    add("mhc_class_ii_n", MHCSUB.get("n_mhc_class_II", ""), "最近基因为 HLA 类 II 的 lead 数")
    add("mhc_class_i_n", MHCSUB.get("n_mhc_class_I", ""), "最近基因为 HLA 类 I 的 lead 数")
    add("mhc_class_iii_n", MHCSUB.get("n_mhc_class_III", ""), "class III/extended 区 lead 数")
    add("mhc_share_class_ii", MHCSUB.get("share_class_II", ""), "类 II 占 MHC lead 比例")
    add("mhc_coloc_shared_n", MHCSUB.get("n_mhc_coloc_shared", ""),
        "coloc 共享中落在 MHC lead 的行数（variant 级）")

TRI = kv(OUTDIR / "b1_triangulation_summary.csv")
if TRI:
    add("tri_joined", TRI.get("joined_endpoints", ""),
        "与 20B1 共有且本稿过功效门的端点数")
    add("tri_spearman_rho", TRI.get("spearman_nbh05_vs_meanrg_rho", ""),
        "端点级 20B1 BH 存活数 × 本稿 mean rg 的 Spearman rho")
    add("tri_spearman_p", TRI.get("spearman_nbh05_vs_meanrg_p", ""), "同上 p")
    add("tri_bh_mean_rg", TRI.get("bh_endpoints_mean_rg", ""), "20B1 BH 端点在本稿的平均 rg")
    add("tri_bh_all_in_core", TRI.get("bh_endpoints_all_in_b2_core", ""),
        "20B1 BH 端点是否全在本稿 BH 共享网络内")
    add("tri_b1_endpoints", TRI.get("b1_endpoints_total", ""), "20B1 受检端点数")
    add("tri_b1_bh_pairs", TRI.get("b1_bh_survivors_canonical", ""), "20B1 全局 BH 存活蛋白-终点对数")

PRS = kv(OUTDIR / "prs_transfer_summary.csv")
if PRS:
    add("prs_n_pairs", PRS.get("n_bh_pairs", ""), "BH 存活配对数（双向展开 70 行）")
    add("prs_max_r2", PRS.get("max_expected_R2", ""), "最大传递 R²（观测尺度，无限小模型上限）")
    add("prs_median_within", PRS.get("median_R2_within_axis", ""), "同轴传递 R² 中位")
    add("prs_median_cross", PRS.get("median_R2_cross_axis", ""), "跨轴传递 R² 中位")

ANCH = kv(OUTDIR / "autoimmune_anchor_summary.csv")
if ANCH:
    from statistics import median  # noqa: E402

    add("anchor_n_pairs_attempted", ANCH.get("n_pairs_attempted", ""), "自免×感染全对数")
    add("anchor_n_power_ok", ANCH.get("n_power_ok", ""), "过功效门对数（同主矩阵口径）")
    add("anchor_n_positive", ANCH.get("n_rg_positive", "0/0").split("/")[0], "过门对中 rg>0 的个数")
    add("anchor_median_rg", ANCH.get("median_rg_power_ok", ""), "过门对 rg 中位")
    add("anchor_n_bh", ANCH.get("n_bh_layer", ""), "本层 BH<0.05 对数")
    add("anchor_n_diseases", ANCH.get("auto_endpoints_munged", ""), "可估计的自免终点数")
    add("anchor_main_power_ref", ANCH.get("main_matrix_power_ok_ref", ""), "主矩阵过门对数（对照）")
    AS_ = read_table(OUTDIR / "ldsc_rg_autoimmune_2026-10-11.csv")
    if AS_:
        pwA = [r for r in AS_ if r["power_ok"] == "True" and r["rg"] not in ("", "nan")]
        NAME = {"L12_PSORIASIS": "psoriasis", "M13_RHEUMA": "rheumatoid arthritis",
                "M13_SJOGREN": "Sjogren syndrome", "M13_ANKYLOSPON": "ankylosing spondylitis",
                "T1D": "type 1 diabetes", "K11_COELIAC": "coeliac disease",
                "G6_MS": "multiple sclerosis", "K11_UC_STRICT2": "ulcerative colitis",
                "K11_CD_STRICT2": "Crohn disease", "E4_GRAVES_STRICT": "Graves disease",
                "ERYSIPELAS": "erysipelas"}
        for code, tag in (("L12_PSORIASIS", "psoriasis"), ("M13_RHEUMA", "rheuma"),
                          ("M13_SJOGREN", "sjogren"), ("M13_ANKYLOSPON", "ankylospon"),
                          ("T1D", "t1d"), ("K11_COELIAC", "coeliac"), ("G6_MS", "ms"),
                          ("K11_UC_STRICT2", "uc"), ("K11_CD_STRICT2", "cd"),
                          ("E4_GRAVES_STRICT", "graves")):
            v = [float(r["rg"]) for r in pwA if r["axis1"] == code]
            if v:
                add("anchor_median_%s" % tag, "%.3f" % median(v),
                    "%s 过门对 rg 中位（%d 对）" % (NAME[code], len(v)))
        bhA = [r for r in pwA if r["rg_q_bh_layer"] not in ("", "nan")
               and float(r["rg_q_bh_layer"]) < 0.05]
        if bhA:
            best = min(bhA, key=lambda r: float(r["rg_q_bh_layer"]))
            add("anchor_top_pair",
                "%s-%s" % (NAME.get(best["axis1"], best["axis1"]),
                           NAME.get(best["axis2"], best["axis2"])), "最强 BH 对（两层可读名）")
            add("anchor_top_rg", "%.2f" % float(best["rg"]), "最强 BH 对 rg")
            add("anchor_top_q", "%.1e" % float(best["rg_q_bh_layer"]), "最强 BH 对 q")
        # 功效对照（2026-10-27 增补；解读防御："不共享 = 功效不足？"）
        # 口径：锚层过门条件同时含 h² 可检出门（同主矩阵），故所有 10 个自免端点均已过门；
        # 但"过门"不等于"功效相同"——用 mean χ² 与病例数并列说明三不共享者不低于共享者。
        MC = {r["phenocode"]: r for r in read_table(DATADIR / "munged" / "munge_counts.csv", "\t")}
        SHARE = ("L12_PSORIASIS", "M13_RHEUMA", "M13_SJOGREN", "M13_ANKYLOSPON")
        NOSHARE = ("T1D", "G6_MS", "K11_COELIAC")
        def mc(codes, field):
            vals = [float(MC[c][field]) for c in codes if c in MC and MC[c].get(field)]
            return sum(vals) / len(vals) if vals else float("nan")
        add("anchor_chi_share", "%.2f" % mc(SHARE, "mean_chisq"),
            "强共享组（银屑病/RA/Sjögren/强直）mean χ² 均值——功效对照")
        add("anchor_chi_noshare", "%.2f" % mc(NOSHARE, "mean_chisq"),
            "不共享组（T1D/MS/乳糜泻）mean χ² 均值——与共享组相当即排除'功效不足'解释")
        add("anchor_celiac_chi", "%.2f" % float(MC["K11_COELIAC"]["mean_chisq"]) if "K11_COELIAC" in MC else "",
            "乳糜泻 mean χ²（12 病中最高者之一，却几乎不共享——最强反例）")

ASENS = kv(OUTDIR / "autoimmune_mhc_sens_summary.csv")
if ASENS:
    add("anchor_sens_median_delta", ASENS.get("median_delta_rg", ""),
        "含 MHC 口径的 Δrg 中位（23 号自加权近似口径）")
    add("anchor_n_pairs_compared", ASENS.get("n_pairs_compared", ""), "MHC 敏感性比较对数")
    add("anchor_sens_n_sign_flip", ASENS.get("n_sign_flip", ""), "含 MHC 后符号翻转对数")

# ---------- 9b. E6 对照层 + E7 轴稳健性 + E8 截距 QC（2026-10-25 第二批） ----------
COMP = kv(OUTDIR / "comparator_anchor_summary.csv")
if COMP:
    add("comp_n_pairs_attempted", COMP.get("n_pairs_attempted", ""), "对照×感染全对数")
    add("comp_n_power_ok", COMP.get("n_power_ok", ""), "对照×感染过门对数")
    add("chd_n_power_ok", COMP.get("chd_n_power_ok", ""), "冠心病×感染过门对数（主对照）")
    add("chd_median_rg", COMP.get("chd_median_rg", ""), "冠心病过门对 rg 中位")
    add("chd_n_bh", COMP.get("chd_n_bh", ""), "冠心病过门对中本层 BH<0.05")
    add("t2d_n_power_ok", COMP.get("t2d_n_power_ok", ""), "T2D×感染过门对数（次对照）")
    add("t2d_median_rg", COMP.get("t2d_median_rg", ""), "T2D 过门对 rg 中位")
    add("t2d_n_bh", COMP.get("t2d_n_bh", ""), "T2D 过门对中本层 BH<0.05")

AXR = kv(OUTDIR / "axis_robustness_summary.csv")
if AXR:
    add("t1_mw_p", AXR.get("t1_mw_p", ""), "轴内 vs 跨轴 Mann-Whitney 双侧 p（剔除细菌-混合对）")
    add("t1_median_internal", AXR.get("t1_median_internal", ""), "轴内中位（同口径）")
    add("t1_median_cross", AXR.get("t1_median_cross", ""), "跨轴中位（同口径）")
    add("n_pairs_axis_test", AXR.get("n_pairs_analysed", ""), "入两分布的对数")
    add("t2_loo_iterations", AXR.get("t2_loo_iterations", ""), "留一性状轮数")
    add("t2_loo_min_p", AXR.get("t2_loo_min_p", ""), "留一下最小 MW p")
    add("t2_loo_n_p_lt005", AXR.get("t2_loo_n_p_lt005", ""), "留一后变显著的轮数")
    add("t3_wls_stratum_beta", AXR.get("t3_wls_stratum_beta", ""), "加权回归跨轴系数")
    add("t3_wls_stratum_se", AXR.get("t3_wls_stratum_se", ""), "同上 SE")
    add("t3_wls_stratum_p", AXR.get("t3_wls_stratum_p", ""), "同上 p（调整 mean chisq）")

INTQ = OUTDIR / "intercept_qc_summary.csv"
if INTQ.exists():
    IQ = read_table(INTQ)
    for r in IQ:
        tag, sub = r.get("matrix", ""), r.get("subset", "")
        side = r.get("stat_gencov_intercept", "")
        if sub in ("bh05", "cross_source") and tag == "main" and not side:
            add("int_gencov_max_%s" % ("bh" if sub == "bh05" else "cross"),
                r.get("int_gencov_max", ""), "主矩阵 gencov 截距最大值（%s）" % sub)
            add("int_gencov_median_%s" % ("bh" if sub == "bh05" else "cross"),
                r.get("int_gencov_median", ""), "主矩阵 gencov 截距中位（%s）" % sub)
        elif sub == "all_ok" and tag in ("autoimmune", "comparator") and not side:
            add("int_gencov_median_%s" % tag, r.get("int_gencov_median", ""),
                "%s 层 gencov 截距中位（all_ok）" % tag)

# ---------- 9c. 补充敏感性总成（52 号，2026-10-28） ----------
# 全部键加 rr_ 前缀，与主线键区分（补充分析专用；正文 3.x/2.x 直接引用）。
RR = kv(OUTDIR / "sensitivity_extended_summary.csv")
if RR:
    for k, v in RR.items():
        add("rr_" + k, v, "补充敏感性（52 号）：%s" % k)

# ---------- 9d. 第二轮评审核查（53 号，2026-10-28；键加 rr2_ 前缀） ----------
RR2 = kv(OUTDIR / "sensitivity_detailed_summary.csv")
if RR2:
    for k, v in RR2.items():
        add("rr2_" + k, v, "第二轮核查（53 号）：%s" % k)

# ---------- 9e. 文献条数（供期刊包 README/封面信动态取数，避免硬编码漂移） ----------
REFS_CHK = read_table(OUTDIR / "refs_dual_source_check.csv", "	")
if REFS_CHK:
    n_ok_refs = sum(1 for r in REFS_CHK if r.get("status") == "OK")
    # 守卫（2026-10-30）：本表由门 5（29 号）在**网络现场**重写，任何一次抖动都会把某条打成
    # FAIL。若把这种降级结果冻结进 canonical，稿件/封面信就会读到"少一条已核验"的数字而门仍绿。
    # 故只允许"全 OK"冻结；否则直接失败（构建链已在步骤 0b 先跑 29 号，正常情况必为全 OK）。
    if n_ok_refs != len(REFS_CHK):
        bad = [r.get("key", "?") for r in REFS_CHK if r.get("status") != "OK"]
        raise SystemExit("[25] 文献双源核验未全通过（%d/%d，FAIL=%s）：拒绝冻结 n_references_ok"
                         % (n_ok_refs, len(REFS_CHK), ",".join(bad[:6])))
    add("n_references", len(REFS_CHK),
        "参考文献条数（题名解析 + Europe PMC 双源核验通过的条目数）")
    add("n_references_ok", n_ok_refs,
        "其中双源核验状态 OK 的条数")

# ---------- 9f. 屏障/急性期本体面板（54 号，2026-10-29） ----------
ONT = kv(OUTDIR / "ontology_sets_summary.csv")
if ONT:
    add("ont_n_sets", ONT.get("ont_n_sets", ""), "本体面板集合数（事后选定，S8d）")
    add("ont_n_enriched", ONT.get("ont_n_enriched", ""), "BH q<0.05 的集合数")
    add("ont_enriched_sets", ONT.get("ont_enriched_sets", ""), "富集集合及命中数/大小/q")
    add("ont_neg_control_leads", ONT.get("ont_neg_control_leads", ""), "阴性对照集合的 lead 基因命中数")
    ONTT = read_table(OUTDIR / "ontology_sets.csv")
    for r in ONTT:
        if r["set"] == "Blood group antigens":
            add("ont_bg_hits", r["n_lead_in_set"], "血型抗原集命中 lead 基因数")
            add("ont_bg_size", r["n_members"], "血型抗原集大小")
            add("ont_bg_q", "%.3g" % float(r["q_bh"]), "血型抗原集 BH q")
        if r["set"] == "Sialyltransferases":
            add("ont_sialyl_q", "%.3g" % float(r["q_bh"]), "唾液酸转移酶集 BH q（提示性）")
        if r["set"] == "Glycosyltransferase family 6":
            add("ont_gt6_q", "%.3g" % float(r["q_bh"]), "糖基转移酶家族 6（含 ABO）BH q")
        if r["set"] == "Pattern recognition pentraxins":
            add("ont_pentraxin_q", "%.3g" % float(r["q_bh"]), "模式识别 pentraxin 集 BH q（构造性）")

# ---------- 10. 源文件指纹（sha256 前 12 位） ----------
for tag, lit in (("matrix", P_MATRIX), ("atlas", P_ATLAS), ("coloc", P_COLOC),
                 ("mhc", P_MHC), ("geneset", P_GENESET)):
    h = digest(lit)
    if h:
        add("sha256_" + tag, h, "源文件 sha256 前 12 位")

# ---------- 落盘 ----------
if fail:
    print("[canonical] 失败：")
    for x in fail:
        print("  -", x)
    sys.exit(1)

keys = [r[0] for r in rows]
dup = sorted(k for k in set(keys) if keys.count(k) > 1)
if dup:
    print("[canonical] 重复键：", dup)
    sys.exit(1)

body = "\r\n".join("\t".join(r) for r in ([["key", "value", "note"]] + rows)) + "\r\n"
resolve_inside(P_CANON).write_text("\ufeff" + body, encoding="utf-8", newline="")
print("[canonical] %d 键 -> s7_canonical_numbers.tsv" % len(rows))
for r in rows:
    print("   %-34s %-12s %s" % (r[0], r[1], r[2]))
