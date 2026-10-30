# -*- coding: utf-8 -*-
""" 53_sensitivity_detailed.py — 第二轮评审的量化核查（2026-10-28）

对应意见：
  M-B 轴比较单一主口径（主分析剔除 mixed；含 mixed 版本作敏感性）+ 三处数字对齐
  M-C 每次敏感性给效应量（Hodges-Lehmann）+ bootstrap 95% CI + 可检测差异；
      并交代加权回归的编码方向（为何原始中位跨轴更高而调整后系数为负）
  M-D 头条 coloc 数字扣除不予采信的 MHC 对
  M-H p–Z 守卫分母（限制前 / panel 限制后 / 最终集三个口径）
  次5 等位链方向：统计回归集中回文 SNP（A/T、C/G）比例，说明参考面板对齐已消除链歧义
  次7 gene-set OR 的 95% CI（Woolf，小计数下注明不稳）
输出：work/output/sensitivity_detailed_summary.csv
用法：python 53_sensitivity_detailed.py
"""
import csv
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATE, OUT  # noqa: E402

MATRIX = OUT / "ldsc_rg" / ("ldsc_matrix_%s.csv" % DATE)
SUM = OUT / "sensitivity_detailed_summary.csv"
MIXED = "PNEUMONIA_MIXED"


def hl_ci(a, b, rng, n_boot=5000):
    """Hodges-Lehmann 位移（跨轴 − 轴内）与 bootstrap 95% CI。"""
    a = np.asarray(a, float); b = np.asarray(b, float)
    diffs = (b[:, None] - a[None, :]).ravel()
    hl = float(np.median(diffs))
    boots = np.empty(n_boot)
    for i in range(n_boot):
        ai = rng.choice(a, len(a), replace=True)
        bi = rng.choice(b, len(b), replace=True)
        boots[i] = np.median((bi[:, None] - ai[None, :]).ravel())
    return hl, float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def strata(m, ax, exclude_mixed=True, any_nonbacterial=False):
    """any_nonbacterial=True 复刻 13 号的口径：凡不是"两性状都是细菌轴"即计入 cross
    （即含 8 对 bacterial-mixed + 1 对 mixed-viral，共 17 对）——仅作口径溯源用。"""
    internal, cross = [], []
    for _, r in m[m["power_ok"] == True].iterrows():  # noqa: E712
        if exclude_mixed and MIXED in (r["axis1"], r["axis2"]):
            continue
        a1 = ax.get(r["axis1"], "").lower(); a2 = ax.get(r["axis2"], "").lower()
        if a1 == "bacterial" and a2 == "bacterial":
            internal.append(float(r["rg"]))
        elif any_nonbacterial:
            cross.append(float(r["rg"]))
        elif "viral" in (a1, a2):
            cross.append(float(r["rg"]))
    return internal, cross


def axis_test(df, label, rows):
    """对给定矩阵子集做轴内 vs 跨轴检验并写入 rows（与 50 号同口径：剔除 mixed 对）。"""
    t1 = pd.read_csv(OUT / "table1.csv")
    ax = t1.set_index("Trait")["Pathogen_axis"].to_dict()
    a, b = strata(df, ax)
    if len(a) < 3 or len(b) < 3:
        rows.append({"key": "%s_n" % label, "value": "insufficient", "note": label})
        return
    p = mannwhitneyu(a, b, alternative="two-sided").pvalue
    rows += [
        {"key": "%s_n_internal" % label, "value": str(len(a)), "note": "轴内对数"},
        {"key": "%s_n_cross" % label, "value": str(len(b)), "note": "跨轴对数"},
        {"key": "%s_median_internal" % label, "value": "%.3f" % np.median(a), "note": "轴内中位"},
        {"key": "%s_median_cross" % label, "value": "%.3f" % np.median(b), "note": "跨轴中位"},
        {"key": "%s_mw_p" % label, "value": "%.3f" % p, "note": "Mann-Whitney 双侧 p"},
    ]


def main():
    rows = []
    m = pd.read_csv(MATRIX)
    ax = pd.read_csv(OUT / "table1.csv").set_index("Trait")["Pathogen_axis"].to_dict()
    rng = np.random.default_rng(20261028)

    # ---------- M-B / M-C：单一主口径 + 效应量 ----------
    cases = [
        ("primary", dict(exclude_mixed=True)),
        ("mixed_inclusive", dict(exclude_mixed=False)),
        ("stratum13_def", dict(exclude_mixed=False, any_nonbacterial=True)),
        ("excl_int_gt_03", dict(exclude_mixed=True, int_max=0.3)),
        ("excl_int_gt_02", dict(exclude_mixed=True, int_max=0.2)),
        ("excl_rg_gt_09", dict(exclude_mixed=True, rg_max=0.9)),
    ]
    for name, kw in cases:
        mm = m.copy()
        if "int_max" in kw:
            mm.loc[mm["int_gencov"].notna() & (mm["int_gencov"] > kw["int_max"]), "power_ok"] = False
        if "rg_max" in kw:
            mm.loc[mm["rg"].notna() & (mm["rg"] > kw["rg_max"]), "power_ok"] = False
        a, b = strata(mm, ax, kw["exclude_mixed"], kw.get("any_nonbacterial", False))
        if len(a) < 3 or len(b) < 3:
            continue
        p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        hl, lo, hi = hl_ci(a, b, rng)
        rows += [
            {"key": "%s_n_internal" % name, "value": str(len(a)), "note": "轴内对数"},
            {"key": "%s_n_cross" % name, "value": str(len(b)), "note": "跨轴对数"},
            {"key": "%s_median_internal" % name, "value": "%.3f" % np.median(a), "note": "轴内中位"},
            {"key": "%s_median_cross" % name, "value": "%.3f" % np.median(b), "note": "跨轴中位"},
            {"key": "%s_mw_p" % name, "value": "%.3f" % p, "note": "Mann-Whitney 双侧 p"},
            {"key": "%s_hl" % name, "value": "%.3f" % hl, "note": "Hodges-Lehmann 位移（跨轴−轴内）"},
            {"key": "%s_hl_lo" % name, "value": "%.3f" % lo, "note": "HL 的 bootstrap 2.5% 分位"},
            {"key": "%s_hl_hi" % name, "value": "%.3f" % hi, "note": "HL 的 bootstrap 97.5% 分位"},
        ]
    # 可检测差异（粗略）：基于 bootstrap 位移分布半宽
    a, b = strata(m, ax)
    _, lo, hi = hl_ci(a, b, rng)
    rows.append({"key": "primary_detectable_shift", "value": "%.3f" % ((hi - lo) / 2),
                 "note": "主分析下 95% CI 半宽（可作为该功效水平下可分辨的最小位移）"})

    # ---------- M-C：加权回归的编码方向 ----------
    mm = m[m["power_ok"] == True].copy()
    mm = mm[~mm["axis1"].eq(MIXED) & ~mm["axis2"].eq(MIXED)]
    mm["cross"] = [1.0 if "viral" in (ax.get(x, "").lower(), ax.get(y, "").lower()) else 0.0
                   for x, y in zip(mm["axis1"], mm["axis2"])]
    mm["chi_mean"] = (mm["mean_chisq_1"] + mm["mean_chisq_2"]) / 2
    X = np.column_stack([np.ones(len(mm)), mm["cross"].values, mm["chi_mean"].values])
    w = 1.0 / np.square(mm["rg_se"].values)
    XtW = X.T * w
    beta = np.linalg.lstsq(XtW @ X, XtW @ mm["rg"].values, rcond=None)[0]
    resid = mm["rg"].values - X @ beta
    sigma2 = float((w * resid ** 2).sum() / (len(mm) - 3))
    cov = sigma2 * np.linalg.inv(XtW @ X)
    se = np.sqrt(np.diag(cov))
    rows += [
        {"key": "wls_beta_cross", "value": "%.3f" % beta[1],
         "note": "编码：cross 虚拟变量=1 表示跨轴；系数为调整 mean χ² 后的跨轴−轴内差"},
        {"key": "wls_se_cross", "value": "%.3f" % se[1], "note": "同上 SE"},
        {"key": "wls_beta_chi", "value": "%.3f" % beta[2], "note": "mean χ² 的系数"},
        {"key": "wls_se_chi", "value": "%.3f" % se[2], "note": "同上 SE"},
        {"key": "wls_corr_chi_rg", "value": "%.3f" % np.corrcoef(mm["chi_mean"], mm["rg"])[0, 1],
         "note": "mean χ² 与 rg 的 Pearson 相关（接近 0：χ² 不是符号变化的成因）"},
        {"key": "wls_se_corr_rg", "value": "%.3f" % np.corrcoef(mm["rg"], mm["rg_se"])[0, 1],
         "note": "rg 与其标准误的相关；显著为正即逆方差权重是信息性的（高 rg 对权重低）——"
                 "这是调整后系数与未加权中位方向相反的原因，故正文不采用该加权模型"},
        {"key": "wls_chi_internal", "value": "%.3f" % mm[mm.cross == 0]["chi_mean"].mean(),
         "note": "轴内对的 mean χ² 均值"},
        {"key": "wls_chi_cross", "value": "%.3f" % mm[mm.cross == 1]["chi_mean"].mean(),
         "note": "跨轴对的 mean χ² 均值"},
    ]
    # χ² 分层内的轴内/跨轴差（不依赖线性假设）
    q = pd.qcut(mm["chi_mean"], 3, labels=["low", "mid", "high"])
    for lab in ("low", "mid", "high"):
        sub = mm[q == lab]
        ai = sub[sub.cross == 0]["rg"]; bi = sub[sub.cross == 1]["rg"]
        rows.append({"key": "chi_stratum_%s" % lab,
                     "value": "%.3f vs %.3f (n=%d/%d)" % (ai.median() if len(ai) else np.nan,
                                                         bi.median() if len(bi) else np.nan,
                                                         len(ai), len(bi)),
                     "note": "χ² 三分位内：轴内 vs 跨轴 rg 中位"})

    # ---------- M-D：头条 coloc 扣除 MHC ----------
    C = pd.read_csv(OUT / "coloc_results.csv")
    L = pd.read_csv(OUT / "locus_atlas.csv")
    REG = {r["SNP"]: r["region_id"] for _, r in L.iterrows()}
    MHC_REG = {r["region_id"] for _, r in L.iterrows() if r["in_mhc"] == True}  # noqa: E712
    best = {}
    for _, r in C.iterrows():
        k = (REG.get(r["locus"], "UNMAPPED"), tuple(sorted([r["trait1"], r["trait2"]])))
        best[k] = max(best.get(k, 0.0), float(r["PP_H4"]))
    nonmhc = {k: v for k, v in best.items() if k[0] not in MHC_REG}
    rows += [
        {"key": "coloc_region_nonmhc_pairs", "value": str(len(nonmhc)),
         "note": "区域级 coloc 对，扣除 MHC 区（原 38）"},
        {"key": "coloc_region_nonmhc_share", "value": str(sum(1 for v in nonmhc.values() if v > 0.8)),
         "note": "区域级 PP_H4>0.8，扣除 MHC（原 12）"},
        {"key": "coloc_mhc_region_included", "value": str(len(best) - len(nonmhc)),
         "note": "被扣除的 MHC 区域级对数（1 对，正文已声明不采信）"},
    ]

    # ---------- M-H：p–Z 守卫的三个分母 ----------
    MC = list(csv.DictReader(open(OUT.parent / "data" / "munged" / "munge_counts.csv",
                                 encoding="utf-8-sig"), delimiter="\t"))
    panel = set()
    from constants import FINNGEN_ENDPOINTS  # noqa: E402
    panel = set(FINNGEN_ENDPOINTS.values())
    tot_all = sum(int(r["total_variants"] or 0) for r in MC)
    tot_panel_merged = sum(int(r["merged_1000g"] or 0) for r in MC if r["phenocode"] in panel)
    tot_all_merged = sum(int(r["merged_1000g"] or 0) for r in MC)
    pz_panel = sum(int(r["pz_mismatch_dropped"] or 0) for r in MC if r["phenocode"] in panel)
    pz_all = sum(int(r["pz_mismatch_dropped"] or 0) for r in MC)
    rows += [
        {"key": "pz_total_before_restriction", "value": str(tot_all),
         "note": "分母①：panel 限制前的输入变异数（全队列求和）"},
        {"key": "pz_total_after_restriction", "value": str(tot_all_merged),
         "note": "分母②：与 1000G 参考面板合并后（守卫实际作用层，全部性状）"},
        {"key": "pz_total_panel_after_restriction", "value": str(tot_panel_merged),
         "note": "分母③：合并后且仅感染面板 25 终点（守卫实际作用层，面板内）"},
        {"key": "pz_dropped_panel", "value": str(pz_panel), "note": "面板内被守卫剔除数（实测 0）"},
        {"key": "pz_dropped_all", "value": str(pz_all), "note": "全部性状被守卫剔除数（实测 57）"},
        {"key": "pz_frac_panel", "value": "%.2e" % (pz_panel / tot_panel_merged),
         "note": "面板内比例（守卫在面板内无作用）"},
        {"key": "pz_frac_all", "value": "%.2e" % (pz_all / tot_all_merged),
         "note": "全部性状比例（≈1.2e-6）"},
    ]

    # ---------- 次5：回文 SNP 比例（面板内一个终点的回归集代理） ----------
    ss = OUT.parent / "data" / "munged" / ("%s.sumstats" % FINNGEN_ENDPOINTS["BACT_NOS"])
    if ss.exists():
        d = pd.read_csv(ss, sep="\t", usecols=["A1", "A2"])
        pal = d.apply(lambda r: {r["A1"], r["A2"]} in ({"A", "T"}, {"C", "G"}), axis=1)
        rows.append({"key": "palindromic_frac", "value": "%.4f" % pal.mean(),
                     "note": "回文（A/T、C/G）SNP 占最终集比例；这些变异已按参考等位编码定向，"
                             "不存在按 MAF 猜链的环节（munge 吸收到 1000G 参考等位）"})

    # ---------- 次7：gene-set OR 的 95% CI（由既有列重建 2x2） ----------
    gs = pd.read_csv(OUT / "gene_set_enrichment.csv")
    n_lead_total = int(pd.read_csv(OUT / "iei_summary.csv").set_index("key")["value"]["lead_genes_n"])
    hg = pd.read_csv(OUT.parent / "data" / "hgnc_complete_set.txt", sep="	", low_memory=False, dtype=str)
    BG = int((hg["locus_group"] == "protein-coding gene").sum())
    hla = gs[gs["gene_set"].astype(str).str.upper().str.startswith("HLA")]
    if len(hla) and {"n_members", "n_lead_in_set"} <= set(gs.columns):
        r0 = hla.iloc[0]
        a = float(r0["n_lead_in_set"]); n_mem = float(r0["n_members"])
        b = n_mem - a; c = n_lead_total - a; d = BG - n_mem - c
        lor = np.log((a * d) / (b * c)); se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        lo, hi = np.exp(lor - 1.96 * se), np.exp(lor + 1.96 * se)
        # 不写回 gene_set_enrichment.csv：该文件由 19 号重生成，若在此加列会造成
        # 端到端门的合法漂移（DRIFT gene_set_enrichment.csv，2026-10-28 实测）。
        # CI 只进 canonical（rr2_hla_or_ci）与 S18。
        rows.append({"key": "hla_or_ci", "value": "%.0f (%.0f-%.0f)" % (r0["odds_ratio"], lo, hi),
                     "note": "HLA 富集 OR 与 Woolf 95%% CI（a=%d/%d lead，背景 %d 蛋白编码基因；CI 宽属小计数常态）"
                             % (int(a), n_lead_total, BG)})
        rows.append({"key": "hla_or_ci_bg", "value": str(BG), "note": "背景基因数（HGNC protein-coding）"})
    else:
        rows.append({"key": "hla_or_ci", "value": "NA", "note": "gene_set_enrichment 缺 n_members/n_lead_in_set"})

    # ---------- M-G：coloc.susie 全局口径汇总（修正后结果；对角口径单列留档） ----------
    SP = OUT / "coloc_susie_results.csv"
    if SP.exists():
        su = pd.read_csv(SP)
        if "PP_H4" in su.columns:
            g = su["PP_H4"].dropna()
            diag = su["PP_H4_diagonal_only"].dropna() if "PP_H4_diagonal_only" in su.columns else pd.Series(dtype=float)
            rows += [
                {"key": "susie_n_pairs", "value": str(len(su)),
                 "note": "区域级去重后的多信号 coloc 对数"},
                {"key": "susie_global_median", "value": "%.3f" % g.median(),
                 "note": "全局归一化 PP.H4 中位（正确口径）"},
                {"key": "susie_global_max", "value": "%.3f" % g.max(), "note": "全局口径最大值"},
                {"key": "susie_global_min", "value": "%.3f" % g.min(), "note": "全局口径最小值"},
                {"key": "susie_regions", "value": str(su["region"].nunique()) if "region" in su.columns else "NA",
                 "note": "涉及区域数"},
                {"key": "susie_global_n_gt08", "value": str(int((g > 0.8).sum())),
                 "note": "全局口径下 PP_H4>0.8 的对数（原稿误报 14/14）"},
                {"key": "susie_diag_median", "value": ("%.3f" % diag.median()) if len(diag) else "NA",
                 "note": "对角行口径中位（留档：行内归一化，不可与 abf 比较）"},
                {"key": "susie_max_combos", "value": str(int(su["n_signal_combinations"].max()))
                 if "n_signal_combinations" in su.columns else "NA",
                 "note": "单对最大信号对数（H4 被稀释的机制）"},
            ]

    # ---------- 次4：liability 尺度 h² 的三档患病率敏感性 ----------
    t1 = pd.read_csv(OUT / "table1.csv")
    if {"Cases", "Controls", "h2_observed"} <= set(t1.columns):
        t1 = t1[(t1["Cases"].astype(str) != "") & (t1["Controls"].astype(str) != "")]
        cc = pd.to_numeric(t1["Cases"], errors="coerce"); ck = pd.to_numeric(t1["Controls"], errors="coerce")
        P = cc / (cc + ck); h2o = pd.to_numeric(t1["h2_observed"], errors="coerce")
        out = []
        for K in (0.02, 0.05, 0.10):
            z = norm.pdf(norm.ppf(1 - K))
            fac = (K ** 2) * ((1 - K) ** 2) / (P * (1 - P) * (z ** 2))
            out.append((K, float((h2o * fac).median())))
        rows += [{"key": "h2_liab_K%02d" % int(K * 100), "value": "%.4f" % v,
                  "note": "liability 尺度 h² 中位（假定人群患病率 K=%.2f，样本患病率作 P）" % K}
                 for K, v in out]
        rows.append({"key": "h2_liab_note", "value": "assumption-based",
                     "note": "三档 K 仅示数量级；K 本身未知，故正文不采用单一转换值"})

    # ---------- 评审 1（统计）：按性状配对分析（8 个细菌经 INFLUENZA 的同一性状比较） ----------
    from scipy.stats import wilcoxon
    mm2 = m[m["power_ok"] == True].copy()
    bact = [t for t in pd.read_csv(OUT / "table1.csv").query("Pathogen_axis=='Bacterial'")["Trait"]]
    cross_pairs = [(r["axis1"], r["axis2"], float(r["rg"])) for _, r in mm2.iterrows()
                   if "INFLUENZA" in (r["axis1"], r["axis2"])
                   and MIXED not in (r["axis1"], r["axis2"])]   # 与主口径一致：剔除 mixed 对
    within = {}
    for _, r in mm2.iterrows():
        a, b, v = r["axis1"], r["axis2"], float(r["rg"])
        if a in bact and b in bact:
            within.setdefault(a, []).append(v); within.setdefault(b, []).append(v)
    diffs = []
    for a, b, v in cross_pairs:
        o = b if a == "INFLUENZA" else a
        if o in within and within[o]:
            diffs.append(float(np.mean(within[o])) - v)
    rows.append({"key": "paired_n", "value": str(len(diffs)),
                 "note": "配对单位 = 细菌性状（其与其余细菌的均值 rg 减去与 INFLUENZA 的 rg）"})
    if len(diffs) >= 5:
        pv = float(wilcoxon(diffs).pvalue)
        rows += [
            {"key": "paired_median_diff", "value": "%.3f" % np.median(diffs),
             "note": "配对中位差（轴内均值 − 跨轴）；正=轴内更高"},
            {"key": "paired_wilcoxon_p", "value": "%.3f" % pv, "note": "Wilcoxon 符号秩双侧 p"},
            {"key": "paired_n_positive", "value": str(int(sum(d > 0 for d in diffs))),
             "note": "配对差为正的性状数"},
        ]
    # 跨轴对是否全部为细菌×INFLUENZA
    rows.append({"key": "cross_all_influenza",
                 "value": "yes" if len(cross_pairs) == len(diffs) and len(cross_pairs) > 0 else "no",
                 "note": "过门跨轴对是否全部以 INFLUENZA 为一侧"})
    # 留一性状的方向稳定性（找出会翻转方向的性状）
    flips = []
    for t in sorted(set(mm2["axis1"]) | set(mm2["axis2"])):
        a2, b2 = strata(mm2[~mm2["axis1"].eq(t) & ~mm2["axis2"].eq(t)], ax)
        if len(a2) >= 3 and len(b2) >= 3 and np.median(a2) > np.median(b2):
            flips.append(t)
    rows.append({"key": "loo_direction_flip", "value": ", ".join(flips) or "none",
                 "note": "留一后轴内中位反超跨轴的性状（方向不稳定的来源）"})

    # ---------- 评审 1/6：更严格功效门的敏感性 ----------
    for tag, cond in (("h2p01", m["h2_1_p"].lt(0.01) & m["h2_2_p"].lt(0.01)),
                      ("h2p001", m["h2_1_p"].lt(0.001) & m["h2_2_p"].lt(0.001)),
                      ("chi105", m["mean_chisq_1"].ge(1.05) & m["mean_chisq_2"].ge(1.05))):
        sub = m.copy()
        sub.loc[~cond | (sub["status"] != "ok"), "power_ok"] = False
        axis_test(sub, "axis_strict_%s" % tag, rows)

    # ---------- 第三轮补：逐病 BH 计数（正文按计数排序，替代人为分档） ----------
    ad = pd.read_csv(OUT / ("ldsc_rg_autoimmune_%s.csv" % DATE))
    adp = ad[ad["power_ok"] == True]
    for code in ("L12_PSORIASIS", "M13_RHEUMA", "M13_SJOGREN", "M13_ANKYLOSPON",
                 "K11_CD_STRICT2", "K11_UC_STRICT2", "E4_GRAVES_STRICT", "T1D",
                 "G6_MS", "K11_COELIAC"):
        sub = adp[adp["axis1"] == code]
        rows.append({"key": "bhcount_%s" % code.lower(),
                     "value": str(int((sub["rg_q_bh_layer"] < 0.05).sum())),
                     "note": "%s 过门对中 BH<0.05 的个数" % code})
    cd = pd.read_csv(OUT / ("ldsc_rg_comparator_%s.csv" % DATE))
    cdp = cd[cd["power_ok"] == True]
    for code in ("I9_CHD", "T2D"):
        sub = cdp[cdp["axis1"] == code]
        rows.append({"key": "bhcount_%s" % code.lower(),
                     "value": str(int((sub["rg_q_bh_layer"] < 0.05).sum())),
                     "note": "%s 过门对中 BH<0.05 的个数（对照层）" % code})

    # ---------- 第三轮评审 M2'：显式功效模拟（替代"CI 半宽=可检出下限"的混淆表述） ----------
    # 设计：以观测分布为总体，跨轴组整体加 δ，重抽样 28 vs 8，Mann-Whitney α=0.05，
    # 1000 次/δ；求 80% 功效对应的 δ。种子固定，结果可复现。
    rng_p = np.random.default_rng(20261028)
    a0, b0 = strata(m, ax)
    grid = np.round(np.arange(0.0, 0.41, 0.05), 2)
    powers = []
    for delta in grid:
        hits = 0
        for _ in range(1000):
            ai = rng_p.choice(a0, len(a0), replace=True)
            bi = rng_p.choice(b0, len(b0), replace=True) + delta
            if mannwhitneyu(ai, bi, alternative="two-sided").pvalue < 0.05:
                hits += 1
        powers.append(hits / 1000)
    rows.append({"key": "power_grid", "value": "; ".join("%.2f:%.2f" % (d, p) for d, p in zip(grid, powers)),
                 "note": "功效模拟：δ 为跨轴组整体位移；每点 1000 次重抽样"})
    p80 = next((d for d, p in zip(grid, powers) if p >= 0.80), None)
    rows.append({"key": "power_shift_80", "value": ("%.2f" % p80) if p80 is not None else ">0.40",
                 "note": "80% 功效对应的跨轴组中位位移（模拟；非 CI 半宽）"})

    # ---------- 评审 2/5：外部层在三个区域的区域信号强度（决定"未复制"是否可解释） ----------
    CD = OUT.parent / "data" / "coloc"
    for lab, loc in (("fut2", "rs492602"), ("abo", "rs495828"), ("crp", "rs3116651")):
        for tr in ("HGI_A2", "HGI_C2"):
            f = CD / ("%s_%s.tsv" % (loc, tr))
            if not f.exists():
                continue
            d = pd.read_csv(f, sep="	")
            if "pval" not in d.columns:
                continue
            rows.append({"key": "%s_%s_minp" % (lab, tr.lower()),
                         "value": "%.1e" % float(d["pval"].min()),
                         "note": "%s 区 %s 在 coloc 窗口内的最小 P（外部层是否有区域信号）" % (lab.upper(), tr)})

    pd.DataFrame(rows).to_csv(SUM, index=False)
    print("[2nd] %d 项 -> %s" % (len(rows), SUM.name))
    for r in rows:
        print("   %-32s %s" % (r["key"], r["value"]))


if __name__ == "__main__":
    main()
