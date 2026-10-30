# -*- coding: utf-8 -*-
""" 52_sensitivity_extended.py — 补充敏感性分析总成（2026-10-28）

对应补充分析条目（编号沿用原编号）：
  M3 样本重叠/共享对照 → 截距校准（含内部对照）+ 剔除高截距对后重估轴检验
  M4 PNEUMONIA_MIXED 定义性重叠 → 剔除 mixed 轴后的位点分类与区域级 coloc
  M8 层间 FDR → 汇总各 FDR 域 + 合并全域 BH 敏感性
  M9 p 值一致性阈值 → 报告实际剔除量占比（阈值依据=数据）
  次6 Table 3 折叠不对称 → 核对跨源对是否本来就每区一对
  次4 rg>1 对清单（表型重叠）
  次7 coloc 先验敏感性 → 跨 p12 三档的支持状态翻转计数
  次9 置换对照 → 带固定种子的 20 次全行置换（补录此前未留种子的缺陷）

只读已有冻结产物与 munge 后的 sumstats；不重跑主矩阵。
输出：work/output/sensitivity_extended_summary.csv + 若干明细 CSV
用法：python 52_sensitivity_extended.py [--skip-perm]
"""
import csv
import importlib.util
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, mannwhitneyu

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATE, OUT  # noqa: E402

MATRIX = OUT / "ldsc_rg" / ("ldsc_matrix_%s.csv" % DATE)
ATLAS = OUT / "locus_atlas.csv"
COLOC = OUT / "coloc_results.csv"
PRIOR = OUT / "coloc_p12_sensitivity.csv"
AUTO = OUT / ("ldsc_rg_autoimmune_%s.csv" % DATE)
COMP = OUT / ("ldsc_rg_comparator_%s.csv" % DATE)
MUNGE = OUT.parent / "data" / "munged" / "munge_counts.csv"
SUM = OUT / "sensitivity_extended_summary.csv"
PERM = OUT / "permutation_null_20.csv"


def kv_summary(rows):
    pd.DataFrame(rows).to_csv(SUM, index=False)


def axis_test(m, label, out_rows):
    """轴内 vs 跨轴 Mann-Whitney（与 50 号同口径：剔除细菌-混合对）。"""
    t1 = pd.read_csv(OUT / "table1.csv")
    AX = t1.set_index("Trait")["Pathogen_axis"].to_dict()
    rec = []
    for _, r in m[m["power_ok"] == True].iterrows():  # noqa: E712
        if "PNEUMONIA_MIXED" in (r["axis1"], r["axis2"]):
            continue
        a1 = AX.get(r["axis1"], "").lower(); a2 = AX.get(r["axis2"], "").lower()
        if a1 == "bacterial" and a2 == "bacterial":
            rec.append(("internal", float(r["rg"])))
        elif "viral" in (a1, a2):
            rec.append(("cross", float(r["rg"])))
    a = [v for k, v in rec if k == "internal"]; b = [v for k, v in rec if k == "cross"]
    if len(a) < 3 or len(b) < 3:
        out_rows.append({"key": "%s_n" % label, "value": "insufficient", "note": label})
        return
    p = mannwhitneyu(a, b, alternative="two-sided").pvalue
    out_rows += [
        {"key": "%s_n_internal" % label, "value": str(len(a)), "note": label},
        {"key": "%s_n_cross" % label, "value": str(len(b)), "note": label},
        {"key": "%s_median_internal" % label, "value": "%.3f" % np.median(a), "note": label},
        {"key": "%s_median_cross" % label, "value": "%.3f" % np.median(b), "note": label},
        {"key": "%s_mw_p" % label, "value": "%.3f" % p, "note": label},
    ]


def main():
    rows = []
    m = pd.read_csv(MATRIX)
    ok = m[m["status"] == "ok"].copy()
    bh = ok[ok["rg_q_bh"].notna() & (ok["rg_q_bh"] < 0.05)]

    # ---------- M3：截距校准（内部对照） ----------
    ctl = {}
    for f, lab in (("ldsc_rg/T2D_I9_CHD.csv", "pos"),
                   ("ldsc_rg/T2D_I9_CHD_shuffled.csv", "perm"),
                   ("ldsc_rg/T2D_I9_CHD_signflip_null.csv", "signflip")):
        p = OUT / f
        if p.exists():
            d = pd.read_csv(p)
            if "gencov_intercept" in d.columns:
                ctl[lab] = float(d["gencov_intercept"].iloc[0])
    rows += [
        {"key": "ctl_gencov_int_positive", "value": ("%.4f" % ctl["pos"]) if "pos" in ctl else "",
         "note": "内部正对照 T2D×CHD 的 cross-trait 截距（已知共享对照池）"},
        {"key": "ctl_gencov_int_permuted", "value": ("%.4f" % ctl["perm"]) if "perm" in ctl else "",
         "note": "全行置换零对照截距（应≈0）"},
        {"key": "ctl_gencov_int_signflip", "value": ("%.4f" % ctl["signflip"]) if "signflip" in ctl else "",
         "note": "符号翻转零分布截距（应≈0）"},
        {"key": "bh_gencov_int_median_ctx", "value": "%.4f" % bh["int_gencov"].median(),
         "note": "显著感染对 gencov 截距中位（对比上两项：低于正对照即重叠较轻）"},
        {"key": "bh_gencov_int_median_3dp", "value": "%.3f" % bh["int_gencov"].median(),
         "note": "同上，正文引用用的三位有效数字（次1：避免 16 位浮点）"},
    ]
    # 剔除高截距对后重估（阈值 0.3 / 0.2 两档；0.3 约为正对照 0.211 的 1.4 倍）
    for thr in (0.3, 0.2):
        sub = m.copy()
        sub.loc[sub["int_gencov"].notna() & (sub["int_gencov"] > thr), "power_ok"] = False
        axis_test(sub, "axis_excl_int_gt_%s" % str(thr).replace(".", ""), rows)
    # 定义性重叠对（rg≈1）剔除后重估
    sub = m.copy()
    sub.loc[sub["rg"].notna() & (sub["rg"] > 0.9), "power_ok"] = False
    axis_test(sub, "axis_excl_rg_gt_09", rows)
    axis_test(m, "axis_all", rows)

    # ---------- M4：剔除 mixed 轴后的位点分类 ----------
    L = pd.read_csv(ATLAS)
    recs = []
    for _, r in L.iterrows():
        axes = [a for a in str(r["axes"]).split(",") if a]
        excl = [a for a in axes if a != "mixed"]
        recs.append({"SNP": r["SNP"], "region_id": r["region_id"], "in_mhc": r["in_mhc"],
                     "axes_all": len(axes), "axes_nomixed": len(excl),
                     "class_nomixed": "pathogen_specific" if len(excl) <= 1 else "cross_axis"})
    R = pd.DataFrame(recs)
    R.to_csv(OUT / "locus_class_no_mixed.csv", index=False)
    rows += [
        {"key": "nomixed_n_specific", "value": str(int((R["class_nomixed"] == "pathogen_specific").sum())),
         "note": "剔除 mixed 轴后：单轴位点数（原 18 特异 + 依赖 mixed 者）"},
        {"key": "nomixed_n_cross", "value": str(int((R["class_nomixed"] == "cross_axis").sum())),
         "note": "剔除 mixed 轴后：跨轴位点数（原 9）"},
        {"key": "nomixed_n_cross_regions", "value": str(R[R["class_nomixed"] == "cross_axis"]["region_id"].nunique()),
         "note": "剔除 mixed 轴后：跨轴区域数（原 4）"},
        {"key": "nomixed_n_cross_regions_nonmhc",
         "value": str(R[(R["class_nomixed"] == "cross_axis") & (R["in_mhc"] != True)]["region_id"].nunique()),  # noqa: E712
         "note": "剔除 mixed 轴后：非 MHC 跨轴区域数（原 3）"},
        {"key": "mixed_dependent_loci",
         "value": str(int(((R["axes_all"] - R["axes_nomixed"]) > 0).sum())),
         "note": "轴覆盖依赖 mixed 的 lead 数（分类会因此改变）"},
    ]
    # 区域级 coloc：剔除 mixed 参与的对
    C = pd.read_csv(COLOC)
    REG = {r["SNP"]: r["region_id"] for _, r in L.iterrows()}
    def regbest(df):
        best = {}
        for _, r in df.iterrows():
            k = (REG.get(r["locus"], "UNMAPPED"), tuple(sorted([r["trait1"], r["trait2"]])))
            best[k] = max(best.get(k, 0.0), float(r["PP_H4"]))
        return best
    allb = regbest(C)
    nomix = C[~C["trait1"].eq("PNEUMONIA_MIXED") & ~C["trait2"].eq("PNEUMONIA_MIXED")]
    nb = regbest(nomix)
    rows += [
        {"key": "coloc_region_pairs_all", "value": str(len(allb)), "note": "区域级 coloc 对数（原 38）"},
        {"key": "coloc_region_share_all", "value": str(sum(1 for v in allb.values() if v > 0.8)),
         "note": "区域级 PP_H4>0.8（原 12）"},
        {"key": "coloc_region_pairs_nomixed", "value": str(len(nb)), "note": "剔除 mixed 参与对后"},
        {"key": "coloc_region_share_nomixed", "value": str(sum(1 for v in nb.values() if v > 0.8)),
         "note": "剔除 mixed 参与对后 PP_H4>0.8"},
    ]
    # 次6：跨源对的区域冗余（应≈0，因每区每对只取一行）
    def is_cross(r):
        return str(r["trait1"]).startswith("HGI") or str(r["trait2"]).startswith("HGI")
    xs = C[C.apply(is_cross, axis=1)]
    red = len(xs) - len(regbest(xs))
    rows.append({"key": "coloc_cross_redundancy", "value": str(red),
                 "note": "跨源对的区域去重删行数（0 = 跨源本来每区每对一行，故 Table 3 折叠不对称）"})

    # ---------- M8：FDR 域汇总 + 合并全域 BH ----------
    a = pd.read_csv(AUTO); c = pd.read_csv(COMP)
    # 2026-10-28 修正（五评委核查）：三层推断域统一为**过功效门**的对，故合并域亦用
    # power_ok 子集（45 + 100 + 20 = 165），与 12/43/49 号的 BH 域一致。
    dom = [("main", m[m["power_ok"] == True]),  # noqa: E712
           ("autoimmune", a[a["power_ok"] == True]),  # noqa: E712
           ("comparator", c[c["power_ok"] == True])]  # noqa: E712
    pooled = []
    for name, d in dom:
        pooled.append(d.assign(domain=name)[["domain", "rg_p", "rg", "axis1", "axis2"]])
        rows.append({"key": "fdr_%s_tests" % name, "value": str(len(d)), "note": "该层过功效门检验数（推断域）"})
    P = pd.concat(pooled, ignore_index=True)
    q = false_discovery_control(P["rg_p"].values, method="bh")
    P["q_global"] = q
    P.to_csv(OUT / "fdr_pooled_all_domains.csv", index=False)
    rows += [
        {"key": "fdr_pooled_tests", "value": str(len(P)), "note": "全域合并检验数（三层过门对之和）"},
        {"key": "fdr_pooled_survivors", "value": str(int((q < 0.05).sum())), "note": "合并全域 BH<0.05"},
        {"key": "fdr_pooled_main_survivors",
         "value": str(int(((P["domain"] == "main") & (q < 0.05)).sum())),
         "note": "合并全域后仍显著的主矩阵对数（对照层内 35）"},
        {"key": "fdr_pooled_auto_survivors",
         "value": str(int(((P["domain"] == "autoimmune") & (q < 0.05)).sum())),
         "note": "合并全域后仍显著的自免锚层对数（对照层内 29）"},
    ]

    # ---------- M9：p 值一致性阈值实际剔除量 ----------
    MC = list(csv.DictReader(open(MUNGE, encoding="utf-8-sig"), delimiter="\t"))
    tot = sum(int(r["total_variants"] or 0) for r in MC)
    pz = sum(int(r["pz_mismatch_dropped"] or 0) for r in MC)
    maf = sum(int(r["maf_out_dropped"] or 0) for r in MC)
    rows += [
        {"key": "qc_total_variants", "value": str(tot), "note": "全部终点输入变异总数（munge 台账求和）"},
        {"key": "qc_pz_dropped", "value": str(pz), "note": "p–Z 一致性 30% 守卫剔除总数"},
        {"key": "qc_pz_dropped_frac", "value": "%.2e" % (pz / tot), "note": "占输入比例（阈值实际影响面）"},
        {"key": "qc_maf_dropped", "value": str(maf), "note": "MAF 越界剔除总数"},
    ]

    # ---------- 次4：rg>1 对 ----------
    over = ok[ok["rg"] > 1.0][["axis1", "axis2", "rg", "rg_q_bh"]]
    over.to_csv(OUT / "rg_over_one.csv", index=False)
    rows.append({"key": "n_rg_gt_1", "value": str(len(over)),
                 "note": "rg>1 的可解释对（绝大多数 q=NA，属表型定义重叠；清单见 rg_over_one.csv）"})

    # ---------- 次7：coloc 先验敏感性 ----------
    # 两级口径并列：**(位点,对)** 是分析单位（同一区域多 lead 会重复测试），
    # **(区域,对)** 才是生物学单位（规则 222 同族）。同时给出主结论区域的三档跨源最大值。
    PR = pd.read_csv(PRIOR)
    loc2ph = defaultdict(dict)
    for (lo, t1, t2), sub in PR.groupby(["locus", "trait1", "trait2"]):
        for _, r in sub.iterrows():
            loc2ph[(lo, t1, t2)][float(r["p12"])] = float(r["PP_H4"])
    tiers = sorted({k for d in loc2ph.values() for k in d})
    flip_loc = sum(1 for d in loc2ph.values() if len({v > 0.8 for v in d.values()}) > 1)
    reg2ph = defaultdict(lambda: defaultdict(float))
    for (lo, t1, t2), d in loc2ph.items():
        rid = REG.get(lo, "UNMAPPED")
        for p12, v in d.items():
            reg2ph[(rid, t1, t2)][p12] = max(reg2ph[(rid, t1, t2)][p12], v)
    flip_reg = sum(1 for d in reg2ph.values() if len({v > 0.8 for v in d.values()}) > 1)
    rows += [
        {"key": "prior_tiers", "value": ",".join("%g" % t for t in tiers),
         "note": "coloc.abf 先验 p12 三档（主结果用 coloc 默认 1e-5）"},
        {"key": "prior_n_sets_locus", "value": str(len(loc2ph)), "note": "位点级覆盖组合数"},
        {"key": "prior_flip_locus", "value": str(flip_loc),
         "note": "位点级 PP_H4>0.8 判定翻转数（同区分 lead 重复测试会放大此数）"},
        {"key": "prior_n_sets_region", "value": str(len(reg2ph)), "note": "区域级覆盖组合数（生物学单位）"},
        {"key": "prior_flip_region", "value": str(flip_reg), "note": "区域级支持状态翻转数"},
    ]
    # 主结论区域：跨源最大值在三档下的稳健性（FUT2 / ABO / CRP）
    for lab, rid in (("fut2", "R21_chr19_49Mb"), ("abo", "R14_chr9_136Mb"),
                     ("crp", "R02_chr1_159Mb")):
        vals = {}
        for (rid2, t1, t2), d in reg2ph.items():
            if rid2 != rid:
                continue
            if str(t1).startswith("HGI") or str(t2).startswith("HGI"):
                for p12, v in d.items():
                    vals[p12] = max(vals.get(p12, 0.0), v)
        for p12 in tiers:
            if p12 in vals:
                rows.append({"key": "%s_cross_max_p12_%g" % (lab, p12),
                             "value": "%.3f" % vals[p12],
                             "note": "%s 区跨源最大 PP_H4（p12=%g）" % (lab.upper(), p12)})

    # ---------- 次9：带种子的 20 次置换零分布 ----------
    if "--skip-perm" not in sys.argv:
        try:
            spec = importlib.util.spec_from_file_location(
                "mr", "/path/to/project/work/code/12_matrix_runner.py")
            mr = importlib.util.module_from_spec(spec); spec.loader.exec_module(mr)
            from ldsc_py3 import regressions as reg  # noqa: F401
            meta, traits, M = mr.load_master({"T2D": "T2D", "I9_CHD": "I9_CHD"})
            t1, t2 = traits["T2D"], traits["I9_CHD"]
            okm, same = mr.align(t1, t2)
            z1 = t1["Z"].values[okm.values]
            z2raw = t2["Z"].values[okm.values]
            z2 = np.where(same.values[okm.values], z2raw, -z2raw)
            N1 = t1["N"].values[okm.values].astype(float).reshape((-1, 1))
            N2 = t2["N"].values[okm.values].astype(float).reshape((-1, 1))
            x = meta["L2"][okm.values].reshape((-1, 1))
            w = meta["W_L2"][okm.values].reshape((-1, 1))
            n = len(z1)
            rng = np.random.default_rng(20261028)   # 固定种子（此前对照未留种子，本轮补录）
            out = []
            for i in range(20):
                zp = rng.permutation(z2)
                try:
                    rg = reg.RG(z1.reshape((n, 1)), zp.reshape((n, 1)), x, w, N1, N2,
                                np.array([[M]]), n_blocks=200)
                    v = float(rg.rg_ratio) if rg._negative_hsq is None else np.nan
                    pv = float(rg.p) if rg._negative_hsq is None else np.nan
                except FloatingPointError:
                    v, pv = np.nan, np.nan
                out.append({"perm": i, "rg": v, "p": pv})
                print("[perm] %2d rg=%.4f p=%.3f" % (i, v, pv))
            pd.DataFrame(out).to_csv(PERM, index=False)
            arr = np.array([o["rg"] for o in out if not np.isnan(o["rg"])])
            rows += [
                {"key": "perm_n", "value": str(len(arr)), "note": "全行置换次数（固定种子 20261028）"},
                {"key": "perm_rg_mean", "value": "%.4f" % arr.mean(), "note": "置换零分布 rg 均值"},
                {"key": "perm_rg_sd", "value": "%.4f" % arr.std(ddof=1), "note": "置换零分布 rg 标准差"},
                {"key": "perm_rg_max_abs", "value": "%.4f" % np.abs(arr).max(),
                 "note": "置换零分布 |rg| 最大（真实对 0.366 应远高于此）"},
            ]
        except Exception as exc:
            rows.append({"key": "perm_status", "value": "failed:%s" % str(exc)[:60], "note": "置换对照"})

    kv_summary(rows)
    print("[review] %d 项 -> %s" % (len(rows), SUM.name))
    for r in rows:
        print("   %-34s %s" % (r["key"], r["value"]))


if __name__ == "__main__":
    main()
