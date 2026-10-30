# -*- coding: utf-8 -*-
""" 21_figures.py — 图件生成（Stage 5，2026-10-11）

图件清单（供稿件主图/补充图）：
  Fig1 figure_rg_heatmap.png            —— 遗传相关矩阵热图（已由 13_stage2_analysis.py 生成）
  Fig2 figure_rg_forest.png             —— 显著配对 rg 森林图（按轴组合分组，带 95% CI 与 q 值标记）
  Fig3 figure_locus_class.png           —— 位点分类条形 + cross-axis 位点轴覆盖矩阵
  Fig4 figure_coloc_regions.png         —— 三区域共定位 PP.H4 条形（CRP / FUT2 / MHC）
用法：python 21_figures.py
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

sys.path.insert(0, "/path/to/project/work/code")
from constants import OUT, FIG  # noqa: E402

m_ok = "#3b6ea5"
m_bad = "#b0b0b0"
AXCOL = {"bacterial": "#3b6ea5", "mixed": "#d08a3e", "viral": "#8e4b8b"}


def axis_of(lab):
    if lab == "PNEUMONIA_MIXED":
        return "mixed"
    if lab in ("INFLUENZA", "VIRAL_WARTS", "COVID19", "HGI_C2", "HGI_A2", "ZOSTER",
               "VIRALPNEUMO", "HERPES_SIMPLEX", "EBV", "VIRAL_CNS", "VARICELLA"):
        return "viral"
    return "bacterial"


def fig_forest():
    df = pd.read_csv(OUT / "ldsc_rg" / "ldsc_matrix_2026-10-11.csv")
    d = df[df["power_ok"] == True].copy()  # noqa: E712
    d["pair"] = d["axis1"] + " × " + d["axis2"]
    d["grp"] = [tuple(sorted([axis_of(a), axis_of(b)])) for a, b in zip(d["axis1"], d["axis2"])]
    d["grpname"] = d["grp"].map(lambda t: " / ".join(t))
    d = d.sort_values(["grpname", "rg"])
    fig, ax = plt.subplots(figsize=(7.4, 12.5))
    y = np.arange(len(d))
    for i, (_, r) in enumerate(d.iterrows()):
        sig = (r["rg_q_bh"] < 0.05)
        col = AXCOL[axis_of(r["axis1"])] if sig else m_bad
        # 形状与颜色双重编码显著性：显著=实心圆，不显著=空心方（色觉障碍下仍可区分）
        ax.errorbar(r["rg"], i, xerr=1.96 * r["rg_se"],
                    fmt="o" if sig else "s", ms=4.6 if sig else 4.0,
                    mfc=col if sig else "none", mec=col, mew=0.9,
                    color=col, ecolor=col, elinewidth=1.0, capsize=2,
                    alpha=0.95 if sig else 0.6)
    ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", mfc="k", mec="k",
                                  label="FDR-significant (q < 0.05)"),
                       plt.Line2D([], [], marker="s", ls="", mfc="none", mec="k",
                                  label="not significant")],
              fontsize=7, frameon=False, loc="lower right")
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.axvline(1, color="k", lw=0.5, ls="--", alpha=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(d["pair"], fontsize=5.6)
    ax.set_xlabel("Genetic correlation (rg) with 95% CI", fontsize=9)
    ax.set_xlim(-0.6, 1.3)
    ax.set_title("Cross-pathogen genetic correlations (bivariate LDSC, FinnGen R12)\n"
                 "power-passed pairs only (n=%d); filled = BH q<0.05 (colour = axis of first trait), "
                 "hollow grey = not significant"
                 % len(d), fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "figure_rg_forest.png", dpi=300)
    plt.close(fig)
    print("[fig] figure_rg_forest.png（%d 对）" % len(d))


def fig_locus_class():
    a = pd.read_csv(OUT / "locus_atlas.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4),
                             gridspec_kw={"width_ratios": [1, 1.7]})
    order = ["broad", "intermediate", "pathogen_specific"]
    cnt = a["locus_class"].value_counts().reindex(order).fillna(0).astype(int)
    cntm = a.groupby(["locus_class", "in_mhc"]).size().unstack(fill_value=0).reindex(order).fillna(0)
    x = np.arange(len(order))
    axes[0].bar(x - 0.2, cntm.get(False, pd.Series(0, index=order)).values, width=0.4,
                label="non-MHC", color="#3b6ea5")
    axes[0].bar(x + 0.2, cntm.get(True, pd.Series(0, index=order)).values, width=0.4,
                label="MHC", color="#c0504d")
    axes[0].set_xticks(x); axes[0].set_xticklabels(order, fontsize=8)
    axes[0].set_ylabel("Number of lead loci", fontsize=8.5)
    axes[0].legend(fontsize=7.5, frameon=False)
    axes[0].set_title("Locus classification (n=%d lead SNPs, %d regions)"
                      % (len(a), a["region_id"].nunique()), fontsize=9)

    cross = a[a["locus_class"] != "pathogen_specific"].copy()
    axes_names = ["bacterial", "mixed", "viral"]
    M = np.zeros((len(cross), len(axes_names)))
    for i, (_, r) in enumerate(cross.iterrows()):
        for k, an in enumerate(axes_names):
            if an in str(r["axes"]).split(","):
                M[i, k] = 1
    M[M == 0] = np.nan
    axes[1].imshow(M, cmap="Blues", vmin=0, vmax=1.6, aspect="auto")
    axes[1].set_xticks(range(len(axes_names))); axes[1].set_xticklabels(axes_names, fontsize=8)
    axes[1].set_yticks(range(len(cross)))
    axes[1].set_yticklabels(["%s (%s)" % (r["SNP"], r["nearby"] if "nearby" in cross.columns
                                          else str(r["nearest_genes"]).split(",")[0])
                             for _, r in cross.iterrows()], fontsize=6.5)
    axes[1].set_title("Pathogen-axis coverage of cross-axis loci", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "figure_locus_class.png", dpi=300)
    plt.close(fig)
    print("[fig] figure_locus_class.png（跨轴位点 %d）" % len(cross))


def fig_coloc():
    c = pd.read_csv(OUT / "coloc_results.csv")
    c["region"] = c["locus"].map({
        "rs3116651": "CRP", "rs876537": "CRP",
        "rs492602": "FUT2", "rs516246": "FUT2",
        "rs495828": "ABO",
        "rs477515": "MHC", "rs1063355": "MHC", "rs17843604": "MHC", "rs2071479": "MHC",
    })
    c = c.dropna(subset=["region"])
    # 跨源配对（含 HGI 端）单列标记，避免把跨源 null 结果与同源结果混读
    c["crosssrc"] = c["trait1"].str.startswith("HGI") | c["trait2"].str.startswith("HGI")
    c["pair"] = c["trait1"] + " × " + c["trait2"]
    c = c.sort_values(["region", "PP_H4"])
    fig, ax = plt.subplots(figsize=(7.2, 8.6))
    rid = {"CRP": "#c0504d", "FUT2": "#3b6ea5", "MHC": "#8064a2", "ABO": "#d08a3e"}
    y = np.arange(len(c))
    for i, (_, r) in enumerate(c.iterrows()):
        col = rid[r["region"]] if r["PP_H4"] > 0.8 else m_bad
        ax.barh(i, r["PP_H4"], color=col, height=0.72,
                hatch="///" if r["crosssrc"] else None,
                edgecolor="k" if r["crosssrc"] else "none", linewidth=0.3)
    ax.axvline(0.8, color="k", ls="--", lw=0.8)
    ax.axvline(0.5, color="k", ls=":", lw=0.7)
    ax.set_yticks(y); ax.set_yticklabels(c["pair"], fontsize=6.2)
    ax.set_xlabel("coloc PP.H4 (shared causal variant)", fontsize=9)
    ax.set_xlim(0, 1.34)
    ax.set_title("Colocalisation at shared regions (coloc.abf; dashed = 0.8 threshold)\n"
                 "colour = region, hatched = cross-source (FinnGen vs COVID HGI)",
                 fontsize=8.6)
    # MHC 组标注：同区多 lead 的重复测试，不计为证据（降低"4 根高柱=4 个证据"的误读）
    mhc_i = [i for i, (_, r) in enumerate(c.iterrows()) if r["region"] == "MHC"]
    if mhc_i:
        y0, y1 = min(mhc_i) - 0.4, max(mhc_i) + 0.4
        ax.annotate("MHC: one pair,\nfour lead variants\n(repeat tests;\nnot evidence)",
                    xy=(1.06, (y0 + y1) / 2), xytext=(1.09, (y0 + y1) / 2),
                    fontsize=6.2, va="center", ha="left",
                    arrowprops=dict(arrowstyle="-", lw=0.7, color="#555555"))
        ax.plot([1.04, 1.04, 1.06], [y0, y1, y1], color="#555555", lw=0.8)
    ax.legend(handles=[Patch(color=v, label=k) for k, v in rid.items()]
              + [Patch(facecolor="white", edgecolor="k", hatch="///", label="cross-source")],
              fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "figure_coloc_regions.png", dpi=300)
    plt.close(fig)
    print("[fig] figure_coloc_regions.png（%d 对）" % len(c))


if __name__ == "__main__":
    fig_forest()
    fig_locus_class()
    fig_coloc()
