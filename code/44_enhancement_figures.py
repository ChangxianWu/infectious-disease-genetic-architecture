# -*- coding: utf-8 -*-
""" 44_enhancement_figures.py — 增强层补充图（2026-10-25）

输出到 work/figures_supp/。
  figS1_mhc_subregion.png —— MHC lead 亚区 × 跨轴覆盖（E2 位置层证据）
  figS2_prs_transfer.png —— BH 对 PRS 传递 R² Top15（E5 量级诚实呈现）
配色沿用 21_figures.py（m_ok 蓝 / 类 III 橙 / 灰）。
用法：python 44_enhancement_figures.py
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "/path/to/project/work/code")
from constants import OUT, FIG  # noqa: E402

SUP = FIG.parent / "figures_supp"
SUP.mkdir(parents=True, exist_ok=True)
C2 = "#3b6ea5"   # HLA class II
C3 = "#d08a3e"   # class III / extended
CROSS = "#8e4b8b"
WITHIN = "#3b6ea5"
GREY = "#b0b0b0"

t1 = pd.read_csv(OUT / "table1.csv")
AX = t1.set_index("Trait")["Pathogen_axis"].to_dict()


def fig_s1():
    d = pd.read_csv(OUT / "mhc_subregion_leads.csv")
    d = d[d["source"] == "lead"].copy()
    d["label"] = d["snp"] + "\n(" + d["classified_gene"] + ")"
    colmap = {"HLA_class_II": C2, "class_III_extended": C3, "HLA_class_I": "#4b9b6e",
              "MHC_other": GREY}
    d["col"] = d["subregion"].map(colmap)
    d = d.sort_values(["n_axes_sig", "subregion"], ascending=[False, True])
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    y = np.arange(len(d))
    ax.barh(y, d["n_axes_sig"], color=d["col"], height=0.55, alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(d["label"], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("Pathogen axes reached (of bacterial / mixed / viral)")
    ax.set_xticks([0, 1, 2, 3])
    for i, (_, r) in enumerate(d.iterrows()):
        if r["n_axes_sig"] >= 2:
            ax.text(r["n_axes_sig"] + 0.05, i, "cross-axis", va="center", fontsize=7.5,
                    color="#444444")
    handles = [plt.Rectangle((0, 0), 1, 1, color=C2, alpha=0.9),
               plt.Rectangle((0, 0), 1, 1, color=C3, alpha=0.9)]
    ax.legend(handles, ["HLA class II (5/7)", "class III/extended (2/7)"],
              loc="lower right", fontsize=8, frameon=False)
    ax.set_title("MHC-proximal lead variants: subregion vs cross-pathogen reach",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(SUP / "figS1_mhc_subregion.png", dpi=300)
    plt.close(fig)
    print("[fig] figS1_mhc_subregion.png")


def fig_s2():
    d = pd.read_csv(OUT / "prs_transfer.csv")
    d["grp"] = [AX.get(s, "") == AX.get(t, "") for s, t in zip(d["prs_from"], d["predicted_endpoint"])]
    d["R2pct"] = d["expected_R2"] * 100
    top = d.sort_values("R2pct", ascending=False).head(15).copy()
    top["label"] = top["prs_from"] + " → " + top["predicted_endpoint"]
    med_w = d[d["grp"]]["R2pct"].median()
    med_c = d[~d["grp"]]["R2pct"].median()
    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    y = np.arange(len(top))
    ax.barh(y, top["R2pct"], color=[WITHIN if g else CROSS for g in top["grp"]],
            height=0.6, alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(top["label"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Expected PRS transfer R² (%) — observed scale")
    ax.axvline(med_w, color=WITHIN, ls="--", lw=1)
    ax.axvline(med_c, color=CROSS, ls=":", lw=1.2)
    ax.text(med_w, -0.8, " within-axis median %.3f%%" % med_w, fontsize=7.5, color=WITHIN)
    ax.text(med_c, len(top) - 0.2, " cross-axis median %.3f%%" % med_c, fontsize=7.5, color=CROSS)
    handles = [plt.Rectangle((0, 0), 1, 1, color=WITHIN, alpha=0.9),
               plt.Rectangle((0, 0), 1, 1, color=CROSS, alpha=0.9)]
    ax.legend(handles, ["same-axis pair", "cross-axis pair"], loc="lower right",
              fontsize=8, frameon=False)
    ax.set_title("Polygenic transfer among BH pairs (top 15 of 70 directed rows)",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(SUP / "figS2_prs_transfer.png", dpi=300)
    plt.close(fig)
    print("[fig] figS2_prs_transfer.png")


def fig_s3():
    """E3 自免锚层：每自免终点对感染终点的过门对中位 rg（含 BH 计数）。"""
    d = pd.read_csv(OUT / "ldsc_rg_autoimmune_2026-10-11.csv")
    pw = d[d["power_ok"] == True]  # noqa: E712
    rows = []
    for a, sub in pw.groupby("axis1"):
        rows.append({"auto": a, "median_rg": sub["rg"].median(),
                     "n_bh": int((sub["rg_q_bh_layer"] < 0.05).sum()),
                     "n_pairs": len(sub), "n_pos": int((sub["rg"] > 0).sum())})
    r = pd.DataFrame(rows).sort_values("median_rg", ascending=True)
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    y = np.arange(len(r))
    cols = [C2 if b > 0 else GREY for b in r["n_bh"]]
    ax.barh(y, r["median_rg"], color=cols, height=0.6, alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(r["auto"], fontsize=8.5)
    ax.set_xlabel("Median rg vs infection endpoints (gate-passing pairs, MHC-excluded)")
    for i, (_, row) in enumerate(r.iterrows()):
        ax.text(row["median_rg"] + 0.006, i, "%d/%d BH, %d/%d pos" %
                (row["n_bh"], row["n_pairs"], row["n_pos"], row["n_pairs"]),
                va="center", fontsize=7.5, color="#444444")
    handles = [plt.Rectangle((0, 0), 1, 1, color=C2, alpha=0.9),
               plt.Rectangle((0, 0), 1, 1, color=GREY, alpha=0.9)]
    ax.legend(handles, ["has BH<0.05 pairs", "no BH pair"], loc="lower right",
              fontsize=8, frameon=False)
    ax.set_title("Autoimmune anchor layer: per-disease sharing with infections",
                 fontsize=10)
    ax.set_xlim(0, r["median_rg"].max() * 1.45)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(SUP / "figS3_autoimmune_anchor.png", dpi=300)
    plt.close(fig)
    print("[fig] figS3_autoimmune_anchor.png")


if __name__ == "__main__":
    fig_s1()
    fig_s2()
    fig_s3()
    print("[fig] 输出目录：%s" % SUP)
