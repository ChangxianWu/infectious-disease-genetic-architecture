# -*- coding: utf-8 -*-
""" 51_fig5_anchors.py — 新主图 Figure 5：疾病锚层全景（2026-10-27）

为何要这张图（投稿成功率论证）：v0.1 的 4 张主图覆盖相关矩阵/森林图/位点分类/共定位，
**最具新颖性的疾病锚层（自免 + 心血管代谢，含"效应生物学分界"这一反直觉结论）没有任何图**。
编辑初筛通常是"标题 + 摘要 + 看一张图"——把最强的新结论做成主图，直接提高送审概率。

图内容（单图双面板，300 dpi）：
  面板 A：12 个疾病（10 自免 + CHD/T2D）× 感染的过门对 rg 分布（点 = 每对 rg，横线 = 中位），
          按中位排序，共享者（≥1 BH）蓝、无 BH 灰；标注 mean χ²（功效对照，回答"不共享=功效不足?"）。
  面板 B：该疾病最强单对的 rg 与 95% CI（用 rg±1.96·se），BH 显著者实心。
配色沿用 21_figures.py 主色（#3b6ea5 蓝 / #b0b0b0 灰 / #d08a3e 橙）。
用法：python 51_fig5_anchors.py
输出：work/figures/figure_anchor_landscape.png（主图槽位；release_verify 的 --expect-figures 随之为 5）
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATE, FIG, OUT  # noqa: E402

import matplotlib.pyplot as _mpl
_CMAP = _mpl.colormaps["Blues"]


def _cmap(n_bh):
    """按 BH 对数连续着色（0 用灰，>0 用 Blues 0.35–0.95）。"""
    if n_bh <= 0:
        return M_GREY
    return _CMAP(0.35 + 0.60 * min(n_bh, 10) / 10.0)


M_OK = "#3b6ea5"
M_GREY = "#b0b0b0"
M_ACC = "#d08a3e"

LABEL = {
    "L12_PSORIASIS": "Psoriasis", "M13_RHEUMA": "Rheumatoid arthritis",
    "M13_SJOGREN": "Sjogren syndrome", "M13_ANKYLOSPON": "Ankylosing spondylitis",
    "K11_CD_STRICT2": "Crohn disease", "K11_UC_STRICT2": "Ulcerative colitis",
    "E4_GRAVES_STRICT": "Graves disease", "T1D": "Type 1 diabetes",
    "G6_MS": "Multiple sclerosis", "K11_COELIAC": "Coeliac disease",
    "I9_CHD": "Coronary heart disease", "T2D": "Type 2 diabetes",
}
AUTO = OUT / ("ldsc_rg_autoimmune_%s.csv" % DATE)
COMP = OUT / ("ldsc_rg_comparator_%s.csv" % DATE)
MUNGE = OUT.parent / "data" / "munged" / "munge_counts.csv"


def load():
    frames = []
    for p in (AUTO, COMP):
        d = pd.read_csv(p)
        d["layer"] = "autoimmune" if p == AUTO else "cardiometabolic"
        frames.append(d)
    d = pd.concat(frames, ignore_index=True)
    d = d[(d["power_ok"] == True) & d["rg"].notna()]  # noqa: E712
    mc = pd.read_csv(MUNGE, sep="\t", encoding="utf-8-sig",
                     usecols=["phenocode", "mean_chisq"], low_memory=False)
    mc["mean_chisq"] = pd.to_numeric(mc["mean_chisq"], errors="coerce")
    CHI = mc.dropna().set_index("phenocode")["mean_chisq"].to_dict()
    agg = []
    for code, sub in d.groupby("axis1"):
        best = sub.loc[sub["rg"].idxmax()]
        agg.append({
            "code": code, "label": LABEL.get(code, code),
            "n": len(sub), "median_rg": sub["rg"].median(),
            "n_bh": int((sub["rg_q_bh_layer"] < 0.05).sum()),
            "best_rg": float(best["rg"]), "best_se": float(best["rg_se"]),
            "best_q": float(best["rg_q_bh_layer"]) if pd.notna(best["rg_q_bh_layer"]) else np.nan,
            "chi": float(CHI.get(code, np.nan)), "layer": sub["layer"].iloc[0],
        })
    return d, pd.DataFrame(agg).sort_values(["n_bh", "median_rg"], ascending=[False, False])


def main():
    d, agg = load()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 6.2),
                                   gridspec_kw={"width_ratios": [1, 1.15]})

    # ---- 面板 A：逐病 rg 分布（条 = 中位，点 = 每对）----
    y = np.arange(len(agg))
    for i, (_, r) in enumerate(agg.iterrows()):
        sub = d[d["axis1"] == r["code"]]
        col = _cmap(r["n_bh"])
        jitter = np.linspace(-0.17, 0.17, len(sub)) if len(sub) > 1 else [0]
        ax1.scatter(sub["rg"], np.full(len(sub), i) + jitter, s=16, color=col, alpha=0.55,
                    linewidths=0, zorder=2)
        ax1.plot([r["median_rg"]] * 2, [i - 0.33, i + 0.33], color=col, lw=2.4, zorder=3)
    ax1.set_yticks(y)
    ax1.set_yticklabels(["%s\n(mean $\\chi^2$ %.2f)" % (r["label"], r["chi"])
                         for _, r in agg.iterrows()], fontsize=8.2)
    ax1.set_xlabel("Genetic correlation with infection endpoints\n(gate-passing pairs; line = median)")
    ax1.axvline(0, color="#888888", lw=0.8, ls=":")
    ax1.set_title("A  Sharing by disease", fontsize=10, loc="left")
    ax1.spines[["top", "right"]].set_visible(False)

    # ---- 面板 B：最强单对 rg + 95% CI ----
    # 着色按**疾病级** BH 状态（与面板 A 同一判据），避免"面板 A 蓝、面板 B 灰"的自相矛盾
    # （Crohn 有 1 对 BH 显著但其 rg 最大的那对不显著——显著性看 q 不看 rg 大小）。
    a2 = agg.sort_values(["n_bh", "best_rg"], ascending=[False, False])
    y2 = np.arange(len(a2))
    for i, (_, r) in enumerate(a2.iterrows()):
        dis_sig = r["n_bh"] > 0
        pair_sig = pd.notna(r["best_q"]) and r["best_q"] < 0.05
        col = _cmap(r["n_bh"])
        ax2.errorbar(r["best_rg"], i, xerr=1.96 * r["best_se"], fmt="o", ms=5, color=col,
                     ecolor=col, elinewidth=1.2, capsize=2.5,
                     alpha=0.95 if dis_sig else 0.65, zorder=3)
        if pair_sig:
            tag = "q=%.1e" % r["best_q"]
        elif dis_sig:
            tag = "largest pair n.s."      # 该病有 BH 对，但 rg 最大的这对不显著
        else:
            tag = "no BH pair"
        ax2.text(1.06, i, tag, fontsize=7.0, va="center",
                 color="#333333" if dis_sig else "#888888")
    ax2.set_yticks(y2)
    ax2.set_yticklabels(a2["label"], fontsize=8.2)
    ax2.set_xlim(-0.1, 1.28)
    ax2.axvline(0, color="#888888", lw=0.8, ls=":")
    ax2.set_xlabel("Strongest pair per disease (rg, 95% CI)")
    ax2.set_title("B  Strongest disease–infection pair", fontsize=10, loc="left")
    ax2.spines[["top", "right"]].set_visible(False)

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", ls="", color=_cmap(v), label="%d BH pair%s" % (v, "" if v == 1 else "s"))
               for v in (0, 1, 2, 4, 6, 10)]
    ax1.legend(handles=handles, loc="lower right", fontsize=7.8, frameon=False)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / "figure_anchor_landscape.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print("[fig5] %s（%d 病；BH 病 %d）" % (out, len(agg), int((agg['n_bh'] > 0).sum())))


if __name__ == "__main__":
    main()
