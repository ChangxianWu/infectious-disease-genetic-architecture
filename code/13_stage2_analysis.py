# -*- coding: utf-8 -*-
""" 13_stage2_analysis.py — 结果判读（Stage 2/4，2026-10-11）

输入：work/output/ldsc_rg/ldsc_matrix_2026-10-11.csv（300 对主矩阵）
      work/data/munged/*.sumstats（含 HGI_C2/HGI_A2 外部层）
输出：
  1. ldsc_matrix_axis_summary.csv —— 轴内/跨轴分层统计（广谱证据的第一步）
  2. hgi_external_pairs.csv —— HGI C2/A2 × 25 FinnGen 端点的外部层配对
  3. sensitivity_no_covid.csv —— 剔除 COVID19 端点后重跑全矩阵（设计文档风险②）
  4. figure_rg_heatmap.png —— 遗传相关矩阵热图（通过功效门的性状）
用法：python 13_stage2_analysis.py
"""
import sys
import time
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, norm

sys.path.insert(0, "/path/to/project/work/code")
sys.path.insert(0, "/path/to/project/work/code/ldsc_py3")
import importlib.util  # noqa: E402
from constants import DATA, OUT, FIG, DATE  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "mr", "/path/to/project/work/code/12_matrix_runner.py")
mr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mr)
mat_runner = mr

AXIS = {
    # 细菌
    "INTESTINAL": "bacterial", "GASTRO_NOS": "bacterial", "PYELONEPHR": "bacterial",
    "ERYSIPELAS": "bacterial", "BACT_NOS": "bacterial", "PNEUMOBACT": "bacterial",
    "SEPSIS": "bacterial", "OSTEOMYEL": "bacterial", "MENINGBACT": "bacterial",
    "ENDOCARD": "bacterial",
    "PNEUMONIA_MIXED": "mixed",
    "COVID19": "viral", "INFLUENZA": "viral", "ZOSTER": "viral", "VIRALPNEUMO": "viral",
    "HERPES_SIMPLEX": "viral", "VIRAL_WARTS": "viral", "EBV": "viral",
    "VIRAL_CNS": "viral", "VARICELLA": "viral",
    "MYCOSES": "fungal", "CANDIDIASIS": "fungal", "DERMATOPHYTOSIS": "fungal",
    "INFESTATIONS": "parasitic", "PROTOZOAL": "parasitic",
}


def annotate(df):
    df = df.copy()
    df["axis_label_1"] = df["axis1"].map(AXIS)
    df["axis_label_2"] = df["axis2"].map(AXIS)
    df["same_axis"] = df["axis_label_1"] == df["axis_label_2"]
    df["bact_bact"] = (df["axis_label_1"] == "bacterial") & (df["axis_label_2"] == "bacterial")
    return df


def main():
    # ---------- 1. 主矩阵轴分层 ----------
    df = annotate(pd.read_csv(OUT / "ldsc_rg" / ("ldsc_matrix_%s.csv" % DATE)))
    ok = df[df["power_ok"] == True].copy()  # noqa: E712
    rows = []
    for name, sub in (("all_power_ok", ok),
                      ("bacterial_internal", ok[ok["bact_bact"]]),
                      ("cross_axis", ok[~ok["same_axis"]]),
                      ("fungal_or_parasitic_any", ok[(ok["axis_label_1"].isin(["fungal", "parasitic"]))
                                                     | (ok["axis_label_2"].isin(["fungal", "parasitic"]))])):
        if len(sub) == 0:
            rows.append({"stratum": name, "n_pairs": 0})
            continue
        rows.append({
            "stratum": name, "n_pairs": len(sub),
            "n_sig_bh05": int((sub["rg_q_bh"] < 0.05).sum()),
            "n_sig_bonf05": int((sub["rg_p_bonf"] < 0.05).sum()),
            "median_rg": float(sub["rg"].median()), "median_abs_rg": float(sub["rg"].abs().median()),
            "max_rg": float(sub["rg"].max()),
        })
    summ = pd.DataFrame(rows)
    summ.to_csv(OUT / "ldsc_matrix_axis_summary.csv", index=False)
    print("[S2] 轴分层:\n%s\n" % summ.to_string(index=False))

    # ---------- 2. HGI 外部层 ----------
    label_to_ph = {lab: ph for lab, ph in mat_runner.FINNGEN_ENDPOINTS.items()
                   if ph not in mat_runner.MTB_BLOCKED
                   and (mat_runner.MUNGED / ("%s.sumstats" % ph)).exists()}
    label_to_ph["HGI_C2"] = "HGI_C2"
    label_to_ph["HGI_A2"] = "HGI_A2"
    print("[S2] 载入主表（含 HGI 外部层）...")
    meta, traits, M = mat_runner.load_master(label_to_ph)
    ext = []
    for hgi in ("HGI_C2", "HGI_A2"):
        for lab, ph in mat_runner.FINNGEN_ENDPOINTS.items():
            if ph in mat_runner.MTB_BLOCKED or lab not in traits:
                continue
            row = mat_runner.run_pair(lab, hgi, traits[lab], traits[hgi], meta, M)
            row["axis1"], row["axis2"] = lab, hgi
            ext.append(row)
            print("[S2-ext %d] %s x %s: rg=%s (%s)" % (len(ext), lab, hgi,
                  ("%.3f" % row["rg"]) if pd.notna(row.get("rg")) else "NA", row["status"]))
    ext_df = pd.DataFrame(ext)
    for i in (1, 2):
        se = ext_df["h2_%d_se" % i]
        ext_df["h2_%d_p" % i] = norm.sf(ext_df["h2_%d" % i] / se)
    ext_df["power_ok"] = (ext_df["status"] == "ok") & (ext_df["h2_1_p"] < 0.05) \
        & (ext_df["h2_2_p"] < 0.05) & (ext_df["mean_chisq_1"] >= 1.02) \
        & (ext_df["mean_chisq_2"] >= 1.02)
    oke = ext_df["power_ok"] == True  # noqa: E712
    ext_df["rg_q_bh"] = np.nan
    if oke.sum():
        ext_df.loc[oke, "rg_q_bh"] = false_discovery_control(
            ext_df.loc[oke, "rg_p"].values, method="bh")
    ext_df.to_csv(OUT / "hgi_external_pairs.csv", index=False)
    print("[S2] HGI 外部层 %d 对（power_ok=%d）" % (len(ext_df), int(oke.sum())))

    # ---------- 3. 敏感性：剔除 COVID19 ----------
    lab_nc = [l for l in label_to_ph if l != "COVID19"]
    sens = []
    for l1, l2 in combinations(sorted(lab_nc), 2):
        if l1.startswith("HGI") or l2.startswith("HGI"):
            continue
        row = mat_runner.run_pair(l1, l2, traits[l1], traits[l2], meta, M)
        sens.append(row)
    sdf = annotate(pd.DataFrame(sens))
    for i in (1, 2):
        se = sdf["h2_%d_se" % i]
        sdf["h2_%d_p" % i] = norm.sf(sdf["h2_%d" % i] / se)
    sdf["power_ok"] = (sdf["status"] == "ok") & (sdf["h2_1_p"] < 0.05) \
        & (sdf["h2_2_p"] < 0.05) & (sdf["mean_chisq_1"] >= 1.02) & (sdf["mean_chisq_2"] >= 1.02)
    oks = sdf["power_ok"] == True  # noqa: E712
    sdf["rg_q_bh"] = np.nan
    sdf["rg_p_bonf"] = np.nan
    if oks.sum():
        sdf.loc[oks, "rg_q_bh"] = false_discovery_control(sdf.loc[oks, "rg_p"].values, method="bh")
        sdf.loc[oks, "rg_p_bonf"] = np.minimum(1.0, sdf.loc[oks, "rg_p"].values * int(oks.sum()))
    sdf.to_csv(OUT / "sensitivity_no_covid.csv", index=False)
    print("[S2] 敏感（剔 COVID）%d 对：power_ok=%d BH<0.05=%d" %
          (len(sdf), int(oks.sum()), int((sdf["rg_q_bh"] < 0.05).sum())))

    # ---------- 4. 热图 ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    good = sorted(set(ok["axis1"]) | set(ok["axis2"]))
    n = len(good)
    Mx = pd.DataFrame(np.nan, index=good, columns=good)
    for _, r in ok.iterrows():
        if r["axis1"] in Mx.index and r["axis2"] in Mx.columns:
            Mx.loc[r["axis1"], r["axis2"]] = r["rg"]
            Mx.loc[r["axis2"], r["axis1"]] = r["rg"]
    order = sorted(good, key=lambda l: (AXIS.get(l, "z"), l))
    Mx = Mx.loc[order, order]
    fig, ax = plt.subplots(figsize=(8.2, 7.0))
    im = ax.imshow(Mx.values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(order, rotation=90, fontsize=7)
    ax.set_yticklabels(order, fontsize=7)
    for i in range(n):
        for j in range(n):
            v = Mx.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=5.6,
                        color="black")
    ax.set_title("Cross-pathogen genetic correlation matrix\n"
                 "(FinnGen R12 infections, EUR, bivariate LDSC; power-passed traits)",
                 fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.8, label="rg")
    fig.tight_layout()
    fig.savefig(FIG / "figure_rg_heatmap.png", dpi=300)
    Mx.to_csv(OUT / "rg_matrix_power_ok.csv")
    print("[S2] 热图 -> %s（%d 性状）" % (FIG / "figure_rg_heatmap.png", n))


if __name__ == "__main__":
    main()
