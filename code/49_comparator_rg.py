# -*- coding: utf-8 -*-
""" 49_comparator_rg.py — E6 锚层特异性对照：非自免疾病 × 感染 rg（2026-10-25/27 改号）

编号说明（2026-10-27）：原编号 46 与并行会话的 46_submission_format_assets.py 撞号；
本脚本顺延为 49 号，输出文件名不变（ldsc_rg_comparator_<DATE>.csv）。

问题：E3 锚层显示自免（银屑病/RA/Sjögren/强直）×感染共享。读者可问：
"这是自免特异，还是任何常见慢病都与感染共享（一般共病伪影）？"

对照设计（复用 43 号全链口径，不另立标准）：
  - I9_CHD（冠心病，纯心血管、与免疫通路无强先验关联）× 25 感染终点 = 25 对
  - T2D（2 型糖尿病，代谢病但带轻度免疫背景）× 25 感染终点 = 25 对（次对照）
  两性状 sumstats 已在 munge 目录（主矩阵验证对照，2026-10-11 已 munge），无需下载。
  功效门 / BH（本层内）/ |rg|>1.2 出界判 nan —— 与 43 号同规。

预期解读（两种都有信息量，预登记于本 docstring）：
  - CHD 对 rg 近零且 BH=0 → 锚层是**自免特异**（与"银屑病/RA=IL-17/TNF 轴"叙事一致）；
  - CHD 对也显著 → 锚层须降级为"常见慢病共病共享"，自免叙事收缩（如实报告）。
  实测为后者（2026-10-27）：CHD 0.246/9BH、T2D 0.303/10BH → 稿件叙事已改为
  "效应生物学分界"（见 3.10 节），本层即该结论的对照证据。
用法：python 49_comparator_rg.py
输出：work/output/ldsc_rg_comparator_2026-10-11.csv + comparator_anchor_summary.csv
"""
import importlib.util
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, norm

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATE, FINNGEN_ENDPOINTS, OUT  # noqa: E402


def load_spec(spec_path, name):
    spec = importlib.util.spec_from_file_location(name, spec_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mr = load_spec("/path/to/project/work/code/12_matrix_runner.py", "matrix_runner_20b2_c")
OUTCSV = OUT / ("ldsc_rg_comparator_%s.csv" % DATE)
SUMCSV = OUT / "comparator_anchor_summary.csv"
COMPARATORS = ["I9_CHD", "T2D"]  # 主对照=CHD；T2D 为次对照（轻度免疫背景，解释留档）


def main():
    inf = {lab: ph for lab, ph in FINNGEN_ENDPOINTS.items()
           if ph not in mr.MTB_BLOCKED and (mr.MUNGED / ("%s.sumstats" % ph)).exists()}
    comp = {ph: ph for ph in COMPARATORS if (mr.MUNGED / ("%s.sumstats" % ph)).exists()}
    missing = [ph for ph in COMPARATORS if ph not in comp]
    if missing:
        raise SystemExit("[comp] 缺 munge 产物：%s（这些是主矩阵验证对照，理应已 munge）" % missing)
    print("[comp] 感染 %d × 对照 %d" % (len(inf), len(comp)))

    all_labels = dict(inf)
    all_labels.update(comp)
    meta, traits, M = mr.load_master(all_labels)
    print("[comp] 共用 SNP 集 %d | M_5_50=%.0f" % (len(meta["SNP"]), M))

    rows = []
    n_total = len(comp) * len(inf)
    done = 0
    for c in sorted(comp):
        for lab2 in sorted(inf):
            row = mr.run_pair(c, lab2, traits[c], traits[lab2], meta, M)
            rows.append(row)
            done += 1
            if done % 10 == 0 or done == n_total:
                print("[comp %d/%d] %s x %s: rg=%s status=%s"
                      % (done, n_total, c, lab2,
                         ("%.3f" % row["rg"]) if pd.notna(row.get("rg")) else "NA", row["status"]))
    df = pd.DataFrame(rows)
    for i in (1, 2):
        se = df["h2_%d_se" % i]
        df["h2_%d_p" % i] = norm.sf(df["h2_%d" % i] / se)
        df["h2_%d_p" % i] = df["h2_%d_p" % i].where(df["h2_%d" % i].notna(), np.nan)
    df["power_ok"] = (df["status"] == "ok") & (df["h2_1_p"] < 0.05) & (df["h2_2_p"] < 0.05) \
        & (df["mean_chisq_1"] >= 1.02) & (df["mean_chisq_2"] >= 1.02)
    # 2026-10-28 统一推断域：过功效门（与主矩阵同规）
    okm = df["power_ok"] == True  # noqa: E712
    df["rg_q_bh_layer"] = np.nan
    if okm.any():
        df.loc[okm, "rg_q_bh_layer"] = false_discovery_control(
            df.loc[okm, "rg_p"].values, method="bh")
    df.to_csv(OUTCSV, index=False)

    pw = df[df["power_ok"]]
    chd = pw[pw["axis1"] == "I9_CHD"]
    t2d = pw[pw["axis1"] == "T2D"]
    summ = pd.DataFrame([
        {"key": "n_pairs_attempted", "value": str(len(df)), "note": "对照×感染全对"},
        {"key": "n_power_ok", "value": str(int(df["power_ok"].sum())), "note": "过功效门"},
        {"key": "chd_n_power_ok", "value": str(len(chd)),
         "note": "CHD×感染过门对（主对照）"},
        {"key": "chd_median_rg", "value": ("%.3f" % chd["rg"].median()) if len(chd) else "NA",
         "note": "CHD 过门对 rg 中位（对照读数：自免强共享组中位 0.29–0.37）"},
        {"key": "chd_n_bh", "value": str(int((chd["rg_q_bh_layer"] < 0.05).sum())) if len(chd) else "0",
         "note": "CHD 过门对中本层 BH<0.05（对照读数：自免银屑病 9/RA 9）"},
        {"key": "t2d_n_power_ok", "value": str(len(t2d)), "note": "T2D×感染过门对（次对照）"},
        {"key": "t2d_median_rg", "value": ("%.3f" % t2d["rg"].median()) if len(t2d) else "NA",
         "note": "T2D 过门对 rg 中位"},
        {"key": "t2d_n_bh", "value": str(int((t2d["rg_q_bh_layer"] < 0.05).sum())) if len(t2d) else "0",
         "note": "T2D 过门对中本层 BH<0.05"},
        {"key": "interpretation_preset", "value": "preregistered",
         "note": "两结局预登记于本脚本 docstring（自免特异 vs 一般共病，均如实报告）"},
    ])
    summ.to_csv(SUMCSV, index=False)
    print("[comp] CHD：过门 %s，中位 %s，BH %s" %
          (len(chd), ("%.3f" % chd["rg"].median()) if len(chd) else "NA",
           int((chd["rg_q_bh_layer"] < 0.05).sum()) if len(chd) else 0))


if __name__ == "__main__":
    main()
