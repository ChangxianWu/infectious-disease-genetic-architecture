# -*- coding: utf-8 -*-
""" 23_mhc_sensitivity.py — MHC 敏感性层（Stage 4b，2026-10-11）

问题：主矩阵用官方 `weights.hm3_noMHC` → MHC（chr6:25–35 Mb）SNP 整体退出回归。
      因此主结果是"非 MHC 架构"，MHC 共享只能靠位点层看到。需给出 MHC **包含**时的
      对照结果，才能量化 MHC 对跨病原共享的贡献。

实现（不引入新数据，只用已有官方资产）：
  - 参考 LD 与权重都用同一份官方 `LDscore.*.l2.ldscore.gz`（含 MHC 的 1000G EUR 面板）；
    权重按 LDSC 权重函数的反条件方差形式自加权（w_ld = 该 SNP 自身的 l2 值），
    等价于"不施加 MHC 排除"的最近似口径。
  - 同一性状对在两种口径下的 rg 并列报告；口径差异即 MHC 贡献的下界估计
    （因为除 MHC 外，两种口径的 SNP 集差异也含 noMHC 面板的其余差异）。

**诚实性声明**：本层不是官方 `--w-ld-chr` 与 `--w-ld-chr(含MHC)` 的严格对照，
而是"用含 MHC 的 l2 自加权"近似；解读上限为"MHC 是否改变结论方向/显著性"，
不得据此给出 MHC 的定量贡献百分比。

用法：python 23_mhc_sensitivity.py
输出：work/output/mhc_sensitivity.csv
"""
import sys
import time
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, norm

sys.path.insert(0, "/path/to/project/work/code")
sys.path.insert(0, "/path/to/project/work/code/ldsc_py3")
import regressions as reg  # noqa: E402
from constants import OUT, GWAS_DATA, DATE  # noqa: E402
import importlib.util  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "mr", "/path/to/project/work/code/12_matrix_runner.py")
mr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mr)

LDSC_REF = GWAS_DATA / "ldsc_ref"
# 复核性状：主矩阵中通过功效门的 10 个（含跨轴对与轴内对，覆盖 CRT/FUT2/ABO 所在性状）
PAIR_LABELS = ["INTESTINAL", "GASTRO_NOS", "PYELONEPHR", "ERYSIPELAS", "BACT_NOS",
               "PNEUMOBACT", "SEPSIS", "PNEUMONIA_MIXED", "INFLUENZA", "OSTEOMYEL"]


def load_ref_full():
    """含 MHC 的参考：l2 与自加权权重同源。"""
    l2s, m550 = [], []
    for i in range(1, 23):
        l2s.append(pd.read_csv(LDSC_REF / "LDscore" / ("LDscore.%d.l2.ldscore.gz" % i),
                               sep="\t", usecols=["SNP", "L2"]))
        m550.append(float(open(LDSC_REF / "LDscore" / ("LDscore.%d.l2.M_5_50" % i)).read().strip()))
    l2 = pd.concat(l2s, ignore_index=True).drop_duplicates("SNP", keep="first")
    return l2, float(np.sum(m550))


def main():
    t0 = time.time()
    l2full, M = load_ref_full()
    print("[mhc] 含 MHC 参考 l2=%d M_5_50=%.0f" % (len(l2full), M))
    label_to_ph = {l: ph for l, ph in mr.FINNGEN_ENDPOINTS.items() if l in PAIR_LABELS}
    meta, traits, _ = mr.load_master(label_to_ph)
    print("[mhc] noMHC 共用索引 %d" % len(meta["SNP"]))

    # 含 MHC 索引：以 l2full 顺序，重取各性状 Z/A1/A2/N
    idx_full = pd.Index(l2full["SNP"])
    full = {}
    for lab, ph in label_to_ph.items():
        d = pd.read_csv(mr.MUNGED / ("%s.sumstats" % ph), sep="\t",
                        usecols=["SNP", "A1", "A2", "N", "Z"]).drop_duplicates("SNP")
        full[lab] = d.set_index("SNP").reindex(idx_full)
    z_l2 = l2full["L2"].values.astype(np.float64)

    rows = []
    pairs = list(combinations(sorted(label_to_ph), 2))
    for k, (l1, l2lab) in enumerate(pairs, 1):
        a, b = full[l1], full[l2lab]
        same = (a["A1"] == b["A1"]) & (a["A2"] == b["A2"])
        swap = (a["A1"] == b["A2"]) & (a["A2"] == b["A1"])
        ok = (same | swap).values
        ok = ok & np.isfinite(a["Z"].values) & np.isfinite(b["Z"].values)
        z1 = a["Z"].values[ok]
        z2raw = b["Z"].values[ok]
        z2 = np.where(same.values[ok], z2raw, -z2raw)
        N1 = a["N"].values[ok].astype(float).reshape((-1, 1))
        N2 = b["N"].values[ok].astype(float).reshape((-1, 1))
        x = z_l2[ok].reshape((-1, 1))
        w = z_l2[ok].reshape((-1, 1))   # 自加权（含 MHC）
        n = len(z1)
        try:
            rg = reg.RG(z1.reshape((n, 1)), z2.reshape((n, 1)), x, w, N1, N2,
                        np.array([[M]]), n_blocks=200)
            status = "ok" if rg._negative_hsq is None else "negative_hsq"
        except FloatingPointError as exc:
            status = "fp_error:%s" % str(exc)[:30]
            rg = None
        # 是否含 MHC SNP
        rows.append({"trait1": l1, "trait2": l2lab, "n_snp_incl_mhc": n,
                     "status": status,
                     "rg_incl_mhc": (rg.rg_ratio if rg and status == "ok" else np.nan),
                     "rg_se_incl_mhc": (rg.rg_se if rg and status == "ok" else np.nan),
                     # p 下限 1e-300：真值可能更小，但 CSV 写精确 0 会被 data_audit 的
                     # "P 值写作 0"断言拦下（该断言为稿件语义而设，数值列同样适用）。
                     # 下限化只压缩浮点下溢区间，不改任何显著性判定（阈值 0.05 远高于此）。
                     "p_incl_mhc": (max(float(rg.p), 1e-300) if rg and status == "ok" else np.nan)})
    df = pd.DataFrame(rows)

    # 与主矩阵（noMHC）并排
    main = pd.read_csv(OUT / "ldsc_rg" / ("ldsc_matrix_%s.csv" % DATE))
    key = lambda r: tuple(sorted([r["axis1"], r["axis2"]]))
    main_by = {}
    for _, r in main.iterrows():
        main_by[tuple(sorted([r["axis1"], r["axis2"]]))] = r
    df["rg_excl_mhc"] = df.apply(
        lambda r: main_by.get(tuple(sorted([r["trait1"], r["trait2"]])), pd.Series()).get("rg", np.nan),
        axis=1)
    df["q_excl_mhc"] = df.apply(
        lambda r: main_by.get(tuple(sorted([r["trait1"], r["trait2"]])), pd.Series()).get("rg_q_bh", np.nan),
        axis=1)
    df["p_excl_mhc"] = df.apply(
        lambda r: main_by.get(tuple(sorted([r["trait1"], r["trait2"]])), pd.Series()).get("rg_p", np.nan),
        axis=1)
    df["pair"] = df["trait1"] + " x " + df["trait2"]
    df["mhc_effect"] = df["rg_incl_mhc"] - df["rg_excl_mhc"]
    # 两侧 p 值统一下限 1e-300：真值可能更小，但 CSV 写精确 0 会被 data_audit 的
    # "P 值写作 0"断言拦下；下限化只压缩浮点下溢区间，不改任何显著性判定（阈值 0.05 远高于此）。
    for _col in ("p_excl_mhc", "p_incl_mhc"):
        df[_col] = pd.to_numeric(df[_col], errors="coerce").clip(lower=1e-300)
    # 含 MHC 侧同样过 BH（与 noMHC 侧口径对齐；规则 40/85：族口径必须一致）
    ok_i = df["status"] == "ok"
    df["q_incl_mhc"] = np.nan
    if ok_i.sum():
        df.loc[ok_i, "q_incl_mhc"] = false_discovery_control(
            df.loc[ok_i, "p_incl_mhc"].values, method="bh")
    # 翻转判定一律**同口径比较**（BH vs BH、p vs p），避免 q 对 p 的伪翻转
    df["flip_bh05"] = ((df["q_excl_mhc"] < 0.05) & (df["q_incl_mhc"] >= 0.05)) | \
                      ((df["q_excl_mhc"] >= 0.05) & (df["q_incl_mhc"] < 0.05))
    df["flip_p05"] = ((df["p_excl_mhc"] < 0.05) & (df["p_incl_mhc"] >= 0.05)) | \
                     ((df["p_excl_mhc"] >= 0.05) & (df["p_incl_mhc"] < 0.05))
    df["conclusion_flip"] = df["flip_bh05"]
    df.to_csv(OUT / "mhc_sensitivity.csv", index=False)
    ok = df[df["status"] == "ok"]
    print("[mhc] %d 对；含 MHC 可估 %d；BH 口径翻转 %d 对 / 名义 p 口径翻转 %d 对（%.1f min）"
          % (len(df), len(ok), int(df["flip_bh05"].sum()), int(df["flip_p05"].sum()),
             (time.time() - t0) / 60))
    print()
    cols = ["pair", "rg_excl_mhc", "q_excl_mhc", "rg_incl_mhc", "q_incl_mhc",
            "mhc_effect", "flip_bh05", "flip_p05"]
    print(ok.sort_values("mhc_effect", key=lambda s: -s.abs())[cols].to_string(index=False))


if __name__ == "__main__":
    main()
