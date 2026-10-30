# -*- coding: utf-8 -*-
""" 32_supp_xlsx.py — 补充表容器（单 xlsx 多 sheet；structure_check 惯例 c）

为何用单 xlsx：structure_check.py:156–169 的"惯例 (c)"——单 xlsx 多 sheet 容器，
sheet 名 `S<N>_<slug>`，**sheet 名即补表实体名**，正文引用按家族号核对，
最不易产生孤儿/悬空引用（20A3/20B1 同规）。

补充表设计（供稿件与本课题主张的完整可核对性）：
  S1  endpoints       25 个 FinnGen 终点 + 2 个外部层：轴、病例数、h²、mean χ²、功效门
  S2  rg_matrix       全部 300 对 bivariate LDSC 结果（含 status/rg_raw 分档与双口径多重检验）
  S3  hgi_external    COVID HGI 外部层 50 对（跨源分层）
  S4  locus_atlas     27 个 lead 位点：分类、覆盖轴、MHC 标注、区归并、最近基因
  S5  coloc_abf       63 对 coloc.abf 四后验 + top SNP + 先验敏感性（三档 p12 并列）
  S6  coloc_susie     14 对多信号共定位（信号结构 + 是否双信号）
  S7  mhc_sens        MHC 敏感性 45 对（含/不含 MHC 的 rg 与 q 并列）
  S8  gene_sets       HGNC 基因集富集 + 逐基因 group 扫描 + Reactome 通路 hub 检查
  S9  qc_munge        逐终点 munge 质控计数（总变异 → 各过滤步 → 最终 SNP）

用法：python 32_supp_xlsx.py
输出：work/attachments/Supplementary_Tables_project.xlsx
"""
import os
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import pandas as pd

BASE = Path(os.path.realpath("/path/to/project"))
ROOT = BASE / "project"
OUT = ROOT / "work/output"
DATA = ROOT / "work/data"
# 产物落两处：work/attachments/ 与 投稿材料包/02_附件与清单/。
# 原因：structure_check --root <课项目录> 是从 **投稿材料包** 底下找 02_附件与清单 的；
# 只写 work/ 会造成"实体集为空 → 假报补表引用悬空"（规则 190 同族的假失败）。
ATT = ROOT / "work/attachments"
PKG_ATT = ROOT / "投稿材料包" / "02_附件与清单"      # 包内实体（门控按此目录识别补表）
ATT.mkdir(parents=True, exist_ok=True)

SHEETS = [
    ("S1_endpoints", OUT / "table1.csv", ","),
    ("S2_rg_matrix", OUT / "ldsc_rg" / "ldsc_matrix_2026-10-11.csv", ","),
    ("S3_hgi_external", OUT / "hgi_external_pairs.csv", ","),
    ("S4_locus_atlas", OUT / "locus_atlas.csv", ","),
    ("S5_coloc_abf", OUT / "coloc_results.csv", ","),
    ("S6_coloc_susie", OUT / "coloc_susie_results.csv", ","),
    ("S7_mhc_sens", OUT / "mhc_sensitivity.csv", ","),
    ("S8_gene_sets", OUT / "gene_set_enrichment.csv", ","),
    ("S9_qc_munge", DATA / "munged" / "munge_counts.csv", "\t"),
]

EXTRA = [
    # 主表之外的补充证据（同一容器内，便于读者一次拿到）
    ("S8b_hgnc_scan", OUT / "gene_group_scan.csv", ","),
    ("S8c_reactome", OUT / "pathway_hub_check.csv", ","),
    # 屏障/急性期本体面板（54 号；事后选定，见脚本诚实性声明）
    ("S8d_ontology", OUT / "ontology_sets.csv", ","),
    ("S5b_coloc_prior", OUT / "coloc_p12_sensitivity.csv", ","),
    ("S1b_axis_summary", OUT / "ldsc_matrix_axis_summary.csv", ","),
    ("S2b_no_covid", OUT / "sensitivity_no_covid.csv", ","),
    # 创新性增强层（E1–E5，2026-10-25；逐层一 sheet，正文按家族号引用）
    ("S10_iei_overlap", OUT / "iei_overlap.csv", ","),
    ("S11_mhc_subregion", OUT / "mhc_subregion_leads.csv", ","),
    ("S12_b1_triangulation", OUT / "b1_triangulation.csv", ","),
    ("S13_prs_transfer", OUT / "prs_transfer.csv", ","),
    ("S14_autoimmune_rg", OUT / "ldsc_rg_autoimmune_2026-10-11.csv", ","),
    ("S15_autoimmune_mhc", OUT / "autoimmune_mhc_sens.csv", ","),
    ("S16_axis_robustness", OUT / "axis_robustness.csv", ","),
    ("S17_comparator_rg", OUT / "ldsc_rg_comparator_2026-10-11.csv", ","),
    # 补充敏感性总成（52 号）：M3/M4/M8/M9 + 次4/6/7/9 全部读数
    ("S18_sensitivity", OUT / "sensitivity_extended_summary.csv", ","),
    # 敏感性读数（rr2_*，85 项）：正文多处引用 S18，故并入同一家族 S18b。
    # 表名用中性词（sensitivity），不使用带轮次的命名——对外交付物不记录补充轮次（规则 230）
    ("S18b_sensitivity_detail", OUT / "sensitivity_detailed_summary.csv", ","),
    # 逐性状 SuSiE 精细定位（2026-10-13）：可信集跨性状一致性 = 变异级"共享"证据
    # S19 系列（并行会话产出的逐性状 fine-mapping）本版**不进补充表**：其 max_pip 在 17/19 行
    # 精确为 1、且几乎每个性状都解析出恰好 5 个信号，与 M-G 同型的退化输出；正文未引用，
    # 结构门会判"未引用的补充实体"。保留源文件于 work/output 供追溯，不进投稿件。
]


def main():
    path = ATT / "Supplementary_Tables_project.xlsx"
    written = []
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, src, sep in SHEETS + EXTRA:
            if not Path(src).exists():
                print("[supp] 缺源文件，跳过 %s ← %s" % (name, src.name))
                continue
            df = pd.read_csv(src, sep=sep, dtype=str)
            # 空单元格统一为 "NA"：读者不必区分"空"与"缺失"（也避免 pandas 把 NaN 写空）
            df = df.fillna("NA")
            df.to_excel(xw, sheet_name=name, index=False)
            written.append((name, len(df), df.shape[1]))
    # 同步到包内实体目录（structure_check 识别补表的唯一位置）
    PKG_ATT.mkdir(parents=True, exist_ok=True)
    import shutil
    shutil.copy2(path, PKG_ATT / path.name)
    print("[supp] %s" % path)
    print("[supp] 已同步 -> %s" % (PKG_ATT / path.name))
    for n, r, c in written:
        print("   %-18s %5d 行 x %2d 列" % (n, r, c))
    # 家族号覆盖自查（structure_check 按 S<N> 家族归并；子表 S5b_ 形式归并到 S5）
    import re as _re
    fam = sorted({int(_re.match(r"S(\d+)", n).group(1)) for n, _, _ in written
                  if _re.match(r"S\d+", n)})
    print("[supp] 家族号: %s（共 %d 个）" % (fam, len(fam)))


if __name__ == "__main__":
    main()
