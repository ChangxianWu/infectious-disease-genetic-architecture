# -*- coding: utf-8 -*-
""" 46_submission_format_assets.py — 投稿格式派生件（HMG L4/L5）

产出三类**派生件**（正文/图/表内容零改动，只改容器格式）：
  1. work/figures_tiff/Figure1..4.tif          主图 TIFF，CMYK，300 dpi，LZW
  2. work/figures_supp_tiff/Figure_S1..S3.tif  补充图同上
  3. work/attachments/Supplementary_Tables_project.xls  补充表可编辑 .xls 版

依据：
  * HMG：图件首选 TIFF、>=300 dpi、印刷彩图 CMYK（非 RGB）；补充材料允许格式清单未列 .xlsx。
  * PNG 保留给 BMC（兜底刊直接可用）与网页预览；TIFF/CMYK 供 HMG 主投。

纪律：
  * 只做格式转换，不重采样、不改像素尺寸（300 dpi 已是源图实际分辨率；不伪造 600 dpi）。
  * RGB→CMYK 为**无 ICC 配置文件的朴素分离**；若刊方要求特定印刷 profile，须用色彩管理
    工具从源图重转（本脚本不声称色彩管理）。
  * .xls 由 .xlsx 单一真源导出，导出后**逐 cell 回读比对**（数值容差 1e-9，字符串精确）；
    任一 sheet/行列不符即失败，不静默。
  * 产物写 manifest（SHA256/尺寸/色彩空间/dpi），供 release 侧复核。

用法：python 46_submission_format_assets.py
"""
import csv
import hashlib
import os
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from PIL import Image

ROOT = Path("/path/to/project") / "project"
FIG = ROOT / "work" / "figures"
FIGSUPP = ROOT / "work" / "figures_supp"
OUT_MAIN = ROOT / "work" / "figures_tiff"
OUT_SUPP = ROOT / "work" / "figures_supp_tiff"
SUPP_XLSX = ROOT / "投稿材料包" / "02_附件与清单" / "Supplementary_Tables_project.xlsx"
SUPP_XLS = ROOT / "work" / "attachments" / "Supplementary_Tables_project.xls"
MANIFEST = ROOT / "work" / "output" / "submission_format_manifest.csv"

# 主图顺序必须与 33_build_journal_packages.py 的投稿槽位一致
# （2026-10-27 增补 Figure 5：疾病锚层全景，由 51_fig5_anchors.py 生成）
MAIN_ORDER = ["figure_rg_heatmap.png", "figure_rg_forest.png",
              "figure_locus_class.png", "figure_coloc_regions.png",
              "figure_anchor_landscape.png"]
SUPP_ORDER = ["figS1_mhc_subregion.png", "figS2_prs_transfer.png",
              "figS3_autoimmune_anchor.png"]
DPI = 300


def sha256(path: Path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def flatten_to_rgb(im: Image.Image) -> Image.Image:
    """把带 alpha 的 PNG 压到白底 RGB；否则直接 RGB。"""
    if im.mode in ("RGBA", "LA") or (im.mode == "P" and "transparency" in im.info):
        im = im.convert("RGBA")
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        bg.alpha_composite(im)
        return bg.convert("RGB")
    return im.convert("RGB")


def convert_one(src: Path, dst: Path) -> dict:
    with Image.open(src) as im:
        src_size = im.size
        src_mode = im.mode
        rgb = flatten_to_rgb(im)
    cmyk = rgb.convert("CMYK")
    dst.parent.mkdir(parents=True, exist_ok=True)
    cmyk.save(dst, format="TIFF", compression="tiff_lzw", dpi=(DPI, DPI))

    # 回读校验：色彩空间/尺寸/dpi/压缩
    with Image.open(dst) as chk:
        mode = chk.mode
        size = chk.size
        dpi = tuple(round(v, 1) for v in (chk.info.get("dpi") or (0, 0)))
        comp = chk.tag_v2.get(259)          # TIFF Compression tag: 5 = LZW
    if mode != "CMYK":
        raise RuntimeError("[FATAL] %s 回读色彩空间为 %s，非 CMYK" % (dst.name, mode))
    if size != src_size:
        raise RuntimeError("[FATAL] %s 尺寸变化 %s -> %s" % (dst.name, src_size, size))
    if abs(dpi[0] - DPI) > 0.5 or abs(dpi[1] - DPI) > 0.5:
        raise RuntimeError("[FATAL] %s dpi 回读 %s，非 %d" % (dst.name, dpi, DPI))
    if comp != 5:
        raise RuntimeError("[FATAL] %s 压缩标记 %s，非 LZW(5)" % (dst.name, comp))
    return {
        "kind": "figure_tiff", "file": str(dst.relative_to(ROOT)), "source": src.name,
        "width": size[0], "height": size[1], "source_mode": src_mode,
        "colorspace": mode, "dpi": "%dx%d" % dpi, "compression": "LZW",
        "bytes": dst.stat().st_size, "sha256": sha256(dst),
    }


def export_xls() -> dict:
    import openpyxl
    import xlrd
    import xlwt

    if not SUPP_XLSX.exists():
        raise RuntimeError("[FATAL] 缺补充表源文件：%s（先跑 32_supp_xlsx.py）" % SUPP_XLSX)
    wb_in = openpyxl.load_workbook(SUPP_XLSX, data_only=True)
    wb_out = xlwt.Workbook(encoding="utf-8")
    bold = xlwt.easyxf("font: bold on")

    cell_count = 0
    for ws in wb_in.worksheets:
        sh = wb_out.add_sheet(ws.title[:31])
        widths = {}
        for r, row in enumerate(ws.iter_rows(values_only=True), start=0):
            for c, val in enumerate(row):
                if val is None:
                    continue
                if isinstance(val, bool):
                    sh.write(r, c, int(val))
                elif isinstance(val, (int, float)):
                    sh.write(r, c, float(val))
                else:
                    s = str(val)
                    if len(s) > 32767:
                        raise RuntimeError("[FATAL] %s!R%dC%d 超过 xls 单元格上限" % (ws.title, r + 1, c + 1))
                    sh.write(r, c, s)
                cell_count += 1
                widths[c] = min(60, max(widths.get(c, 8), len(str(val)) + 2))
        for c, w in widths.items():
            sh.col(c).width = w * 256

    SUPP_XLS.parent.mkdir(parents=True, exist_ok=True)
    wb_out.save(str(SUPP_XLS))

    # 逐 cell 回读比对（sheet 名、维度、值）
    xls = xlrd.open_workbook(str(SUPP_XLS), formatting_info=False)
    in_names = wb_in.sheetnames
    out_names = xls.sheet_names()
    if in_names != out_names:
        raise RuntimeError("[FATAL] xls sheet 名不一致：%s vs %s" % (in_names, out_names))
    mismatches = []
    for name in in_names:
        a = wb_in[name]
        b = xls.sheet_by_name(name)
        if a.max_row != b.nrows or a.max_column != b.ncols:
            mismatches.append("%s 维度 %dx%d vs %dx%d" %
                              (name, a.max_row, a.max_column, b.nrows, b.ncols))
            continue
        for r, row in enumerate(a.iter_rows(values_only=True)):
            for c, v in enumerate(row):
                got = b.cell_value(r, c)
                if v is None:
                    if got not in ("", None):
                        mismatches.append("%s!R%dC%d None vs %r" % (name, r + 1, c + 1, got))
                elif isinstance(v, bool):
                    if float(got) != float(int(v)):
                        mismatches.append("%s!R%dC%d bool %r vs %r" % (name, r + 1, c + 1, v, got))
                elif isinstance(v, (int, float)):
                    try:
                        if abs(float(got) - float(v)) > 1e-9 * max(1.0, abs(float(v))):
                            mismatches.append("%s!R%dC%d num %r vs %r" % (name, r + 1, c + 1, v, got))
                    except Exception:
                        mismatches.append("%s!R%dC%d num %r vs non-number %r" % (name, r + 1, c + 1, v, got))
                else:
                    if str(got) != str(v):
                        mismatches.append("%s!R%dC%d str %r vs %r" % (name, r + 1, c + 1, v, got))
    if mismatches:
        raise RuntimeError("[FATAL] xls 回读不一致 %d 处（前 5）：%s" %
                           (len(mismatches), " | ".join(mismatches[:5])))
    return {
        "kind": "supplementary_xls", "file": str(SUPP_XLS.relative_to(ROOT)),
        "source": SUPP_XLSX.name, "sheets": len(in_names), "cells": cell_count,
        "bytes": SUPP_XLS.stat().st_size, "sha256": sha256(SUPP_XLS),
    }


def main():
    OUT_MAIN.mkdir(parents=True, exist_ok=True)
    OUT_SUPP.mkdir(parents=True, exist_ok=True)
    rows = []

    if len(MAIN_ORDER) != 5:
        raise RuntimeError("[FATAL] 主图槽位应 5 个（Figure 1–4 + Figure 5 疾病锚层全景）")
    for i, name in enumerate(MAIN_ORDER, start=1):
        src = FIG / name
        if not src.exists():
            raise RuntimeError("[FATAL] 缺主图：%s" % src)
        rows.append(convert_one(src, OUT_MAIN / ("Figure%d.tif" % i)))

    for i, name in enumerate(SUPP_ORDER, start=1):
        src = FIGSUPP / name
        if not src.exists():
            raise RuntimeError("[FATAL] 缺补充图：%s（先跑 44_enhancement_figures.py）" % src)
        rows.append(convert_one(src, OUT_SUPP / ("Figure_S%d.tif" % i)))

    rows.append(export_xls())

    with open(MANIFEST, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["kind", "file", "source", "width", "height",
                                           "source_mode", "colorspace", "dpi",
                                           "compression", "sheets", "cells", "bytes", "sha256"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print("[fmt] TIFF/CMYK 主图 %d + 补充图 %d（300 dpi, LZW，无重采样）" % (len(MAIN_ORDER), len(SUPP_ORDER)))
    print("[fmt] 补充表 .xls：%d sheets / %d cells，逐 cell 回读一致"
          % (rows[-1]["sheets"], rows[-1]["cells"]))
    print("[fmt] manifest -> %s" % MANIFEST.relative_to(ROOT))
    print("SUBMISSION_FORMAT_ASSETS_OK")


if __name__ == "__main__":
    main()
