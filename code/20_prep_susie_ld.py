# -*- coding: utf-8 -*-
""" 20_prep_susie_ld.py — 为 coloc.susie 生成 LD 缓存（plink r-square 矩阵）

（Stage 3b，2026-10-11）

口径：LD 一律来自官方 LDSC 1000G EUR QC bfiles（与主矩阵、LD score、clumping 同源参考）；
      保留顺序先由 plink --write-snplist 落盘（该文件保持 bim 行序），再按同一顺序取 LD，
      避免"顺序错位静默产生错误 LD"（本文件 .ldsnp 即 SNP 行序真源）。
      plink 输出 .ld 无 SNP 名列，行序必须与 .ldsnp 严格一致。

输入安全：LOCI 为脚本内硬编码常量，逐个经 ^rs[0-9]+$ 白名单正则校验；
          染色体来自本地注释表（int() 强制转换）；所有 plink 调用均为
          字面量参数列表 + shell=False（不拼接任何字符串）。

用法：python 20_prep_susie_ld.py
输出：work/cache/susie_ld/<LOCUS>.ld + <LOCUS>.ldsnp
"""
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATA, OUT, GWAS_DATA  # noqa: E402

REF = GWAS_DATA / "ldsc_ref" / "1000G_EUR_Phase3_plink"
PLINK = GWAS_DATA / "tools" / "plink.exe"
CDIR = DATA / "coloc"
LDIR = DATA.parent / "cache" / "susie_ld"
LDIR.mkdir(parents=True, exist_ok=True)
SAFE = re.compile(r"^rs[0-9]+$")
# 需 susie 的位点：非 MHC 跨轴区（MHC 区 LD 极强、susie 独立信号假设不成立，单列）
LOCI = ("rs3116651", "rs876537", "rs492602", "rs516246", "rs495828")


def main():
    ann = pd.read_csv(OUT / "lead_gene_annotation.tsv", sep="\t")
    for loc in LOCI:
        if not SAFE.match(loc):
            print("[ld] 位点名未过白名单 %r，跳过" % loc)
            continue
        row = ann[ann["SNP"] == loc]
        c = int(row.iloc[0]["CHR_hg38"])
        win = next(CDIR.glob(loc + "_*.tsv"))
        d = pd.read_csv(win, sep="\t")
        snplist = LDIR.joinpath(loc + ".snps")
        snplist.write_text("\n".join(d["SNP"].astype(str)) + "\n", encoding="utf-8")

        bfile = str(REF / ("1000G.EUR.QC.%d" % c))
        opref = str(LDIR.joinpath(loc + "_order"))
        subprocess.run([str(PLINK), "--bfile", bfile, "--extract", str(snplist),
                        "--write-snplist", "--out", opref],
                       capture_output=True, text=True, shell=False, timeout=900)
        ofile = Path(opref + ".snplist")
        if not ofile.exists():
            print("[ld] %s 取顺序失败" % loc)
            continue
        order = [s for s in ofile.read_text().split() if s]

        pref = str(LDIR.joinpath(loc + "_plink"))
        p = subprocess.run([str(PLINK), "--bfile", bfile, "--extract", str(ofile),
                            "-r", "square", "--out", pref],
                           capture_output=True, text=True, shell=False, timeout=1800)
        ldfile = Path(pref + ".ld")
        if not ldfile.exists():
            print("[ld] %s 失败：%s" % (loc, (p.stderr or p.stdout)[-200:]))
            continue
        LDIR.joinpath(loc + ".ldsnp").write_text("\n".join(order) + "\n", encoding="utf-8")
        LDIR.joinpath(loc + ".ld").write_bytes(ldfile.read_bytes())

        n = len(order)
        ldm = pd.read_csv(ldfile, sep=r"\s+", header=None)
        shape_ok = ldm.shape == (n, n)
        diag_ok = bool((abs(ldm.values.diagonal() - 1) < 1e-3).all())
        print("[ld] %s chr%d: %d SNP | 形状 %s | 尺寸%s | 对角=1: %s"
              % (loc, c, n, ldm.shape, "一致" if shape_ok else "**不一致**", diag_ok))


if __name__ == "__main__":
    main()
