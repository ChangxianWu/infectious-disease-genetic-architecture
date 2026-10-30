# -*- coding: utf-8 -*-
""" 14_clump_leads.py — 显著配对候选 lead SNP 提取（plink 1.9 clumping）
（Stage 2b，2026-10-11）

口径：p1 = 5e-8（全基因组显著）；r2 = 0.01；窗口 500 kb；参考面板 = 官方 LDSC
      1000G EUR QC bfiles（与主矩阵/LD score 同一参考，避免参考体系混用）。
实现：.sumstats 只存 CHISQ → 还原 p = chi2.sf(CHISQ, 1)；逐染色体 clump 后合并去重。
安全：trait 名一律经白名单正则校验（^[A-Za-z0-9_]+$），subprocess 全程列表参数 +
      shell=False，不拼接任何 shell 字符串。
输出：work/output/clump/<trait>.leads + leads_annot.tsv + leads_summary.csv
用法：python 14_clump_leads.py [TRAIT1 TRAIT2 ...]
"""
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

sys.path.insert(0, "/path/to/project/work/code")
from constants import DATA, OUT, GWAS_DATA, FINNGEN_ENDPOINTS  # noqa: E402

REF = GWAS_DATA / "ldsc_ref" / "1000G_EUR_Phase3_plink"
PLINK = GWAS_DATA / "tools" / "plink.exe"
MUNGED = DATA / "munged"
WORK = DATA.parent / "cache" / "clump"
OUTD = OUT / "clump"
for _d in (WORK, OUTD):
    _d.mkdir(parents=True, exist_ok=True)

SAFE_NAME = re.compile(r"^[A-Za-z0-9_]+$")

DEFAULT_TRAITS = ["INTESTINAL", "GASTRO_NOS", "PYELONEPHR", "ERYSIPELAS", "BACT_NOS",
                  "PNEUMOBACT", "SEPSIS", "PNEUMONIA_MIXED", "INFLUENZA", "OSTEOMYEL"]


def to_plink_input(sumstats: Path, dst: Path, p_thresh=5e-8):
    """.sumstats → SNP A1 A2 P（仅显著 SNP，减少 plink 负担）。"""
    df = pd.read_csv(sumstats, sep="\t", usecols=["SNP", "A1", "A2", "CHISQ"])
    df["P"] = chi2.sf(df["CHISQ"].values, 1)
    df = df[df["P"] < p_thresh]
    df[["SNP", "A1", "A2", "P"]].to_csv(dst, sep="\t", index=False)
    return len(df)


def clump_trait(label: str, pheno: str):
    if not SAFE_NAME.match(label) or not SAFE_NAME.match(pheno):
        print("[clump] 非法名称，跳过：%r / %r" % (label, pheno))
        return None
    ss = MUNGED / ("%s.sumstats" % pheno)
    if not ss.exists():
        print("[clump] %s 无 sumstats，跳过" % label)
        return None
    pin = WORK / ("%s.forclump" % label)
    n_sig = to_plink_input(ss, pin)
    if n_sig == 0:
        print("[clump] %s 无全基因组显著 SNP" % label)
        (OUTD / ("%s.leads" % label)).write_text("", encoding="utf-8")
        return {"trait": label, "pheno": pheno, "n_sig_snps": 0, "n_leads": 0, "min_p": np.nan}

    leads, rows = [], []
    for c in range(1, 23):
        bfile = REF / ("1000G.EUR.QC.%d" % c)
        pref = WORK / ("%s.chr%d" % (label, c))
        args = [str(PLINK), "--bfile", str(bfile),
                "--clump", str(pin), "--clump-p1", "5e-8", "--clump-r2", "0.01",
                "--clump-kb", "500", "--clump-field", "P", "--clump-snp-field", "SNP",
                "--out", str(pref)]
        subprocess.run(args, capture_output=True, text=True, shell=False, timeout=900)
        cf = Path(str(pref) + ".clumped")
        if cf.exists():
            d = pd.read_csv(cf, sep=r"\s+")
            if len(d):
                d = d.rename(columns={d.columns[2]: "SNP"})
                for _, r in d.iterrows():
                    leads.append(r["SNP"])
                    rows.append({"trait": label, "CHR": r.get("CHR", c), "SNP": r["SNP"],
                                 "BP": r.get("BP", np.nan), "P": r.get("P", np.nan)})
    leads = sorted(set(s for s in leads if isinstance(s, str) and s.startswith("rs")))
    (OUTD / ("%s.leads" % label)).write_text("\n".join(leads) + ("\n" if leads else ""),
                                             encoding="utf-8")
    if rows:
        pd.DataFrame(rows).to_csv(OUTD / ("%s.leads_annot.tsv" % label), sep="\t", index=False)
    summ = {"trait": label, "pheno": pheno, "n_sig_snps": n_sig, "n_leads": len(leads),
            "min_p": float(pd.read_csv(pin, sep="\t")["P"].min())}
    print("[clump] %s: 显著 SNP %d → lead %d" % (label, n_sig, len(leads)))
    return summ


def main():
    traits = sys.argv[1:] or DEFAULT_TRAITS
    out = []
    for t in traits:
        ph = FINNGEN_ENDPOINTS.get(t, t)
        s = clump_trait(t, ph)
        if s:
            out.append(s)
    if out:
        df = pd.DataFrame(out)
        df.to_csv(OUTD / "leads_summary.csv", index=False)
        print()
        print(df.to_string(index=False))


if __name__ == "__main__":
    main()
