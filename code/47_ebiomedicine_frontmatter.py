# -*- coding: utf-8 -*-
""" 47_ebiomedicine_frontmatter.py — eBioMedicine 前置件（Research in context + Highlights）

依据：eBioMedicine 官方 Information for Authors 要求 Original Research 提供
  * Research in context（Evidence before this study / Added value / Implications）
  * 2–6 条 Highlights


纪律：
  * 所有数字从 work/output/s7_canonical_numbers.tsv 取，禁手抄（规则 11/33/51/55）；
  * 只写已在稿件/补充表中落地的结论；E3 自免锚层用**全链完成值**（100/91/29、中位 0.224），
    不再用增强方案里的 "pending full E3" 草案措辞；
  * 写出后回读：无残留 token，canonical 关键数字全部在场；
  * docx 的 zip 时间戳归一化到 1980-01-01，保证重复构建字节一致（本文件不进门 6 哈希，
    但保持确定性便于归档核对）。

用法：python 47_ebiomedicine_frontmatter.py
产出：work/attachments/eBioMedicine_Research_in_context_and_Highlights.docx
      work/attachments/eBioMedicine_Research_in_context_and_Highlights.md（同一内容的可读源）
"""
import hashlib
import io
import os
import sys
import zipfile
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from docx import Document
from docx.shared import Pt

ROOT = Path("/path/to/project") / "project"
CANON = ROOT / "work" / "output" / "s7_canonical_numbers.tsv"
OUT_DIR = ROOT / "work" / "attachments"
DOCX = OUT_DIR / "eBioMedicine_Research_in_context_and_Highlights.docx"
MD = OUT_DIR / "eBioMedicine_Research_in_context_and_Highlights.md"
FIXED_DATE = (1980, 1, 1, 0, 0, 0)


def load_canon():
    c = {}
    with open(CANON, encoding="utf-8-sig", newline="") as fh:
        for line in fh:
            f = line.rstrip("\r\n").split("\t")
            if len(f) >= 2 and f[0].strip():
                c[f[0].strip()] = f[1].strip()
    return c


def f2(c, k):
    return "%.2f" % float(c[k])


def f3(c, k):
    return "%.3f" % float(c[k])


def build_text(c):
    # eBioMedicine 的 Research in context 是"2 句既往 / 1 句增量 / 1 句意义"的固定短面板
    # （2026-10-25 官网核实），不是正文摘要；因此这里只保留可核查的框架性数字。
    evidence = (
        "Host genetic correlations among infection endpoints, and between infections and autoimmune "
        "diseases, have been reported piecemeal. No harmonised cross-pathogen atlas had mapped which loci "
        "act across pathogen types, or tested whether that architecture reflects monogenic "
        "immunodeficiency, major-histocompatibility-complex (MHC) subregions, or autoimmune genetics."
    )
    added = (
        "Across %s infection-endpoint pairs (%s passing the power gate and %s surviving false-discovery "
        "correction), sharing concentrated on mucosal-barrier and acute-phase loci and on HLA class II "
        "(%s of %s MHC-proximal leads), was disjoint from monogenic immunodeficiency (%s of %s lead genes), "
        "and overlapped autoimmune disease in a disease-specific pattern (%s of %s pairs powered; %s "
        "positive; %s false-discovery-significant) that persisted after retaining the MHC (median "
        "delta-rg change %s)."
    ) % (
        c["n_pairs_attempted"], c["n_pairs_power_ok"], c["n_pairs_bh05"],
        c["mhc_class_ii_n"], c["mhc_leads_n"],
        c["iei_lead_overlap"], c["iei_lead_n"],
        c["anchor_n_power_ok"], c["anchor_n_pairs_attempted"],
        c["anchor_n_positive"], c["anchor_n_bh"], f3(c, "anchor_sens_median_delta"),
    )
    implications = (
        "Cross-pathogen susceptibility is therefore a barrier and antigen-presentation phenomenon rather "
        "than a type I interferon phenomenon, with polygenic transfer below 1%% of observed-scale "
        "variance, so the atlas informs pathway biology and study design rather than risk prediction."
    ) % ()
    highlights = [
        "Across %s infection-endpoint pairs, %s survived FDR and all were positive; the "
        "bacterial-versus-viral boundary does not predict sharing."
        % (c["n_pairs_attempted"], c["n_pairs_bh05"]),
        "Shared loci concentrate on mucosal-barrier, acute-phase and HLA class II biology, not type I "
        "interferon.",
        "The shared architecture is disjoint from monogenic immunodeficiency genes (%s/%s lead genes)."
        % (c["iei_lead_overlap"], c["iei_lead_n"]),
        "Autoimmune overlap is broad but disease-specific (%s FDR-significant pairs; null for type 1 "
        "diabetes, coeliac disease and multiple sclerosis)." % c["anchor_n_bh"],
        "Cross-source replication separates FUT2 and ABO (PP.H4 %s and %s) from CRP (%s)."
        % (f3(c, "coloc_fut2_cross_max_h4"), f3(c, "coloc_abo_cross_max_h4"),
           f3(c, "coloc_crp_cross_max_h4")),
    ]
    return evidence, added, implications, highlights


def normalize_zip(path: Path):
    """把 docx 内 zip 条目时间戳归一化，保证重复构建字节一致。"""
    with zipfile.ZipFile(path) as zin:
        items = [(i.filename, zin.read(i.filename)) for i in zin.infolist()]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zout:
        for name, data in sorted(items, key=lambda x: x[0]):
            info = zipfile.ZipInfo(name, date_time=FIXED_DATE)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o600 << 16
            zout.writestr(info, data)
    path.write_bytes(buf.getvalue())


def main():
    c = load_canon()
    for k in ["n_endpoints_munged", "n_pairs_attempted", "n_pairs_power_ok", "n_pairs_bh05",
              "anchor_n_power_ok", "anchor_n_positive", "anchor_n_bh", "iei_lead_overlap",
              "mhc_class_ii_n", "mhc_leads_n"]:
        if k not in c:
            raise RuntimeError("[FATAL] canonical 缺键：%s" % k)

    evidence, added, implications, highlights = build_text(c)

    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(11)

    doc.add_heading("Research in context", level=1)
    doc.add_paragraph("Shared host genetic architecture of susceptibility to multiple pathogens: "
                      "a cross-pathogen atlas of genetic correlation, pleiotropy and colocalisation")
    doc.add_heading("Evidence before this study", level=2)
    doc.add_paragraph(evidence)
    doc.add_heading("Added value", level=2)
    doc.add_paragraph(added)
    doc.add_heading("Implications of the available evidence", level=2)
    doc.add_paragraph(implications)
    doc.add_heading("Highlights", level=1)
    for h in highlights:
        doc.add_paragraph(h, style="List Bullet")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    doc.save(str(DOCX))
    normalize_zip(DOCX)

    md_lines = ["# Research in context", "",
                "**Shared host genetic architecture of susceptibility to multiple pathogens: "
                "a cross-pathogen atlas of genetic correlation, pleiotropy and colocalisation**", "",
                "## Evidence before this study", "", evidence, "",
                "## Added value", "", added, "",
                "## Implications of the available evidence", "", implications, "",
                "## Highlights", ""] + ["- " + h for h in highlights] + [""]
    MD.write_text("\n".join(md_lines), encoding="utf-8", newline="")

    # 回读校验：无残留 token；关键数字在场
    with zipfile.ZipFile(DOCX) as z:
        xml = z.read("word/document.xml").decode("utf-8", "replace")
    for bad in ["{{", "}}", "TODO", "TBD", "undefined", "NaN", "pending full E3"]:
        if bad in xml:
            raise RuntimeError("[FATAL] 回读发现残留 token：%s" % bad)
    prose = " ".join([evidence, added, implications] + highlights)
    must = [c["n_pairs_attempted"], c["n_pairs_power_ok"], c["n_pairs_bh05"],
            c["anchor_n_power_ok"], c["anchor_n_positive"],
            c["anchor_n_bh"], c["iei_lead_overlap"], c["mhc_leads_n"],
            c["mhc_class_ii_n"], f3(c, "anchor_sens_median_delta")]
    missing = [m for m in must if m not in prose]
    if missing:
        raise RuntimeError("[FATAL] 前置件缺关键数字：%s" % ", ".join(missing))
    words = len(prose.split())

    print("[ebm] Research in context + %d Highlights -> %s" % (len(highlights), DOCX.name))
    print("[ebm] 正文 %d 词；关键 canonical 数字 %d 项全部在场；回读无残留 token" % (words, len(must)))
    print("[ebm] sha256 = %s" % hashlib.sha256(DOCX.read_bytes()).hexdigest()[:16])
    print("EBM_FRONTMATTER_OK")


if __name__ == "__main__":
    main()
