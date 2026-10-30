# -*- coding: utf-8 -*-
"""project constants.py — Python 侧口径单一来源。
跨病原感染易感性共享遗传架构 atlas：LDSC 遗传相关 + 共定位 + 广谱/特异位点分类。
"""
from pathlib import Path

COURSE = "project"
ROOT = Path("/path/to/project")
WORK = ROOT / "work"
OUT = WORK / "output"
FIG = WORK / "figures"
DATA = WORK / "data"
CACHE = WORK / "cache"
for _d in (OUT, FIG, DATA, CACHE):
    _d.mkdir(parents=True, exist_ok=True)

# ---- 共享数据仓 ----
GWAS_DATA = Path("/path/to/project/gwas_data")
FINNGEN = GWAS_DATA / "finngen_r12"
FINNGEN_URL_BASE = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
THOUSANDG = GWAS_DATA / "1000g_eur"
PLINK = GWAS_DATA / "tools" / "plink.exe"

ANCESTRY = "EUR"

# ---- 跨病原感染终点面板（与 course_config.R FINNGEN_ENDPOINTS 一致；立项卡 §三 定案 2026-10-11）----
# case 数只从官方 manifest（work/data/finngen_R12_pheno_n.tsv）读取，禁止手抄。
FINNGEN_ENDPOINTS = {
    # 细菌轴（10）
    "INTESTINAL": "AB1_INTESTINAL_INFECTIONS",
    "GASTRO_NOS": "AB1_GASTROENTERITIS_NOS",
    "PYELONEPHR": "N14_PYELONEPHR",
    "ERYSIPELAS": "AB1_ERYSIPELAS",
    "BACT_NOS": "AB1_BACTINF_NOS",
    "PNEUMOBACT": "J10_PNEUMOBACT",
    "SEPSIS": "AB1_OTHER_SEPSIS",
    "OSTEOMYEL": "M13_OSTEOMYELITIS",
    "MENINGBACT": "G6_MENINGBACT",
    "ENDOCARD": "I9_ENDOCARD",
    # 混合轴（1）
    "PNEUMONIA_MIXED": "J10_PNEUMONIA",
    # 病毒轴（9）
    "COVID19": "U22_COVID19_CONFIRMED",
    "INFLUENZA": "J10_INFLUENZA",
    "ZOSTER": "AB1_ZOSTER",
    "VIRALPNEUMO": "J10_VIRALPNEUMO",
    "HERPES_SIMPLEX": "AB1_HERPES_SIMPLEX",
    "VIRAL_WARTS": "AB1_VIRAL_WARTS",
    "EBV": "AB1_EBV",
    "VIRAL_CNS": "AB1_VIRAL_CNS",
    "VARICELLA": "AB1_VARICELLA",
    # 真菌轴（3）
    "MYCOSES": "AB1_MYCOSES",
    "CANDIDIASIS": "AB1_CANDIDIASIS",
    "DERMATOPHYTOSIS": "AB1_DERMATOPHYTOSIS",
    # 寄生虫轴（2）
    "INFESTATIONS": "AB1_PEDICULOSIS_ACARIASIS_OTHERINFEST",
    "PROTOZOAL": "AB1_PROTOZOAL",
}
# Mtb 系本地有文件但入主面板禁用（方向边界：Mtb 归 06/06A/TB 系）
MTB_ENDPOINTS_BLOCKED = {"AB1_TUBERCULOSIS", "AB1_RESP_TUBERCU_CONF",
                         "AB1_RESP_TUBERCU_NONCONF", "AB1TUBERCU_MILIARY",
                         "AB1_TUBERCU_OTH_ORG", "AB1_SEQULAE_TUBERCU", "TBC_RESP"}

# ---- 自免锚层终点（E3 增强层，2026-10-25；与感染面板同库同版本 FinnGen R12）----
# 单一真源：43_autoimmune_runner.py（跑批）与 25_canonical_numbers.py（munge 分类守卫）
# 都从这里取，禁止两处各写一份列表。
AUTOIMMUNE_ENDPOINTS = ["T1D", "M13_RHEUMA", "G6_MS", "K11_COELIAC", "K11_UC_STRICT2",
                        "K11_CD_STRICT2", "L12_PSORIASIS", "M13_ANKYLOSPON",
                        "M13_SJOGREN", "E4_GRAVES_STRICT"]

# ---- 方法常量 ----
LDSC_R2_WIN_CM = 1.0
MAF_MIN = 0.01
GENOME_P = 5e-8
GW_P = 1e-6
NCHISQ_FLOOR = 30  # LDSC 功效门：mean chi2 太低(χ²均值)的 trait 剔除（设计文档风险①）

# 日期锚（所有落盘文件名带日期）
DATE = "2026-10-11"
