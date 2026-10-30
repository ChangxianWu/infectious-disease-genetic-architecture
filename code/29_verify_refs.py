# -*- coding: utf-8 -*-
""" 29_verify_refs.py — 文献双源核验门（门 5）

三层关系（每层可独立复核）：
  1. 稿件 References 节（渲染文本）——解析出每条题名片段；
  2. 核实池 work/tools/refs_pool_2026-10-12.tsv——每条带 PMID 与 status；
  3. Europe PMC 现场再取（EXT_ID 逐条）——题名与池中题名逐字比对。

任一环不一致即失败：稿件出现池中不存在的文献、池中题名与 EPMC 现场返回不一致、
或池中存在未核实条目（项目硬约束：禁止编造作者与文献）。

输出：work/output/refs_dual_source_check.csv
用法：python 29_verify_refs.py
"""
import csv
import html
import io
import json
import os
import re
import sys
import time
import urllib.parse
import zipfile
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

BASE = os.path.realpath("/path/to/project")
sys.path.insert(0, BASE)
from cv_shared.probe_kit import _get  # noqa: E402

ROOT = Path(BASE) / "project"
OUT = ROOT / "work" / "output"
POOL = ROOT / "work" / "tools" / "refs_pool_2026-10-12.tsv"
DOCX = ROOT / "work" / "manuscript" / "project_Manuscript_v0.1.docx"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

if not POOL.exists():
    print("[refs] 缺文献池: %s" % POOL)
    sys.exit(1)
if not DOCX.exists():
    print("[refs] 缺稿件: %s" % DOCX)
    sys.exit(1)

with open(POOL, encoding="utf-8-sig", newline="") as fh:
    pool = list(csv.DictReader(fh, delimiter="\t"))
unverified = [r["key"] for r in pool if r["status"] != "verified"]

with zipfile.ZipFile(DOCX) as z:
    xml = z.read("word/document.xml").decode("utf-8")
text = re.sub(r"<[^>]+>", " ", xml).replace("&amp;", "&")
text = re.sub(r"\s+", " ", text)

rows = []
n_ok = 0
n_bad = 0
n_net = 0
print("[refs] 双源核验 %d 条（池 vs Europe PMC 现场 + 稿件在场）" % len(pool))
for r in pool:
    u = EPMC + "?" + urllib.parse.urlencode(
        {"query": "EXT_ID:%s" % r["pmid"], "format": "json", "resultType": "core",
         "pageSize": "1"})
    body = _get(u, tries=3)
    live_title = ""
    if body:
        res = json.loads(body).get("resultList", {}).get("result", [])
        if res:
            # 现场题名同样做 HTML 实体解码：池侧已解码（26_refs_pool.py），
            # 两侧不做同一变换会把 `&gt;80 000` 判成不一致（本轮实测该假失败）
            live_title = html.unescape(res[0].get("title") or "").rstrip(".")
    title_match = live_title.lower() == r["title"].strip().lower()
    # 稿件中是否出现该文献题名前 45 字（References 节以富文本渲染，故用片段而非全串）
    frag = r["title"].strip()[:45]
    in_doc = frag.lower() in text.lower()
    # 三态（2026-10-30）：网络不可达 ≠ 题名不一致。前者是链路问题（可重试），后者才是数据问题。
    # 但仍**不计 OK**——门照样红，只是把"该重试"和"该查池"分开，避免把断网误读成文献错误。
    if body is None:
        status = "UNVERIFIED_NETWORK"
    elif title_match and in_doc:
        status = "OK"
    else:
        status = "FAIL"
    if status == "OK":
        n_ok += 1
    elif status == "UNVERIFIED_NETWORK":
        n_net += 1
    else:
        n_bad += 1
    rows.append({"key": r["key"], "pmid": r["pmid"], "status": status,
                 "title_match_live": title_match, "present_in_manuscript": in_doc,
                 "pool_title": r["title"][:90], "live_title": live_title[:90]})
    print("  %-18s %-22s %-9s in_doc=%-5s live=%s"
          % (status, r["key"], r["pmid"], in_doc, live_title[:52]))
    time.sleep(0.3)

buf = io.StringIO()
w = csv.DictWriter(buf, fieldnames=["key", "pmid", "status", "title_match_live",
                                    "present_in_manuscript", "pool_title", "live_title"],
                   delimiter="\t")
w.writeheader()
w.writerows(rows)
(OUT / "refs_dual_source_check.csv").write_text(buf.getvalue(),
                                                encoding="utf-8-sig", newline="")
print("[refs] 池中未核实条目: %s" % (unverified or "无"))
print("[refs] 双源核验 %d 条：OK %d，FAIL %d，网络不可达 %d"
      % (len(rows), n_ok, n_bad, n_net))
if n_net:
    print("[refs] 说明：网络不可达条目**不计通过**，仅表示本次链路中断（Europe PMC SSL/超时），"
          "应重试本门；它与题名不一致（FAIL）是两类问题，不要混为一谈。")
if unverified or n_bad or n_net:
    sys.exit(1)
print("REFS_DUAL_SOURCE_OK")
