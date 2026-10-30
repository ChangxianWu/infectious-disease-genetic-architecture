# -*- coding: utf-8 -*-
"""移植 bulik/ldsc 官方 py2 源码 → py3（project Stage 1，规则 96：忠实移植不改算法）。
机械变换清单（只做语法级，不动任何数值逻辑）：
  1. xrange( → range(
  2. .iteritems() → .items()、.itervalues() → .values()、.iterkeys() → .keys()
  3. print >> fh, x / print 'x' 语句 → py3 写法（ldsc.py / munge_sumstats.py 内）
  4. izip → zip（如出现）
  5. basestring → str
不动：任何回归/权重/jackknife 数学；np.matrix 仍可用（numpy 2.x 保留）。
源：work/cache/ldsc_src/（bulik/ldsc master，2026-10-11 经 net_guard 取回）。
"""
import re
from pathlib import Path

SRC = Path("/path/to/project/work/cache/ldsc_src")
DST = Path("/path/to/project/work/code/ldsc_py3")
DST.mkdir(parents=True, exist_ok=True)

FILES = ["parse.py", "sumstats.py", "regressions.py", "irwls.py", "jackknife.py",
         "ldsc.py", "munge_sumstats.py"]

def py2_to_py3(text, fname):
    n_changes = 0

    def sub_count(pat, rep, t):
        nonlocal n_changes
        new, k = re.subn(pat, rep, t)
        n_changes += k
        return new

    t = text
    t = sub_count(r"\bxrange\(", "range(", t)
    t = sub_count(r"\.iteritems\(\)", ".items()", t)
    t = sub_count(r"\.itervalues\(\)", ".values()", t)
    t = sub_count(r"\.iterkeys\(\)", ".keys()", t)
    t = sub_count(r"\bizip\(", "zip(", t)
    t = sub_count(r"\bbasestring\b", "str", t)

    # print 语句（仅这几处模式，逐个替换）
    #   print >>fout, x  →  print(x, file=fout)
    def print_binfile(m):
        return "print(%s, file=%s)" % (m.group(2).strip(), m.group(1).strip())
    t, k = re.subn(r"print >>\s*([A-Za-z_][\w.\[\]]*),\s*(.+)", print_binfile, t)
    n_changes += k
    #   print 'x', y  /  print "x"  （无括号 print；保守：仅当行首是 print 后跟非括号）
    lines = t.split("\n")
    for i, ln in enumerate(lines):
        m = re.match(r"^(\s*)print\s+(?![(>])(.+)$", ln)
        if m:
            indent, rest = m.group(1), m.group(2)
            # 处理尾随逗号
            rest = rest.rstrip()
            if rest.endswith(","):
                rest = rest[:-1]
                lines[i] = "%sprint(%s, end=' ')" % (indent, rest)
            else:
                lines[i] = "%sprint(%s)" % (indent, rest)
            n_changes += 1
    t = "\n".join(lines)
    return t, n_changes

report = []
for f in FILES:
    src = SRC / f
    if not src.exists():
        report.append("%-20s SKIP（源缺失）" % f)
        continue
    t, n = py2_to_py3(src.read_text(encoding="utf-8"), f)
    (DST / f).write_text(t, encoding="utf-8", newline="")
    report.append("%-20s %3d 处变换" % (f, n))

buf = "移植报告（bulik/ldsc py2 → py3，2026-10-11）\n" + "\n".join(report) + "\n"
(DST / "_PORT_NOTES.md").write_text(buf, encoding="utf-8", newline="")
print(buf)
