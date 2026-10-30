移植报告（bulik/ldsc py2 → py3，2026-10-11）
parse.py               1 处变换
sumstats.py            4 处变换
regressions.py         1 处变换
irwls.py               2 处变换
jackknife.py           4 处变换
ldsc.py               12 处变换
munge_sumstats.py      2 处变换
- numpy>=2 兼容修复（不改算法）：regressions.py 三处 1x1 矩阵 float() 改 .item()（2026-10-11）
