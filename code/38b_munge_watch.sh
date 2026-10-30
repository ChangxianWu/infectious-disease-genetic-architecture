#!/usr/bin/env bash
# ============================================================
#  38b_munge_watch.sh — 自免终点落盘即 munge 的看门循环
# ------------------------------------------------------------
# 与 38 号下载并行：每个 finngen_R12_<ph>.gz 完整落地（无 .p 分片残留）即调用
# 10_munge_finngen.py CLI munge；43 号的 munge_autoimmune 对已有 sumstats 自动跳过。
# 10 个全部 munge 完成后退出。日志：work/data/_munge_autoimmune.log
# ============================================================
set -u
ROOT="${HORZONTAL_project:-/path/to/project}"
ROOT="${HORZ_ROOT:-$ROOT}"
RAW="$ROOT/gwas_data/finngen_r12"
MUNGED="$ROOT/project/work/data/munged"
MAN="$ROOT/project/work/data/finngen_R12_pheno_n.tsv"
CODE="10_munge_finngen.py"
AENDPOINTS="T1D M13_RHEUMA G6_MS K11_COELIAC K11_UC_STRICT2 K11_CD_STRICT2 L12_PSORIASIS M13_ANKYLOSPON M13_SJOGREN E4_GRAVES_STRICT"

n_munged() {
  local c=0
  for ph in $AENDPOINTS; do
    [ -s "$MUNGED/$ph.sumstats" ] && c=$((c+1))
  done
  echo $c
}

while [ "$(n_munged)" -lt 10 ]; do
  for ph in $AENDPOINTS; do
    [ -s "$MUNGED/$ph.sumstats" ] && continue
    [ -s "$RAW/finngen_R12_${ph}.gz" ] || continue
    # 有 .p 分片残留 = 未完成拼接，跳过等下一轮
    ls "$RAW/finngen_R12_${ph}.gz".p[0-9]* >/dev/null 2>&1 && continue
    nc=$(awk -F'\t' -v p="$ph" '$2==p{print $4}' "$MAN")
    nk=$(awk -F'\t' -v p="$ph" '$2==p{print $5}' "$MAN")
    if [ -z "${nc:-}" ] || [ -z "${nk:-}" ]; then echo "[watch] manifest 无 $ph"; continue; fi
    echo "[watch] $(date +%H:%M:%S) munge $ph (case=$nc ctrl=$nk)"
    python "$ROOT/project/work/code/$CODE" "$ph" "$nc" "$nk" 2>&1 | tail -1
  done
  sleep 45
done
echo "[watch] 全部 10 个自免终点 munge 完成 $(date +%H:%M:%S)"
