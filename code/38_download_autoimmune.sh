#!/usr/bin/env bash
# ============================================================
# project E3 自身免疫锚层：FinnGen R12 自免终点分片并行下载
# ------------------------------------------------------------
# 复用 00_download_assets.sh 的 dl_sharded 分片下载器（20B1/project 实测纪律）：
# 单连接 120–180 KB/s → 24 MB 分片 × 并发；Content-Length 校验；断点续传。
# 10 端点 × ~0.8 GB ≈ 8 GB；磁盘余 167 GB（2026-10-25 实测）。
# 端点选择（立项卡 E3 层，2026-10-25 定案，病例数取自 finngen_R12_pheno_n.tsv）：
#   T1D 4721 / M13_RHEUMA 16314 / G6_MS 2926 / K11_COELIAC 5130 /
#   K11_UC_STRICT2 7220 / K11_CD_STRICT2 2489 / L12_PSORIASIS 12760 /
#   M13_ANKYLOSPON 3838 / M13_SJOGREN 3309 / E4_GRAVES_STRICT 3962
#   （L12_LUPUS 850 例 lambda=0.98 不入——功效不足，如实记录于方案文档）
# 用法: bash project/work/code/38_download_autoimmune.sh [workers] [shard_mb]
# ============================================================
set -u
ROOT="${HORZ_ROOT:-/path/to/project}"
W="${1:-10}"
SHMB="${2:-24}"

BASE="https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
RAW="$ROOT/gwas_data/finngen_r12"
mkdir -p "$RAW"
AENDPOINTS="T1D M13_RHEUMA G6_MS K11_COELIAC K11_UC_STRICT2 K11_CD_STRICT2 L12_PSORIASIS M13_ANKYLOSPON M13_SJOGREN E4_GRAVES_STRICT"

MAN="$RAW/_20b2_autoimmune_manifest.tsv"
[ -f "$MAN" ] || printf 'phenocode\tfile\tsize\tn_shards\turl\n' > "$MAN"

dl_sharded() {
  local url="$1" dst="$2" label="$3"
  if [ -s "$dst" ]; then
    local have; have=$(curl -sI --max-time 40 "$url" | tr -d '\r' | awk '/[Cc]ontent-[Ll]ength/{print $2; exit}')
    local cur; cur=$(stat -c %s "$dst" 2>/dev/null || echo 0)
    if [ "${have:-0}" = "$cur" ]; then echo "[跳过] $label (${cur} B)"; return 0; fi
  fi
  local size; size=$(curl -sIL --max-time 60 "$url" | tr -d '\r' | awk '/[Cc]ontent-[Ll]ength/{v=$2} END{print v}')
  if [ -z "${size:-}" ] || [ "$size" -lt 1000 ]; then echo "[失败] $label 无 Content-Length"; return 1; fi

  local chunk=$(( SHMB * 1024 * 1024 ))
  local n=$(( (size + chunk - 1) / chunk ))
  echo "[下载] $label  $(( size / 1048576 )) MB  ${n} 片"

  local i p
  pids=()
  for i in $(seq 0 $(( n - 1 ))); do
    local pf="$dst.p$(printf '%04d' $i)"
    local want=$(( chunk )); local off=$(( i * chunk ))
    if [ $(( off + want )) -gt "$size" ]; then want=$(( size - off )); fi
    if [ -s "$pf" ] && [ "$(stat -c %s "$pf")" = "$want" ]; then continue; fi
    (
      for att in 1 2 3 4 5; do
        curl -s --max-time 900 --retry 2 -r "$off-$(( off + want - 1 ))" -o "$pf" "$url" && \
          [ "$(stat -c %s "$pf" 2>/dev/null || echo 0)" = "$want" ] && exit 0
        sleep $(( att * 5 ))
      done
      echo "[分片失败] $label p$i"; exit 1
    ) &
    pids+=($!)
    while [ "$(jobs -rp | wc -l)" -ge "$W" ]; do sleep 2; done
  done
  wait
  local ok=1
  for i in $(seq 0 $(( n - 1 ))); do
    local pf="$dst.p$(printf '%04d' $i)"
    local want=$(( chunk )); local off=$(( i * chunk ))
    if [ $(( off + want )) -gt "$size" ]; then want=$(( size - off )); fi
    if [ ! -s "$pf" ] || [ "$(stat -c %s "$pf")" != "$want" ]; then echo "[缺片] $label p$i"; ok=0; fi
  done
  [ "$ok" = "1" ] || { echo "[失败] $label 分片不全"; return 1; }
  cat $(for i in $(seq 0 $(( n - 1 ))); do echo "$dst.p$(printf '%04d' $i)"; done) > "$dst.part"
  if [ "$(stat -c %s "$dst.part")" = "$size" ]; then
    mv "$dst.part" "$dst"
    rm -f "$dst".p[0-9]*
    echo "[完成] $label  $size B"
    printf '%s\t%s\t%s\t%s\t%s\n' "$label" "$(basename "$dst")" "$size" "$n" "$url" >> "$MAN"
    return 0
  else
    echo "[失败] $label 拼接字节数不符"; return 1
  fi
}

FAIL=0
for ph in $AENDPOINTS; do
  dl_sharded "$BASE/finngen_R12_${ph}.gz" "$RAW/finngen_R12_${ph}.gz" "$ph" || FAIL=1
done

echo "===== 全部队列处理完毕 ====="
if [ "$FAIL" = "0" ]; then echo "[全部完成]"; else echo "[部分失败] 详见日志"; exit 1; fi
