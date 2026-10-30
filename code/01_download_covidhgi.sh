#!/usr/bin/env bash
# ============================================================
# project Stage 0.5b：追加下载 —— J10_INFLUENZA（面板补入）+ COVID HGI freeze7 EUR
# 等主队列（00_download_assets.sh）结束后运行；同款分片纪律（规则 69）。
# 用法: bash project/work/code/01_download_covidhgi.sh [workers] [shard_mb]
# ============================================================
set -u
ROOT="${HORZ_ROOT:-/path/to/project}"
W="${1:-10}"
SHMB="${2:-24}"
BASE_FG="https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
RAW="$ROOT/gwas_data/finngen_r12"
CDIR="$ROOT/project/work/data/covidhgi"
mkdir -p "$CDIR"
HGI="https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats"

source /dev/stdin <<'EOF'
dl_sharded() {
  local url="$1" dst="$2" label="$3"
  if [ -s "$dst" ]; then
    local have; have=$(curl -sIL --max-time 60 "$url" | tr -d '\r' | awk '/[Cc]ontent-[Ll]ength/{v=$2} END{print v}')
    local cur; cur=$(stat -c %s "$dst" 2>/dev/null || echo 0)
    if [ "${have:-0}" = "$cur" ]; then echo "[跳过] $label (${cur} B)"; return 0; fi
  fi
  local size; size=$(curl -sIL --max-time 60 "$url" | tr -d '\r' | awk '/[Cc]ontent-[Ll]ength/{v=$2} END{print v}')
  if [ -z "${size:-}" ] || [ "$size" -lt 1000 ]; then echo "[失败] $label 无 Content-Length"; return 1; fi
  local chunk=$(( SHMB * 1024 * 1024 ))
  local n=$(( (size + chunk - 1) / chunk ))
  echo "[下载] $label  $(( size / 1048576 )) MB  ${n} 片"
  local i
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
    mv "$dst.part" "$dst"; rm -f "$dst".p[0-9]*
    echo "[完成] $label  $size B"
  else
    echo "[失败] $label 拼接大小不符"; return 1
  fi
}
EOF

echo "===== 追加：FinnGen J10_INFLUENZA ====="
dl_sharded "$BASE_FG/finngen_R12_J10_INFLUENZA.gz" "$RAW/finngen_R12_J10_INFLUENZA.gz" "J10_INFLUENZA"

echo "===== 追加：COVID HGI freeze7 EUR（GRCh38 主版）====="
dl_sharded "$HGI/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz" "$CDIR/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz" "HGI_C2_eur"
dl_sharded "$HGI/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz" "$CDIR/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz" "HGI_A2_eur"
echo "[追加完成]"
