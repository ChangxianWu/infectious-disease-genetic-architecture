#!/usr/bin/env bash
# ============================================================
# project Stage 0.5：LDSC 官方资产（Zenodo 7768714）+ FinnGen R12 感染终点 分片并行下载
# ------------------------------------------------------------
# 复用 20B1 实测纪律（2026-10-10）：单连接经代理 120–180 KB/s，0.8 GB 需 1.5 h；
# GCS/Zenodo 均支持 HTTP Range，分片并发。每片独立 .pNN，全部完成后按序拼接，
# 核对总字节 == Content-Length 后改名（规则 69）；断点续传：大小正确的片跳过。
# 用法: bash project/work/code/00_download_assets.sh [workers] [shard_mb]
# ============================================================
set -u
ROOT="${HORZ_ROOT:-/path/to/project}"
W="${1:-10}"
SHMB="${2:-24}"

# ---- 队列 1：LDSC 官方最小资产（Zenodo 7768714，Alkes group README 指路）----
ZEN="https://zenodo.org/records/7768714/files"
LDIR="$ROOT/gwas_data/ldsc_ref"
mkdir -p "$LDIR"
ZFILES="1000G_Phase3_ldscores.tgz 1000G_Phase3_weights_hm3_no_MHC.tgz 1000G_Phase3_plinkfiles.tgz 1000G_Phase3_frq.tgz"

# ---- 队列 2：FinnGen R12 感染终点（20 个新文件；本地已有 6 个感染文件 + TB 系不入面板）----
BASE="https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
RAW="$ROOT/gwas_data/finngen_r12"
mkdir -p "$RAW"
FENDPOINTS="AB1_BACTINF_NOS AB1_ERYSIPELAS J10_PNEUMOBACT N14_PYELONEPHR
G6_MENINGBACT I9_ENDOCARD M13_OSTEOMYELITIS AB1_HERPES_SIMPLEX AB1_EBV
AB1_VIRAL_WARTS AB1_VARICELLA AB1_VIRAL_CNS J10_VIRALPNEUMO AB1_MYCOSES
AB1_CANDIDIASIS AB1_DERMATOPHYTOSIS AB1_PROTOZOAL AB1_PEDICULOSIS_ACARIASIS_OTHERINFEST
U22_COVID19_CONFIRMED AB1_GASTROENTERITIS_NOS"

MAN="$RAW/_20b2_manifest.tsv"
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
  # 校验拼接
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
  else
    echo "[失败] $label 拼接大小不符"; return 1
  fi
}

echo "===== 队列 1：LDSC Zenodo 资产 ====="
for f in $ZFILES; do
  dl_sharded "$ZEN/$f?download=1" "$LDIR/$f" "$f"
done

echo "===== 队列 2：FinnGen R12 感染终点 ====="
for ph in $FENDPOINTS; do
  f="finngen_R12_${ph}.gz"
  url="$BASE/$f"
  dl_sharded "$url" "$RAW/$f" "$ph"
  size=$(curl -sI --max-time 40 "$url" | tr -d '\r' | awk '/[Cc]ontent-[Ll]ength/{print $2; exit}')
  n=$(( (size + SHMB * 1024 * 1024 - 1) / (SHMB * 1024 * 1024) ))
  printf '%s\t%s\t%s\t%s\t%s\n' "$ph" "$f" "$size" "$n" "$url" >> "$MAN"
done
echo "[全部完成]"
