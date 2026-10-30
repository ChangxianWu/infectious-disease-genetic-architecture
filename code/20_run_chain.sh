#!/usr/bin/env bash
# ============================================================
# project 接力链：等主下载队列完成 → 追加下载（J10_INFLUENZA + COVID HGI）
#   → munge 全部 25 终点 → 主矩阵 300 对（bivariate LDSC）
# 全程日志：work/data/_chain_2026-10-11.log
# 用法：bash project/work/code/20_run_chain.sh
# ============================================================
set -u
ROOT="${HORZ_ROOT:-/path/to/project}"
LOG="$ROOT/project/work/data/_chain_2026-10-11.log"
DLOG="$ROOT/project/work/data/_download_00.log"

echo "[chain $(date +%T)] 等待主下载队列（_download_00.log 出现 [全部完成]）..." | tee -a "$LOG"
for i in $(seq 1 240); do
  if grep -q "全部完成" "$DLOG" 2>/dev/null; then
    echo "[chain $(date +%T)] 主队列完成（等待 $((i)) 分钟内）" | tee -a "$LOG"
    break
  fi
  sleep 60
  if [ "$i" = "240" ]; then
    echo "[chain $(date +%T)] 等待超时（4h），主队列未完成——人工检查 _download_00.log" | tee -a "$LOG"
    exit 1
  fi
done

echo "[chain $(date +%T)] 追加下载：J10_INFLUENZA + COVID HGI freeze7 EUR" | tee -a "$LOG"
bash "$ROOT/project/work/code/01_download_covidhgi.sh" 10 24 2>&1 | tee -a "$LOG"

echo "[chain $(date +%T)] munge 全部 25 终点 + 主矩阵 300 对" | tee -a "$LOG"
cd "$ROOT"
python "project/work/code/12_matrix_runner.py" 2>&1 | tee -a "$LOG"
echo "[chain $(date +%T)] 接力链结束" | tee -a "$LOG"
