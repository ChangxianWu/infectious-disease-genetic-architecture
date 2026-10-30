#  19d_coloc_susie_final.R — 多信号共定位（coloc.susie，±100 kb 窗口，有界）
# ------------------------------------------------------------------
# 口径与限制（必读）：
#   * LD = 官方 LDSC 1000G EUR QC（plink -r square 导出，与主矩阵/LD score 同源）
#   * 窗口收窄到 ±100 kb：coloc.susie 的代价随信号数×信号数增长（实测 541-SNP/8 信号
#     窗口单对 2.2 min；±500 kb 的 CRP 窗口 2,907 SNP/10+ 信号曾跑 5 h 未收敛）——
#     故 susie 层只作**多信号对照**，主口径仍是 ±500 kb 的 coloc.abf（台账 §五）。
#   * 结果表 summary 是 (hit1 × hit2) 的多行表，**每行四假设之和为 1**（行内归一化）。
#     因此 **不能把对角行（hit1==hit2）的 PP.H4 当作区域级后验**：对角行在该行内部等于 1，
#     只是"若这两个信号就是同一信号，则共享"的条件后验。区域级后验必须**跨全部信号对
#     求和后再归一化**：PP.H4_global = Σ_rows PP.H4 / Σ_rows Σ_h PP.Hh。
#     第二轮评审 M-G 实测（FUT2 GASTRO_NOS×INTESTINAL）：9×6=54 行，Σ H3=53、Σ H4=1
#     → 对角口径 1.000、正确口径 **0.019**。故本版**同时输出两个口径**，正文只用全局口径。
#     注：K_loc 个 lead 同属一个区域时同一对会被重复测试（规则 222），输出按
#     (region, trait pair) 去重，并记录重复行数。
# 输出：work/output/coloc_susie_results.csv（逐对：信号数 + 双口径 H4 + 方法）
ROOT <- "/path/to/project"
CDIR <- file.path(ROOT, "work/data/coloc"); LDIR <- file.path(ROOT, "work/cache/susie_ld")
OUT  <- file.path(ROOT, "work/output")
suppressPackageStartupMessages({ library(coloc); library(susieR); library(data.table) })
Ntab <- fread(file.path(CDIR, "trait_N.tsv"))
Nof <- function(t) { v <- Ntab$N[Ntab$trait == t]; if (length(v)) as.numeric(v[1]) else 1e5 }
HALF <- 100000
LOCI <- c("rs492602", "rs516246", "rs495828", "rs3116651", "rs876537")
# lead → 区域（用于同区多 lead 去重；与 locus_atlas 的 region_id 一致）
ATLAS <- fread(file.path(OUT, "locus_atlas.csv"))
REGOF <- setNames(ATLAS$region_id, ATLAS$SNP)
res <- list(); t0 <- Sys.time()
for (LOCUS in LOCI) {
  ld_f <- file.path(LDIR, paste0(LOCUS, ".ld")); snp_f <- file.path(LDIR, paste0(LOCUS, ".ldsnp"))
  if (!file.exists(ld_f)) { cat("no LD", LOCUS, "\n"); next }
  snps <- readLines(snp_f); LD <- as.matrix(fread(ld_f, header = FALSE))
  colnames(LD) <- snps; rownames(LD) <- snps
  pat <- paste0("^", LOCUS, "_[A-Za-z0-9_]+[.]tsv$")
  fs <- list.files(CDIR, pattern = pat, full.names = TRUE)
  fs <- fs[!grepl("trait_N", fs)]
  if (length(fs) < 2) { cat("skip", LOCUS, "\n"); next }
  tr_of <- function(f) sub("[.]tsv$", "", sub(paste0("^", LOCUS, "_"), "", basename(f)))
  dd <- lapply(fs, function(f) as.data.frame(fread(f))); names(dd) <- sapply(fs, tr_of)
  focus <- dd[[1]]$pos[match(LOCUS, dd[[1]]$SNP)][1]
  keep <- snps[snps %in% dd[[1]]$SNP & abs(dd[[1]]$pos[match(snps, dd[[1]]$SNP)] - focus) <= HALF]
  cat(sprintf("[susie] %s: 窗口 %d SNP / %d 性状\n", LOCUS, length(keep), length(dd)))
  if (length(dd) >= 2 && length(keep) >= 50) for (i in seq_len(length(dd) - 1)) for (j in (i + 1):length(dd)) {
    a <- dd[[i]]; b <- dd[[j]]
    com <- Reduce(intersect, list(keep, a$SNP, b$SNP))
    if (length(com) < 50) next
    ca <- a[match(com, a$SNP), ]; cb <- b[match(com, b$SNP), ]
    LDs <- LD[com, com, drop = FALSE]
    D1 <- list(beta = ca$beta, varbeta = ca$sebeta^2, snp = ca$SNP, position = ca$pos,
               MAF = ca$MAF, N = Nof(names(dd)[i]), type = "cc", LD = LDs)
    D2 <- list(beta = cb$beta, varbeta = cb$sebeta^2, snp = cb$SNP, position = cb$pos,
               MAF = cb$MAF, N = Nof(names(dd)[j]), type = "cc", LD = LDs)
    m <- "susie"; ns1 <- ns2 <- NA_integer_; h4 <- NA_real_
    r <- tryCatch(suppressWarnings(suppressMessages(coloc.susie(D1, D2))),
                  error = function(e) { cat("  fallback:", conditionMessage(e), "\n"); NULL })
    if (is.null(r)) {
      m <- "abf"
      r2 <- tryCatch(suppressWarnings(coloc.abf(D1, D2)), error = function(e) NULL)
      if (!is.null(r2)) { h4 <- as.numeric(r2$summary[["PP.H4.abf"]]); h4_diag <- NA_real_; ncomb <- NA_integer_ }
    } else {
      s <- as.data.frame(r$summary)
      ns1 <- max(s$hit1, na.rm = TRUE); ns2 <- max(s$hit2, na.rm = TRUE)
      hh <- sprintf("PP.H%d.abf", 0:4)
      tot <- sum(sapply(hh, function(h) sum(s[[h]], na.rm = TRUE)))   # 跨行归一化分母
      diag <- s[s$hit1 == s$hit2, ]
      h4_diag <- if (nrow(diag)) max(diag$PP.H4.abf, na.rm = TRUE) else NA_real_
      h4 <- sum(s[["PP.H4.abf"]], na.rm = TRUE) / tot                 # 区域级（正确口径）
      ncomb <- nrow(s)
    }
    if (!is.finite(h4)) { cat("  NA H4:", m, names(dd)[i], names(dd)[j], "\n"); next }
    res[[length(res) + 1]] <- data.frame(locus = LOCUS,
      region = ifelse(is.na(REGOF[[LOCUS]]), LOCUS, REGOF[[LOCUS]]),
      trait1 = names(dd)[i], trait2 = names(dd)[j],
      method = m, window_kb = HALF / 1000, n_snp = length(com), n_signal_t1 = ns1,
      n_signal_t2 = ns2, n_signal_combinations = ncomb,
      PP_H4 = h4, PP_H4_diagonal_only = h4_diag, stringsAsFactors = FALSE)
    cat(sprintf("  [%s] %s x %s: H4_global=%.4f (diag-only=%.3f, %s combos) n=%d 信号=%s/%s (%.1f min)\n", m,
        names(dd)[i], names(dd)[j], h4, h4_diag, ncomb, length(com), ns1, ns2,
        as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
if (length(res)) {
  out <- do.call(rbind, res)
  out$verdict <- ifelse(out$PP_H4 > 0.8, "share", ifelse(out$PP_H4 > 0.5, "suggest", "not_support"))
  # 同区多 lead 的重复测试去重（规则 222）：一个 (区域, 对) 只保留第一行，记录删除数
  before <- nrow(out)
  out <- out[!duplicated(paste(out$region, out$trait1, out$trait2)), ]
  cat(sprintf("[susie] 区域级去重：%d -> %d 行（删除 %d 次重复测试）\n", before, nrow(out), before - nrow(out)))
  fwrite(out, file.path(OUT, "coloc_susie_results.csv"))
  cat(sprintf("\n[susie] 完成 %d 对，用时 %.1f min\n", nrow(out), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  print(out)
} else cat("[susie] 无结果\n")
