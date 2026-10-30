# ============================================================
#  18_coloc.R — 共享位点共定位（R 4.4.2；coloc）
# 口径（规则 71/96；与 20B1/TB2 沿用同一判定带）：
#   输入 = work/data/coloc/<LOCUS>_<TRAIT>.tsv（FinnGen 官方 beta/sebeta 与 af_alt 派生 MAF）
#   样本量 = work/data/coloc/trait_N.tsv（官方 manifest，case+control）
#   病例比例 s = ncase/(ncase+ncontrol)（case-control 型必须给 s，否则放置 NA 会告警）
#   方法 = coloc.abf（单因果变异先验）；PP.H4>0.8 支持共享 / 0.5–0.8 提示 / <0.5 不支持
#   敏感性：p12 网格 1e-6 / 1e-5 / 1e-4 三档并列（p12 为先验共享概率）
# 用法：Rscript 18_coloc.R [LOCUS]
# 输出：work/output/coloc_results.csv、coloc_p12_sensitivity.csv
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
ONLY <- if (length(args) >= 1) args[1] else NULL

ROOT <- "/path/to/project"
CDIR <- file.path(ROOT, "work/data/coloc")
OUT  <- file.path(ROOT, "work/output")
suppressPackageStartupMessages({ library(coloc); library(data.table) })

Ntab <- fread(file.path(CDIR, "trait_N.tsv"))
Nof <- function(t) { v <- Ntab$N[Ntab$trait == t]; if (length(v)) as.numeric(v[1]) else 1e5 }

files <- list.files(CDIR, pattern = "\\.tsv$", full.names = TRUE)
files <- files[!grepl("trait_N\\.tsv$", files)]
if (!is.null(ONLY)) files <- files[grepl(paste0("/", ONLY, "_"), files)]
snp_of <- function(f) sub("_.*$", "", basename(f))
trait_of <- function(f) sub("\\.tsv$", "", sub("^[^_]+_", "", basename(f)))
groups <- split(files, sapply(files, snp_of))
message(sprintf("[coloc] %d 个位点、%d 个窗口文件", length(groups), length(files)))

res <- list(); sens <- list()
for (snp in names(groups)) {
  fs <- groups[[snp]]
  if (length(fs) < 2) { message("[coloc] ", snp, " 仅 ", length(fs), " 个性状，跳过"); next }
  dat <- lapply(fs, function(f) {
    d <- fread(f)
    d <- d[nchar(SNP) > 0 & is.finite(beta) & is.finite(sebeta) & sebeta > 0]
    list(d = d, trait = trait_of(f))
  })
  names(dat) <- sapply(dat, function(x) x$trait)
  if (length(unique(names(dat))) < length(dat)) next
  for (i in seq_len(length(dat) - 1)) {
    for (j in (i + 1):length(dat)) {
      a <- dat[[i]]; b <- dat[[j]]
      ca <- a$d; cb <- b$d
      common <- intersect(ca$SNP, cb$SNP)
      if (length(common) < 100) {
        message("[coloc] ", snp, " ", a$trait, " x ", b$trait, " 交集 ", length(common), " SNP，跳过")
        next
      }
      ca <- ca[match(common, ca$SNP)]; cb <- cb[match(common, cb$SNP)]
      mk <- function(x, tr) {
        # 若窗口自带逐 SNP 样本量（HGI 端由 N 列注入）则用向量，否则用 manifest 标量
        Nv <- if ("N" %in% names(x)) as.numeric(x$N) else Nof(tr)
        list(beta = x$beta, varbeta = x$sebeta^2, snp = x$SNP,
             position = x$pos, MAF = x$MAF, N = Nv, type = "cc")
      }
      D1 <- mk(ca, a$trait); D2 <- mk(cb, b$trait)
      r <- tryCatch(suppressWarnings(coloc.abf(D1, D2)), error = function(e) NULL)
      if (is.null(r)) { message("[coloc] ", snp, " ", a$trait, "x", b$trait, " 失败"); next }
      s <- r$summary
      res[[length(res) + 1]] <- data.frame(
        locus = snp, trait1 = a$trait, trait2 = b$trait,
        n_snp = s[["nsnps"]], PP_H0 = s[["PP.H0.abf"]], PP_H1 = s[["PP.H1.abf"]],
        PP_H2 = s[["PP.H2.abf"]], PP_H3 = s[["PP.H3.abf"]], PP_H4 = s[["PP.H4.abf"]],
        top_snp = paste(r$results$snp[order(r$results$SNP.PP.H4, decreasing = TRUE)][1:3],
                        collapse = ";"),
        stringsAsFactors = FALSE)
      for (p12 in c(1e-6, 1e-5, 1e-4)) {
        r2 <- tryCatch(suppressWarnings(coloc.abf(D1, D2, p12 = p12)), error = function(e) NULL)
        if (!is.null(r2)) sens[[length(sens) + 1]] <- data.frame(
          locus = snp, trait1 = a$trait, trait2 = b$trait, p12 = p12,
          PP_H4 = r2$summary[["PP.H4.abf"]], n_snp = r2$summary[["nsnps"]])
      }
      message(sprintf("[coloc] %s %s x %s: PP.H4=%.3f (n=%d)",
                      snp, a$trait, b$trait, s[["PP.H4.abf"]], s[["nsnps"]]))
    }
  }
}
if (length(res)) {
  out <- do.call(rbind, res)
  out$verdict <- ifelse(out$PP_H4 > 0.8, "share",
                 ifelse(out$PP_H4 > 0.5, "suggest", "not_support"))
  fwrite(out, file.path(OUT, "coloc_results.csv"))
  message("[coloc] 落盘 ", nrow(out), " 对 -> coloc_results.csv")
  print(table(out$verdict))
  if (length(sens)) fwrite(do.call(rbind, sens), file.path(OUT, "coloc_p12_sensitivity.csv"))
} else {
  message("[coloc] 无结果")
}
