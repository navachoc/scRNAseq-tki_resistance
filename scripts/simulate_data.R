# simulate_data.R
# Fictitious 10x-style data to run the whole pipeline without real data.
# One group = one parental and one resistant sample, written in the format of
# Cell Ranger's filtered_feature_bc_matrix (barcodes, features, matrix).
#
# What the simulation contains (so that the figures look like real cell lines
# rather than separated blobs):
#   - a continuous cell-cycle program (S and G2/M genes follow a phase),
#   - 6 transcriptional states of unequal size that blend into each other,
#   - a "resistance" program (up- and down-regulated genes) whose strength is
#     high in most resistant cells, absent in most parental cells and present
#     in a small parental subpopulation (concentrated in two states),
#   - a sample-specific technical offset per gene (something for Harmony to fix),
#   - library-size variation and negative-binomial counts.
#
# Run (from the repo root):
#   GROUP <- "MET"; source("scripts/simulate_data.R")     # or "EGFR"
# Output: {GROUP}/data/{parental}/filtered_feature_bc_matrix/  and
#         {GROUP}/data/{resistant}/filtered_feature_bc_matrix/  and
#         {GROUP}/data/signatures/{GROUP}_resistance_signature_conservative.csv
# Sample names come from config/pipeline_config.R (SAMPLES).

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(org.Hs.eg.db)
})

if (!exists("GROUP")) GROUP <- "MET"
source("config/colors.R")
source("config/pipeline_config.R")

set.seed(if (GROUP == "MET") 11 else 23)
N_PAR <- if (GROUP == "MET") 4900 else 3800     # cells per sample
N_RES <- if (GROUP == "MET") 6000 else 4300
N_GENES <- 4000

# ------------------------------------------------------------------------------
# 1. Genes and baseline expression
# ------------------------------------------------------------------------------
symbols <- unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys(org.Hs.eg.db, "ENTREZID"),
                                                "SYMBOL", "ENTREZID")$SYMBOL))
symbols <- symbols[grepl("^[A-Z][A-Z0-9-]+$", symbols) & !grepl("^(LOC|LINC|MIR|SNOR|MT-)", symbols)]
cycle   <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
keep    <- unique(c(cycle, "MET", "EGFR", "NF2", "NRG1", "MDM2", "FGFR3", "TACC3", "YAP1", "AXL", "VIM"))
genes   <- unique(c(keep, sample(setdiff(symbols, keep), N_GENES)))[1:N_GENES]
mito    <- c("MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4", "MT-CYB")
genes   <- c(genes, mito)
ng      <- length(genes)
is_mito <- genes %in% mito
base    <- exp(rnorm(ng, -1.3, 1.4)); base[is_mito] <- 5

# ------------------------------------------------------------------------------
# 2. Programs (log-fold changes per gene)
# ------------------------------------------------------------------------------
free <- setdiff(seq_len(ng), c(which(is_mito), match(cycle, genes)))
K <- 6
state_lfc <- matrix(0, ng, K)
for (k in 1:K) { g <- sample(free, 90); state_lfc[g, k] <- runif(90, 0.6, 1.2) }
res_up   <- sample(free, 150); res_dn <- sample(setdiff(free, res_up), 100)
s_idx    <- match(cc.genes.updated.2019$s.genes, genes);   s_idx   <- s_idx[!is.na(s_idx)]
g2m_idx  <- match(cc.genes.updated.2019$g2m.genes, genes); g2m_idx <- g2m_idx[!is.na(g2m_idx)]
state_p  <- c(.30, .22, .17, .14, .10, .07)

simulate_sample <- function(n, res_mean, pre_frac) {
  # latent variables
  phase <- runif(n, 0, 2 * pi)
  main  <- sample(1:K, n, TRUE, state_p)
  w     <- matrix(0, n, K); w[cbind(seq_len(n), main)] <- 1
  nb    <- ((main + sample(c(-1, 1), n, TRUE) - 1) %% K) + 1       # neighbouring state: soft boundaries
  mixw  <- rbeta(n, 1.2, 4); w <- w * (1 - mixw); w[cbind(seq_len(n), nb)] <- w[cbind(seq_len(n), nb)] + mixw
  r <- if (is.null(pre_frac)) rbeta(n, res_mean * 6, (1 - res_mean) * 6) else {
    pre <- runif(n) < pre_frac * ifelse(main %in% c(2, 5), 3.5, 0.5)   # pre-existing cells concentrate in two states
    ifelse(pre, rbeta(n, 5, 3), rbeta(n, 1, 25))
  }
  # log-mean
  lm <- matrix(log(base), ng, n)
  lm <- lm + state_lfc %*% t(w)
  lm[s_idx, ]   <- lm[s_idx, ]   + 0.9 * rep(pmax(cos(phase - 1.2), 0), each = length(s_idx))
  lm[g2m_idx, ] <- lm[g2m_idx, ] + 0.9 * rep(pmax(cos(phase - 3.2), 0), each = length(g2m_idx))
  lm[res_up, ]  <- lm[res_up, ]  + 0.8 * rep(r, each = length(res_up))
  lm[res_dn, ]  <- lm[res_dn, ]  - 0.5 * rep(r, each = length(res_dn))
  list(lm = lm, r = r)
}

write_10x <- function(counts, sample_name, cells) {
  d <- file.path(GROUP, "data", sample_name, "filtered_feature_bc_matrix")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  mm <- file.path(d, "matrix.mtx"); writeMM(counts, mm); R.utils::gzip(mm, overwrite = TRUE)
  bc <- paste0(vapply(seq_len(ncol(counts)), function(i) paste(sample(c("A","C","G","T"), 16, TRUE), collapse = ""), ""), "-1")
  con <- gzfile(file.path(d, "barcodes.tsv.gz"), "w"); writeLines(make.unique(bc), con); close(con)
  con <- gzfile(file.path(d, "features.tsv.gz"), "w")
  writeLines(paste(paste0("ENSG", sprintf("%011d", seq_len(ng))), genes, "Gene Expression", sep = "\t"), con); close(con)
}

for (s in c("par", "res")) {
  n  <- if (s == "par") N_PAR else N_RES
  sim <- if (s == "par") simulate_sample(n, NA, 0.10) else simulate_sample(n, 0.85, NULL)
  batch <- rnorm(ng, 0, 0.12)                                   # technical offset of this sample
  depth <- exp(rnorm(n, 0, 0.4))
  mu <- exp(sim$lm + batch) * rep(depth, each = ng) * 1.4
  counts <- Matrix(matrix(rnbinom(length(mu), mu = mu, size = 3), ng, n), sparse = TRUE)
  rm(mu)
  nm <- SAMPLES[[GROUP]][if (s == "par") 1 else 2]
  write_10x(counts, nm, n)
  cat(sprintf("  %s: %d cells x %d genes written\n", nm, n, ng))
}

# ------------------------------------------------------------------------------
# 3. Signature (what a bulk RNA-seq comparison would hand to sc_05)
# ------------------------------------------------------------------------------
sig_dir <- file.path(GROUP, "data", "signatures"); dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)
sig <- data.frame(Gene = c(genes[res_up], genes[res_dn]),
                  LFC_Res_vs_Par = c(rep(2.5, length(res_up)), rep(-2.5, length(res_dn))),
                  Classification = c(rep("Upregulated", length(res_up)), rep("Downregulated", length(res_dn))))
write.csv(sig, file.path(sig_dir, sprintf("%s_resistance_signature_conservative.csv", GROUP)), row.names = FALSE)
cat(sprintf("Signature: %d genes -> %s\n", nrow(sig), sig_dir))
