# sc_01_qc.R
# QC, initial filtering and exploration
# scRNA-seq pipeline · NSCLC TKI resistance
#
# Run from the project root:
#   (open R inside the repo; the root is set by the .here file)
#   source("scripts/sc_01_qc.R")
#   o bien: Rscript scripts/sc_01_qc.R
#
# Output:
#   {GROUP}/results/figures/qc_01_violin_raw.png
#   {GROUP}/results/figures/qc_02_scatter_raw.png
#   {GROUP}/results/figures/qc_03_violin_filtered.png
#   {GROUP}/results/tables/qc_stats_pre_filter.csv
#   {GROUP}/results/tables/qc_stats_post_filter.csv
#   {GROUP}/results/reports/sc_01_qc_report.md
#   {GROUP}/results/reports/sc_01_qc_report.html
#   {GROUP}/results/{sampleA}_filtered.rds
#   {GROUP}/results/{sampleB}_filtered.rds

# ==============================================================================
# LIBRARIES
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# ==============================================================================
# CONFIG Y RUTAS
# ==============================================================================
setwd(here::here())   # repo root (set by the .here marker)
source("config/colors.R")
# ==============================================================================
# CONFIGURATION — edit only this section before running
# ==============================================================================

GROUP <- "MET"   # MET | EGFR

# ==============================================================================
source("config/pipeline_config.R")
set.seed(42)

sample_a <- SAMPLE_PARENTAL
sample_b <- SAMPLE_RESISTANT
cond_a   <- COND_PARENTAL
cond_b   <- COND_RESISTANT

PROJECT_ROOT <- "."
data_dir     <- file.path(PROJECT_ROOT, GROUP, "data")
output_dir   <- file.path(PROJECT_ROOT, GROUP, "results")
fig_dir      <- file.path(output_dir, "figures")
pfig <- function(f) file.path(fig_dir, paste0("sc_01_", f))
tab_dir      <- file.path(output_dir, "tables")
rep_dir      <- file.path(output_dir, "reports")

for (d in c(fig_dir, tab_dir, rep_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
dir_a <- file.path(data_dir, sample_a, "filtered_feature_bc_matrix")
dir_b <- file.path(data_dir, sample_b, "filtered_feature_bc_matrix")

if (!dir.exists(dir_a)) stop("Not found: ", dir_a)
if (!dir.exists(dir_b)) stop("Not found: ", dir_b)

cat(sprintf("Loading %s...\n", sample_a))
obj_a <- CreateSeuratObject(counts = Read10X(data.dir = dir_a),
                            project = sample_a, min.cells = 3, min.features = 200)
obj_a$sample    <- sample_a
obj_a$condition <- cond_a
cat(sprintf("  %s: %d cells, %d genes\n", sample_a, ncol(obj_a), nrow(obj_a)))

cat(sprintf("Loading %s...\n", sample_b))
obj_b <- CreateSeuratObject(counts = Read10X(data.dir = dir_b),
                            project = sample_b, min.cells = 3, min.features = 200)
obj_b$sample    <- sample_b
obj_b$condition <- cond_b
cat(sprintf("  %s: %d cells, %d genes\n", sample_b, ncol(obj_b), nrow(obj_b)))

# ==============================================================================
# 2. QC METRICS
# ==============================================================================
obj_a[["percent.mt"]]   <- PercentageFeatureSet(obj_a, pattern = "^MT-")
obj_b[["percent.mt"]]   <- PercentageFeatureSet(obj_b, pattern = "^MT-")
obj_a[["percent.ribo"]] <- PercentageFeatureSet(obj_a, pattern = "^RP[SL]")
obj_b[["percent.ribo"]] <- PercentageFeatureSet(obj_b, pattern = "^RP[SL]")

# Merge for pre-filter visualisation only
obj_raw <- merge(obj_a, y = obj_b, add.cell.ids = c(sample_a, sample_b))

# ==============================================================================
# 3. QC EXPLORATION — BEFORE FILTERING
# ==============================================================================
cat("\nGenerating pre-filter QC plots...\n")

qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")

p_violin_raw <- VlnPlot(
  obj_raw, features = qc_metrics, group.by = "sample",
  cols = sample_colors, pt.size = 0.05, ncol = 2
) & theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title  = element_text(size = 11))

ggsave(pfig("qc_01_violin_raw.png"),
       p_violin_raw, width = 10, height = 8, dpi = 800)
cat("  Saved: qc_01_violin_raw.png\n")

p_s1 <- FeatureScatter(obj_raw, "nCount_RNA",   "nFeature_RNA",
                        group.by = "sample", cols = sample_colors) + NoLegend()
p_s2 <- FeatureScatter(obj_raw, "nCount_RNA",   "percent.mt",
                        group.by = "sample", cols = sample_colors) + NoLegend()
p_s3 <- FeatureScatter(obj_raw, "nFeature_RNA", "percent.mt",
                        group.by = "sample", cols = sample_colors)

ggsave(pfig("qc_02_scatter_raw.png"),
       p_s1 / (p_s2 | p_s3), width = 14, height = 10, dpi = 800)
cat("  Saved: qc_02_scatter_raw.png\n")

# ==============================================================================
# 4. PRE-FILTER STATISTICS
# ==============================================================================
qc_stats_pre <- obj_raw@meta.data %>%
  group_by(sample) %>%
  summarise(
    n_cells       = n(),
    median_genes  = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_mito   = round(median(percent.mt), 2),
    pct_mito_gt30 = round(mean(percent.mt > 30) * 100, 1),
    median_ribo   = round(median(percent.ribo), 2),
    .groups = "drop"
  )

cat("\n--- QC pre-filter ---\n")
print(as.data.frame(qc_stats_pre))
write.csv(qc_stats_pre, file.path(tab_dir, "qc_stats_pre_filter.csv"), row.names = FALSE)

# ==============================================================================
# 5. FILTRADO
# Thresholds defined in config/pipeline_config.R → QC_PARAMS[[GROUP]]
# ⚠ REVIEW qc_01_violin_raw.png AND qc_02_scatter_raw.png BEFORE ADJUSTING
# ==============================================================================
MIN_FEATURES <- QC_PARAMS[[GROUP]]$MIN_FEATURES
MAX_FEATURES <- QC_PARAMS[[GROUP]]$MAX_FEATURES
MIN_COUNTS   <- QC_PARAMS[[GROUP]]$MIN_COUNTS
MAX_COUNTS   <- QC_PARAMS[[GROUP]]$MAX_COUNTS
MAX_MITO_PCT <- QC_PARAMS[[GROUP]]$MAX_MITO_PCT

cat(sprintf("\nApplying filters: genes [%d–%d], counts [%d–%d], mito <%.0f%%\n",
    MIN_FEATURES, MAX_FEATURES, MIN_COUNTS, MAX_COUNTS, MAX_MITO_PCT))

filter_cells <- function(obj) {
  subset(obj,
    nFeature_RNA > MIN_FEATURES & nFeature_RNA < MAX_FEATURES &
    nCount_RNA   > MIN_COUNTS   & nCount_RNA   < MAX_COUNTS   &
    percent.mt   < MAX_MITO_PCT)
}

obj_a_filt <- filter_cells(obj_a)
obj_b_filt <- filter_cells(obj_b)

cat("Cells retained:\n")
cat(sprintf("  %s: %d → %d  (%.1f%% retenidas)\n",
    sample_a, ncol(obj_a), ncol(obj_a_filt), ncol(obj_a_filt) / ncol(obj_a) * 100))
cat(sprintf("  %s: %d → %d  (%.1f%% retenidas)\n",
    sample_b, ncol(obj_b), ncol(obj_b_filt), ncol(obj_b_filt) / ncol(obj_b) * 100))

# ==============================================================================
# 6. PLOTS POST-FILTRADO
# ==============================================================================
cat("\nGenerating post-filter QC plots...\n")

obj_filt <- merge(obj_a_filt, y = obj_b_filt, add.cell.ids = c(sample_a, sample_b))

p_violin_filt <- VlnPlot(
  obj_filt, features = qc_metrics, group.by = "sample",
  cols = sample_colors, pt.size = 0.05, ncol = 2
) & theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title  = element_text(size = 11))

ggsave(pfig("qc_03_violin_filtered.png"),
       p_violin_filt, width = 10, height = 8, dpi = 800)
cat("  Saved: qc_03_violin_filtered.png\n")

# ==============================================================================
# 7. POST-FILTER STATISTICS
# ==============================================================================
qc_stats_post <- obj_filt@meta.data %>%
  group_by(sample) %>%
  summarise(
    n_cells       = n(),
    median_genes  = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_mito   = round(median(percent.mt), 2),
    .groups = "drop"
  )

cat("\n--- QC post-filter ---\n")
print(as.data.frame(qc_stats_post))
write.csv(qc_stats_post, file.path(tab_dir, "qc_stats_post_filter.csv"), row.names = FALSE)

# ==============================================================================
# 8. SAVE FILTERED OBJECTS
# ==============================================================================
saveRDS(obj_a_filt, file.path(output_dir, sprintf("%s_filtered.rds", sample_a)))
saveRDS(obj_b_filt, file.path(output_dir, sprintf("%s_filtered.rds", sample_b)))
cat(sprintf("\nObjects saved in %s\n", output_dir))

# ==============================================================================
# 9. GENERATE .md REPORT
# ==============================================================================
pct_ret_a <- round(ncol(obj_a_filt) / ncol(obj_a) * 100, 1)
pct_ret_b <- round(ncol(obj_b_filt) / ncol(obj_b) * 100, 1)

report_lines <- c(
  "================================================================",
  sprintf("QC REPORT — %s: %s vs %s (scRNA-seq)", GROUP, sample_a, sample_b),
  sprintf("Fecha: %s", Sys.time()),
  "================================================================",
  "",
  "RAW DATA (pre-filter):",
  sprintf("  %s  |  cells: %d  |  median genes: %.0f  |  median counts: %.0f  |  median mito: %.2f%%  |  %% cells mito>30: %.1f%%",
    sample_a,
    qc_stats_pre$n_cells[qc_stats_pre$sample == sample_a],
    qc_stats_pre$median_genes[qc_stats_pre$sample == sample_a],
    qc_stats_pre$median_counts[qc_stats_pre$sample == sample_a],
    qc_stats_pre$median_mito[qc_stats_pre$sample == sample_a],
    qc_stats_pre$pct_mito_gt30[qc_stats_pre$sample == sample_a]),
  sprintf("  %s  |  cells: %d  |  median genes: %.0f  |  median counts: %.0f  |  median mito: %.2f%%  |  %% cells mito>30: %.1f%%",
    sample_b,
    qc_stats_pre$n_cells[qc_stats_pre$sample == sample_b],
    qc_stats_pre$median_genes[qc_stats_pre$sample == sample_b],
    qc_stats_pre$median_counts[qc_stats_pre$sample == sample_b],
    qc_stats_pre$median_mito[qc_stats_pre$sample == sample_b],
    qc_stats_pre$pct_mito_gt30[qc_stats_pre$sample == sample_b]),
  "",
  "THRESHOLDS APPLIED:",
  sprintf("  MIN_FEATURES = %d  |  MAX_FEATURES = %d", MIN_FEATURES, MAX_FEATURES),
  sprintf("  MIN_COUNTS   = %d  |  MAX_COUNTS   = %d", MIN_COUNTS, MAX_COUNTS),
  sprintf("  MAX_MITO_PCT = %.0f%%", MAX_MITO_PCT),
  "",
  "METHODOLOGICAL NOTE — THRESHOLDS PER DATASET:",
  "  QC thresholds are set separately for each group (MET/EGFR) based on visual",
  "  inspection of the nFeature, nCount and percent.mt distributions.",
  "  The MET and EGFR models are biologically distinct cell lines with different",
  "  transcriptional complexity and sequencing depth; applying uniform thresholds",
  "  across datasets would be arbitrary. What matters for analytical validity is that both",
  "  samples within each comparison use identical thresholds, which keeps the",
  "  integration and differential analysis unbiased within each model.",
  "",
  "FILTERED DATA (post-filter):",
  sprintf("  %s  |  cells: %d → %d  (%.1f%% retained)  |  median genes: %.0f  |  median counts: %.0f  |  median mito: %.2f%%",
    sample_a, ncol(obj_a), ncol(obj_a_filt), pct_ret_a,
    qc_stats_post$median_genes[qc_stats_post$sample == sample_a],
    qc_stats_post$median_counts[qc_stats_post$sample == sample_a],
    qc_stats_post$median_mito[qc_stats_post$sample == sample_a]),
  sprintf("  %s  |  cells: %d → %d  (%.1f%% retained)  |  median genes: %.0f  |  median counts: %.0f  |  median mito: %.2f%%",
    sample_b, ncol(obj_b), ncol(obj_b_filt), pct_ret_b,
    qc_stats_post$median_genes[qc_stats_post$sample == sample_b],
    qc_stats_post$median_counts[qc_stats_post$sample == sample_b],
    qc_stats_post$median_mito[qc_stats_post$sample == sample_b]),
  "",
  "FIGURES GENERATED:",
  "  • qc_01_violin_raw.png      — distributions pre-filter",
  "  • qc_02_scatter_raw.png     — correlaciones nCount/nFeature/mito",
  "  • qc_03_violin_filtered.png — distributions post-filter",
  "",
  "SAVED OBJECTS:",
  sprintf("  • %s/results/%s_filtered.rds  (%d cells)", GROUP, sample_a, ncol(obj_a_filt)),
  sprintf("  • %s/results/%s_filtered.rds  (%d cells)", GROUP, sample_b, ncol(obj_b_filt)),
  "",
  "NEXT STEP: sc_02_normalize_integrate.R"
)

writeLines(report_lines, file.path(rep_dir, "sc_01_qc_report.md"))
cat(sprintf("  Report .md saved: %s/sc_01_qc_report.md\n", rep_dir))

cat("\n=== sc_01_qc.R COMPLETED ===\n")
