# sc_02_normalize_integrate.R
# Normalisation (SCTransform), integration (Harmony), PCA and UMAP
# scRNA-seq pipeline · NSCLC TKI resistance
#
# Run from the project root:
#   (open R inside the repo; the root is set by the .here file)
#   source("scripts/sc_02_normalize_integrate.R")
#   o bien: Rscript scripts/sc_02_normalize_integrate.R
#
# Input:
#   {GROUP}/results/{sampleA}_filtered.rds   (sc_01)
#   {GROUP}/results/{sampleB}_filtered.rds   (sc_01)
#
# Output:
#   {GROUP}/results/figures/norm_01_elbow.png
#   {GROUP}/results/figures/norm_02_qc_features.png
#   {GROUP}/results/figures/norm_03_before_after_integration.png
#   {GROUP}/results/reports/sc_02_normalize_integrate_report.md
#   {GROUP}/results/reports/sc_02_normalize_integrate_report.html
#   {GROUP}/results/{GROUP}_integrated.rds

# ==============================================================================
# LIBRARIES
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(glmGamPoi)
  library(harmony)
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

# Safety lock — this script ALWAYS recomputes SCTransform + Harmony from
# scratch and overwrites {GROUP}_integrated.rds; it has no "figures only" mode (the
# LogNormalize panel depends on an intermediate object that is not checkpointed).
# With RERUN_VIZ_ONLY <- TRUE, if the checkpoint already exists, the script stops
# without touching it. Set FALSE explicitly only when you really want to
# re-normalise / re-integrate from scratch (e.g. changed QC or parameters).
RERUN_VIZ_ONLY <- FALSE

# ==============================================================================
source("config/pipeline_config.R")
set.seed(42)

sample_a <- SAMPLE_PARENTAL
sample_b <- SAMPLE_RESISTANT

PROJECT_ROOT <- "."
output_dir   <- file.path(PROJECT_ROOT, GROUP, "results")
fig_dir      <- file.path(output_dir, "figures")
pfig <- function(f) file.path(fig_dir, paste0("sc_02_", f))
tab_dir      <- file.path(output_dir, "tables")
rep_dir      <- file.path(output_dir, "reports")

for (d in c(fig_dir, tab_dir, rep_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

integrated_path <- file.path(output_dir, sprintf("%s_integrated.rds", GROUP))
if (RERUN_VIZ_ONLY && file.exists(integrated_path)) {
  stop(
    "[RERUN_VIZ_ONLY] ", integrated_path, " already exists — script stopped so as not to ",
    "overwrite it (SCTransform/Harmony are not deterministic across package reinstalls). ",
    "If you really want to recompute from scratch, set RERUN_VIZ_ONLY <- FALSE."
  )
}

N_FEATURES <- 3000   # variable genes for PCA
N_PCS      <- 50     # PCs to compute
# DIMS_USE and N_NEIGHBORS come from pipeline_config.R (change them there if the ElbowPlot suggests it)

REGRESS_CELL_CYCLE <- TRUE   # TRUE: regress out the cell-cycle effect in SCTransform (recommended if G2/M dominates the UMAP)

# ==============================================================================
# 1. LOAD FILTERED OBJECTS (sc_01)
# ==============================================================================
cat("Loading filtered objects from sc_01...\n")

path_a <- file.path(output_dir, sprintf("%s_filtered.rds", sample_a))
path_b <- file.path(output_dir, sprintf("%s_filtered.rds", sample_b))

if (!file.exists(path_a)) stop("Not found: ", path_a, " — run sc_01_qc.R first")
if (!file.exists(path_b)) stop("Not found: ", path_b, " — run sc_01_qc.R first")

obj_a <- readRDS(path_a)
obj_b <- readRDS(path_b)

cat(sprintf("  %s: %d cells\n", sample_a, ncol(obj_a)))
cat(sprintf("  %s: %d cells\n", sample_b, ncol(obj_b)))

# ==============================================================================
# 1.5 CELL CYCLE SCORING
# ==============================================================================
# Adds S.Score, G2M.Score and Phase to each cell.
# Uses Seurat's canonical gene lists (Tirosh et al. 2015, converted to human symbols).
# If REGRESS_CELL_CYCLE = TRUE, SCTransform will regress S.Score and G2M.Score in step 2.
cat("\nComputing cell-cycle scores...\n")

s_genes   <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes

# Scoring on normalised RNA (quick LogNormalize for this step only)
obj_a <- NormalizeData(obj_a, verbose = FALSE)
obj_b <- NormalizeData(obj_b, verbose = FALSE)
obj_a <- CellCycleScoring(obj_a, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
obj_b <- CellCycleScoring(obj_b, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

cat(sprintf("  %s — fases: %s\n", sample_a, paste(table(obj_a$Phase), collapse = " | ")))
cat(sprintf("  %s — fases: %s\n", sample_b, paste(table(obj_b$Phase), collapse = " | ")))
if (REGRESS_CELL_CYCLE) cat("  REGRESS_CELL_CYCLE=TRUE: S.Score and G2M.Score will be regressed in SCTransform\n")

# ==============================================================================
# 1.6 PRE-SCTransform UMAP (LogNormalize — baseline for comparison)
# ==============================================================================
# Shows the structure without advanced correction to compare the effect of SCTransform + Harmony.
cat("\nGenerating pre-SCTransform UMAP (LogNormalize baseline)...\n")

raw_merged <- merge(obj_a, y = obj_b, add.cell.ids = c(sample_a, sample_b))
raw_merged <- NormalizeData(raw_merged, normalization.method = "LogNormalize",
                             scale.factor = 10000, verbose = FALSE)
raw_merged <- FindVariableFeatures(raw_merged, selection.method = "vst",
                                    nfeatures = N_FEATURES, verbose = FALSE)
raw_merged <- ScaleData(raw_merged, vars.to.regress = "percent.mt", verbose = FALSE)
raw_merged <- RunPCA(raw_merged, npcs = N_PCS, verbose = FALSE)
raw_merged <- RunUMAP(raw_merged, reduction = "pca", dims = DIMS_USE,
                       reduction.name = "umap.raw", verbose = FALSE)

p_raw_sample <- DimPlot(raw_merged, reduction = "umap.raw", group.by = "sample",
                         cols = sample_colors, pt.size = 0.3) +
  labs(title = "LogNormalize · no integration · by sample",
       x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "right")

p_raw_condition <- DimPlot(raw_merged, reduction = "umap.raw", group.by = "condition",
                            cols = condition_colors, pt.size = 0.3) +
  labs(title = "LogNormalize · no integration · by condition",
       x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "right")

rm(raw_merged); gc()

# ==============================================================================
# 2. SCTransform (per sample, before merging)
# ==============================================================================
# SCTransform run separately per sample: better parameter estimation per condition.
# glmGamPoi: faster implementation of the model fit.

options(future.globals.maxSize = 8000 * 1024^2)

cat("\nNormalising with SCTransform...\n")

vars_regress <- if (REGRESS_CELL_CYCLE) c("percent.mt", "S.Score", "G2M.Score") else "percent.mt"

cat(sprintf("  %s...\n", sample_a))
obj_a <- SCTransform(obj_a, vars.to.regress = vars_regress,
                     method = "glmGamPoi", verbose = FALSE)

cat(sprintf("  %s...\n", sample_b))
obj_b <- SCTransform(obj_b, vars.to.regress = vars_regress,
                     method = "glmGamPoi", verbose = FALSE)

cat("  SCTransform completed.\n")

# ==============================================================================
# 3. FEATURE SELECTION AND MERGE
# ==============================================================================
cat("\nSelecting variable genes and merging objects...\n")

features  <- SelectIntegrationFeatures(object.list = list(obj_a, obj_b),
                                       nfeatures   = N_FEATURES)
list_sct  <- PrepSCTIntegration(object.list     = list(obj_a, obj_b),
                                 anchor.features = features,
                                 verbose         = FALSE)
merged    <- merge(list_sct[[1]], y = list_sct[[2]],
                   add.cell.ids = c(sample_a, sample_b),
                   merge.data   = TRUE)

VariableFeatures(merged) <- features

cat(sprintf("  Merged object: %d cells, %d variable genes\n",
            ncol(merged), length(VariableFeatures(merged))))

# ==============================================================================
# 4. PCA Y ELBOW PLOT
# ==============================================================================
cat("\nComputing PCA...\n")

merged <- RunPCA(merged, npcs = N_PCS, verbose = FALSE)

pca_stdev <- merged[["pca"]]@stdev
pct_var   <- round(pca_stdev^2 / sum(pca_stdev^2) * 100, 2)

cat(sprintf("  Variance explained PC1–5: %.1f%%, %.1f%%, %.1f%%, %.1f%%, %.1f%%\n",
            pct_var[1], pct_var[2], pct_var[3], pct_var[4], pct_var[5]))

p_elbow <- ElbowPlot(merged, ndims = N_PCS) +
  geom_vline(xintercept = max(DIMS_USE), linetype = "dashed", color = "grey50") +
  labs(title = sprintf("ElbowPlot — PCA (dashed line = PC%d used for UMAP/Harmony)", max(DIMS_USE))) +
  white_bg

ggsave(pfig("norm_01_elbow.png"),
       p_elbow, width = 8, height = 5, dpi = 800)
cat("  Saved: norm_01_elbow.png\n")

# ==============================================================================
# 5. PRE-INTEGRATION UMAP (raw PCA)
# ==============================================================================
cat("\nComputing pre-integration UMAP (PCA)...\n")

merged <- RunUMAP(merged, reduction = "pca", dims = DIMS_USE,
                  reduction.name = "umap.pca", verbose = FALSE)

p_pre_sample <- DimPlot(merged, reduction = "umap.pca", group.by = "sample",
                         cols = sample_colors, pt.size = 0.3) +
  labs(title = "Pre-Harmony · by sample", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "right")

p_pre_condition <- DimPlot(merged, reduction = "umap.pca", group.by = "condition",
                            cols = condition_colors, pt.size = 0.3) +
  labs(title = "Pre-Harmony · by condition", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "right")

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
qc_gradient  <- c("#f0f0f0", "#c0392b")   # gris → rojo

p_pre_features <- lapply(qc_features, function(f)
  FeaturePlot(merged, reduction = "umap.pca", features = f,
              pt.size = 0.3, cols = qc_gradient) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(plot.title = element_text(size = 10))
)
p_pre_features[[1]] <- p_pre_features[[1]] + ggtitle("Pre-Harmony  |  nFeature_RNA")

# ==============================================================================
# 6. HARMONY (technical integration by sample)
# ==============================================================================
# group.by.vars = "sample": corrects technical variation between 10x runs.
# Do NOT integrate by "condition": preserve parental/resistant differences.

cat("\nRunning Harmony...\n")

merged <- RunHarmony(merged, group.by.vars = "sample", reduction = "pca",
                     reduction.save = "harmony", verbose = FALSE)
cat("  Harmony completed.\n")

# ==============================================================================
# 7. POST-INTEGRATION UMAP (Harmony)
# ==============================================================================
cat("Computing post-integration UMAP (Harmony)...\n")

merged <- RunUMAP(merged, reduction = "harmony", dims = DIMS_USE,
                  reduction.name = "umap.harmony", verbose = FALSE)

p_post_sample <- DimPlot(merged, reduction = "umap.harmony", group.by = "sample",
                          cols = sample_colors, pt.size = 0.3) +
  labs(title = "Post-Harmony · by sample", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "right")

p_post_condition <- DimPlot(merged, reduction = "umap.harmony", group.by = "condition",
                             cols = condition_colors, pt.size = 0.3) +
  labs(title = "Post-Harmony · by condition", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "right")

p_post_phase <- DimPlot(merged, reduction = "umap.harmony", group.by = "Phase",
                         pt.size = 0.3) +
  labs(title = "Post-Harmony · by cell cycle phase", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "right")

p_post_features <- lapply(qc_features, function(f)
  FeaturePlot(merged, reduction = "umap.harmony", features = f,
              pt.size = 0.3, cols = qc_gradient) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(plot.title = element_text(size = 10))
)
p_post_features[[1]] <- p_post_features[[1]] + ggtitle("Post-Harmony  |  nFeature_RNA")

# ==============================================================================
# 8. QC FEATURES PRE/POST HARMONY (figure norm_02)
# ==============================================================================
cat("Generating QC features pre/post Harmony...\n")

p_pre_row <- wrap_plots(p_pre_features,  ncol = 3) +
  plot_annotation(title = "Pre-Harmony",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))

p_post_row <- wrap_plots(p_post_features, ncol = 3) +
  plot_annotation(title = "Post-Harmony",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))

p_qc <- (p_pre_row / p_post_row) +
  plot_annotation(
    title    = sprintf("QC metrics — pre vs post Harmony integration (%s)", GROUP),
    subtitle = "Uniform distribution after integration indicates depth-driven structure removed",
    theme    = theme(plot.title    = element_text(size = 12, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(pfig("norm_02_qc_features.png"),
       p_qc, width = 15, height = 10, dpi = 800)
cat("  Saved: norm_02_qc_features.png\n")

# ==============================================================================
# 9. PANEL COMPARATIVO PRE/POST HARMONY (figure norm_03)
# ==============================================================================
cat("Generating pre/post integration comparison panel...\n")

p_comp <- (
  (p_raw_sample     + labs(title = "LogNormalize  |  by sample")) |
  (p_pre_sample     + labs(title = "Pre-Harmony   |  by sample")) |
  (p_post_sample    + labs(title = "Post-Harmony  |  by sample"))
) / (
  (p_raw_condition  + labs(title = "LogNormalize  |  by condition")) |
  (p_pre_condition  + labs(title = "Pre-Harmony   |  by condition")) |
  (p_post_condition + labs(title = "Post-Harmony  |  by condition"))
) +
  plot_annotation(
    title    = sprintf("Harmony integration effect — %s: %s vs %s", GROUP, sample_a, sample_b),
    subtitle = "Top: batch correction by sample · Bottom: parental/resistant biological separation",
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 10, color = "grey40"))
  )

ggsave(pfig("norm_03_before_after_integration.png"),
       p_comp, width = 21, height = 12, dpi = 800)
cat("  Saved: norm_03_before_after_integration.png\n")

# Cell cycle phase UMAP (norm_04)
phase_tbl <- table(merged$Phase)
p_cc <- p_post_phase +
  labs(subtitle = sprintf("G1: %d  S: %d  G2M: %d  |  REGRESS_CELL_CYCLE = %s",
                          phase_tbl["G1"], phase_tbl["S"], phase_tbl["G2M"],
                          REGRESS_CELL_CYCLE)) +
  white_bg
ggsave(pfig("norm_04_cell_cycle_phase.png"),
       p_cc, width = 8, height = 6, dpi = 800)
cat("  Saved: norm_04_cell_cycle_phase.png\n")

# ==============================================================================
# 9. SAVE INTEGRATED OBJECT
# ==============================================================================
out_rds <- file.path(output_dir, sprintf("%s_integrated.rds", GROUP))
saveRDS(merged, out_rds)
cat(sprintf("\nObject saved: %s\n", out_rds))
cat(sprintf("  Reducciones disponibles: %s\n",
            paste(names(merged@reductions), collapse = ", ")))

# Export to .h5ad for interoperability with Python/scanpy
if (requireNamespace("SeuratDisk", quietly = TRUE)) {
  library(SeuratDisk)
  out_h5seurat <- file.path(output_dir, sprintf("%s_integrated.h5seurat", GROUP))
  out_h5ad     <- file.path(output_dir, sprintf("%s_integrated.h5ad", GROUP))
  tryCatch({
    SaveH5Seurat(merged, filename = out_h5seurat, overwrite = TRUE, verbose = FALSE)
    Convert(out_h5seurat, dest = "h5ad", overwrite = TRUE, verbose = FALSE)
    file.remove(out_h5seurat)
    cat(sprintf("  h5ad saved: %s\n", out_h5ad))
  }, error = function(e) {
    cat(sprintf("  WARNING: could not export h5ad — %s\n", conditionMessage(e)))
  })
} else {
  cat("  WARNING: SeuratDisk not installed — skipping h5ad export\n")
  cat("         Install with: remotes::install_github('mojaveazure/seurat-disk')\n")
}

# ==============================================================================
# 10. GENERATE .md REPORT
# ==============================================================================
n_total    <- ncol(merged)
n_a        <- sum(merged$sample == sample_a)
n_b        <- sum(merged$sample == sample_b)
n_var_feat <- length(VariableFeatures(merged))

report_lines <- c(
  "================================================================",
  sprintf("NORMALIZE/INTEGRATE REPORT — %s: %s vs %s (scRNA-seq)", GROUP, sample_a, sample_b),
  sprintf("Fecha: %s", Sys.time()),
  "================================================================",
  "",
  "INTEGRATED OBJECT:",
  sprintf("  Total cells:          %d", n_total),
  sprintf("  %s (parental):   %d", sample_a, n_a),
  sprintf("  %s  (resistant): %d", sample_b, n_b),
  sprintf("  Variable genes:       %d  (SelectIntegrationFeatures)", n_var_feat),
  "",
  "CELL CYCLE SCORING:",
  sprintf("  S genes: %d  |  G2M genes: %d", length(s_genes), length(g2m_genes)),
  sprintf("  Distribution (merged) — G1: %d  S: %d  G2M: %d",
          sum(merged$Phase == "G1"), sum(merged$Phase == "S"), sum(merged$Phase == "G2M")),
  sprintf("  Regression applied: %s", ifelse(REGRESS_CELL_CYCLE, "YES (S.Score + G2M.Score)", "NO")),
  "",
  "NORMALISATION:",
  "  Method:    SCTransform (glmGamPoi)",
  sprintf("  Regression: %s", paste(vars_regress, collapse = " + ")),
  "  Run separately per sample before merging",
  "",
  sprintf("PCA  (npcs = %d):", N_PCS),
  sprintf("  PC1: %.2f%%  PC2: %.2f%%  PC3: %.2f%%  PC4: %.2f%%  PC5: %.2f%%",
          pct_var[1], pct_var[2], pct_var[3], pct_var[4], pct_var[5]),
  sprintf("  Suma varianza PC1–%d: %.1f%%", max(DIMS_USE), sum(pct_var[DIMS_USE])),
  "",
  "INTEGRATION:",
  "  Method:        Harmony",
  "  group.by.vars: sample  (corrects the technical batch between 10x runs)",
  "  condition NOT integrated: preserves parental/resistant differences",
  sprintf("  Dims usadas:   1:%d", max(DIMS_USE)),
  "",
  "SAVED REDUCTIONS:",
  "  pca          → basis for the pre-integration UMAP",
  "  harmony      → basis for the post-integration UMAP and clustering (sc_03)",
  "  umap.pca     → UMAP without batch correction",
  "  umap.harmony → post-Harmony UMAP (use by default)",
  "",
  "FIGURES GENERATED:",
  "  • norm_01_elbow.png                    — ElbowPlot, PC selection",
  "  • norm_02_qc_features.png              — QC metrics pre/post Harmony (2x3 grid)",
  "  • norm_03_before_after_integration.png — Comparative panel LogNormalize/Pre/Post-Harmony (3x2)",
  "  • norm_04_cell_cycle_phase.png         — UMAP coloured by cell-cycle phase",
  "",
  sprintf("SAVED OBJECT: %s/results/%s_integrated.rds  (%d cells)", GROUP, GROUP, n_total),
  sprintf("              %s/results/%s_integrated.h5ad  (Python/scanpy — if SeuratDisk is available)", GROUP, GROUP),
  "",
  "NEXT STEP: sc_03_clustering.R",
  sprintf("  → FindNeighbors(reduction = 'harmony', dims = 1:%d)", max(DIMS_USE)),
  "  → FindClusters (resoluciones 0.2, 0.4, 0.6, 0.8, 1.0)",
  "  → DimPlot by cluster and condition"
)

writeLines(report_lines, file.path(rep_dir, "sc_02_normalize_integrate_report.md"))
cat(sprintf("  Report .md saved: %s/sc_02_normalize_integrate_report.md\n", rep_dir))

cat("\n=== sc_02_normalize_integrate.R COMPLETED ===\n")
