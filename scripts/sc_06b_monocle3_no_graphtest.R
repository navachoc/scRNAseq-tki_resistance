# ==============================================================================
# sc_06b_monocle3_no_graphtest.R — Trajectories and pseudotime (Monocle3), WITHOUT graph_test
# Variant of sc_06_monocle3.R that skips Section 6 (graph_test, Moran's I),
# the slow step (~40-50 min, single core over ~18,800 genes). Useful to
# rebuild the CDS and the trajectory UMAP quickly after re-clustering,
# without waiting for the table/figure of genes ranked by autocorrelation.
#
# To also get traj_graph_test_results.csv and the full traj_06 (with the
#   top-genes-by-Moran's-I panel), run sc_06_monocle3.R with
# RERUN_VIZ_ONLY <- FALSE.
#
# MET: MET1 vs MET1R  |  Adaptable a EGFR: EGFR1 vs EGFR1_R
# Agosto 2026
# ==============================================================================
#
# Biological question:
#   Do parental cells with a high resistance signature (pre-existing)
#   occupy intermediate or advanced pseudotime positions towards the resistant state?
#   If so, pseudotime links the pre-existing state to full resistance.
#
# Prerequisitos:
#   - monocle3:       BiocManager::install("monocle3")
#   - SeuratWrappers: remotes::install_github("satijalab/seurat-wrappers")
#   - Input:          MET/results/MET_signature.rds  (output sc_05)
#
# Run (from the repo root):
#   source("scripts/sc_06b_monocle3_no_graphtest.R")
#   Rscript scripts/sc_06b_monocle3_no_graphtest.R
#
# Outputs (MET):
#   figures/ : traj_01_graph_umap.png
#              traj_02_pseudotime_split.png
#              traj_03_pseudotime_all.png
#              traj_04_pseudotime_by_sample.png
#              traj_05_preexistence_pseudotime.png
#              traj_06_genes_pseudotime.png  (driver-gene panel only, no graph_test)
#   tables/  : traj_pseudotime_per_cell.csv
#   reports/ : sc_06_monocle3_report.md
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tibble)
  if (requireNamespace("ggridges", quietly = TRUE)) library(ggridges)
})

setwd(here::here())   # repo root (set by the .here marker)

# ==============================================================================
# CONFIGURATION — edit only this section before running
# ==============================================================================

GROUP          <- "MET"     # MET | EGFR
root_method    <- "condition" # "cluster" | "condition" | "signature"
root_cluster   <- "0"       # only if root_method == "cluster"

# ==============================================================================
source("config/colors.R")
source("config/pipeline_config.R")   # define SAMPLE_PARENTAL/RESISTANT, COND_*, sample_colors

# Method to define the pseudotime root cells:
#   "cluster"   — parental cells of the most enriched cluster in MET1
#   "condition" — all parental cells (MET1)
#   "signature" — parental cells with the lowest UCell (percentile root_score_pct)
root_score_pct <- 10    # UCell percentile; only if root_method == "signature"

# Key genes to display along pseudotime — overwritten after source(pipeline_config.R)
GENES_PSEUDOTIME <- GENES_RESISTANCE_PANEL[[GROUP]]
# ==============================================================================

set.seed(42)

output_dir <- file.path(".", GROUP, "results")
fig_dir    <- file.path(output_dir, "figures")
tab_dir    <- file.path(output_dir, "tables")
rep_dir    <- file.path(output_dir, "reports")
pfig <- function(f) file.path(fig_dir, paste0("sc_06_", f))
for (d in c(fig_dir, tab_dir, rep_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. LOAD SEURAT OBJECT (output of sc_05)
# ==============================================================================

cat("\n[Sec 1] Loading Seurat object...\n")
rds_in <- file.path(output_dir, sprintf("%s_signature.rds", GROUP))
if (!file.exists(rds_in)) {
  rds_in <- file.path(output_dir, sprintf("%s_clustered.rds", GROUP))
  cat("  Using MET_clustered.rds (sc_05 not available)\n")
}
if (!file.exists(rds_in)) stop("Missing: ", rds_in, " — run sc_03 or sc_05 first")

seurat_obj <- readRDS(rds_in)
Idents(seurat_obj) <- "seurat_clusters"
cat(sprintf("  Cells: %d | Clusters: %s\n",
            ncol(seurat_obj),
            paste(levels(seurat_obj$seurat_clusters), collapse = ", ")))

# Check that the score columns exist (from sc_05)
ucell_col          <- sprintf("%s_resistance_up_UCell", GROUP)
ucell_combined_col <- sprintf("%s_UCell_combined", GROUP)
has_score          <- ucell_col %in% colnames(seurat_obj@meta.data)
has_combined       <- ucell_combined_col %in% colnames(seurat_obj@meta.data)
if (!has_score)    cat("  Note: UCell scores not found — traj_05 not available\n")
if (!has_combined) cat("  Note: UCell_combined not found — NA will be exported in the table\n")

# ==============================================================================
# 2. SEURAT → MONOCLE3 CONVERSION
# ==============================================================================

cat("\n[Sec 2] Converting Seurat → CellDataSet (Monocle3)...\n")

# as.cell_data_set transfers UMAP, PCA, clusters and metadata from Seurat
cds <- as.cell_data_set(seurat_obj)

# plot_cells(..., genes=) requires rowData$gene_short_name
rowData(cds)$gene_short_name <- rownames(cds)

# Make sure 'condition' and 'sample' are in colData
colData(cds)$condition <- seurat_obj$condition[rownames(colData(cds))]
colData(cds)$sample    <- seurat_obj$sample[rownames(colData(cds))]

# Transfer the Seurat UMAP (umap.harmony) to Monocle3 as the "UMAP" reduction
# Needed because as.cell_data_set does not recognise reductions with non-standard names
reducedDims(cds)[["UMAP"]] <- Embeddings(seurat_obj, "umap.harmony")

# Transfer Seurat clusters to Monocle3
cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3)

cat(sprintf("  CellDataSet: %d cells, %d genes\n",
            ncol(cds), nrow(cds)))

# ==============================================================================
# 3. LEARN THE TRAJECTORY
# ==============================================================================

cat("\n[Sec 3] Learning the trajectory graph...\n")

# use_partition = FALSE: forces a single connected graph (more interpretable with 2 conditions)
# learn_graph is the most computationally expensive step of this section (~2-5 min with 11k cells)
cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)

cat("  Grafo aprendido\n")

# Checkpoint: save the CDS after learn_graph
cds_path <- file.path(output_dir, sprintf("%s_cds_trajectory.rds", GROUP))
saveRDS(cds, cds_path)
cat(sprintf("  CDS saved: %s\n", basename(cds_path)))

# ==============================================================================
# 4. ORDER CELLS IN PSEUDOTIME
# ==============================================================================

cat("\n[Sec 4] Ordering cells in pseudotime...\n")
cat(sprintf("  root_method = '%s'\n", root_method))

# --- Option A: by cluster ---
if (root_method == "cluster") {
  root_cells <- colnames(seurat_obj)[seurat_obj$seurat_clusters == root_cluster &
                                       seurat_obj$condition == COND_PARENTAL]
  if (length(root_cells) == 0) {
    cat("  Warning: no parental cells in root_cluster — using the whole cluster\n")
    root_cells <- colnames(seurat_obj)[seurat_obj$seurat_clusters == root_cluster]
  }
  root_label <- sprintf("cluster %s (%s)", root_cluster, COND_PARENTAL)

# --- Option B: all parental cells ---
} else if (root_method == "condition") {
  root_cells <- colnames(seurat_obj)[seurat_obj$condition == COND_PARENTAL]
  root_label <- sprintf("all %s cells (n=%d)", COND_PARENTAL, length(root_cells))

# --- Option D: parental cells with the lowest resistance signature ---
} else if (root_method == "signature") {
  if (!has_score) stop("root_method='signature' requires UCell scores from sc_05")
  parental_idx  <- seurat_obj$condition == COND_PARENTAL
  scores_par    <- seurat_obj@meta.data[[ucell_col]][parental_idx]
  threshold     <- quantile(scores_par, root_score_pct / 100, na.rm = TRUE)
  root_cells    <- colnames(seurat_obj)[parental_idx &
                     seurat_obj@meta.data[[ucell_col]] <= threshold]
  root_label    <- sprintf("parental with UCell ≤ p%d (n=%d)", root_score_pct, length(root_cells))

} else {
  stop("root_method must be 'cluster', 'condition' or 'signature'")
}

cat(sprintf("  Root cells: %d  (%s)\n", length(root_cells), root_label))
cds <- order_cells(cds, root_cells = root_cells, reduction_method = "UMAP")

pseudotime_vals <- pseudotime(cds)
cat(sprintf("  Pseudotime: min=%.2f, max=%.2f, NA=%d\n",
            min(pseudotime_vals, na.rm = TRUE),
            max(pseudotime_vals, na.rm = TRUE),
            sum(is.infinite(pseudotime_vals))))

# Transfer pseudotime to the Seurat object for Seurat-based visualisations
seurat_obj$pseudotime <- ifelse(is.infinite(pseudotime_vals[colnames(seurat_obj)]),
                                 NA_real_,
                                 pseudotime_vals[colnames(seurat_obj)])

# ==============================================================================
# 5. VISUALIZACIONES MONOCLE3
# ==============================================================================

cat("\n[Sec 5] Visualizaciones...\n")

# --- traj_01: UMAP with trajectory graph (coloured by Seurat cluster) ---
p1 <- plot_cells(
  cds,
  color_cells_by          = "seurat_clusters",
  label_groups_by_cluster = FALSE,
  label_leaves            = FALSE,
  label_branch_points     = FALSE,
  graph_label_size        = 0,
  cell_size               = 0.5,
  trajectory_graph_segment_size = 1.0,
  show_trajectory_graph   = TRUE
) +
  scale_color_manual(values = cluster_colors) +
  labs(
    title    = sprintf("Monocle3 trajectory — %s: %s vs %s", GROUP, COND_PARENTAL, COND_RESISTANT),
    subtitle = sprintf("Graph over Harmony UMAP · root = %s", root_label)
  ) +
  theme_classic(base_size = 12) + white_bg

ggsave(pfig("traj_01_graph_umap.png"),
       p1, width = 8, height = 7, dpi = 300)
cat("  traj_01_graph_umap.png\n")

# --- traj_02: Pseudotime split by condition ---
colData(cds)$pseudotime_sqrt <- sqrt(pseudotime(cds))
colData(cds)$pseudotime_sqrt[is.infinite(colData(cds)$pseudotime_sqrt)] <- NA

p2 <- plot_cells(
  cds,
  color_cells_by        = "pseudotime",
  label_cell_groups     = FALSE,
  label_leaves          = FALSE,
  label_branch_points   = FALSE,
  graph_label_size      = 0,
  cell_size             = 0.5,
  trajectory_graph_segment_size = 0.8,
  show_trajectory_graph = TRUE
) +
  scale_color_viridis_c(option = "viridis", trans = "sqrt",
                        name = "Pseudotime\n(√ scale)") +
  facet_wrap(~sample) +
  labs(
    title    = "Pseudotime split by sample",
    subtitle = sprintf("%s vs %s — gradient yellow→purple = start→end of trajectory",
                       SAMPLE_PARENTAL, SAMPLE_RESISTANT)
  ) +
  theme_classic(base_size = 12) + white_bg

ggsave(pfig("traj_02_pseudotime_split.png"),
       p2, width = 14, height = 6, dpi = 300)
cat("  traj_02_pseudotime_split.png\n")

# --- traj_03: Pseudotime, all cells ---
p3 <- plot_cells(
  cds,
  color_cells_by        = "pseudotime",
  label_cell_groups     = FALSE,
  label_leaves          = FALSE,
  label_branch_points   = FALSE,
  graph_label_size      = 0,
  cell_size             = 0.5,
  trajectory_graph_segment_size = 0.8,
  show_trajectory_graph = TRUE
) +
  scale_color_viridis_c(option = "viridis", trans = "sqrt",
                        name = "Pseudotime\n(√ scale)") +
  labs(
    title    = "Pseudotime — all cells",
    subtitle = "Gradiente viridis (escala √): amarillo = inicio (parental), violeta = fin (resistant)"
  ) +
  theme_classic(base_size = 12) + white_bg

ggsave(pfig("traj_03_pseudotime_all.png"),
       p3, width = 8, height = 7, dpi = 300)
cat("  traj_03_pseudotime_all.png\n")

pt_df <- data.frame(
  barcode    = colnames(seurat_obj),
  pseudotime = seurat_obj$pseudotime,
  cluster    = as.character(seurat_obj$seurat_clusters),
  condition  = seurat_obj$condition,
  sample     = seurat_obj$sample
) |> filter(!is.na(pseudotime))

# --- traj_04: Pseudotime distribution per sample ---
if (requireNamespace("ggridges", quietly = TRUE)) {
  p4 <- ggplot(pt_df, aes(x = pseudotime, y = sample, fill = sample)) +
    ggridges::geom_density_ridges(
      alpha = 0.75, scale = 1.2,
      quantile_lines = TRUE, quantiles = 2,
      bandwidth = 0.25
    ) +
    scale_fill_manual(values = sample_colors) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    labs(
      title    = sprintf("Pseudotime distribution — %s vs %s", SAMPLE_PARENTAL, SAMPLE_RESISTANT),
      subtitle = "Central line = median · resistant cells shifted towards higher pseudotime",
      x = "Pseudotime", y = NULL
    ) +
    theme_classic(base_size = 13) + white_bg + NoLegend()
} else {
  p4 <- ggplot(pt_df, aes(x = pseudotime, fill = sample, color = sample)) +
    geom_density(alpha = 0.5, linewidth = 0.8) +
    scale_fill_manual(values  = sample_colors) +
    scale_color_manual(values = sample_colors) +
    labs(
      title    = sprintf("Pseudotime distribution — %s vs %s", SAMPLE_PARENTAL, SAMPLE_RESISTANT),
      subtitle = "A resistant distribution shifted towards higher pseudotime indicates the trajectory direction",
      x = "Pseudotime", y = "Density"
    ) +
    theme_classic(base_size = 13) + white_bg
}

ggsave(pfig("traj_04_pseudotime_by_sample.png"),
       p4, width = 10, height = 6, dpi = 300)
cat("  traj_04_pseudotime_by_sample.png\n")

# --- traj_05: 3 grupos — Parental low score / Parental high score / Resistente ---
if (has_score && "high_score_parental" %in% colnames(seurat_obj@meta.data)) {
  pt_df$high_score <- seurat_obj$high_score_parental[match(pt_df$barcode, colnames(seurat_obj))]

  pt_3grp <- pt_df |>
    mutate(pt_group = case_when(
      condition == COND_RESISTANT                             ~ "Resistant",
      condition == COND_PARENTAL & grepl("high score", high_score, ignore.case = TRUE) ~ "Parental high score",
      condition == COND_PARENTAL & grepl("low score",  high_score, ignore.case = TRUE) ~ "Parental low score",
      TRUE ~ NA_character_
    )) |>
    filter(!is.na(pt_group), !is.na(pseudotime))

  pt_3grp$pt_group <- factor(pt_3grp$pt_group,
                           levels = c("Parental low score", "Parental high score", "Resistant"))

  color_3grp <- c(
    "Parental low score" = "#BDC3C7",
    "Parental high score" = "#E17055",
    "Resistant"          = "#2D3436"
  )

  n_3grp <- table(pt_3grp$pt_group)

  wt_ab <- wilcox.test(
    x     = pt_3grp$pseudotime[pt_3grp$pt_group == "Parental high score"],
    y     = pt_3grp$pseudotime[pt_3grp$pt_group == "Parental low score"],
    exact = FALSE
  )
  wt_ar <- wilcox.test(
    x     = pt_3grp$pseudotime[pt_3grp$pt_group == "Parental high score"],
    y     = pt_3grp$pseudotime[pt_3grp$pt_group == "Resistant"],
    exact = FALSE
  )
  p7_pval_ab <- wt_ab$p.value
  p7_pval_ar <- wt_ar$p.value

  p_ab <- formatC(wt_ab$p.value, format = "e", digits = 2)
  p_ar <- formatC(wt_ar$p.value, format = "e", digits = 2)

  y_max <- max(pt_3grp$pseudotime, na.rm = TRUE)

  pt_high <- pt_3grp |> filter(pt_group == "Parental high score")

  p7 <- ggplot(pt_3grp, aes(x = pt_group, y = pseudotime, fill = pt_group)) +
    geom_violin(trim = TRUE, alpha = 0.80) +
    geom_jitter(
      data  = pt_high,
      aes(x = pt_group, y = pseudotime),
      width = 0.08, size = 0.6, alpha = 0.25,
      color = "#c0392b", inherit.aes = FALSE
    ) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.75) +
    scale_fill_manual(values = color_3grp) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    annotate("segment", x = 1, xend = 2,
             y = y_max * 1.02, yend = y_max * 1.02, linewidth = 0.6) +
    annotate("text", x = 1.5, y = y_max * 1.07,
             label = sprintf("p = %s", p_ab), size = 3.8) +
    annotate("segment", x = 2, xend = 3,
             y = y_max * 1.13, yend = y_max * 1.13, linewidth = 0.6) +
    annotate("text", x = 2.5, y = y_max * 1.18,
             label = sprintf("p = %s", p_ar), size = 3.8) +
    labs(
      title    = "Pre-existence evidence — Pseudotime per group",
      subtitle = sprintf("Parental low score (n=%d) · Parental high score (n=%d) · Resistant (n=%d)",
                         as.integer(n_3grp["Parental low score"]),
                         as.integer(n_3grp["Parental high score"]),
                         as.integer(n_3grp["Resistant"])),
      x = NULL, y = "Pseudotime"
    ) +
    theme_classic(base_size = 13) + white_bg + NoLegend()

  ggsave(pfig("traj_05_preexistence_pseudotime.png"),
         p7, width = 9, height = 7, dpi = 300)
  cat("  traj_05_preexistence_pseudotime.png\n")
} else {
  cat("  traj_05: sc_05 scores not available — figure skipped\n")
  p7_pval_ab <- NA_real_
  p7_pval_ar <- NA_real_
}

# ==============================================================================
# 6. DRIVER-GENE FIGURE ALONG THE TRAJECTORY (no graph_test)
# ==============================================================================

cat("\n[Sec 6] Driver genes along the trajectory (no graph_test)...\n")

genes_present <- intersect(GENES_PSEUDOTIME, rownames(seurat_obj))

if (length(genes_present) >= 2) {
  p6a <- plot_cells(
    cds,
    genes                 = genes_present,
    show_trajectory_graph = FALSE,
    label_cell_groups     = FALSE,
    cell_size             = 0.4
  ) +
    scale_color_viridis_c(option = "magma", begin = 0.1, end = 0.95,
                          name = "Expr\n(% Max)") +
    labs(title = sprintf("Driver genes (%s)", paste(GENES_PSEUDOTIME, collapse = " / "))) +
    theme_classic(base_size = 9) + white_bg

  n_rows_driver <- ceiling(length(genes_present) / 4)

  ggsave(pfig("traj_06_genes_pseudotime.png"),
         p6a, width = 16, height = n_rows_driver * 4, dpi = 300, bg = "white")
  cat("  traj_06_genes_pseudotime.png (driver panel only, no graph_test top genes)\n")
} else {
  cat("  traj_06: driver genes not found — figure skipped\n")
}

# ==============================================================================
# 7. TABLAS
# ==============================================================================

cat("\n[Sec 7] Saving tables...\n")

pt_export <- data.frame(
  barcode          = colnames(seurat_obj),
  sample           = seurat_obj$sample,
  condition        = seurat_obj$condition,
  cluster          = as.character(seurat_obj$seurat_clusters),
  pseudotime       = seurat_obj$pseudotime,
  ucell_score      = if (has_score)    seurat_obj@meta.data[[ucell_col]]          else NA_real_,
  ucell_combined   = if (has_combined) seurat_obj@meta.data[[ucell_combined_col]] else NA_real_
)
write.csv(pt_export,
          file.path(tab_dir, "traj_pseudotime_per_cell.csv"),
          row.names = FALSE)
cat("  traj_pseudotime_per_cell.csv\n")
cat("  (traj_graph_test_results.csv NOT generated — use sc_06_monocle3.R for that)\n")

# ==============================================================================
# 8. SAVE OBJECT
# ==============================================================================

cat("\n[Sec 8] Saving object with pseudotime...\n")
out_rds <- file.path(output_dir, sprintf("%s_trajectory.rds", GROUP))
saveRDS(seurat_obj, out_rds)
cat(sprintf("  %s_trajectory.rds saved\n", GROUP))

# ==============================================================================
# 9. REPORT
# ==============================================================================

cat("\n[Sec 9] Generating report...\n")

n_cells_pt     <- sum(!is.na(seurat_obj$pseudotime))
n_cells_pt_par <- sum(!is.na(seurat_obj$pseudotime) & seurat_obj$condition == COND_PARENTAL)
n_cells_pt_res <- sum(!is.na(seurat_obj$pseudotime) & seurat_obj$condition == COND_RESISTANT)
pt_med_par     <- median(pt_df$pseudotime[pt_df$condition == COND_PARENTAL],    na.rm = TRUE)
pt_med_res     <- median(pt_df$pseudotime[pt_df$condition == COND_RESISTANT],   na.rm = TRUE)

report_lines <- c(
  sprintf("# sc_06b_monocle3_no_graphtest — %s: %s vs %s", GROUP, COND_PARENTAL, COND_RESISTANT),
  "",
  sprintf("**Fecha:** %s", Sys.Date()),
  "",
  "> Report generated without graph_test (Moran's I) — see sc_06_monocle3.R for that part.",
  "",
  "## Trajectory",
  sprintf("- Cells with valid pseudotime: %d / %d", n_cells_pt, ncol(seurat_obj)),
  sprintf("- %s cells with valid pseudotime: %d", COND_PARENTAL,  n_cells_pt_par),
  sprintf("- %s cells with valid pseudotime: %d", COND_RESISTANT, n_cells_pt_res),
  sprintf("- Pseudotime mediana %s: %.2f", COND_PARENTAL,  pt_med_par),
  sprintf("- Pseudotime mediana %s: %.2f", COND_RESISTANT, pt_med_res),
  sprintf("- Root: %s", root_label),
  "",
  "## Pre-existence analysis",
  if (!is.na(p7_pval_ab))
    sprintf("- Wilcoxon Parental high vs Parental low score: p = %.2e", p7_pval_ab)
  else
    "- Wilcoxon high vs low score: not computed (sc_05 scores not available)",
  if (!is.na(p7_pval_ar))
    sprintf("- Wilcoxon Parental high score vs Resistant: p = %.2e", p7_pval_ar)
  else
    "- Wilcoxon high score vs Resistant: not computed",
  "",
  "## Figures generated",
  "| # | File | Description |",
  "|---|---------|-------------|",
  "| 1 | traj_01_graph_umap.png | UMAP with trajectory graph (clusters) |",
  "| 2 | traj_02_pseudotime_split.png | Pseudotime split by sample |",
  "| 3 | traj_03_pseudotime_all.png | Pseudotime, all cells |",
  "| 4 | traj_04_pseudotime_by_sample.png | Pseudotime distribution per sample (ridge/density) |",
  "| 5 | traj_05_preexistence_pseudotime.png | 3 groups: Parental low/high score + Resistant (Wilcoxon) |",
  "| 6 | traj_06_genes_pseudotime.png | Driver genes only (no graph_test top genes) |",
  "",
  "## Next step",
  "- If traj_graph_test_results.csv or the top-genes-by-Moran's-I panel is needed, run sc_06_monocle3.R with RERUN_VIZ_ONLY <- FALSE",
  "- To run the other group set `GROUP <- 'EGFR'` (needs its own data and signature)"
)

writeLines(report_lines,
           file.path(rep_dir, "sc_06_monocle3_report.md"))
cat("  sc_06_monocle3_report.md\n")

cat("\n=== sc_06b_monocle3_no_graphtest.R completed ===\n")
