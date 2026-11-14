# ==============================================================================
# sc_05_signature.R — Projection of the resistance signature (UCell)
# MET: MET1 vs MET1R  |  EGFR: EGFR1 vs EGFR1_R  |  change GROUP in CONFIG
# Mayo 2026
# ==============================================================================
#
# Run (from the repo root):
#   source("scripts/sc_05_signature.R")
#   Rscript scripts/sc_05_signature.R
#
# Outputs (MET):
#   figures/ : sig_01_umap_score_UCell.png
#              sig_02_violin_UCell.png
#              sig_04_gene_coverage.png
#              sig_05_high_score_umap.png
#              sig_06_preexistence_barplot.png
#              sig_07_signature_heatmap_cluster.png
#              sig_08_signature_heatmap_cells.png
#              sig_08b_signature_heatmap_sample.png
#              sig_09_featureplot_key_genes.png
#              sig_10_dotplot_key_genes.png
#              sig_11_ucell_landscape_2d.png
#              sig_12_density_combined.png
#   tables/  : sig_scores_per_cell.csv
#              sig_preexistence_summary.csv
#              sig_gene_coverage.csv
#   reports/ : sc_05_signature_report.md
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggridges)
  library(pheatmap)
  library(tibble)
})

setwd(here::here())   # repo root (set by the .here marker)

# ==============================================================================
# CONFIGURATION — edit only this section before running
# ==============================================================================

GROUP          <- "MET"   # MET | EGFR
HIGH_SCORE_PCTL <- 0.90   # percentile used to define "high-score cells"
N_SIG_HEATMAP   <- 100    # top genes for the signature heatmap (60 up + 40 down)
RERUN_VIZ_ONLY  <- FALSE  # TRUE = load the checkpoint and regenerate figures only

# ==============================================================================

# Load config first — SAMPLE_PARENTAL and the other variables must exist
# before validating file paths
source("config/colors.R")
source("config/pipeline_config.R")

# Load the signature CSV (gene_lists.R refers to it, so it must be loaded before that file)
sig_dir      <- file.path(GROUP, "data", "signatures")
sig_csv      <- file.path(sig_dir, sprintf("%s_resistance_signature_conservative.csv", GROUP))
sig_csv_full <- file.path(sig_dir, sprintf("%s_resistance_signature.csv", GROUP))
if (!file.exists(sig_csv)) stop(
  "Signature not found: ", sig_csv,
  " — copy it from the bulk RNA-seq results (Gene and Classification columns, see data/README)"
)
sig_df      <- read.csv(sig_csv)
sig_df_full <- if (file.exists(sig_csv_full)) read.csv(sig_csv_full) else NULL
cat(sprintf("Conservative signature: %d genes (%s)\n", nrow(sig_df), basename(sig_csv)))
if (!is.null(sig_df_full))
  cat(sprintf("Full signature:         %d genes (%s)\n", nrow(sig_df_full), basename(sig_csv_full)))

# Signature variables with dynamic names — gene_lists.R picks them up via get0()
sig_key <- tolower(GROUP)
assign(paste0(sig_key, "_resistance_signature_721"), sig_df$Gene)
assign(paste0(sig_key, "_resistance_signature_up"),  sig_df$Gene[sig_df$Classification == "Upregulated"])
assign(paste0(sig_key, "_resistance_signature_dn"),  sig_df$Gene[sig_df$Classification == "Downregulated"])

source("config/gene_lists.R")
set.seed(42)

output_dir <- file.path(".", GROUP, "results")
fig_dir    <- file.path(output_dir, "figures")
tab_dir    <- file.path(output_dir, "tables")
rep_dir    <- file.path(output_dir, "reports")
pfig <- function(f) file.path(fig_dir, paste0("sc_05_", f))
for (d in c(fig_dir, tab_dir, rep_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


# Score column names — consistent between full run and VIZ_ONLY
ucell_up_key <- paste0(GROUP, "_resistance_up")
ucell_dn_key <- paste0(GROUP, "_resistance_dn")
UCELL_COL    <- paste0(ucell_up_key, "_UCell")
UCELL_DN_COL <- paste0(ucell_dn_key, "_UCell")

# ==============================================================================
# 1/3/4. LOAD OBJECT AND COMPUTE SCORES, OR LOAD FROM CHECKPOINT
# ==============================================================================

if (!RERUN_VIZ_ONLY) {

  cat("\n[Sec 1] Loading object...\n")
  rds_file <- file.path(output_dir, sprintf("%s_clustered.rds", GROUP))
  if (!file.exists(rds_file)) stop("Missing: ", rds_file, " — run sc_03 first")
  merged <- readRDS(rds_file)
  gc()
  Idents(merged) <- "seurat_clusters"
  CLUSTERS_ALL <- levels(merged$seurat_clusters)
  cat(sprintf("  Cells: %d | Clusters: %s\n",
              ncol(merged), paste(CLUSTERS_ALL, collapse = ", ")))

  # ── 2. Filter signature genes ──────────────────────────────────────────────
  cat("\n[Sec 2] Filtering signature genes in the object...\n")
  genes_up_all     <- get(paste0(tolower(GROUP), "_resistance_signature_up"))
  genes_dn_all     <- get(paste0(tolower(GROUP), "_resistance_signature_dn"))
  genes_up_present <- intersect(genes_up_all, rownames(merged))
  genes_dn_present <- intersect(genes_dn_all, rownames(merged))

  cat(sprintf("  Upregulated:   %d / %d capturados (%.1f%%)\n",
              length(genes_up_present), length(genes_up_all),
              100 * length(genes_up_present) / max(length(genes_up_all), 1)))
  cat(sprintf("  Downregulated: %d / %d capturados (%.1f%%)\n",
              length(genes_dn_present), length(genes_dn_all),
              100 * length(genes_dn_present) / max(length(genes_dn_all), 1)))

  if (length(genes_up_present) < 10)
    warning("Fewer than 10 upregulated genes captured — check the matrix")

  # ── 3. UCell ───────────────────────────────────────────────────────────────
  cat("\n[Sec 3] UCell...\n")
  gc()
  ucell_features <- list(genes_up_present)
  names(ucell_features) <- ucell_up_key
  if (length(genes_dn_present) >= 5)
    ucell_features[[ucell_dn_key]] <- genes_dn_present

  DefaultAssay(merged) <- "RNA"
  merged <- AddModuleScore_UCell(merged, features = ucell_features, name = "_UCell")
  DefaultAssay(merged) <- "SCT"
  cat(sprintf("  UCell scores: %s\n", paste(names(ucell_features), collapse = ", ")))

  # Combined UCell score (up − down): captures full bidirectional signature
  # Pre-existing parental cells should have high UCell_up AND low UCell_dn
  if (all(c(UCELL_COL, UCELL_DN_COL) %in% colnames(merged@meta.data))) {
    merged[[paste0(GROUP, "_UCell_combined")]] <-
      merged@meta.data[[UCELL_COL]] - merged@meta.data[[UCELL_DN_COL]]
    cat(sprintf("  UCell combined (up - down) → %s_UCell_combined\n", GROUP))
  }

} else {

  sig_rds <- file.path(output_dir, sprintf("%s_signature.rds", GROUP))
  if (!file.exists(sig_rds)) stop(
    "Checkpoint not found: ", sig_rds,
    " — run with RERUN_VIZ_ONLY=FALSE first"
  )
  merged <- readRDS(sig_rds)
  Idents(merged) <- "seurat_clusters"
  CLUSTERS_ALL <- levels(merged$seurat_clusters)
  cat(sprintf("[RERUN_VIZ_ONLY] %s_signature.rds loaded: %d cells\n",
              GROUP, ncol(merged)))

  # Rebuild genes_up_all/dn_all for coverage (sig_05)
  genes_up_all     <- get(paste0(tolower(GROUP), "_resistance_signature_up"))
  genes_dn_all     <- get(paste0(tolower(GROUP), "_resistance_signature_dn"))
  genes_up_present <- intersect(genes_up_all, rownames(merged))
  genes_dn_present <- intersect(genes_dn_all, rownames(merged))

  # Recompute combined UCell if it is not in the checkpoint
  ucell_combined_col <- paste0(GROUP, "_UCell_combined")
  if (!ucell_combined_col %in% colnames(merged@meta.data) &&
      all(c(UCELL_COL, UCELL_DN_COL) %in% colnames(merged@meta.data))) {
    merged[[ucell_combined_col]] <- merged@meta.data[[UCELL_COL]] - merged@meta.data[[UCELL_DN_COL]]
    cat(sprintf("  UCell combined recalculado → %s\n", ucell_combined_col))
  }

}

# Global variable for combined UCell (available in both if/else branches)
UCELL_COMBINED_COL <- paste0(GROUP, "_UCell_combined")
if (!UCELL_COMBINED_COL %in% colnames(merged@meta.data)) UCELL_COMBINED_COL <- NULL

# ==============================================================================
# 5. VISUALIZACIONES
# ==============================================================================

cat("\n[Sec 5] Visualizaciones...\n")

# Signature genes ordered by LFC (used by heatmaps and featureplots)
sig_up_ordered <- sig_df$Gene[sig_df$Classification == "Upregulated"]
sig_up_ordered <- sig_up_ordered[order(
  sig_df$LFC_Res_vs_Par[sig_df$Classification == "Upregulated"],
  decreasing = TRUE
)]
sig_dn_ordered <- sig_df$Gene[sig_df$Classification == "Downregulated"]
sig_dn_ordered <- sig_dn_ordered[order(
  sig_df$LFC_Res_vs_Par[sig_df$Classification == "Downregulated"],
  decreasing = FALSE
)]

# Score dataframe and thresholds (shared across figures)
score_df <- data.frame(
  barcode   = colnames(merged),
  condition = merged$condition,
  sample    = merged$sample,
  cluster   = as.character(merged$seurat_clusters),
  UCell     = merged@meta.data[[UCELL_COL]],
  UCell_dn  = if (UCELL_DN_COL %in% colnames(merged@meta.data))
                merged@meta.data[[UCELL_DN_COL]] else NA_real_
)
if (!is.null(UCELL_COMBINED_COL))
  score_df$UCell_combined <- merged@meta.data[[UCELL_COMBINED_COL]]

threshold_ucell          <- quantile(score_df$UCell[score_df$condition == COND_PARENTAL], HIGH_SCORE_PCTL)
threshold_ucell_combined <- if (!is.null(UCELL_COMBINED_COL) && "UCell_combined" %in% names(score_df))
  quantile(score_df$UCell_combined[score_df$condition == COND_PARENTAL], HIGH_SCORE_PCTL) else NA_real_

# --- sig_01: FeaturePlot UCell — split by sample (MET1 / MET1R) ---
p_ucell_fp <- FeaturePlot(
  merged, features = UCELL_COL,
  reduction = "umap.harmony",
  split.by  = "sample",
  cols      = score_gradient,
  pt.size   = 0.3, order = TRUE
) &
  labs(title = sprintf("%s resistance signature — UCell score\n%s  |  %s",
                       GROUP, SAMPLE_PARENTAL, SAMPLE_RESISTANT)) &
  theme(plot.title = element_text(face = "bold")) & white_bg

ggsave(pfig("sig_01_umap_score_UCell.png"),
       p_ucell_fp, width = 12, height = 5, dpi = 800)
cat("  sig_01_umap_score_UCell.png\n")

# --- sig_02: VlnPlot UCell by condition ---
p_v1 <- VlnPlot(
  merged, features = UCELL_COL,
  group.by = "condition", cols = condition_colors, pt.size = 0
) + labs(title = "UCell resistance score", subtitle = "Rank-based · dropout-robust",
         y = "Score") + NoLegend() + white_bg

ggsave(pfig("sig_02_violin_UCell.png"),
       p_v1, width = 7, height = 6, dpi = 800)
cat("  sig_02_violin_UCell.png\n")

# --- sig_04: Gene coverage bulk → sc ---
cov_df <- data.frame(
  Direction = c("Upregulated", "Downregulated"),
  Total     = c(length(genes_up_all),     length(genes_dn_all)),
  Captured  = c(length(genes_up_present), length(genes_dn_present))
) |>
  dplyr::mutate(
    Pct   = round(Captured / Total * 100, 1),
    Label = sprintf("%d / %d\n(%.1f%%)", Captured, Total, Pct)
  )

p_cov <- ggplot(cov_df, aes(x = Direction, y = Pct, fill = Direction)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = Label), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = coverage_colors) +
  scale_y_continuous(limits = c(0, 115)) +
  labs(
    title    = sprintf("Gene coverage %s signature — bulk RNA-seq → scRNA-seq", GROUP),
    subtitle = "Conservative signature genes detected in 10x matrix",
    x = NULL, y = "% genes detected"
  ) +
  theme_classic(base_size = 13) + white_bg + NoLegend()

ggsave(pfig("sig_04_gene_coverage.png"),
       p_cov, width = 7, height = 6, dpi = 800)
cat("  sig_04_gene_coverage.png\n")

# --- sig_05: UMAP preexistence — English labels, distinct colors ---
hs_high   <- "Parental high score"
hs_low    <- "Parental low score"
hs_resist <- "Resistant"

# Primary threshold: combined (up - down) if available, else UCell_up
primary_score_col <- if (!is.null(UCELL_COMBINED_COL) && !is.na(threshold_ucell_combined))
  "UCell_combined" else "UCell"
primary_threshold <- if (primary_score_col == "UCell_combined") threshold_ucell_combined else threshold_ucell
primary_label     <- if (primary_score_col == "UCell_combined") "UCell combined (up−down)" else "UCell up"

score_df$high_score_parental <- dplyr::case_when(
  score_df$condition == COND_PARENTAL & score_df[[primary_score_col]] > primary_threshold ~ hs_high,
  score_df$condition == COND_PARENTAL                                                      ~ hs_low,
  TRUE                                                                                     ~ hs_resist
)
merged$high_score_parental <- score_df$high_score_parental

hs_colors <- c(
  "Parental high score" = "#013c58",   # dark navy
  "Parental low score"  = "#b2bec3",   # silver grey
  "Resistant"           = "#e67e22"    # burnt orange
)

p_hs <- DimPlot(
  merged, reduction = "umap.harmony",
  group.by = "high_score_parental",
  cols = hs_colors, pt.size = 0.5, order = TRUE
) +
  labs(
    title    = sprintf("Preexistence: %s cells with high resistance score", COND_PARENTAL),
    subtitle = sprintf("%s > p%d (%s) = %.5f",
                       primary_label, round(HIGH_SCORE_PCTL * 100), COND_PARENTAL, primary_threshold)
  ) +
  theme_classic(base_size = 12) + white_bg +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(pfig("sig_05_high_score_umap.png"),
       p_hs, width = 10, height = 7, dpi = 800)
cat("  sig_05_high_score_umap.png\n")

# --- sig_06: Barplot % high-score cells per cluster ---
pct_high_cluster <- score_df %>%
  dplyr::group_by(cluster, sample) %>%
  dplyr::summarise(
    n_total  = dplyr::n(),
    n_high   = sum(.data[[primary_score_col]] > primary_threshold),
    pct_high = 100 * n_high / n_total,
    .groups  = "drop"
  ) %>%
  dplyr::mutate(cluster = factor(cluster, levels = sort(unique(as.integer(cluster)))))

p_bar <- ggplot(pct_high_cluster,
                aes(x = cluster, y = pct_high, fill = sample)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", pct_high)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = sprintf("Preexistence by cluster — %% high-score cells (%s > p%d)",
                       primary_label, round(HIGH_SCORE_PCTL * 100)),
    subtitle = sprintf("Threshold p%d (%s) = %.5f",
                       round(HIGH_SCORE_PCTL * 100), COND_PARENTAL, primary_threshold),
    x = "Cluster", y = "% high-score cells",
    fill = NULL
  ) +
  theme_classic(base_size = 13) + white_bg +
  theme(legend.position = "top")

ggsave(pfig("sig_06_preexistence_barplot.png"),
       p_bar, width = 10, height = 6, dpi = 800)
cat("  sig_06_preexistence_barplot.png\n")

# --- sig_09: FeaturePlot key resistance genes ---
# Biologically relevant genes: driver, resistance mechanisms, downstream pathway
# targets and cross-talk partners, as defined in GENES_KEY_SIGNATURE
# (config/pipeline_config.R)
key_genes_fp <- GENES_KEY_SIGNATURE[[GROUP]]
key_genes_fp <- key_genes_fp[key_genes_fp %in% rownames(merged)]

if (length(key_genes_fp) > 0) {
  nc_fp <- min(5, length(key_genes_fp))
  nr_fp <- ceiling(length(key_genes_fp) / nc_fp)
  p_key <- FeaturePlot(
    merged, features = key_genes_fp,
    reduction = "umap.harmony",
    cols      = expr_gradient,
    pt.size   = 0.2, ncol = nc_fp, order = TRUE
  ) & white_bg
  ggsave(pfig("sig_09_featureplot_key_genes.png"),
         p_key, width = nc_fp * 5, height = nr_fp * 5, dpi = 800)
  cat("  sig_09_featureplot_key_genes.png\n")
}

# --- sig_10: DotPlot key resistance genes by cluster ---
if (length(key_genes_fp) > 0) {
  p_dot <- DotPlot(
    merged, features = key_genes_fp,
    group.by  = "seurat_clusters",
    assay     = "SCT",
    cols      = expr_gradient,
    dot.scale = 6
  ) +
    labs(title = sprintf("Key resistance genes — %s vs %s", SAMPLE_PARENTAL, SAMPLE_RESISTANT),
         x = NULL, y = "Cluster") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
    white_bg
  ggsave(pfig("sig_10_dotplot_key_genes.png"),
         p_dot, width = 12, height = 7, dpi = 800)
  cat("  sig_10_dotplot_key_genes.png\n")
}

# --- sig_11: 2D landscape UCell_up vs UCell_down — trajectory of states ---
if (!is.null(UCELL_COMBINED_COL) && !all(is.na(score_df$UCell_dn))) {

  # Sample down to max 5000 pts per condition for readability
  set.seed(42)
  score_df_sub12 <- score_df %>%
    dplyr::group_by(condition) %>%
    dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(5000L, nrow(.x)))) %>%
    dplyr::ungroup()

  # Highlight parental high-combined cells
  score_df_sub12$group <- dplyr::case_when(
    score_df_sub12$condition == COND_PARENTAL & score_df_sub12$UCell_combined > threshold_ucell_combined
      ~ "Parental high (preexistent)",
    score_df_sub12$condition == COND_PARENTAL ~ "Parental low",
    TRUE ~ "Resistant"
  )
  group_colors12 <- c(
    "Parental high (preexistent)" = "#013c58",
    "Parental low"                = "#b2bec3",
    "Resistant"                   = "#e67e22"
  )

  p12 <- ggplot(score_df_sub12, aes(x = UCell, y = UCell_dn, color = group)) +
    geom_point(
      data  = dplyr::filter(score_df_sub12, group != "Parental high (preexistent)"),
      alpha = 0.15, size = 0.5
    ) +
    geom_point(
      data  = dplyr::filter(score_df_sub12, group == "Parental high (preexistent)"),
      alpha = 0.6, size = 0.8
    ) +
    scale_color_manual(values = group_colors12) +
    labs(
      title    = sprintf("%s signature — 2D UCell landscape", GROUP),
      subtitle = "UP score (X) vs DOWN score (Y) · Parental high = UCell_combined > p90",
      x        = sprintf("UCell score — UP genes (%d detected)", length(genes_up_present)),
      y        = sprintf("UCell score — DOWN genes (%d detected)", length(genes_dn_present)),
      color    = NULL
    ) +
    theme_classic(base_size = 13) + white_bg +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8)))

  ggsave(pfig("sig_11_ucell_landscape_2d.png"),
         p12, width = 10, height = 7, dpi = 800)
  cat("  sig_11_ucell_landscape_2d.png\n")
}

# --- sig_12: Density bidireccional (score combinado) ---
if (!is.null(UCELL_COMBINED_COL) && "UCell_combined" %in% names(score_df)) {
  p13 <- ggplot(score_df, aes(x = UCell_combined, y = condition, fill = condition)) +
    geom_density_ridges(alpha = 0.75, scale = 1.2, rel_min_height = 0.01) +
    geom_vline(xintercept = threshold_ucell_combined, linetype = "dashed",
               color = "grey30", linewidth = 0.8) +
    scale_fill_manual(values = condition_colors) +
    labs(
      title    = sprintf("%s resistance score — bidireccional (UCell up − down)", GROUP),
      subtitle = sprintf("Combined score · dashed = p%d (%s): %.5f",
                         round(HIGH_SCORE_PCTL * 100), COND_PARENTAL, threshold_ucell_combined),
      x = "UCell Combined Score (up − down)", y = NULL
    ) +
    theme_classic(base_size = 13) + white_bg + NoLegend()

  ggsave(pfig("sig_12_density_combined.png"),
         p13, width = 10, height = 6, dpi = 800)
  cat("  sig_12_density_combined.png\n")
}

# ==============================================================================
# 5b. HEATMAPS
# ==============================================================================

cat("\n[Sec 5b] Heatmaps...\n")

cluster_comp  <- table(merged$seurat_clusters, merged$sample)
pct_a         <- round(cluster_comp[, SAMPLE_PARENTAL]  / rowSums(cluster_comp) * 100, 1)
pct_b         <- round(cluster_comp[, SAMPLE_RESISTANT] / rowSums(cluster_comp) * 100, 1)
n_cells_cl    <- rowSums(cluster_comp)
ucell_by_cl   <- tapply(merged@meta.data[[UCELL_COL]], merged$seurat_clusters, mean)
ucell_by_cl   <- round(ucell_by_cl[names(pct_a)], 4)

ann_col_base  <- data.frame(
  N_cells        = as.integer(n_cells_cl),
  Mean_UCell     = as.numeric(ucell_by_cl),
  Pct_Parental   = as.numeric(pct_a),
  Pct_Resistant  = as.numeric(pct_b),
  row.names      = names(pct_a)
)

n_up_heat <- round(N_SIG_HEATMAP * 0.6)
n_dn_heat <- N_SIG_HEATMAP - n_up_heat

top_up_present <- head(intersect(sig_up_ordered, rownames(merged)), n_up_heat)
top_dn_present <- head(intersect(sig_dn_ordered, rownames(merged)), n_dn_heat)

# Force driver genes into the gene list regardless of rank
for (g in GENES_FORCE_HEATMAP[[GROUP]]) {
  if (g %in% rownames(merged) && !g %in% c(top_up_present, top_dn_present)) {
    if (g %in% sig_up_ordered)      top_up_present <- c(top_up_present, g)
    else if (g %in% sig_dn_ordered) top_dn_present <- c(top_dn_present, g)
    else                            top_up_present <- c(top_up_present, g)
  }
}
top_sig_present <- c(top_up_present, top_dn_present)

row_ann <- data.frame(
  Direction = c(rep("Upregulated", length(top_up_present)),
                rep("Downregulated", length(top_dn_present))),
  row.names = top_sig_present
)
row_ann_colors <- list(Direction = c(Upregulated = "#E74C3C", Downregulated = "#3498DB"))

# --- sig_07: pheatmap top genes × cluster averages ---
if (length(top_sig_present) >= 5) {
  avg_mat <- AverageExpression(
    merged,
    features = top_sig_present,
    group.by = "seurat_clusters",
    assays   = "SCT"
  )$SCT

  avg_scaled <- t(scale(t(as.matrix(avg_mat))))
  avg_scaled[!is.finite(avg_scaled)] <- 0

  colnames(avg_scaled) <- sub("^g", "", colnames(avg_scaled))
  ann_col8 <- ann_col_base[colnames(avg_scaled), , drop = FALSE]

  sc_max <- max(abs(avg_scaled), na.rm = TRUE)
  if (sc_max == 0) sc_max <- 1
  sc_breaks <- seq(-sc_max, sc_max, length.out = 101)

  p8 <- pheatmap(
    avg_scaled,
    color             = heatmap_zscore,
    breaks            = sc_breaks,
    annotation_col    = ann_col8,
    annotation_row    = row_ann,
    annotation_colors = row_ann_colors,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    show_rownames     = TRUE,
    fontsize_row      = 8,
    fontsize_col      = 10,
    main = sprintf("Top %d signature genes (up=%d, down=%d) × clusters — z-score SCT",
                   length(top_sig_present), length(top_up_present), length(top_dn_present)),
    silent = TRUE
  )
  ggsave(pfig("sig_07_signature_heatmap_cluster.png"),
         plot = p8$gtable, width = 14, height = 26, dpi = 800, bg = "white")
  cat("  sig_07_signature_heatmap_cluster.png\n")

  # --- sig_08: DoHeatmap individual cells (~100 per cluster) ---
  set.seed(42)
  cells_sub <- merged@meta.data %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::slice_sample(n = 100, replace = FALSE) %>%
    dplyr::pull(barcode)

  p9 <- DoHeatmap(
    merged[, cells_sub],
    features     = top_sig_present,
    group.by     = "seurat_clusters",
    group.colors = cluster_colors,
    disp.min     = -2.5,
    disp.max     =  2.5,
    slot         = "scale.data",
    assay        = "SCT"
  ) +
    scale_fill_gradientn(
      colors   = c("#3498DB", "white", "#E74C3C"),
      na.value = "white",
      name     = "z-score"
    ) +
    labs(title = sprintf("%s signature — top %d genes · ~100 cells/cluster",
                         GROUP, length(top_sig_present))) +
    white_bg +
    theme(axis.text.y = element_text(size = 6))

  ggsave(pfig("sig_08_signature_heatmap_cells.png"),
         p9, width = 16, height = 22, dpi = 800, bg = "white")
  cat("  sig_08_signature_heatmap_cells.png\n")

  # --- sig_08b: DoHeatmap individual cells — grouped by sample (MET1 vs MET1R) ---
  set.seed(42)
  cells_sub_s <- merged@meta.data %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::group_by(sample) %>%
    dplyr::slice_sample(n = 500, replace = FALSE) %>%
    dplyr::pull(barcode)

  p9b <- DoHeatmap(
    merged[, cells_sub_s],
    features     = top_sig_present,
    group.by     = "sample",
    group.colors = sample_colors,
    disp.min     = -2.5,
    disp.max     =  2.5,
    slot         = "scale.data",
    assay        = "SCT"
  ) +
    scale_fill_gradientn(
      colors   = c("#3498DB", "white", "#E74C3C"),
      na.value = "white",
      name     = "z-score"
    ) +
    labs(title = sprintf("%s signature — top %d genes · ~500 cells/sample",
                         GROUP, length(top_sig_present))) +
    white_bg +
    theme(axis.text.y = element_text(size = 6))

  ggsave(pfig("sig_08b_signature_heatmap_sample.png"),
         p9b, width = 12, height = 22, dpi = 800, bg = "white")
  cat("  sig_08b_signature_heatmap_sample.png\n")

} else {
  cat("  Warning: <5 genes present in signature — sig_07/08 skipped\n")
}

# ==============================================================================
# 6. STATISTICAL TEST OF PRE-EXISTENCE
# ==============================================================================

cat("\n[Sec 6] Statistical test of pre-existence...\n")

wt_ucell <- wilcox.test(UCell ~ condition, data = score_df, exact = FALSE)

n_high_parental   <- sum(score_df$condition == COND_PARENTAL  & score_df[[primary_score_col]] > primary_threshold)
n_total_parental  <- sum(score_df$condition == COND_PARENTAL)
n_high_resistant  <- sum(score_df$condition == COND_RESISTANT & score_df[[primary_score_col]] > primary_threshold)
n_total_resistant <- sum(score_df$condition == COND_RESISTANT)

cat(sprintf("  UCell Wilcoxon p = %.2e\n", wt_ucell$p.value))
cat(sprintf("  High score %s:  %d / %d (%.1f%%)\n",
            COND_PARENTAL,  n_high_parental,  n_total_parental,
            100 * n_high_parental  / n_total_parental))
cat(sprintf("  High score %s: %d / %d (%.1f%%)\n",
            COND_RESISTANT, n_high_resistant, n_total_resistant,
            100 * n_high_resistant / n_total_resistant))

# ==============================================================================
# 7. TABLAS
# ==============================================================================

cat("\n[Sec 7] Saving tables...\n")

score_df$is_high_UCell    <- score_df$UCell > threshold_ucell
if ("UCell_combined" %in% names(score_df) && !is.na(threshold_ucell_combined))
  score_df$is_high_combined <- score_df$UCell_combined > threshold_ucell_combined
write.csv(score_df,
          file.path(tab_dir, "sig_scores_per_cell.csv"),
          row.names = FALSE)
cat("  sig_scores_per_cell.csv\n")

preexist_summary <- data.frame(
  Metric = c(
    "N cells parental", "N cells resistant",
    "N high-score parental", "N high-score resistant",
    "Pct high-score parental (%)", "Pct high-score resistant (%)",
    "UCell Wilcoxon p-value",
    "UCell threshold (pctl90)"
  ),
  Value = c(
    n_total_parental, n_total_resistant,
    n_high_parental,  n_high_resistant,
    round(100 * n_high_parental  / n_total_parental,  2),
    round(100 * n_high_resistant / n_total_resistant, 2),
    formatC(wt_ucell$p.value, format = "e", digits = 3),
    round(threshold_ucell, 5)
  )
)
write.csv(preexist_summary,
          file.path(tab_dir, "sig_preexistence_summary.csv"),
          row.names = FALSE)
cat("  sig_preexistence_summary.csv\n")

gene_coverage <- data.frame(
  Direction    = c("Upregulated", "Downregulated"),
  Total_bulk   = c(length(genes_up_all),     length(genes_dn_all)),
  Captured_sc  = c(length(genes_up_present), length(genes_dn_present)),
  Pct_captured = round(c(
    100 * length(genes_up_present) / max(length(genes_up_all), 1),
    100 * length(genes_dn_present) / max(length(genes_dn_all), 1)
  ), 1)
)
write.csv(gene_coverage,
          file.path(tab_dir, "sig_gene_coverage.csv"),
          row.names = FALSE)
cat("  sig_gene_coverage.csv\n")

# ==============================================================================
# 8. SAVE OBJECT
# ==============================================================================

cat("\n[Sec 8] Saving object with scores...\n")
out_rds <- file.path(output_dir, sprintf("%s_signature.rds", GROUP))
saveRDS(merged, out_rds)
cat(sprintf("  %s_signature.rds saved\n", GROUP))

# ==============================================================================
# 9. REPORT
# ==============================================================================

cat("\n[Sec 9] Generating report...\n")

report_lines <- c(
  sprintf("# sc_05_signature — %s: %s vs %s", GROUP, COND_PARENTAL, COND_RESISTANT),
  "",
  sprintf("**Fecha:** %s", Sys.Date()),
  "",
  "## Resistance signature",
  sprintf("- File: `%s`", basename(sig_csv)),
  sprintf("- Total genes: %d (up: %d, down: %d)",
          nrow(sig_df), length(genes_up_all), length(genes_dn_all)),
  sprintf("- Captured in scRNA-seq: up %d/%d (%.1f%%), down %d/%d (%.1f%%)",
          length(genes_up_present), length(genes_up_all),
          100 * length(genes_up_present) / max(length(genes_up_all), 1),
          length(genes_dn_present), length(genes_dn_all),
          100 * length(genes_dn_present) / max(length(genes_dn_all), 1)),
  "",
  "## Pre-existence results",
  sprintf("- %s cells: %d total | %d high score (%.1f%%, %s > p%d = %.5f)",
          COND_PARENTAL, n_total_parental, n_high_parental,
          100 * n_high_parental / n_total_parental,
          primary_label, round(HIGH_SCORE_PCTL * 100), primary_threshold),
  sprintf("- %s cells: %d total | %d high score (%.1f%%)",
          COND_RESISTANT, n_total_resistant, n_high_resistant,
          100 * n_high_resistant / n_total_resistant),
  sprintf("- Wilcoxon UCell: p = %.2e", wt_ucell$p.value),
  "",
  "## Figures generated",
  "| # | File | Description |",
  "|---|---------|-------------|",
  "| 1  | sig_01_umap_score_UCell.png | UMAP UCell score — split by sample |",
  "| 2  | sig_02_violin_UCell.png | VlnPlot UCell score by condition |",
  "| 4  | sig_04_gene_coverage.png | Gene coverage bulk → sc |",
  "| 5  | sig_05_high_score_umap.png | UMAP preexistence (Parental high/low · Resistant) |",
  "| 6  | sig_06_preexistence_barplot.png | % high-score cells per cluster |",
  "| 7  | sig_07_signature_heatmap_cluster.png | Top signature genes × clusters (z-score) |",
  "| 8  | sig_08_signature_heatmap_cells.png | DoHeatmap individual cells — per cluster |",
  "| 8b | sig_08b_signature_heatmap_sample.png | DoHeatmap individual cells — per sample (MET1 vs MET1R) |",
  "| 9  | sig_09_featureplot_key_genes.png | FeaturePlot of the key genes + top signature genes |",
  "| 10 | sig_10_dotplot_key_genes.png | DotPlot key genes by cluster |",
  "| 11 | sig_11_ucell_landscape_2d.png | 2D landscape UCell_up vs UCell_down — trajectory of states |",
  "| 12 | sig_12_density_combined.png | Density score combinado (up − down) |",
  "",
  "## Next step",
  "- To run the other group set `GROUP <- 'EGFR'` (needs its own data and signature)"
)

writeLines(report_lines,
           file.path(rep_dir, "sc_05_signature_report.md"))
cat("  sc_05_signature_report.md\n")

cat("\n=== sc_05_signature.R completed ===\n")
