# ==============================================================================
# sc_03_clustering.R
# Louvain clustering on the Harmony space, multiple resolutions
# scRNA-seq pipeline · NSCLC TKI resistance
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(clustree)
})

setwd(here::here())   # repo root (set by the .here marker)
source("config/colors.R")
set.seed(42)

# ==============================================================================
# CONFIGURATION — edit only this section before running
# ==============================================================================

GROUP <- "MET"   # MET | EGFR

# ==============================================================================
source("config/pipeline_config.R")

sample_a <- SAMPLE_PARENTAL
sample_b <- SAMPLE_RESISTANT

# DIMS_USE and N_NEIGHBORS come from pipeline_config.R
RESOLUTIONS <- c(0.2, 0.4, 0.6, 0.8, 1.0)
DEFAULT_RES <- 0.4

# Set TRUE to skip FindNeighbors/FindClusters and regenerate only the figures
# from the already saved *_clustered.rds checkpoint.
RERUN_VIZ_ONLY <- FALSE

PROJECT_ROOT <- "."
output_dir   <- file.path(PROJECT_ROOT, GROUP, "results")
fig_dir      <- file.path(output_dir, "figures")
pfig <- function(f) file.path(fig_dir, paste0("sc_03_", f))
tab_dir      <- file.path(output_dir, "tables")
rep_dir      <- file.path(output_dir, "reports")

for (d in c(fig_dir, tab_dir, rep_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# 1/2/3. CLUSTERING OR LOAD FROM CHECKPOINT
# ==============================================================================

if (!RERUN_VIZ_ONLY) {

  cat(sprintf("Loading %s_integrated.rds...\n", GROUP))
  merged <- readRDS(file.path(output_dir, sprintf("%s_integrated.rds", GROUP)))
  cat(sprintf("  Cells: %d (%s=%d, %s=%d)\n",
              ncol(merged),
              sample_a, sum(merged$sample == sample_a),
              sample_b, sum(merged$sample == sample_b)))

  cat(sprintf("\nFindNeighbors (reduction=harmony, dims=1:%d, k=%d)...\n",
              max(DIMS_USE), N_NEIGHBORS))
  merged <- FindNeighbors(merged,
                          reduction  = "harmony",
                          dims       = DIMS_USE,
                          k.param    = N_NEIGHBORS,
                          graph.name = c("harmony_nn", "harmony_snn"),
                          verbose    = FALSE)
  cat("  Grafo harmony_snn construido\n")

  cat("FindClusters (resoluciones:", paste(RESOLUTIONS, collapse = ", "), ")...\n")
  merged <- FindClusters(merged,
                         graph.name = "harmony_snn",
                         resolution = RESOLUTIONS,
                         verbose    = FALSE)

  for (res in RESOLUTIONS) {
    col  <- paste0("harmony_snn_res.", res)
    n_cl <- nlevels(merged[[col, drop = TRUE]])
    cat(sprintf("  res %.1f → %d clusters\n", res, n_cl))
  }

} else {

  clust_rds <- file.path(output_dir, sprintf("%s_clustered.rds", GROUP))
  if (!file.exists(clust_rds)) stop(
    "Checkpoint not found: ", clust_rds,
    " — run with RERUN_VIZ_ONLY=FALSE first"
  )
  merged <- readRDS(clust_rds)
  Idents(merged) <- "seurat_clusters"
  cat(sprintf("[RERUN_VIZ_ONLY] %s_clustered.rds loaded: %d cells\n",
              GROUP, ncol(merged)))

}

# ==============================================================================
# 4. FIGURE 1 — Resolution grid
# ==============================================================================

cat("\nGenerating clust_01_resolution_grid.png...\n")

plot_list <- lapply(RESOLUTIONS, function(res) {
  col_name <- paste0("harmony_snn_res.", res)
  n_cl     <- nlevels(merged[[col_name, drop = TRUE]])
  DimPlot(merged, reduction = "umap.harmony", group.by = col_name,
          label = TRUE, repel = TRUE, cols = cluster_colors, pt.size = 0.3) +
    NoLegend() +
    ggtitle(sprintf("res = %.1f  (%d clusters)", res, n_cl)) +
    theme(plot.title = element_text(size = 11))
})

ggsave(pfig("clust_01_resolution_grid.png"),
       wrap_plots(plot_list, ncol = 3), width = 20, height = 12, dpi = 800)
cat("  Saved: clust_01_resolution_grid.png\n")

# ==============================================================================
# 5. FIGURE 2 — Clustree
# ==============================================================================

cat("Generating clust_02_clustree.png...\n")

ggsave(pfig("clust_02_clustree.png"),
       clustree(merged, prefix = "harmony_snn_res."),
       width = 10, height = 12, dpi = 800)
cat("  Saved: clust_02_clustree.png\n")

# ==============================================================================
# 6. SET DEFAULT RESOLUTION
# ==============================================================================

default_col <- paste0("harmony_snn_res.", DEFAULT_RES)

if (!RERUN_VIZ_ONLY) {
  Idents(merged)         <- default_col
  merged$seurat_clusters <- merged[[default_col, drop = TRUE]]
}

n_clusters_default <- nlevels(merged$seurat_clusters)
cat(sprintf("\nActive resolution: %.1f → %d clusters\n", DEFAULT_RES, n_clusters_default))

# ==============================================================================
# 7. FIGURE 3 — UMAP clusters + condition
# ==============================================================================

cat("Generating clust_03_umap_clusters_condition.png...\n")

p_clust <- DimPlot(merged, reduction = "umap.harmony", group.by = "seurat_clusters",
                   label = TRUE, repel = TRUE, cols = cluster_colors, pt.size = 0.3) +
  labs(title = sprintf("Clusters (res = %.1f)", DEFAULT_RES)) +
  theme(legend.position = "right")

p_cond <- DimPlot(merged, reduction = "umap.harmony", group.by = "condition",
                  cols = condition_colors, pt.size = 0.3) +
  labs(title = "Condition") +
  theme(legend.position = "right")

ggsave(pfig("clust_03_umap_clusters_condition.png"),
       p_clust | p_cond, width = 16, height = 7, dpi = 800)
cat("  Saved: clust_03_umap_clusters_condition.png\n")

# ==============================================================================
# 8. FIGURE 4 — Composition per cluster (sample % barplot)
# ==============================================================================

cat("Generating clust_04_composition_barplot.png...\n")

comp_df_ord <- merged@meta.data |>
  group_by(seurat_clusters, sample) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(seurat_clusters) |>
  mutate(pct   = n / sum(n) * 100,
         total = sum(n)) |>
  ungroup() |>
  mutate(seurat_clusters = factor(seurat_clusters,
                                  levels = as.character(sort(as.integer(levels(seurat_clusters))))))

p_comp_ord <- ggplot(comp_df_ord, aes(x = seurat_clusters, y = pct, fill = sample)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Cluster", y = "% cells", fill = "Sample",
       title = sprintf("Composition per cluster — n cells (res = %.1f)", DEFAULT_RES)) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

ggsave(pfig("clust_04_composition_barplot.png"),
       p_comp_ord, width = 10, height = 5, dpi = 800)
cat("  Saved: clust_04_composition_barplot.png\n")

# ==============================================================================
# 8c. FIGURE 5 — FeaturePlot of key genes (resistance)
# ==============================================================================

cat("Generating clust_05_featureplot_resistance_genes.png...\n")

genes_feature <- GENES_RESISTANCE_PANEL[[GROUP]]
genes_feature <- genes_feature[genes_feature %in% rownames(merged)]

ggsave(pfig("clust_05_featureplot_resistance_genes.png"),
       FeaturePlot(merged, features = genes_feature, reduction = "umap.harmony",
                   cols = expr_gradient, pt.size = 0.2, ncol = 4, order = TRUE),
       width = 20, height = ceiling(length(genes_feature) / 4) * 5, dpi = 800)
cat(sprintf("  Saved: clust_05_featureplot_resistance_genes.png (%s)\n",
            paste(genes_feature, collapse = ", ")))

# ==============================================================================
# 8d. FIGURE 6 — DotPlot of markers per cluster
# ==============================================================================

cat("Generating clust_06_dotplot.png...\n")

dot_genes <- GENES_RESISTANCE_PANEL[[GROUP]]
dot_genes <- dot_genes[dot_genes %in% rownames(merged)]
Idents(merged) <- "seurat_clusters"

ggsave(pfig("clust_06_dotplot.png"),
       DotPlot(merged, features = dot_genes, assay = "SCT",
               cols = expr_gradient, dot.scale = 6) +
         labs(title = sprintf("Expression of resistance genes per cluster — %s", GROUP),
              x = NULL, y = "Cluster") +
         theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
         white_bg,
       width = 14, height = 6, dpi = 800)
cat("  Saved: clust_06_dotplot.png\n")

# ==============================================================================
# 9. TABLE — Numeric composition (dynamic per sample)
# ==============================================================================

cat(sprintf("Generating cluster_composition_res%.1f.csv...\n", DEFAULT_RES))

comp_table <- merged@meta.data |>
  group_by(seurat_clusters, sample) |>
  summarise(n = n(), .groups = "drop") |>
  pivot_wider(names_from = sample, values_from = n, values_fill = 0) |>
  rename(n_a = all_of(sample_a), n_b = all_of(sample_b)) |>
  mutate(
    Total                           = n_a + n_b,
    !!sprintf("pct_%s", sample_a) := round(n_a / Total * 100, 1),
    !!sprintf("pct_%s", sample_b) := round(n_b / Total * 100, 1)
  ) |>
  rename(!!sample_a := n_a, !!sample_b := n_b) |>
  arrange(seurat_clusters)

write.csv(comp_table,
          file.path(tab_dir, sprintf("cluster_composition_res%.1f.csv", DEFAULT_RES)),
          row.names = FALSE)
cat(sprintf("  Saved: cluster_composition_res%.1f.csv\n", DEFAULT_RES))

# ==============================================================================
# 10. SAVE OBJECT
# ==============================================================================

if (!RERUN_VIZ_ONLY) {
  cat("\nPreparing SCT models for FindMarkers...\n")
  merged <- PrepSCTFindMarkers(merged)

  out_rds <- file.path(output_dir, sprintf("%s_clustered.rds", GROUP))
  cat(sprintf("\nSaving %s_clustered.rds...\n", GROUP))
  saveRDS(merged, out_rds)
  cat(sprintf("  %s_integrated.rds — NOT modified\n", GROUP))
  cat(sprintf("  %s_clustered.rds  — %d cells, seurat_clusters (res %.1f)\n",
              GROUP, ncol(merged), DEFAULT_RES))
} else {
  cat("\n[RERUN_VIZ_ONLY] Object not re-saved (no changes in clustering)\n")
}

# ==============================================================================
# 11. REPORT .md
# ==============================================================================

cat("\nGenerating sc_03_clustering_report.md...\n")

res_summary <- lapply(RESOLUTIONS, function(res) {
  col  <- paste0("harmony_snn_res.", res)
  n_cl <- nlevels(merged[[col, drop = TRUE]])
  data.frame(Resolucion = res, N_clusters = n_cl)
}) |> bind_rows()

pct_a_col <- sprintf("pct_%s", sample_a)
pct_b_col <- sprintf("pct_%s", sample_b)
comp_rows <- apply(comp_table, 1, function(r) {
  sprintf("| %s | %s | %s | %s | %s%% | %s%% |",
          r["seurat_clusters"],
          r[sample_a], r[sample_b], r["Total"],
          r[pct_a_col], r[pct_b_col])
})

report_lines <- c(
  sprintf("# sc_03_clustering — %s: %s vs %s", GROUP, sample_a, sample_b),
  "",
  sprintf("**Fecha:** %s", Sys.Date()),
  sprintf("**Cells:** %d (%s=%d, %s=%d)",
          ncol(merged),
          sample_a, sum(merged$sample == sample_a),
          sample_b, sum(merged$sample == sample_b)),
  "",
  "## Parameters",
  sprintf("- Reduction: harmony (dims 1:%d)", max(DIMS_USE)),
  sprintf("- k vecinos: %d", N_NEIGHBORS),
  sprintf("- Resoluciones: %s", paste(RESOLUTIONS, collapse = ", ")),
  sprintf("- **Default resolution: %.1f (%d clusters)**", DEFAULT_RES, n_clusters_default),
  "",
  "## Clusters per resolution",
  "",
  "| Resolution | N clusters |",
  "|------------|------------|",
  paste(apply(res_summary, 1, function(r)
    sprintf("| %.1f | %s |", as.numeric(r[1]), r[2])), collapse = "\n"),
  "",
  sprintf("## Composition per cluster (res %.1f)", DEFAULT_RES),
  "",
  sprintf("| Cluster | %s | %s | Total | %% %s | %% %s |",
          sample_a, sample_b, sample_a, sample_b),
  "|---------|------|-----|-------|--------|-------|",
  paste(comp_rows, collapse = "\n"),
  "",
  "## Figures",
  "- `clust_01_resolution_grid.png` — UMAPs at multiple resolutions",
  "- `clust_02_clustree.png` — cluster tree across resolutions",
  "- `clust_03_umap_clusters_condition.png` — clusters + condition",
  "- `clust_04_composition_barplot.png` — barplot % MET1/MET1R per cluster (absolute n)",
  "",
  "## Next step: sc_04_markers_all.R",
  sprintf("- Input: `%s/results/%s_clustered.rds`", GROUP, GROUP),
  "- FindMarkers per cluster or per condition within a cluster"
)

writeLines(report_lines, file.path(rep_dir, "sc_03_clustering_report.md"))
cat("  sc_03_clustering_report.md\n")

cat("\n=== sc_03_clustering.R completed ===\n")
