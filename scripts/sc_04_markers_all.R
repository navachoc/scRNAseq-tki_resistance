# ==============================================================================
# sc_04_markers_all.R — Unified differential expression & enrichment
# MET: MET1 vs MET1R  |  Adaptable to EGFR: EGFR1 vs EGFR1_R
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(presto)
  library(ggrepel)
  library(tibble)
  library(pheatmap)
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

LOGFC_THRESHOLD    <- 0.25
MIN_PCT            <- 0.10
PADJ_THRESHOLD     <- 0.05
MIN_CELLS_PREEXIST <- 30   # min parental cells per cluster

# GENES_DTP defined in config/pipeline_config.R per GROUP
GENES_DTP <- GENES_DTP[[GROUP]]

# Set TRUE to skip FindMarkers/enrichment and regenerate only the figures
# from the already saved *_de_results.rds checkpoint.
RERUN_VIZ_ONLY <- FALSE
# ==============================================================================

output_dir <- file.path(".", GROUP, "results")
fig_dir    <- file.path(output_dir, "figures")
tab_dir    <- file.path(output_dir, "tables")
rep_dir    <- file.path(output_dir, "reports")
pfig <- function(f) file.path(fig_dir, paste0("sc_04_", f))
for (d in c(fig_dir, tab_dir, rep_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# HELPER FUNCTIONS (defined once)
# ==============================================================================

run_enrichment <- function(gene_vec, universe_entrez = NULL) {
  tryCatch({
    entrez_df <- bitr(gene_vec, fromType = "SYMBOL",
                      toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
    cat(sprintf("      bitr: %d/%d genes mapped\n", nrow(entrez_df), length(gene_vec)))
    if (nrow(entrez_df) < 5) return(list(go = NULL, kegg = NULL))
    ego <- enrichGO(gene          = entrez_df$ENTREZID,
                    universe      = universe_entrez,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.25,
                    minGSSize     = 5,
                    readable      = TRUE)
    ekegg <- enrichKEGG(gene         = entrez_df$ENTREZID,
                        universe     = universe_entrez,
                        organism     = "hsa",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.25,
                        minGSSize    = 5)
    cat(sprintf("      GO-BP: %d sig. terms | KEGG: %d sig. terms\n",
                if (!is.null(ego))   nrow(as.data.frame(ego))   else 0L,
                if (!is.null(ekegg)) nrow(as.data.frame(ekegg)) else 0L))
    list(go = ego, kegg = ekegg)
  }, error = function(e) {
    cat(sprintf("    Warning enrichment: %s\n", e$message))
    list(go = NULL, kegg = NULL)
  })
}

save_enrich_tables <- function(enrich, prefix, dir) {
  go_df   <- if (!is.null(enrich$go))   as.data.frame(enrich$go)
             else data.frame(Note = "No significant terms")
  kegg_df <- if (!is.null(enrich$kegg)) as.data.frame(enrich$kegg)
             else data.frame(Note = "No significant terms")
  write.csv(go_df,   file.path(dir, paste0(prefix, "_go.csv")),   row.names = FALSE)
  write.csv(kegg_df, file.path(dir, paste0(prefix, "_kegg.csv")), row.names = FALSE)
}

make_enrich_barplot <- function(enrich_obj, title, fill_color, n = 15) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "No significant enriched terms", size = 5) +
             theme_void() +
             theme(panel.background = element_rect(fill = "white", color = NA)) +
             labs(title = title))
  }
  df <- as.data.frame(enrich_obj) |>
    dplyr::slice_min(order_by = p.adjust, n = n) |>
    dplyr::mutate(Description = factor(Description,
                    levels = Description[order(p.adjust, decreasing = TRUE)]))
  ggplot(df, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = fill_color, high = "grey80",
                        name = "p adj",
                        guide = guide_colorbar(reverse = TRUE)) +
    labs(title = title, x = "Number of genes", y = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.y       = element_text(size = 9),
          panel.background  = element_rect(fill = "white", color = NA),
          plot.background   = element_rect(fill = "white", color = NA))
}

make_volcano <- function(vol_df, label_up, label_down,
                          logfc_thr, padj_thr,
                          color_up   = "#ef745c",
                          color_down = "#12CBC4",
                          title      = "",
                          n_label    = 20) {
  vol_df <- vol_df |>
    dplyr::mutate(
      neg_log10_padj = -log10(p_val_adj + 1e-300),
      significance   = dplyr::case_when(
        p_val_adj < padj_thr & avg_log2FC >  logfc_thr ~ label_up,
        p_val_adj < padj_thr & avg_log2FC < -logfc_thr ~ label_down,
        TRUE ~ "NS"
      )
    )
  top_label <- vol_df |>
    dplyr::filter(significance != "NS") |>
    dplyr::slice_max(order_by = abs(avg_log2FC), n = n_label)
  ggplot(vol_df, aes(x = avg_log2FC, y = neg_log10_padj, color = significance)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_vline(xintercept = c(-logfc_thr, logfc_thr),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(padj_thr),
               linetype = "dashed", color = "grey50") +
    ggrepel::geom_text_repel(data = top_label, aes(label = gene),
                              size = 3, max.overlaps = 20) +
    scale_color_manual(values = setNames(
      c(color_up, color_down, "grey70"),
      c(label_up, label_down, "NS")
    )) +
    labs(title = title, x = "log2 Fold Change", y = "-log10(p adj)", color = NULL) +
    theme_classic(base_size = 13) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA))
}



# ==============================================================================
# 1. LOAD OBJECT
# ==============================================================================

rds_file <- file.path(output_dir, sprintf("%s_clustered.rds", GROUP))
cat(sprintf("Loading %s...\n", basename(rds_file)))
merged <- readRDS(rds_file)
Idents(merged) <- "seurat_clusters"
CLUSTERS_ALL <- levels(merged$seurat_clusters)
  cat(sprintf("  Cells: %d | Clusters: %s\n",
            ncol(merged), paste(CLUSTERS_ALL, collapse = ", ")))

# Enrichment universe: all detected genes → Entrez IDs
cat("  Building Entrez universe...\n")
all_genes_detected <- rownames(merged)
universe_entrez_df <- suppressMessages(
  bitr(all_genes_detected, fromType = "SYMBOL",
       toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
)
UNIVERSE_ENTREZ <- universe_entrez_df$ENTREZID
cat(sprintf("  Universe: %d genes → %d Entrez IDs\n",
            length(all_genes_detected), length(UNIVERSE_ENTREZ)))

# ==============================================================================
# 2/3/4. DIFFERENTIAL EXPRESSION OR LOAD FROM CHECKPOINT
# ==============================================================================

de_checkpoint <- file.path(output_dir, sprintf("%s_de_results.rds", GROUP))

if (!RERUN_VIZ_ONLY) {

  # ── 2. FindAllMarkers ──────────────────────────────────────────────────────
  cat("\n[Sec 2] FindAllMarkers...\n")
  all_markers <- FindAllMarkers(
    merged, assay = "SCT", slot = "data",
    logfc.threshold = LOGFC_THRESHOLD, min.pct = MIN_PCT,
    test.use = "wilcox", only.pos = TRUE, verbose = FALSE
  ) |> dplyr::filter(p_val_adj < PADJ_THRESHOLD)

  write.csv(all_markers,
            file.path(tab_dir, "mall_markers_all_clusters.csv"),
            row.names = FALSE)
  cat(sprintf("  mall_markers_all_clusters.csv (%d genes)\n", nrow(all_markers)))

  # ── 3. FindMarkers loop per cluster ───────────────────────────────────────
  cat("\n[Sec 3] FindMarkers loop per cluster...\n")
  results_list <- list()

  for (cl in CLUSTERS_ALL) {
    cat(sprintf("\n  [Cluster %s] FindMarkers vs resto...\n", cl))

    vol_df <- FindMarkers(
      merged, ident.1 = cl,
      assay = "SCT", slot = "data",
      logfc.threshold = 0, min.pct = 0.05,
      test.use = "wilcox", verbose = FALSE
    ) |> tibble::rownames_to_column("gene")

    markers_cl <- vol_df |>
      dplyr::filter(p_val_adj < PADJ_THRESHOLD,
                    abs(avg_log2FC) >= LOGFC_THRESHOLD) |>
      dplyr::arrange(desc(avg_log2FC))

    write.csv(markers_cl,
              file.path(tab_dir, sprintf("mall_markers_cluster%s_vs_rest.csv", cl)),
              row.names = FALSE)
    cat(sprintf("    Sig: %d (up: %d, down: %d)\n",
                nrow(markers_cl),
                sum(markers_cl$avg_log2FC > 0),
                sum(markers_cl$avg_log2FC < 0)))

    results_list[[cl]] <- list(vol_df = vol_df, markers = markers_cl)
  }

  # ── 4. Bucle preexistencia ─────────────────────────────────────────────────
  cat(sprintf("\n[Sec 4] Bucle preexistencia %s...\n", SAMPLE_PARENTAL))

  parental_cells  <- colnames(merged)[merged$sample == SAMPLE_PARENTAL]
  merged_parental <- subset(merged, cells = parental_cells)
  Idents(merged_parental) <- "seurat_clusters"
  preexist_results <- list()

  for (cl in CLUSTERS_ALL) {
    n_cl   <- sum(merged_parental$seurat_clusters == cl)
    n_rest <- sum(merged_parental$seurat_clusters != cl)

    if (n_cl < MIN_CELLS_PREEXIST || n_rest < MIN_CELLS_PREEXIST) {
      cat(sprintf("  [Cluster %s] Skip — %d cells %s (min %d)\n",
                  cl, n_cl, SAMPLE_PARENTAL, MIN_CELLS_PREEXIST))
      next
    }

    cat(sprintf("  [Cluster %s] %s: %d vs resto %d...\n",
                cl, SAMPLE_PARENTAL, n_cl, n_rest))

    markers_pre <- FindMarkers(
      merged_parental, ident.1 = cl,
      assay = "SCT", slot = "data",
      logfc.threshold = LOGFC_THRESHOLD, min.pct = MIN_PCT,
      test.use = "wilcox", verbose = FALSE
    ) |>
      tibble::rownames_to_column("gene") |>
      dplyr::filter(p_val_adj < PADJ_THRESHOLD) |>
      dplyr::arrange(desc(avg_log2FC))

    write.csv(markers_pre,
              file.path(tab_dir, sprintf("mall_preexist_cluster%s_vs_%srest.csv",
                                          cl, SAMPLE_PARENTAL)),
              row.names = FALSE)
    preexist_results[[cl]] <- markers_pre
    cat(sprintf("    Sig: %d genes\n", nrow(markers_pre)))
  }

  # ── Save DE checkpoint ─────────────────────────────────────────────────────
  saveRDS(
    list(all_markers      = all_markers,
         results_list     = results_list,
         preexist_results = preexist_results),
    de_checkpoint
  )
  cat(sprintf("\n  DE checkpoint saved: %s\n", basename(de_checkpoint)))

} else {

  if (!file.exists(de_checkpoint)) stop(
    "Checkpoint not found: ", de_checkpoint,
    " — run with RERUN_VIZ_ONLY=FALSE first"
  )
  de_data          <- readRDS(de_checkpoint)
  all_markers      <- de_data$all_markers
  results_list     <- de_data$results_list
  preexist_results <- de_data$preexist_results
  cat(sprintf("[RERUN_VIZ_ONLY] DE checkpoint loaded: %d genes all_markers, %d clusters\n",
              nrow(all_markers), length(results_list)))

}

# ==============================================================================
# 5. FIGURES — HEATMAPS (always)
# ==============================================================================

cat("\n[Sec 5] Heatmaps top10 markers...\n")

top10_genes <- all_markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(order_by = avg_log2FC, n = 10) |>
  dplyr::ungroup() |>
  dplyr::pull(gene) |>
  unique()

avail_cols <- cluster_colors[CLUSTERS_ALL]

# mall_01a — heatmap per cluster (original DoHeatmap)
p_heat <- DoHeatmap(
  merged, features = top10_genes, assay = "SCT", slot = "scale.data",
  group.colors = avail_cols
) +
  scale_fill_gradientn(colors = heatmap_zscore, name = "Z-score",
                       na.value = "white") +
  white_bg

ggsave(pfig("mall_01a_heatmap_top10_by_cluster.png"),
       p_heat, width = 20, height = 14, dpi = 800)
cat("  mall_01a_heatmap_top10_by_cluster.png\n")

# mall_01b — same gene set, grouped by sample (MET1 vs MET1R)
p_heat_sample <- DoHeatmap(
  merged, features = top10_genes, assay = "SCT", slot = "scale.data",
  group.by = "sample", group.colors = sample_colors
) +
  scale_fill_gradientn(colors = heatmap_zscore, name = "Z-score",
                       na.value = "white") +
  white_bg

ggsave(pfig("mall_01b_heatmap_top10_by_sample.png"),
       p_heat_sample, width = 20, height = 14, dpi = 800)
cat("  mall_01b_heatmap_top10_by_sample.png\n")

# mall_01c — heatmap per cluster + sample annotation (pheatmap)
genes_ok  <- top10_genes[top10_genes %in% rownames(merged[["SCT"]]@scale.data)]
mat       <- merged[["SCT"]]@scale.data[genes_ok, ]
mat       <- pmin(pmax(mat, -2.5), 2.5)

ann_col <- data.frame(
  Cluster = factor(merged$seurat_clusters, levels = CLUSTERS_ALL),
  Sample  = factor(merged$sample),
  row.names = colnames(merged)
)
cell_order <- order(ann_col$Cluster, ann_col$Sample)
mat     <- mat[, cell_order]
ann_col <- ann_col[cell_order, , drop = FALSE]

ann_colors <- list(
  Cluster = setNames(cluster_colors[CLUSTERS_ALL], CLUSTERS_ALL),
  Sample  = sample_colors
)

pheatmap(
  mat,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  show_colnames     = FALSE,
  color             = heatmap_zscore,
  breaks            = seq(-2.5, 2.5, length.out = length(heatmap_zscore) + 1),
  border_color      = NA,
  fontsize_row      = 8,
  annotation_names_col = TRUE,
  main = sprintf("Top 10 markers per cluster — %s vs %s", SAMPLE_PARENTAL, SAMPLE_RESISTANT),
  filename = pfig("mall_01c_heatmap_top10_cluster_sample.png"),
  width = 20, height = 14
)
cat("  mall_01c_heatmap_top10_cluster_sample.png\n")

# ==============================================================================
# 6. FIGURES — VOLCANOES + ENRICHMENT PER CLUSTER (always)
# ==============================================================================

cat("\n[Sec 6] Volcanoes + enrichment per cluster...\n")

for (cl in CLUSTERS_ALL) {
  vol_df     <- results_list[[cl]]$vol_df
  markers_cl <- results_list[[cl]]$markers

  if (is.null(vol_df)) {
    cat(sprintf("  [Cluster %s] No data — skipped\n", cl))
    next
  }

  p_vol <- make_volcano(
    vol_df,
    label_up   = sprintf("Up cluster %s", cl),
    label_down = sprintf("Down cluster %s", cl),
    logfc_thr  = LOGFC_THRESHOLD,
    padj_thr   = PADJ_THRESHOLD,
    title = sprintf("Cluster %s vs rest", cl)
  )
  ggsave(pfig(sprintf("mall_02_%s_volcano_cluster%s.png", cl, cl)),
         p_vol, width = 10, height = 7, dpi = 800)
  cat(sprintf("  mall_02_%s_volcano_cluster%s.png\n", cl, cl))

  genes_up <- markers_cl |> dplyr::filter(avg_log2FC > 0) |> dplyr::pull(gene)
  if (length(genes_up) >= 5) {
    enrich_up <- run_enrichment(genes_up, universe_entrez = UNIVERSE_ENTREZ)
    save_enrich_tables(enrich_up, sprintf("mall_enrich_cluster%s_up", cl), tab_dir)
    p_up <- make_enrich_barplot(
      enrich_up$go,
      title = sprintf("GO-BP up — cluster %s", cl), fill_color = "#ef745c"
    ) /
    make_enrich_barplot(
      enrich_up$kegg,
      title = sprintf("KEGG up — cluster %s", cl), fill_color = "#ef745c"
    )
    ggsave(pfig(sprintf("mall_03_%s_enrichment_up_cluster%s.png", cl, cl)),
           p_up, width = 12, height = 14, dpi = 800)
  }

  genes_down <- markers_cl |> dplyr::filter(avg_log2FC < 0) |> dplyr::pull(gene)
  if (length(genes_down) >= 5) {
    enrich_down <- run_enrichment(genes_down, universe_entrez = UNIVERSE_ENTREZ)
    save_enrich_tables(enrich_down, sprintf("mall_enrich_cluster%s_down", cl), tab_dir)
    p_dn <- make_enrich_barplot(
      enrich_down$go,
      title = sprintf("GO-BP down — cluster %s", cl), fill_color = "#12CBC4"
    ) /
    make_enrich_barplot(
      enrich_down$kegg,
      title = sprintf("KEGG down — cluster %s", cl), fill_color = "#12CBC4"
    )
    ggsave(pfig(sprintf("mall_04_%s_enrichment_down_cluster%s.png", cl, cl)),
           p_dn, width = 12, height = 14, dpi = 800)
  }

  cat(sprintf("  Cluster %s completed.\n", cl))
}

# ==============================================================================
# 7. FIGURES — PRE-EXISTENCE ENRICHMENT (always)
# ==============================================================================

cat(sprintf("\n[Sec 7] Enrichment preexistencia %s...\n", SAMPLE_PARENTAL))

for (cl in CLUSTERS_ALL) {
  pre <- preexist_results[[cl]]
  if (is.null(pre)) next

  genes_pre_up <- pre |> dplyr::filter(avg_log2FC > 0) |> dplyr::pull(gene)
  if (length(genes_pre_up) >= 5) {
    enrich_pre <- run_enrichment(genes_pre_up, universe_entrez = UNIVERSE_ENTREZ)
    save_enrich_tables(enrich_pre,
                       sprintf("mall_enrich_preexist_cluster%s_up", cl), tab_dir)
    p_pre <- make_enrich_barplot(
      enrich_pre$go,
      title = sprintf("GO-BP — %s cl%s preexistencia", SAMPLE_PARENTAL, cl),
      fill_color = "#6C5CE7"
    ) /
    make_enrich_barplot(
      enrich_pre$kegg,
      title = sprintf("KEGG — %s cl%s preexistencia", SAMPLE_PARENTAL, cl),
      fill_color = "#6C5CE7"
    )
    ggsave(pfig(sprintf("mall_05_%s_preexist_enrich_cluster%s.png", cl, cl)),
           p_pre, width = 12, height = 14, dpi = 800)
    cat(sprintf("  mall_05_%s_preexist_enrich_cluster%s.png\n", cl, cl))
  }
}

# ==============================================================================
# 8. ANNOTATED UMAP + PRE-EXISTENCE SUMMARY TABLE (always)
# ==============================================================================

cat("\n[Sec 8] UMAP annotated with transcriptional states...\n")

# Annotations defined in pipeline_config.R (CLUSTER_LABELS) — fill in after sc_03
cluster_labels <- CLUSTER_LABELS[[GROUP]]

# Clusters without an annotation keep their numeric id as label
all_clusters <- levels(factor(merged$seurat_clusters))
cluster_labels_full <- setNames(all_clusters, all_clusters)
cluster_labels_full[names(cluster_labels)] <- unname(cluster_labels)

merged$cluster_label <- dplyr::recode(
  as.character(merged$seurat_clusters),
  !!!cluster_labels_full
)

# Map cluster_colors (by number) to the labels; recycle the palette if there are more clusters
label_cols <- cluster_colors[names(cluster_labels_full)]
label_cols[is.na(label_cols)] <- rep_len(unname(cluster_colors), sum(is.na(label_cols)))
label_colors <- setNames(unname(label_cols), unname(cluster_labels_full))

p_umap_ann <- DimPlot(
  merged, reduction = "umap.harmony",
  group.by = "cluster_label", label = TRUE, repel = TRUE,
  label.size = 3.5, pt.size = 0.4,
  cols = label_colors
) +
  labs(title = sprintf("%s — Transcriptional state (%s vs %s)",
                       GROUP, SAMPLE_PARENTAL, SAMPLE_RESISTANT)) +
  theme_classic(base_size = 12) + white_bg +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(pfig("mall_06_umap_annotated.png"),
       p_umap_ann, width = 11, height = 7, dpi = 800)
cat("  mall_06_umap_annotated.png\n")

# Pre-existence summary table
cat("\n[Sec 8b] Pre-existence summary table...\n")

preexist_summary <- dplyr::bind_rows(lapply(CLUSTERS_ALL, function(cl) {
  n_cl   <- sum(merged$seurat_clusters == cl   & merged$sample == SAMPLE_PARENTAL)
  n_rest <- sum(merged$seurat_clusters != cl   & merged$sample == SAMPLE_PARENTAL)
  pre    <- preexist_results[[cl]]
  n_up   <- if (!is.null(pre)) sum(pre$avg_log2FC > 0) else NA_integer_
  n_down <- if (!is.null(pre)) sum(pre$avg_log2FC < 0) else NA_integer_
  dtp_hit <- if (!is.null(pre)) {
    paste(intersect(pre$gene[pre$avg_log2FC > 0], GENES_DTP), collapse = ";")
  } else NA_character_
  data.frame(
    cluster       = cl,
    label         = unname(cluster_labels_full[as.character(cl)]),
    n_parental    = n_cl,
    n_rest        = n_rest,
    analyzed      = !is.null(pre) && n_cl >= MIN_CELLS_PREEXIST,
    DEG_up        = n_up,
    DEG_down      = n_down,
    DTP_markers   = dtp_hit,
    stringsAsFactors = FALSE
  )
}))

write.csv(preexist_summary,
          file.path(tab_dir, "mall_preexist_summary.csv"),
          row.names = FALSE)
cat("  mall_preexist_summary.csv\n")
print(preexist_summary)

# ==============================================================================
# 8. REPORT
# ==============================================================================

cat("\n[Sec 9] Generating report...\n")

report_lines <- c(
  sprintf("# sc_04_markers_all — %s: %s vs %s", GROUP, SAMPLE_RESISTANT, SAMPLE_PARENTAL),
  "",
  sprintf("**Fecha:** %s", Sys.Date()),
  "",
  "## Analyses performed",
  "1. FindAllMarkers → heatmap per cluster (mall_01a) + heatmap per sample (mall_01b) + heatmap cluster + cell origin (mall_01c)",
  "2. Cluster loop: volcano per cluster vs rest (mall_02) + GO/KEGG enrichment up/down (mall_03, mall_04)",
  sprintf("3. Pre-existence loop %s: FindMarkers per cluster → GO/KEGG enrichment up (mall_05)",
          SAMPLE_PARENTAL),
  "4. UMAP annotated with transcriptional states (mall_06)",
  "",
  "## Summary per cluster"
)

for (cl in CLUSTERS_ALL) {
  if (is.null(results_list[[cl]])) next
  res <- results_list[[cl]]$markers
  report_lines <- c(
    report_lines, "",
    sprintf("### Cluster %s", cl),
    sprintf("- Sig genes vs rest: %d (up: %d, down: %d)",
            nrow(res), sum(res$avg_log2FC > 0), sum(res$avg_log2FC < 0)),
    sprintf("- Top 5 up: %s",
            paste(head(res$gene[res$avg_log2FC > 0], 5), collapse = ", "))
  )
  pre <- preexist_results[[cl]]
  if (!is.null(pre)) {
    report_lines <- c(
      report_lines,
      sprintf("- Pre-existence %s: %d sig genes (up: %d)",
              SAMPLE_PARENTAL, nrow(pre), sum(pre$avg_log2FC > 0)),
      sprintf("- Top 5 up preexistencia: %s",
              paste(head(pre$gene[pre$avg_log2FC > 0], 5), collapse = ", "))
    )
  }
}

report_lines <- c(
  report_lines, "",
  "## Next step: sc_05_signature.R (AddModuleScore)"
)

writeLines(report_lines,
           file.path(rep_dir, "sc_04_markers_all_report.md"))
cat("  sc_04_markers_all_report.md\n")

cat("\n=== sc_04_markers_all.R completed ===\n")
