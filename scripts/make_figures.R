# make_figures.R
# Clean, minimalist version of the four main single-cell figures of a group
# (parental vs resistant), drawn from the checkpoints of the pipeline:
#   1. figure_umaps_composition.png   UMAPs (sample / state) + composition barplot
#   2. figure_umap_high_score.png     UMAP of the pre-existence classification
#   3. figure_preexistence_barplot.png  % of high-score cells per state and sample
#   4. figure_pseudotime_groups.png   pseudotime per group with Wilcoxon tests
# No titles, no legend titles and no sample or model names: the two samples
# are called "Parental" and "Resistant". Group identity comes from the
# `condition` column, so nothing in this script depends on sample names.
#
# Run (from the repo root; GROUP and RESULTS_ROOT are the only things to set):
#   GROUP <- "MET"                       # "MET" or "EGFR"
#   RESULTS_ROOT <- "."                  # folder that holds config/ and <GROUP>/results/ (the original project)
#   source("scripts/make_figures.R")
#
# Input : {RESULTS_ROOT}/{GROUP}/results/{GROUP}_trajectory.rds   (output of sc_06)
#         {RESULTS_ROOT}/config/colors.R and pipeline_config.R    (state labels, colours)
# Output: results/{GROUP}/ in this repo (or OUT_DIR if you set it before running)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

if (!exists("GROUP"))        GROUP <- "MET"
if (!exists("RESULTS_ROOT")) RESULTS_ROOT <- "."

in_rds  <- file.path(RESULTS_ROOT, GROUP, "results", sprintf("%s_trajectory.rds", GROUP))
out_dir <- if (exists("OUT_DIR")) OUT_DIR else file.path("results", GROUP)   # default: results/<GROUP>/ of this repo
if (!file.exists(in_rds)) stop(
  "Not found: ", in_rds, "
",
  "Set RESULTS_ROOT to the project folder that holds ", GROUP, "/results/ and config/, e.g.
",
  "  RESULTS_ROOT <- \"C:/path/to/scRNA-seq\"
",
  "(RESULTS_ROOT is now: ", normalizePath(RESULTS_ROOT, mustWork = FALSE), ")")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Colours and state labels come from the project config
old_wd <- setwd(RESULTS_ROOT)
source("config/colors.R")
source("config/pipeline_config.R")
setwd(old_wd)

# ==============================================================================
# 1. DATA
# ==============================================================================
obj <- readRDS(in_rds)
DefaultAssay(obj) <- "SCT"
cl_col <- "harmony_snn_res.0.4"

cl_chr  <- as.character(obj@meta.data[[cl_col]])
# State names: generic ("State 1", "State 2", ...) by default; set USE_STATE_NAMES <- TRUE
# before running to use the annotation in CLUSTER_LABELS (config/pipeline_config.R) instead
if (!exists("USE_STATE_NAMES")) USE_STATE_NAMES <- FALSE
if (USE_STATE_NAMES) {
  labels <- CLUSTER_LABELS[[GROUP]]
} else {
  ids    <- as.character(sort(unique(suppressWarnings(as.integer(cl_chr)))))
  labels <- setNames(paste("State", seq_along(ids)), ids)
}
obj$state <- factor(ifelse(cl_chr %in% names(labels), labels[cl_chr], cl_chr))

# Two samples called Parental / Resistant, taken from `condition`
cond <- as.character(obj$condition)
obj$group2 <- factor(ifelse(cond == COND_PARENTAL, "Parental", "Resistant"),
                     levels = c("Parental", "Resistant"))
grp_cols <- setNames(unname(sample_colors)[1:2], c("Parental", "Resistant"))

cl_ids <- levels(factor(cl_chr))
cl_cols <- cluster_colors[cl_ids]
cl_cols[is.na(cl_cols)] <- rep_len(unname(cluster_colors), sum(is.na(cl_cols)))
names(cl_cols) <- cl_ids
state_cols <- setNames(unname(cl_cols), ifelse(cl_ids %in% names(labels), labels[cl_ids], cl_ids))
state_cols <- state_cols[!duplicated(names(state_cols))]

# State order: from most to least parental, as in the thesis figures
tab <- table(obj$state, obj$group2)
state_order <- names(sort(tab[, "Parental"] / rowSums(tab), decreasing = TRUE))
obj$state_ordered <- factor(as.character(obj$state), levels = state_order)

cat(sprintf("Cells: %d | Parental: %d | Resistant: %d | states: %d\n",
            ncol(obj), sum(obj$group2 == "Parental"), sum(obj$group2 == "Resistant"), nlevels(obj$state)))

# ==============================================================================
# 2. STYLE
# ==============================================================================
base <- theme_classic(base_size = 13) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

umap_theme <- base + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                           aspect.ratio = 1, legend.position = "right")

umap_by <- function(var, cols) {
  DimPlot(obj, reduction = "umap.harmony", group.by = var, cols = cols, pt.size = 0.3) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(override.aes = list(size = 3))) + umap_theme
}

# ==============================================================================
# 3. FIGURE 1 — TRIPLE UMAP + COMPOSITION
# ==============================================================================
if (USE_STATE_NAMES) {
  p_a <- umap_by("group2", grp_cols) | umap_by(cl_col, cl_cols) | umap_by("state", state_cols)
  tags <- c("A", "", "", "B")
} else {
  # generic states are the clusters renamed, so a separate cluster panel would repeat the state panel
  p_a <- umap_by("group2", grp_cols) | umap_by("state", state_cols)
  tags <- c("A", "", "B")
}

comp <- obj@meta.data |>
  count(state_ordered, group2) |>
  group_by(state_ordered) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup()

p_b <- ggplot(comp, aes(state_ordered, pct, fill = group2)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 3.2, colour = "white", fontface = "bold") +
  scale_fill_manual(values = grp_cols) +
  labs(x = NULL, y = "Cells (%)") +
  base + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"))

p_fig <- p_a / p_b + plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = list(tags))
ggsave(file.path(out_dir, "figure_umaps_composition.png"), p_fig,
       width = 19, height = 13, dpi = 300, bg = "white")

# ==============================================================================
# 4. FIGURE 2 — UMAP OF THE PRE-EXISTENCE CLASSIFICATION
# ==============================================================================
hs <- as.character(obj$high_score_parental)
obj$hs_group <- factor(ifelse(obj$group2 == "Resistant", "Resistant", hs),
                       levels = c("Parental high score", "Parental low score", "Resistant"))
hs_cols <- c("Parental high score" = "#013c58", "Parental low score" = "#b2bec3", "Resistant" = "#e67e22")

p_hs <- DimPlot(obj, reduction = "umap.harmony", group.by = "hs_group", cols = hs_cols,
                pt.size = 0.4, order = TRUE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 3))) + umap_theme
ggsave(file.path(out_dir, "figure_umap_high_score.png"), p_hs,
       width = 9, height = 6, dpi = 300, bg = "white")

# ==============================================================================
# 5. FIGURE 3 — % HIGH-SCORE CELLS PER STATE AND SAMPLE
# ==============================================================================
is_high <- obj$hs_group != "Parental low score"   # high-score parental cells and all resistant cells

pct_high <- obj@meta.data |>
  mutate(is_high = is_high) |>
  group_by(state_ordered, group2) |>
  summarise(n_total = n(), n_high = sum(is_high), pct = 100 * n_high / n_total, .groups = "drop")

p_pre <- ggplot(pct_high, aes(state_ordered, pct, fill = group2)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3) +
  scale_fill_manual(values = grp_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "% high-score cells") +
  base + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
ggsave(file.path(out_dir, "figure_preexistence_barplot.png"), p_pre,
       width = 9, height = 6, dpi = 300, bg = "white")

# ==============================================================================
# 6. FIGURE 4 — PSEUDOTIME PER GROUP (WILCOXON)
# ==============================================================================
pt_df <- data.frame(pseudotime = obj$pseudotime, group = obj$hs_group) |>
  filter(!is.na(pseudotime), is.finite(pseudotime))
pt_df$group <- factor(pt_df$group, levels = c("Parental low score", "Parental high score", "Resistant"))
pt_cols <- c("Parental low score" = "#b2bec3", "Parental high score" = "#e67e22", "Resistant" = "#013c58")

wt <- function(a, b) formatC(wilcox.test(pt_df$pseudotime[pt_df$group == a],
                                         pt_df$pseudotime[pt_df$group == b], exact = FALSE)$p.value,
                             format = "e", digits = 2)
p_ab <- wt("Parental high score", "Parental low score")
p_ar <- wt("Parental high score", "Resistant")
y_max <- max(pt_df$pseudotime)

p_pt <- ggplot(pt_df, aes(group, pseudotime, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.80, scale = "width", width = 0.9) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", alpha = 0.75) +
  scale_fill_manual(values = pt_cols) +
  # square-root axis: pseudotime is right-skewed (most cells near the root), so a
  # linear axis flattens the violins; sqrt keeps 0 and spreads the bulk
  scale_y_sqrt(breaks = c(0, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 50),
               expand = expansion(mult = c(0.02, 0.12))) +
  annotate("segment", x = 1, xend = 2, y = y_max * 1.08, yend = y_max * 1.08, linewidth = 0.6) +
  annotate("text", x = 1.5, y = y_max * 1.20, label = sprintf("p = %s", p_ab), size = 3.8) +
  annotate("segment", x = 2, xend = 3, y = y_max * 1.40, yend = y_max * 1.40, linewidth = 0.6) +
  annotate("text", x = 2.5, y = y_max * 1.55, label = sprintf("p = %s", p_ar), size = 3.8) +
  labs(x = NULL, y = "Pseudotime (sqrt scale)") +
  base + theme(legend.position = "none")
ggsave(file.path(out_dir, "figure_pseudotime_groups.png"), p_pt,
       width = 9, height = 7, dpi = 300, bg = "white")

cat(sprintf("\nFigures saved in %s\n", out_dir))
