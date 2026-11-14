# colors.R
# Colour palettes — scRNA-seq NSCLC TKI resistance pipeline
# System aligned with RNAseq_bulk/config/colors.R
# scRNA-seq pipeline · NSCLC TKI resistance

library(RColorBrewer)

# ==============================================================================
# CONDITION / SAMPLE
# ==============================================================================

# Biological condition — aligned with RNAseq_bulk/config/colors.R
condition_colors <- c(
  "Parental"  = "#e6a817",   # dorado
  "Resistant" = "#013c58"    # azul marino
)

# Technical sample identity — turquoise/salmon
# Distinguishes the sample (10x run) from the biological condition
met_sample_colors <- c(
  "MET1" = "#12CBC4",   # parental — turquesa
  "MET1R"  = "#ef745c"    # resistant — salmon
)

egfr_sample_colors <- c(
  "EGFR1"   = "#12CBC4",   # parental — turquesa
  "EGFR1_R" = "#ef745c"    # resistant in vivo — salmon
)

# ==============================================================================
# CLUSTERS (colours for clusters 0-8, maximum perceptual separation)
# Colours are assigned by cluster id, independent of any biological annotation.
# ==============================================================================

cluster_colors <- c(
  "0"  = "#3498db",   # blue
  "1"  = "#27ae60",   # green
  "2"  = "#f39c12",   # amber
  "3"  = "#e74c3c",   # red
  "4"  = "#95a5a6",   # grey
  "5"  = "#8e44ad",   # purple
  "6"  = "#fd79a8",   # pink
  # Extra slots for datasets with more clusters
  "7"  = "#00b894",
  "8"  = "#00cec9",
  "9"  = "#6c5ce7",
  "10" = "#fdcb6e",
  "11" = "#e17055",
  "12" = "#636e72",
  "13" = "#b2bec3",
  "14" = "#2d3436",
  "15" = "#74b9ff",
  "16" = "#a29bfe",
  "17" = "#55efc4",
  "18" = "#fab1a0",
  "19" = "#dfe6e9"
)

# ==============================================================================
# GENE EXPRESSION (gradients)
# ==============================================================================

# scRNA-seq's own colour identity: teal / magenta, distinct from
# RNAseq_bulk (negro-gris-dorado), WES (morado) y scDNA-seq (azul-naranja/rojo,
# black-blue-red). Two variants: sequential (magnitude, no negatives, e.g.
# raw expression or score) and diverging (above/below the mean, e.g.
# z-score), both anchored on the same colour pair so that they read
# as the same family.
GRADIENT_SCRNA_SEQ  <- c("lightgrey", "#c2185b")
GRADIENT_SCRNA_DIV  <- c("#16a085", "white", "#c2185b")

# FeaturePlot and DotPlot of raw expression: grey -> magenta (previously grey -> dark red)
expr_gradient <- GRADIENT_SCRNA_SEQ

# Heatmap z-score (bidireccional): verde azulado -> blanco -> magenta
heatmap_zscore   <- colorRampPalette(GRADIENT_SCRNA_DIV)(100)
heatmap_standard <- colorRampPalette(GRADIENT_SCRNA_DIV)(100)
heatmap_distances <-
  colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# ==============================================================================
# RESISTANCE SIGNATURE SCORE
# FeaturePlot scores UCell/AMS y heatmap pathways
# Grey -> magenta (previously pearl grey -> burnt orange -> dark red)
# ==============================================================================

score_gradient <- colorRampPalette(GRADIENT_SCRNA_SEQ)(100)

# ==============================================================================
# PRE-EXISTENCE — colours for high_score_parental and traj_07 groups
# Concept: parental low signature → parental high signature → resistant
# ==============================================================================

preexistence_colors <- c(
  "Parental low score" = "#b2bec3",   # silver grey  — no signature
  "Parental high score" = "#e67e22",   # burnt orange — precursor
  "Resistant"          = "#013c58"    # azul marino     — = Resistant
)

high_score_high_color <- "#e67e22"   # parental high score
high_score_low_color  <- "#b2bec3"   # parental low score

# ==============================================================================
# ENRICHMENT BARPLOTS (mall_08 up, mall_09 down, mall_10 preexistencia)
# ==============================================================================

enrich_up_color       <- "#ef745c"   # salmon orange — up genes
enrich_down_color     <- "#12CBC4"   # turquoise     — down genes
enrich_preexist_color <- "#6C5CE7"   # violeta        — preexistencia up

# ==============================================================================
# GENE COVERAGE (sig_08)
# ==============================================================================

coverage_colors <- c(
  "Upregulated"   = "#e67e22",   # naranja quemado (= high_score_high_color)
  "Downregulated" = "#2980b9"    # azul acero
)

# ==============================================================================
# ANOTACIONES pheatmap
# ==============================================================================

annotation_colors_met <- list(
  Sample    = met_sample_colors,
  Condition = condition_colors
)

annotation_colors_egfr <- list(
  Sample    = egfr_sample_colors,
  Condition = condition_colors
)

# ==============================================================================
# GENERAL CATEGORICAL PALETTE
# ==============================================================================

categorical_palette <- c(
  "#e74c3c", "#3498db", "#27ae60", "#f39c12", "#8e44ad",
  "#1abc9c", "#34495e", "#e67e22", "#95a5a6", "#16a085"
)

# ==============================================================================
# BASE THEME (white background, no grid)
# ==============================================================================

white_bg <- theme(
  panel.background  = element_rect(fill = "white", colour = NA),
  plot.background   = element_rect(fill = "white", colour = NA),
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank()
)

cat("colors.R loaded\n")
