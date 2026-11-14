# pipeline_config.R
# Central configuration of the scRNA-seq pipeline
# Change GROUP here or override it before source()
# scRNA-seq pipeline · NSCLC TKI resistance

if (!exists("GROUP")) GROUP <- "MET"   # fallback — GROUP must be set in the CONFIGURATION block of each script

# Dimensionality reduction parameters — shared by sc_02 and sc_03
# Change here if the ElbowPlot suggests a different elbow
if (!exists("DIMS_USE"))    DIMS_USE    <- 1:15
if (!exists("N_NEIGHBORS")) N_NEIGHBORS <- 20

SAMPLES <- list(
  MET  = c("MET1",  "MET1R"),
  EGFR = c("EGFR1",  "EGFR1_R")
)

CONDITIONS <- list(
  MET  = c("Parental", "Resistant"),
  EGFR = c("Parental", "Resistant")
)

GENES_DRIVER <- list(
  MET  = c("MET"),
  EGFR = c("EGFR")
)

SAMPLE_COLORS <- list(
  MET  = met_sample_colors,
  EGFR = egfr_sample_colors
)

GENES_RESISTANCE_PANEL <- list(
  MET  = c("NF2", "WWTR1", "YAP1", "CTGF", "CYR61", "ANKRD1",
           "TEAD1", "TEAD2", "TEAD3", "TEAD4", "NRG1", "MET"),
  EGFR = c("EGFR", "MDM2", "FGFR3", "TACC3", "TP53", "KRAS",
           "ERBB2", "ERBB3", "CDK4", "CCND1")
)

# Transcriptional-state annotation of clusters: fill in after sc_04 (marker genes and GO/KEGG).
# Names must match the cluster ids of the chosen resolution. Unlisted clusters keep
# their numeric id. Example:
#   MET = c("0" = "S-phase", "1" = "M-phase", "2" = "Hypoxic")
CLUSTER_LABELS <- list(
  MET  = character(0),
  EGFR = character(0)
)

GENES_KEY_SIGNATURE <- list(
  MET  = c("NF2", "NRG1", "EGFR", "MET", "ANKRD1", "CTGF", "ERBB4"),
  EGFR = c("EGFR", "MDM2", "FGFR3", "TACC3", "MET", "KRAS", "TP53",
           "NF2", "ERBB2", "ERBB3")
)

GENES_FORCE_HEATMAP <- list(
  MET  = c("NF2", "NRG1"),
  EGFR = c("EGFR", "MDM2")
)

# Drug-tolerant persister (DTP) / EMT markers used in sc_04 to annotate the pre-existence table
# (literature-based; edit to suit the model)
GENES_DTP <- list(
  MET  = c("AXL", "VIM", "ZEB1", "ALDH1A1", "CD44", "ABCB1", "NR2F1", "NGFR"),
  EGFR = c("AXL", "VIM", "ZEB1", "ALDH1A1", "CD44", "ABCB1", "NR2F1", "NGFR")
)

# QC thresholds per dataset — adjusted to the distributions observed in sc_01
# ⚠ Edit here after reviewing qc_01_violin_raw.png and qc_02_scatter_raw.png
QC_PARAMS <- list(
  MET = list(
    MIN_FEATURES = 200,
    MAX_FEATURES = 7000,   # permissive: high-complexity libraries
    MIN_COUNTS   = 500,
    MAX_COUNTS   = 50000,
    MAX_MITO_PCT = 30
  ),
  EGFR = list(
    MIN_FEATURES = 500,
    MAX_FEATURES = 8000,   # higher transcriptional complexity
    MIN_COUNTS   = 1000,
    MAX_COUNTS   = 60000,  # doublets visible above 60k in the scatter
    MAX_MITO_PCT = 20      # dead cells clearly separated below this threshold
  )
)

if (!GROUP %in% names(SAMPLES)) {
  stop(sprintf("GROUP '%s' not recognised. Options: %s",
               GROUP, paste(names(SAMPLES), collapse = ", ")))
}

SAMPLE_PARENTAL  <- SAMPLES[[GROUP]][1]
SAMPLE_RESISTANT <- SAMPLES[[GROUP]][2]
COND_PARENTAL    <- CONDITIONS[[GROUP]][1]
COND_RESISTANT   <- CONDITIONS[[GROUP]][2]
sample_colors    <- SAMPLE_COLORS[[GROUP]]

cat(sprintf("pipeline_config.R: GROUP=%s (%s vs %s)\n",
            GROUP, SAMPLE_PARENTAL, SAMPLE_RESISTANT))
