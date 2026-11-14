# gene_lists.R
# Gene lists used by the scRNA-seq pipeline
# scRNA-seq pipeline · NSCLC TKI resistance

# ==============================================================================
# MET — MET-amplified model
# ==============================================================================

# MET signalling and bypass
met_core <- c(
  "MET", "EGFR", "KRAS", "BRAF", "AXL", "ERBB2", "ERBB3", "IGF1R",
  "MAPK1", "AKT1", "MTOR", "PIK3CA", "TP53", "CDK4", "CDK6", "CCND1", "MYC"
)

met_bypass <- c(
  "EGFR", "KRAS", "BRAF", "AXL", "RET", "FGFR1", "FGFR2", "ERBB2", "IGF1R"
)

# AKT/ERK/mTOR — reactivation common to all resistant populations
met_pathway_common <- c(
  "AKT1", "AKT2", "AKT3", "MAPK1", "MAPK3", "MTOR",
  "RPS6KB1", "EIF4EBP1", "MAP2K1", "MAP2K2"
)

# NF2 inactivation + YAP/TAZ — resistant populations with NF2 loss
met_nf2_yap <- c(
  "NF2", "WWTR1", "YAP1", "CTGF", "CYR61", "ANKRD1",
  "TEAD1", "TEAD2", "TEAD3", "TEAD4"
)

# EGFR/HER bypass — resistant populations with EGFR/HER bypass
met_egfr_bypass <- c(
  "EGFR", "ERBB2", "ERBB3", "ERBB4", "HBEGF", "AREG", "EREG", "NRG1"
)

# Lista maestra MET resistance core
met_resistance_core <- unique(c(
  met_core, met_bypass,
  "NF2", "YAP1", "WWTR1", "CTGF", "CYR61",
  "STAT3", "JAK1", "JAK2",
  "BCL2", "MCL1", "BCL2L1",
  "VIM", "CDH1", "CDH2", "ZEB1", "ZEB2", "SNAI1", "SNAI2",
  "CCND1", "CDK4", "CDK6", "RB1", "CDKN2A"
))

# ==============================================================================
# EGFR — EGFR-mutant model
# ==============================================================================

egfr_core <- c(
  "EGFR", "ERBB2", "ERBB3", "ERBB4",
  "KRAS", "NRAS", "BRAF", "RAF1",
  "MAP2K1", "MAP2K2", "MAPK1", "MAPK3",
  "PIK3CA", "AKT1", "AKT2", "AKT3", "MTOR", "PTEN"
)

# Additional alterations acquired in vivo
egfr_acquired_alterations <- c(
  "MDM2", "TP53",           # MDM2 amplification
  "FGFR3", "TACC3",         # FGFR3-TACC3 fusion
  "FGFR1", "FGFR2", "FGFR4" # FGFR pathway
)

egfr_resistance_bypass <- c(
  "MET", "AXL", "IGF1R", "ERBB2", "FGFR1", "FGFR2", "FGFR3",
  "KRAS", "BRAF", "PIK3CA"
)

# ==============================================================================
# GENERAL PATHWAYS
# ==============================================================================

mapk_erk_pathway <- c(
  "KRAS", "HRAS", "NRAS", "BRAF", "RAF1",
  "MAP2K1", "MAP2K2", "MAPK1", "MAPK3",
  "DUSP6", "SPRY2", "FOS", "JUN"
)

pi3k_akt_pathway <- c(
  "PIK3CA", "PIK3CB", "AKT1", "AKT2", "AKT3",
  "MTOR", "PTEN", "PDPK1", "IRS1", "IGF1R",
  "TSC1", "TSC2", "FOXO1", "FOXO3"
)

hippo_pathway <- c(
  "NF2", "LATS1", "LATS2", "MST1", "MST2",
  "YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4",
  "CTGF", "CYR61", "ANKRD1"
)

wnt_pathway <- c(
  "WNT5A", "WNT3A", "FZD1", "FZD7", "LRP5", "LRP6",
  "CTNNB1", "APC", "AXIN1", "GSK3B",
  "TCF7L2", "MYC", "CCND1", "SNAI1"
)

jak_stat_pathway <- c(
  "JAK1", "JAK2", "JAK3",
  "STAT1", "STAT3", "STAT5A", "STAT5B",
  "IL6", "IL6R"
)

# ==============================================================================
# EMT
# ==============================================================================

emt_tfs <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")

mesenchymal_markers <- c(
  "VIM", "FN1", "CDH2", "CDH11",
  "MMP2", "MMP9", "MMP14", "SPARC",
  "ITGB1", "ITGA5", "TNC"
)

epithelial_markers <- c(
  "CDH1", "OCLN", "CLDN1", "CLDN4", "EPCAM", "KRT8", "KRT18"
)

emt_signature <- unique(c(emt_tfs, mesenchymal_markers, epithelial_markers))

# ==============================================================================
# APOPTOSIS Y CICLO CELULAR
# ==============================================================================

apoptosis_genes <- c(
  "BCL2", "BCL2L1", "MCL1", "BCL2L11", "BAX", "BAK1",
  "CASP3", "CASP8", "CASP9", "TP53", "MDM2", "CDKN1A", "CDKN2A"
)

cell_cycle <- c(
  "CCND1", "CCNE1", "CCNA2", "CCNB1",
  "CDK2", "CDK4", "CDK6",
  "CDKN1A", "CDKN1B", "CDKN2A",
  "RB1", "E2F1", "MKI67"
)

stemness_genes <- c(
  "NANOG", "SOX2", "KLF4", "MYC",
  "ALDH1A1", "CD44", "PROM1", "BMI1", "LIN28A"
)

# ==============================================================================
# LISTA EXPORTADA
# ==============================================================================

all_gene_lists <- list(
  # MET — bulk RNA-seq signature (defined in sc_05 if GROUP == "MET")
  met_resistance_signature_721 = get0("met_resistance_signature_721"),
  met_resistance_signature_up  = get0("met_resistance_signature_up"),
  met_resistance_signature_dn  = get0("met_resistance_signature_dn"),
  # EGFR — bulk RNA-seq signature (defined in sc_05 if GROUP == "EGFR")
  egfr_resistance_signature_721 = get0("egfr_resistance_signature_721"),
  egfr_resistance_signature_up  = get0("egfr_resistance_signature_up"),
  egfr_resistance_signature_dn  = get0("egfr_resistance_signature_dn"),
  # MET
  met_core            = met_core,
  met_bypass          = met_bypass,
  met_resistance_core = met_resistance_core,
  met_pathway_common      = met_pathway_common,
  met_nf2_yap      = met_nf2_yap,
  met_egfr_bypass = met_egfr_bypass,
  # EGFR
  egfr_core              = egfr_core,
  egfr_resistance_bypass = egfr_resistance_bypass,
  egfr_acquired_alterations      = egfr_acquired_alterations,
  # Pathways
  mapk_erk_pathway    = mapk_erk_pathway,
  pi3k_akt_pathway    = pi3k_akt_pathway,
  hippo_pathway       = hippo_pathway,
  wnt_pathway         = wnt_pathway,
  jak_stat_pathway    = jak_stat_pathway,
  # Fenotipo
  emt_signature       = emt_signature,
  emt_tfs             = emt_tfs,
  mesenchymal_markers = mesenchymal_markers,
  epithelial_markers  = epithelial_markers,
  apoptosis_genes     = apoptosis_genes,
  cell_cycle          = cell_cycle,
  stemness_genes      = stemness_genes
)

cat("gene_lists.R loaded\n")
cat("  met_resistance_core:", length(met_resistance_core), "genes\n")
cat("  egfr_core:          ", length(egfr_core), "genes\n")
