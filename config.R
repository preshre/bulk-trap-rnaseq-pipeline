# ============================================================================
# config.R - Central configuration for RNA-seq analysis pipeline
# ============================================================================

# === Data paths ===
BASE_DIRS <- c(
  "/gpfs/scratch/lpunepallera/Trap-Seq/g-drive-data/CamK2a.TRAP_PTC_BO_Trained/star_salmon",
  "/gpfs/scratch/lpunepallera/Trap-Seq/g-drive-data/CamK2a.iPKRTRAP_PTC_BO_Trained/star_salmon",
  "/gpfs/scratch/lpunepallera/Trap-Seq/g-drive-data/CamK2a.TRAP_PTC_BO_Trained-4h/star_salmon"
)

TX2GENE_FILENAME <- "tx2gene.tsv"

# === Output directory ===
# Generate timestamped output directory for each run
# Results stored in 'results' folder alongside scripts
SCRIPT_DIR <- dirname(sys.frame(1)$ofile)
if (is.null(SCRIPT_DIR) || SCRIPT_DIR == "") {
  SCRIPT_DIR <- getwd()
}
RESULTS_BASE <- file.path(SCRIPT_DIR, "results")
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(RESULTS_BASE, paste0("results_", TIMESTAMP))

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", OUTPUT_DIR)

# === DEG thresholds (consistent across all analyses) ===
LOG2FC_THRESHOLD <- 0.58
PADJ_CUTOFF      <- 0.05
PVALUE_CUTOFF    <- 0.05

# === Sample definitions ===
# vTRAP BO vs Trained (PT1h)
SAMPLES_vTRAP_BO_Trained <- list(
  samples = c(
    "CamK2a-vTRAP-BO-PT1h-IP1", "CamK2a-vTRAP-BO-PT1h-IP2", "CamK2a-vTRAP-BO-PT1h-IP3",
    "CamK2a-vTRAP-Trained-PT1h-IP1", "CamK2a-vTRAP-Trained-PT1h-IP3", "CamK2a-vTRAP-Trained-PT1h-IP4"
  ),
  conditions = c(rep("vTRAP_BO", 3), rep("vTRAP_Trained", 3)),
  named_vector = c(
    "CamK2a-vTRAP-BO-PT1h-IP1" = "vTRAP_BO",
    "CamK2a-vTRAP-BO-PT1h-IP2" = "vTRAP_BO",
    "CamK2a-vTRAP-BO-PT1h-IP3" = "vTRAP_BO",
    "CamK2a-vTRAP-Trained-PT1h-IP1" = "vTRAP_Trained",
    "CamK2a-vTRAP-Trained-PT1h-IP3" = "vTRAP_Trained",
    "CamK2a-vTRAP-Trained-PT1h-IP4" = "vTRAP_Trained"
  ),
  label = "vTRAP_BO_vTRAP_Trained",
  ref_level = "vTRAP_BO"
)

# vTRAP BO vs Trained (PT4h)
SAMPLES_vTRAP_PT4h <- list(
  samples = c(
    "CamK2a-vTRAP-BO-PT4h-IP1", "CamK2a-vTRAP-BO-PT4h-IP2", "CamK2a-vTRAP-BO-PT4h-IP3",
    "CamK2a-vTRAP-Trained-PT4h-IP1", "CamK2a-vTRAP-Trained-PT4h-IP2", "CamK2a-vTRAP-Trained-PT4h-IP3"
  ),
  conditions = c(rep("vTRAP_BO_PT4h", 3), rep("vTRAP_Trained_PT4h", 3)),
  named_vector = c(
    "CamK2a-vTRAP-BO-PT4h-IP1" = "vTRAP_BO_PT4h",
    "CamK2a-vTRAP-BO-PT4h-IP2" = "vTRAP_BO_PT4h",
    "CamK2a-vTRAP-BO-PT4h-IP3" = "vTRAP_BO_PT4h",
    "CamK2a-vTRAP-Trained-PT4h-IP1" = "vTRAP_Trained_PT4h",
    "CamK2a-vTRAP-Trained-PT4h-IP2" = "vTRAP_Trained_PT4h",
    "CamK2a-vTRAP-Trained-PT4h-IP3" = "vTRAP_Trained_PT4h"
  ),
  label = "vTRAP_BO_PT4h_vTRAP_Trained_PT4h",
  ref_level = "vTRAP_BO_PT4h"
)

# iPKRTRAP BO vs Trained
SAMPLES_iPKRTRAP_BO_Trained <- list(
  samples = c(
    "CamK2a-iPKRTRAP-ASV-BO-IP1", "CamK2a-iPKRTRAP-ASV-BO-IP2", "CamK2a-iPKRTRAP-ASV-BO-IP3",
    "CamK2a-iPKRTRAP-ASV-Trained-IP1", "CamK2a-iPKRTRAP-ASV-Trained-IP2", "CamK2a-iPKRTRAP-ASV-Trained-IP3"
  ),
  conditions = c(rep("iPKRTRAP_BO", 3), rep("iPKRTRAP_Trained", 3)),
  named_vector = c(
    "CamK2a-iPKRTRAP-ASV-BO-IP1" = "iPKRTRAP_BO",
    "CamK2a-iPKRTRAP-ASV-BO-IP2" = "iPKRTRAP_BO",
    "CamK2a-iPKRTRAP-ASV-BO-IP3" = "iPKRTRAP_BO",
    "CamK2a-iPKRTRAP-ASV-Trained-IP1" = "iPKRTRAP_Trained",
    "CamK2a-iPKRTRAP-ASV-Trained-IP2" = "iPKRTRAP_Trained",
    "CamK2a-iPKRTRAP-ASV-Trained-IP3" = "iPKRTRAP_Trained"
  ),
  label = "iPKRTRAP_BO_iPKRTRAP_Trained",
  ref_level = "iPKRTRAP_BO"
)

# vTRAP BO vs iPKRTRAP BO
SAMPLES_vTRAP_BO_iPKRTRAP_BO <- list(
  samples = c(
    "CamK2a-vTRAP-BO-PT1h-IP1", "CamK2a-vTRAP-BO-PT1h-IP2", "CamK2a-vTRAP-BO-PT1h-IP3",
    "CamK2a-iPKRTRAP-ASV-BO-IP1", "CamK2a-iPKRTRAP-ASV-BO-IP2", "CamK2a-iPKRTRAP-ASV-BO-IP3"
  ),
  conditions = c(rep("vTRAP_BO", 3), rep("iPKRTRAP_BO", 3)),
  named_vector = c(
    "CamK2a-vTRAP-BO-PT1h-IP1" = "vTRAP_BO",
    "CamK2a-vTRAP-BO-PT1h-IP2" = "vTRAP_BO",
    "CamK2a-vTRAP-BO-PT1h-IP3" = "vTRAP_BO",
    "CamK2a-iPKRTRAP-ASV-BO-IP1" = "iPKRTRAP_BO",
    "CamK2a-iPKRTRAP-ASV-BO-IP2" = "iPKRTRAP_BO",
    "CamK2a-iPKRTRAP-ASV-BO-IP3" = "iPKRTRAP_BO"
  ),
  label = "vTRAP_BO_iPKRTRAP_BO",
  ref_level = "vTRAP_BO"
)

# vTRAP Trained vs iPKRTRAP Trained
SAMPLES_vTRAP_Trained_iPKRTRAP_Trained <- list(
  samples = c(
    "CamK2a-vTRAP-Trained-PT1h-IP1", "CamK2a-vTRAP-Trained-PT1h-IP3", "CamK2a-vTRAP-Trained-PT1h-IP4",
    "CamK2a-iPKRTRAP-ASV-Trained-IP1", "CamK2a-iPKRTRAP-ASV-Trained-IP2", "CamK2a-iPKRTRAP-ASV-Trained-IP3"
  ),
  conditions = c(rep("vTRAP_Trained", 3), rep("iPKRTRAP_Trained", 3)),
  named_vector = c(
    "CamK2a-vTRAP-Trained-PT1h-IP1" = "vTRAP_Trained",
    "CamK2a-vTRAP-Trained-PT1h-IP3" = "vTRAP_Trained",
    "CamK2a-vTRAP-Trained-PT1h-IP4" = "vTRAP_Trained",
    "CamK2a-iPKRTRAP-ASV-Trained-IP1" = "iPKRTRAP_Trained",
    "CamK2a-iPKRTRAP-ASV-Trained-IP2" = "iPKRTRAP_Trained",
    "CamK2a-iPKRTRAP-ASV-Trained-IP3" = "iPKRTRAP_Trained"
  ),
  label = "vTRAP_Trained_iPKRTRAP_Trained",
  ref_level = "vTRAP_Trained"
)

# === Comparisons to Run ===
# Add/remove comparisons from this list to control pipeline execution
ALL_COMPARISONS <- list(
  SAMPLES_vTRAP_BO_Trained,
  SAMPLES_vTRAP_PT4h,
  SAMPLES_iPKRTRAP_BO_Trained,
  SAMPLES_vTRAP_BO_iPKRTRAP_BO,
  SAMPLES_vTRAP_Trained_iPKRTRAP_Trained
)

# === Gene Marker Files (for volcano plots) ===
GENE_DATA_DIR <- "../gene_data"

# === Helper Function ===
get_pval_suffix <- function(use_pvalue) {
  if (use_pvalue) "_pvalue" else "_padj"
}

message("✅ Config loaded: ", length(ALL_COMPARISONS), " comparisons defined")
