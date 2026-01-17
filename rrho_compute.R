# ============================================================================
# rrho_compute.R - RRHO2 Analysis Computation
# ============================================================================
#
# Computes Rank-Rank Hypergeometric Overlap (RRHO) between two ranked gene lists.
# Generates heatmaps, Venn diagrams, and overlap gene lists.
#
# ============================================================================

run_rrho2_analysis <- function(vec1_path, vec2_path, output_dir, label1, label2, threshold = -1.3) {
  # Load libraries
  if (!requireNamespace("RRHO2", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
    remotes::install_github("RRHO2/RRHO2", build_opts = c("--no-resave-data", "--no-manual"))
  }
  library(RRHO2)
  library(ggplot2)

  # Load vectors
  vec1 <- readRDS(vec1_path)
  vec2 <- readRDS(vec2_path)

  # Convert to data.frame
  df1 <- na.omit(data.frame(Genes = names(vec1), Score = as.numeric(vec1)))
  df2 <- na.omit(data.frame(Genes = names(vec2), Score = as.numeric(vec2)))

  # Keep common genes in both lists
  common_genes <- intersect(df1$Genes, df2$Genes)
  df1 <- df1[df1$Genes %in% common_genes, ]
  df2 <- df2[df2$Genes %in% common_genes, ]

  # Create output dirs
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  plots_dir <- file.path(output_dir, "plots")
  csv_dir   <- file.path(output_dir, "csv")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(csv_dir,   showWarnings = FALSE, recursive = TRUE)

  # Run RRHO2
  rrho <- RRHO2_initialize(df1, df2, labels = c(label1, label2), log10.ind = TRUE)

  # Save the object
  saveRDS(rrho, file = file.path(output_dir, "RRHO2_result.rds"))

  # Generate plots
  pdf(file.path(plots_dir, "RRHO2_heatmap.pdf"), width = 8, height = 6)
  RRHO2_heatmap(rrho)
  dev.off()

  for (tp in c("dd", "uu", "du", "ud")) {
    pdf(file.path(plots_dir, paste0("RRHO2_vennDiagram_", tp, ".pdf")))
    RRHO2_vennDiagram(rrho, type = tp)
    dev.off()
  }

  # CSV exports
  .to_df <- function(x) {
    if (is.null(x)) return(data.frame(Genes = character(), stringsAsFactors = FALSE))
    if (is.vector(x)) return(data.frame(Genes = as.character(x), stringsAsFactors = FALSE))
    if (is.data.frame(x)) {
      if (ncol(x) == 0) return(data.frame(Genes = character(), stringsAsFactors = FALSE))
      if (!"Genes" %in% colnames(x)) {
        colnames(x)[1] <- "Genes"
      }
      x$Genes <- as.character(x$Genes)
      return(x[, "Genes", drop = FALSE])
    }
    # Fallback
    data.frame(Genes = as.character(x), stringsAsFactors = FALSE)
  }

  .write_csv <- function(df, path) {
    write.csv(df, file = path, row.names = FALSE)
  }

  all_overlap <- list()

  export_one_type <- function(type_key) {
    slot_name <- paste0("genelist_", type_key)
    obj <- rrho[[slot_name]]
    if (is.null(obj)) return(invisible())

    g1   <- .to_df(obj[[paste0("gene_list1_", type_key)]])
    g2   <- .to_df(obj[[paste0("gene_list2_", type_key)]])
    over <- .to_df(obj[[paste0("gene_list_overlap_", type_key)]])

    .write_csv(g1,   file.path(csv_dir, paste0("genes_list1_", type_key, ".csv")))
    .write_csv(g2,   file.path(csv_dir, paste0("genes_list2_", type_key, ".csv")))
    .write_csv(over, file.path(csv_dir, paste0("genes_overlap_", type_key, ".csv")))

    if (nrow(over) > 0) {
      over$quadrant <- rep(type_key, nrow(over))
      over$label1   <- rep(label1,   nrow(over))
      over$label2   <- rep(label2,   nrow(over))
      all_overlap[[type_key]] <- over
    } else {
      message("No overlap genes for quadrant: ", type_key)
    }
  }

  for (tp in c("dd", "uu", "du", "ud")) export_one_type(tp)

  if (length(all_overlap)) {
    combo <- do.call(rbind, all_overlap)
    combo <- combo[, c("Genes", "quadrant", "label1", "label2"), drop = FALSE]
    write.csv(combo, file.path(csv_dir, "genes_overlap_all_quadrants.csv"), row.names = FALSE)
  } else {
    write.csv(data.frame(Genes=character(), quadrant=character(), label1=character(), label2=character()),
              file.path(csv_dir, "genes_overlap_all_quadrants.csv"),
              row.names = FALSE)
  }


  message("RRHO2 analysis complete. Outputs in: ", normalizePath(output_dir))
}