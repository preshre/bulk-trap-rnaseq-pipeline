# ============================================================================
# generate_DE_plots.R - MA Plot Generation
# ============================================================================
#
# Generates MA plots from DESeq2 differential expression analysis.
# Supports both adjusted p-value and raw p-value cutoffs.
#
# ============================================================================

run_deseq_analysis_with_ma_plot <- function(
  samples,
  conditions,
  base_dirs,
  tx2gene_filename,
  output_dir,
  comparison_label,
  lfc_cutoff    = 0.58,
  padj_cutoff   = 0.05,
  use_pvalue    = FALSE,
  pvalue_cutoff = 0.05
) {
  suppressMessages({
    library(tximport)
    library(readr)
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
  })

  message("▶ Starting DESeq2 for: ", comparison_label)

  # Output directory
  plot_dir <- file.path(output_dir)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  # Locate quant.sf files
  files <- sapply(samples, function(sample) {
    for (base in base_dirs) {
      f <- file.path(base, sample, "quant.sf")
      if (file.exists(f)) return(f)
    }
    stop(paste("❌ quant.sf not found for sample:", sample))
  })
  names(files) <- samples

  # Load tx2gene
  tx2gene <- NULL
  for (base in base_dirs) {
    tx2gene_path <- file.path(base, tx2gene_filename)
    if (file.exists(tx2gene_path)) {
      tx2gene <- read_tsv(tx2gene_path, col_names = TRUE, show_col_types = FALSE)
      break
    }
  }
  if (is.null(tx2gene)) stop("❌ tx2gene.tsv not found in any base directory.")

  # Import & DESeq2
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  coldata <- data.frame(row.names = samples, condition = factor(conditions, levels = unique(conditions)))
  dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

  dds <- DESeq(dds)
  res <- results(dds)

  # Save DE results
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  csv_path <- file.path(output_dir, paste0("DE_results_", comparison_label, ".csv"))
  write.csv(res_df, csv_path, row.names = FALSE)
  message("✅ DE results saved to: ", csv_path)

  # Choose significance column & cutoff
  pcol   <- if (use_pvalue) "pvalue" else "padj"
  cutoff <- if (use_pvalue) pvalue_cutoff else padj_cutoff

  res_df[[pcol]][is.na(res_df[[pcol]])] <- 1
  res_df$significance <- ifelse(res_df[[pcol]] < cutoff, "Significant", "Not Significant")

  top_up <- subset(res_df, significance == "Significant" & log2FoldChange >  0.56)
  top_up <- top_up[order(top_up[[pcol]]), ][seq_len(min(5, nrow(top_up))), ]

  top_down <- subset(res_df, significance == "Significant" & log2FoldChange < -0.56)
  top_down <- top_down[order(top_down[[pcol]]), ][seq_len(min(5, nrow(top_down))), ]

  top_genes <- rbind(top_up, top_down)

  # MA Plot
  ma_plot <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
    scale_x_log10() +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "blue", size = 0.5) +
    theme_minimal() +
    labs(
      title = paste0("Enhanced MA Plot: ", comparison_label, " (", pcol, " < ", cutoff, ")"),
      x = "Mean of Normalized Counts (log10 scale)",
      y = "Log2 Fold Change",
      color = "Significance"
    )

  ggsave(
    filename = file.path(output_dir, paste0("MA_plot_", comparison_label, ".pdf")),
    plot = ma_plot, width = 8, height = 6
  )

  message("✅ Enhanced MA plot saved.\n")
}
