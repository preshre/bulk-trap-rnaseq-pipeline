# ============================================================================
# deg_analysis.R - Differential Expression Gene Analysis with Heatmap
# ============================================================================
#
# Performs DESeq2 analysis and generates:
#   - DEG lists (up/down regulated)
#   - Top 2000 genes by significance
#   - Heatmap of top DEGs
#
# ============================================================================

run_deg_analysis <- function(
  samples,
  conditions,
  base_dirs,
  tx2gene_filename,
  output_dir,
  label,
  use_pvalue = FALSE,
  pvalue_cutoff = 0.05,
  padj_cutoff = 0.05,
  log2fc_threshold = 0.58,
  baseMean_min = 10,
  top_n_heatmap = 30
) {
  library(tximport)
  library(readr)
  library(DESeq2)
  library(pheatmap)
  library(grid)

  message("▶ Starting DEG analysis for: ", label)

  # Create output directory
  output_path <- file.path(output_dir, label, "DEG")
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

  # Locate quant.sf files
  files <- sapply(samples, function(sample) {
    for (base in base_dirs) {
      f <- file.path(base, sample, "quant.sf")
      if (file.exists(f)) return(f)
    }
    stop("❌ quant.sf not found for sample: ", sample)
  })
  names(files) <- samples

  # Load tx2gene
  tx2gene <- NULL
  for (base in base_dirs) {
    tx2gene_path <- file.path(base, tx2gene_filename)
    if (file.exists(tx2gene_path)) {
      tx2gene <- read_tsv(tx2gene_path, col_names = c("TXNAME", "GENEID"), show_col_types = FALSE)
      break
    }
  }
  if (is.null(tx2gene)) stop("❌ tx2gene.tsv not found")

  # Import with tximport
  names(conditions) <- samples
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

  # Set up conditions
  col_data <- data.frame(condition = factor(conditions), row.names = names(conditions))

  # Filter low-count genes
  filter_condition <- rowSums(txi$counts >= 20) >= (ncol(txi$counts) / 2)
  txi$counts <- txi$counts[filter_condition, ]
  txi$abundance <- txi$abundance[filter_condition, ]
  txi$length <- txi$length[filter_condition, ]

  # DESeq2 analysis
  dds <- DESeqDataSetFromTximport(txi, colData = col_data, design = ~ condition)
  dds <- DESeq(dds)
  levs <- levels(col_data$condition)
  
  res <- results(dds, contrast = c("condition", levs[1], levs[2]))
  message("  Contrast: ", levs[1], " vs ", levs[2])

  # Determine cutoff based on use_pvalue
  sig_col <- if (use_pvalue) "pvalue" else "padj"
  cutoff  <- if (use_pvalue) pvalue_cutoff else padj_cutoff

  # Export top 2000 genes by significance
  res_df <- as.data.frame(res)
  res_df <- res_df[order(res_df[[sig_col]]), ]
  res_top2000 <- head(res_df, 2000)
  write.csv(res_top2000, file.path(output_path, paste0("Top2000_genes_by_", sig_col, ".csv")))

  # Filter DEGs
  DEGs <- res[which(res[[sig_col]] < cutoff & abs(res$log2FoldChange) > log2fc_threshold & res$baseMean > baseMean_min), ]
  DEGs <- DEGs[order(DEGs$log2FoldChange, decreasing = TRUE), ]

  up_lev1 <- subset(DEGs, log2FoldChange > 0)
  up_lev2 <- subset(DEGs, log2FoldChange < 0)

  write.csv(DEGs, file.path(output_path, "DEGS.csv"))
  write.csv(up_lev1, file.path(output_path, paste0("DEGS_upregulated_", levs[1], ".csv")))
  write.csv(up_lev2, file.path(output_path, paste0("DEGS_upregulated_", levs[2], ".csv")))

  message("  Total DEGs: ", nrow(DEGs), " (", nrow(up_lev1), " up in ", levs[1], ", ", nrow(up_lev2), " up in ", levs[2], ")")

  # Heatmap of top DEGs
  if (nrow(DEGs) == 0) {
    message("  ⚠ No DEGs found - skipping heatmap")
    return(invisible(list(dds = dds, res = res, DEGs = DEGs)))
  }

  # Select top genes for heatmap
  n_per_direction <- floor(top_n_heatmap / 2)
  if (nrow(DEGs) > top_n_heatmap) {
    top_up <- head(up_lev1, n_per_direction)
    top_down <- head(up_lev2, n_per_direction)
    top_combined <- rbind(top_down, top_up)
    write.csv(top_combined, file.path(output_path, paste0("DEGS_TOP_", top_n_heatmap, ".csv")))
  } else {
    top_combined <- rbind(up_lev1, up_lev2)
  }

  # Generate heatmap
  vsd <- vst(dds, blind = FALSE)
  top_expr <- assay(vsd)[rownames(top_combined), ]

  heatmap_plot <- pheatmap(top_expr,
                           scale = "row",
                           clustering_distance_rows = "euclidean",
                           clustering_distance_cols = "euclidean",
                           clustering_method = "complete",
                           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                           main = paste("Top DEGs Heatmap:", label),
                           fontsize_row = 6,
                           fontsize_col = 8,
                           angle_col = 45,
                           cellwidth = 18,
                           cellheight = 6,
                           legend = TRUE,
                           border_color = "gray",
                           silent = TRUE)

  # Save heatmap
  save_pheatmap_pdf <- function(x, filename, width = 6, height = 10) {
    pdf(filename, width = width, height = height)
    grid.newpage()
    grid.draw(x$gtable)
    dev.off()
  }

  save_pheatmap_pdf(heatmap_plot, file.path(output_path, paste0("Heatmap_Top", top_n_heatmap, ".pdf")))

  message("✅ DEG analysis complete for: ", label)

  invisible(list(dds = dds, res = res, DEGs = DEGs))
}
