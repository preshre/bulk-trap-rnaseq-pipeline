# ============================================================================
# pca_volcano_plot.R - PCA and Volcano Plot Generation
# ============================================================================
#
# Generates PCA plots and volcano plots from DESeq2 analysis.
# Features:
#   - Removes glial markers from volcano plots
#   - Supports custom gene labeling
#   - Supports y-axis breaks for high -log10(p) values
#
# ============================================================================

suppressMessages({
  library(tximport)
  library(readr)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(readxl)
  library(svglite)
  library(ggbreak)
  library(matrixStats)
})

# Load marker lists (lowercased)
glial_markers      <- tolower(unique(read_excel("../gene_data/Glial marker genes.xlsx")$Gene))
exc_neuron_markers <- tolower(unique(read_excel("../gene_data/Excitatory neuron marker genes.xlsx")$Gene))
neuron_markers     <- tolower(unique(read_excel("../gene_data/Neuron marker genes.xlsx")$Gene))
plasticity_genes   <- tolower(unique(read_excel("../gene_data/PFC_Plasticity_Genes.xlsx")$Gene))

run_analysis <- function(
  samples,
  conditions,
  label,
  base_dirs,
  output_dir,
  use_pvalue = FALSE,
  pvalue_cutoff = 0.05,
  log2fc_threshold = 0.58,
  tx2gene_filename = "tx2gene.tsv",
  custom_up_genes = NULL,
  custom_down_genes = NULL,
  use_break = FALSE,
  break_first = 18,
  break_last  = 31
) {
  # Output directory
  plot_dir <- file.path(output_dir, label)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  # Find quant.sf files
  files <- sapply(samples, function(sample) {
    for (base in base_dirs) {
      f <- file.path(base, sample, "quant.sf")
      if (file.exists(f)) return(f)
    }
    stop(paste("❌ quant.sf not found for sample:", sample))
  }, USE.NAMES = FALSE)
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

  # Import quantifications
  cat("\n=== Verifying sample files for:", label, "===\n")
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    cat("Missing quant.sf files:\n"); print(missing_files)
  } else {
    cat("All quant.sf files found.\n")
  }

  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

  # DESeq2 setup
  levs <- unique(conditions)
  colData <- data.frame(
    row.names = colnames(txi$counts),
    condition = factor(conditions, levels = levs)
  )

  dds <- DESeqDataSetFromTximport(txi, colData = colData, design = ~ condition)

  # Filter: keep genes with counts >= 20 in at least half the samples
  keep <- rowSums(counts(dds) >= 20) >= (ncol(dds) / 2)
  dds  <- dds[keep, ]

  dds <- DESeq(dds)
  vsd <- vst(dds, blind = TRUE)

  res <- results(dds)

  # PCA plot
  top_var_genes <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), 500)
  pca <- prcomp(t(assay(vsd)[top_var_genes, ]))
  percent_var <- round(100 * (pca$sdev)^2 / sum(pca$sdev^2), 1)

  pca_data <- data.frame(
    Sample    = rownames(pca$x),
    PC1       = pca$x[, 1],
    PC2       = pca$x[, 2],
    condition = colData$condition
  )

  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = Sample), size = 3) +
    labs(
      title = paste("PCA of", label, "Samples"),
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance")
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )

  ggsave(file.path(plot_dir, "PCA_plot.pdf"), plot = pca_plot, width = 8, height = 6)

  # Volcano plot
  res_df <- as.data.frame(res)
  p_column <- if (use_pvalue) "pvalue" else "padj"

  res_df$significance <- "Not Significant"
  res_df[[p_column]][is.na(res_df[[p_column]])] <- 1

  res_df$significance[
    res_df[[p_column]] < pvalue_cutoff & res_df$log2FoldChange >  log2fc_threshold
  ] <- "Upregulated"

  res_df$significance[
    res_df[[p_column]] < pvalue_cutoff & res_df$log2FoldChange < -log2fc_threshold
  ] <- "Downregulated"

  # Remove glial markers
  rownames(res_df) <- tolower(rownames(res_df))
  res_df <- res_df[!rownames(res_df) %in% glial_markers, ]

  # Marker checks
  cat("Genes after glial removal:", nrow(res_df), "\n\n")
  cat("Excitatory neuron markers present:", length(intersect(rownames(res_df), exc_neuron_markers)), "\n")
  cat("Neuron markers present:", length(intersect(rownames(res_df), neuron_markers)), "\n")
  cat("Plasticity genes present:", length(intersect(rownames(res_df), plasticity_genes)), "\n")

  # Top genes for labeling
  if (!is.null(custom_up_genes) || !is.null(custom_down_genes)) {
    custom_genes <- tolower(unique(c(custom_up_genes, custom_down_genes)))
    top_genes <- res_df[rownames(res_df) %in% custom_genes, , drop = FALSE]
  } else {
    top_up   <- res_df[res_df$significance == "Upregulated", , drop = FALSE]
    top_up   <- top_up[order(top_up[[p_column]]), ][seq_len(min(10, nrow(top_up))), , drop = FALSE]
    top_down <- res_df[res_df$significance == "Downregulated", , drop = FALSE]
    top_down <- top_down[order(top_down[[p_column]]), ][seq_len(min(10, nrow(top_down))), , drop = FALSE]
    top_genes <- rbind(top_up, top_down)
  }

  # Volcano axis limits
  xlim <- c(-max(abs(res_df$log2FoldChange), na.rm = TRUE),
             max(abs(res_df$log2FoldChange), na.rm = TRUE))

  y_max_raw <- max(-log10(res_df[[p_column]]), na.rm = TRUE)
  # nice top tick (1 decimal)
  y_max_tick <- ceiling(y_max_raw / 10) * 10
  ylim  <- c(0, y_max_tick)

  if (use_break) {
    pre_ticks <- pretty(c(0, break_first), n = 5)
    pre_ticks <- pre_ticks[pre_ticks <= break_first & pre_ticks >= 0]
    y_breaks <- unique(c(pre_ticks, break_first, break_last, y_max_tick))
  } else {
    y_breaks <- pretty(c(0, y_max_tick), n = 6)
  }

  base_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(.data[[p_column]]), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(
      "Upregulated"     = "red",
      "Downregulated"   = "blue",
      "Not Significant" = "gray",
      "NA"              = "black"
    )) +
    geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), size = 3, max.overlaps = 50) +
    geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = "darkgray") +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "darkgray") +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_classic() +
    labs(
      title = paste0("Volcano Plot (", label, " - ", ifelse(use_pvalue, "p", "padj"), ")"),
      x = "Log2 Fold Change",
      y = paste0("-Log10(", ifelse(use_pvalue, "p-value", "adjusted p-value"), ")")
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 12, margin = margin(b = 10)),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right  = element_blank()
    )

  if (use_break) {
    volcano_plot <- base_volcano +
      ggbreak::scale_y_break(breaks = c(break_first, break_last)) +
      scale_y_continuous(breaks = y_breaks, limits = ylim, expand = expansion(mult = c(0.02, 0.05)))
  } else {
    volcano_plot <- base_volcano +
      scale_y_continuous(breaks = y_breaks, limits = ylim, expand = expansion(mult = c(0.02, 0.05)))
  }

  ggsave(file.path(plot_dir, "volcano_plot.pdf"),
         plot = volcano_plot, width = 8, height = 6, device = cairo_pdf)

  return(res_df)
}

