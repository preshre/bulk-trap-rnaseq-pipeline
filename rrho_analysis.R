# ============================================================================
# rrho_analysis.R - RRHO Vector Generation from DESeq2
# ============================================================================
#
# Generates RRHO ranking vectors from DESeq2 differential expression analysis.
# The ranking metric is: sign(log2FoldChange) * -log10(padj)
#
# ============================================================================

run_deseq_to_rrho_vector <- function(
  sample_names,
  conditions,
  base_dirs,
  tx2gene_filename,
  label,
  output_dir,
  ref_level = NULL,
  min_count = 20,
  min_samples_prop = 0.5
) {
  suppressMessages({
    library(tximport)
    library(DESeq2)
    library(readr)
  })

  message("▶ Starting DESeq2 analysis for: ", label)

  # Output directory
  plot_dir <- file.path(output_dir, label, "RRHO")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  # Locate quant.sf files
  files <- sapply(sample_names, function(sample) {
    for (base in base_dirs) {
      f <- file.path(base, sample, "quant.sf")
      if (file.exists(f)) return(f)
    }
    stop(paste("❌ quant.sf not found for sample:", sample))
  })
  names(files) <- sample_names

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
  if (!all(c("TXNAME","GENEID") %in% names(tx2gene))) {
    names(tx2gene)[1:2] <- c("TXNAME","GENEID")
  }

  # Import with tximport
  txi <- tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene,
    countsFromAbundance = "lengthScaledTPM"
  )

  # Setup colData
  stopifnot(length(sample_names) == length(conditions))
  coldata <- data.frame(condition = factor(conditions), row.names = sample_names)
  if (!is.null(ref_level)) {
    if (!ref_level %in% levels(coldata$condition)) {
      stop("ref_level '", ref_level, "' is not in conditions: ", paste(levels(coldata$condition), collapse=", "))
    }
    coldata$condition <- relevel(coldata$condition, ref = ref_level)
  }
  stopifnot(all(rownames(coldata) == colnames(txi$counts)))

  # Gene filtering
  min_samples <- ceiling(ncol(txi$counts) * min_samples_prop)
  keep <- rowSums(txi$counts >= min_count) >= min_samples
  txi$counts    <- txi$counts[keep, , drop = FALSE]
  txi$abundance <- txi$abundance[keep, , drop = FALSE]
  txi$length    <- txi$length[keep, , drop = FALSE]

  # DESeq2 analysis
  dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)

  # Save DE results
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  write.csv(res_df, file.path(plot_dir, paste0(label, "_DESeq2_results.csv")), row.names = FALSE)

  # Compute rank metric: sign(LFC) * -log10(padj)
  rank_metric <- function(log2fc, pval, padj, gene_ids) {
    eps <- 1e-300
    use_p <- ifelse(is.na(padj), pval, padj)
    use_p[is.na(use_p)] <- 1
    use_p <- pmax(use_p, eps)
    score <- sign(log2fc) * (-log10(use_p))
    names(score) <- gene_ids
    ok <- is.finite(score)
    score <- score[ok]
    score <- tapply(score, names(score), function(v) v[which.max(abs(v))])
    score <- unlist(score, use.names = TRUE)
    sort(score, decreasing = TRUE)
  }

  vec <- rank_metric(res$log2FoldChange, res$pvalue, res$padj, rownames(res))

  # Save RDS + a human-readable CSV
  saveRDS(vec, file = file.path(plot_dir, paste0(label, "_RRHO_vector.rds")))
  rank_tbl <- data.frame(
    gene  = names(vec),
    score = as.numeric(vec),
    rank  = seq_along(vec),
    row.names = NULL
  )
  write.csv(rank_tbl, file.path(plot_dir, paste0(label, "_RRHO_ranked_genes.csv")), row.names = FALSE)

  message("✅ Completed: ", label, "  (", length(vec), " genes in RRHO vector)")
  invisible(vec)
}
