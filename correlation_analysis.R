# ============================================================================
# correlation_analysis.R - Replicate Correlation Analysis
# ============================================================================
#
# Generates Pearson correlation heatmaps for RNA-seq replicate samples.
# Outputs include:
#   - Global correlation matrix across all samples
#   - Within-group correlation matrices for each condition
#   - Pairwise correlation summaries
#
# ============================================================================

run_replicate_correlations <- function(
  base_dirs,
  samples_named_vector,
  output_dir,
  label,
  tx2gene_filename = "tx2gene.tsv",
  transform_method = c("vst", "rlog"),
  blind = TRUE
) {
  suppressMessages({
    library(tximport)
    library(readr)
    library(DESeq2)
    library(pheatmap)
    library(ggplot2)
    library(stats)
    library(utils)
  })

  transform_method <- match.arg(transform_method)
  message("▶ Starting Pearson correlation analysis for label: ", label)

  # Locate quant.sf files
  files <- sapply(names(samples_named_vector), function(sample) {
    for (base in base_dirs) {
      f <- file.path(base, sample, "quant.sf")
      if (file.exists(f)) return(f)
    }
    stop(paste("❌ quant.sf not found for sample:", sample))
  })
  names(files) <- names(samples_named_vector)

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

  # tximport
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

  # DESeq2 setup & transform
  col_data <- data.frame(
    row.names = names(samples_named_vector),
    condition = factor(unname(samples_named_vector))
  )
  dds <- DESeqDataSetFromTximport(txi, colData = col_data, design = ~ condition)

  keep <- rowSums(counts(dds)) > 0
  dds <- dds[keep, ]

  if (transform_method == "vst") {
    dt <- vst(dds, blind = blind)
  } else {
    dt <- rlog(dds, blind = blind)
  }
  mat <- assay(dt)

  # Output directory
  out_dir <- file.path(output_dir, label, "ReplicateCorrelations")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Color palette and breaks
  blue_pal <- colorRampPalette(c("#9ecae1","#4292c6","#08519c"))
  mk_breaks <- function(m) {
    rng <- range(m, na.rm = TRUE)
    if (isTRUE(all.equal(rng[1], rng[2]))) rng <- rng + c(-1e-6, 1e-6)
    seq(rng[1], rng[2], length.out = 101)
  }


  # Global Pearson correlation matrix
  cor_all <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  write.csv(cor_all, file.path(out_dir, "correlation_all_samples.csv"))

  ann <- data.frame(Group = col_data$condition)
  rownames(ann) <- rownames(col_data)
  pdf(file.path(out_dir, "correlation_heatmap_all_samples.pdf"), width = 8, height = 7)
  pheatmap(
    cor_all,
    # annotation_col = ann,
    main = paste0("Pearson Correlation - All Samples"),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    color  = blue_pal(100),
    breaks = mk_breaks(cor_all)
  )
  dev.off()

  # Within-group correlations
  groups <- levels(col_data$condition)
  summary_list <- list()

  for (g in groups) {
    grp_samples <- rownames(col_data)[col_data$condition == g]
    if (length(grp_samples) < 2) {
      message("⚠ Group '", g, "' has <2 samples; skipping within-group correlation.")
      next
    }

    cor_grp <- cor(mat[, grp_samples, drop = FALSE],
                   method = "pearson",
                   use = "pairwise.complete.obs")

    write.csv(cor_grp, file.path(out_dir, paste0("correlation_", g, ".csv")))

    pdf(file.path(out_dir, paste0("correlation_heatmap_", g, ".pdf")), width = 6, height = 5)
    pheatmap(
      cor_grp,
      main = paste0("Within-Group Pearson - ", g),
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      color  = blue_pal(100),
      breaks = mk_breaks(cor_grp)
    )
    dev.off()

    if (ncol(cor_grp) >= 2) {
      upper_idx <- which(upper.tri(cor_grp))
      r_vals <- cor_grp[upper_idx]
      summary_list[[g]] <- data.frame(
        group = g,
        n_reps = length(grp_samples),
        n_pairs = length(r_vals),
        r_mean = mean(r_vals),
        r_median = median(r_vals),
        r_min = min(r_vals),
        r_max = max(r_vals),
        stringsAsFactors = FALSE
      )

      combs <- t(combn(grp_samples, 2))
      pair_df <- data.frame(
        sample_i = combs[,1],
        sample_j = combs[,2],
        r = apply(combs, 1, function(ix) cor_grp[ix[1], ix[2]]),
        row.names = NULL
      )
      write.csv(pair_df,
                file.path(out_dir, paste0("pairwise_r_", g, ".csv")),
                row.names = FALSE)
    }
  }

  if (length(summary_list) > 0) {
    summary_df <- do.call(rbind, summary_list)
    write.csv(summary_df, file.path(out_dir, "within_group_correlation_summary.csv"), row.names = FALSE)
  } else {
    summary_df <- data.frame()
  }

  message("✅ Pearson correlation analysis completed for label: ", label)
  invisible(list(
    transform = transform_method,
    global_correlation = cor_all,
    within_group_summary = summary_df
  ))
}
