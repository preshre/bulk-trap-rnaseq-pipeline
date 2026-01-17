# ============================================================================
# go_analysis.R - GO Enrichment Analysis
# ============================================================================
#
# Performs Gene Ontology enrichment analysis on DEGs.
# Generates bubble plots for BP, MF, and CC ontologies showing
# top enriched terms for upregulated and downregulated genes.
#
# ============================================================================

run_go_analysis <- function(
  base_dirs,
  samples_named_vector,
  output_dir,
  label,
  tx2gene_filename = "tx2gene.tsv",
  log2fc_threshold = 0.58,
  pvalue_cutoff = 0.05,
  use_pvalue = FALSE,
  top_per_direction = 3,
  ref_level = NULL
) {
  suppressMessages({
    library(tximport)
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(DESeq2)
    library(readr)
    library(ggplot2)
    library(dplyr)
  })

  message("▶ Starting GO analysis for label: ", label)

  # Locate inputs
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
      tx2gene <- read_tsv(tx2gene_path, col_names = TRUE)
      break
    }
  }
  if (is.null(tx2gene)) stop("❌ tx2gene.tsv not found in any base directory.")

  # Import & DESeq2
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

  col_data <- data.frame(
    row.names = names(samples_named_vector),
    condition = factor(samples_named_vector)
  )
  
  # Determine reference and alt levels
  lvl <- levels(col_data$condition)
  if (length(lvl) != 2) stop("❌ Expected exactly 2 groups")
  
  if (!is.null(ref_level)) {
    if (!ref_level %in% lvl) stop("❌ ref_level not in levels")
    alt_level <- setdiff(lvl, ref_level)
  } else {
    ref_level <- lvl[1]
    alt_level <- lvl[2]
  }
  col_data$condition <- factor(col_data$condition, levels = c(ref_level, alt_level))

  dds <- DESeqDataSetFromTximport(txi, colData = col_data, design = ~ condition)
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", alt_level, ref_level))

  res_df <- as.data.frame(res)
  res_df <- res_df[!is.na(res_df$pvalue), , drop = FALSE]

  # DEG filters
  thr_fun <- function(df) {
    if (use_pvalue) {
      df$pvalue < pvalue_cutoff & abs(df$log2FoldChange) > log2fc_threshold
    } else {
      df$padj   < pvalue_cutoff & abs(df$log2FoldChange) > log2fc_threshold
    }
  }

  sig <- res_df[thr_fun(res_df), , drop = FALSE]
  up_in_alt_syms <- rownames(sig[sig$log2FoldChange > 0, , drop = FALSE])  # Up in alt (tested)
  up_in_ref_syms <- rownames(sig[sig$log2FoldChange < 0, , drop = FALSE])  # Up in ref (baseline)

  message("• DEGs up in ", alt_level, ": ", length(up_in_alt_syms))
  message("• DEGs up in ", ref_level, ": ", length(up_in_ref_syms))

  if (length(up_in_alt_syms) + length(up_in_ref_syms) < 5) {
    warning("⚠ Too few DEGs after filtering. Skipping GO enrichment.")
    return(NULL)
  }

  # Convert to ENTREZ IDs
  to_entrez <- function(syms) {
    unique(bitr(syms, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID)
  }
  up_alt_eg <- to_entrez(up_in_alt_syms)
  up_ref_eg <- to_entrez(up_in_ref_syms)

  labeled_output_dir <- file.path(output_dir, label, "GO")
  dir.create(labeled_output_dir, recursive = TRUE, showWarnings = FALSE)

  # GO enrichment
  go_ontologies <- c("BP", "MF", "CC")
  go_up_ref <- go_up_alt <- list()

  enrich_one <- function(genes, ont) {
    if (length(genes) < 5) return(NULL)
    enrichGO(
      gene           = genes,
      OrgDb          = org.Mm.eg.db,
      keyType        = "ENTREZID",
      ont            = ont,
      pvalueCutoff   = 0.05,
      pAdjustMethod  = ifelse(use_pvalue, "none", "bonferroni"),
      qvalueCutoff   = ifelse(use_pvalue, 1.0, 0.05),
      readable       = TRUE
    )
  }

  for (ont in go_ontologies) {
    go_up_ref[[ont]] <- enrich_one(up_ref_eg, ont)
    go_up_alt[[ont]] <- enrich_one(up_alt_eg, ont)

    # Save per-direction tables
    for (lst_name in c("ref", "alt")) {
      lst <- if (lst_name == "ref") go_up_ref[[ont]] else go_up_alt[[ont]]
      if (!is.null(lst) && nrow(lst) > 0) {
        nm <- if (lst_name == "ref") paste0("UpIn_", gsub("[/ ]", "_", ref_level))
              else paste0("UpIn_", gsub("[/ ]", "_", alt_level))
        write.csv(
          as.data.frame(lst),
          file = file.path(labeled_output_dir, paste0("GO_", ont, "_", nm, "_results.csv")),
          row.names = FALSE
        )
      }
    }
  }

  # Build per-ontology plots
  pick_top_n <- function(res, n_each) {
    if (is.null(res) || nrow(res) == 0) return(NULL)
    as.data.frame(res) |>
      arrange(p.adjust, desc(Count)) |>
      slice_head(n = n_each)
  }

  for (ont in go_ontologies) {
    top_ref <- pick_top_n(go_up_ref[[ont]], top_per_direction)
    top_alt <- pick_top_n(go_up_alt[[ont]], top_per_direction)

    df <- bind_rows(
      if (!is.null(top_ref)) mutate(top_ref, Direction = paste0("Up in ", ref_level)),
      if (!is.null(top_alt)) mutate(top_alt, Direction = paste0("Up in ", alt_level))
    )

    if (is.null(df) || nrow(df) == 0) {
      message("⚠ No enriched terms to plot for ", ont)
      next
    }

    # Keep only columns we need and order y by p.adjust (smallest at top)
    df <- df |>
      transmute(Term = Description, Direction, Count = Count, p.adjust = p.adjust) |>
      arrange(p.adjust, desc(Count))

    df$Term      <- factor(df$Term, levels = rev(unique(df$Term)))  # top at top of axis
    df$Direction <- factor(df$Direction, levels = c(paste0("Up in ", ref_level),
                                                    paste0("Up in ", alt_level)))

    p <- ggplot(df, aes(x = Direction, y = Term, size = Count, color = p.adjust)) +
      geom_point() +
      theme_minimal(base_size = 12) +
      labs(
        title = paste0(ont, " — Top ", top_per_direction, " per Direction (", label, ")"),
        x = NULL,
        y = paste0("GO ", ont, " term")
      ) +
      scale_color_continuous(name = "p.adjust", trans = "reverse", guide = guide_colorbar(reverse = TRUE)) +
      scale_size_continuous(name = "Count")

    out_pdf <- file.path(labeled_output_dir, paste0("GO_bubbles_", ont, "_", label, ".pdf"))
    ggsave(out_pdf, plot = p, width = 7, height = 5)
    message("✅ Wrote: ", out_pdf)
  }

  message("✅ GO analysis completed for label: ", label)
  invisible(list(up_ref = go_up_ref, up_alt = go_up_alt,
                 ref_level = ref_level, alt_level = alt_level))
}
