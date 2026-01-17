# ============================================================================
# run_pipeline.R - Master Pipeline for RNA-seq Analysis
# ============================================================================
#
# This pipeline generates:
#   1. Correlation analysis (replicate correlations)
#   2. PCA plots
#   3. Volcano plots
#   4. MA plots
#   5. GO enrichment analysis
#   6. DEG analysis with heatmap
#   7. RRHO analysis (optional)
#
# Usage:
#   Rscript run_pipeline.R                    # Run with padj (default)
#   Rscript run_pipeline.R --use-pvalue       # Run with raw p-values
#   Rscript run_pipeline.R --both             # Run both padj and pvalue
#
# ============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
RUN_PVALUE <- "--use-pvalue" %in% args || "--both" %in% args
RUN_PADJ   <- !"--use-pvalue" %in% args || "--both" %in% args

message("\n", paste(rep("=", 70), collapse = ""))
message("RNA-seq Analysis Pipeline")
message(paste(rep("=", 70), collapse = ""))
message("Run with padj:   ", RUN_PADJ)
message("Run with pvalue: ", RUN_PVALUE)
message(paste(rep("=", 70), collapse = ""), "\n")

# Load configuration
source("config.R")

# ============================================================================
# 1. CORRELATION ANALYSIS
# ============================================================================
run_all_correlations <- function() {
  message("\n", paste(rep("-", 50), collapse = ""))
  message("STEP 1: Running Correlation Analysis")
  message(paste(rep("-", 50), collapse = ""))
  
  source("correlation_analysis.R")
  
  for (comp in ALL_COMPARISONS) {
    message("\n▶ Correlation: ", comp$label)
    tryCatch({
      run_replicate_correlations(
        base_dirs = BASE_DIRS,
        samples_named_vector = comp$named_vector,
        output_dir = OUTPUT_DIR,
        label = comp$label,
        tx2gene_filename = TX2GENE_FILENAME,
        transform_method = "vst"
      )
    }, error = function(e) {
      message("❌ Error in correlation analysis for ", comp$label, ": ", e$message)
    })
  }
}

# ============================================================================
# 2. PCA & VOLCANO PLOTS
# ============================================================================
run_all_pca_volcano <- function(use_pvalue) {
  pval_suffix <- if (use_pvalue) "_pvalue" else ""
  
  message("\n", paste(rep("-", 50), collapse = ""))
  message("STEP 2: Running PCA & Volcano Plots (", ifelse(use_pvalue, "p-value", "padj"), ")")
  message(paste(rep("-", 50), collapse = ""))
  
  source("pca_volcano_plot.R")
  
  for (comp in ALL_COMPARISONS) {
    label_out <- paste0(comp$label, pval_suffix)
    message("\n▶ PCA/Volcano: ", label_out)
    
    tryCatch({
      run_analysis(
        samples = comp$samples,
        conditions = comp$conditions,
        label = label_out,
        base_dirs = BASE_DIRS,
        output_dir = OUTPUT_DIR,
        use_pvalue = use_pvalue,
        pvalue_cutoff = if (use_pvalue) PVALUE_CUTOFF else PADJ_CUTOFF,
        log2fc_threshold = LOG2FC_THRESHOLD,
        tx2gene_filename = TX2GENE_FILENAME
      )
    }, error = function(e) {
      message("❌ Error in PCA/Volcano for ", label_out, ": ", e$message)
    })
  }
}

# ============================================================================
# 3. MA PLOTS
# ============================================================================
run_all_ma_plots <- function(use_pvalue) {
  pval_suffix <- if (use_pvalue) "_pvalue" else ""
  
  message("\n", paste(rep("-", 50), collapse = ""))
  message("STEP 3: Running MA Plots (", ifelse(use_pvalue, "p-value", "padj"), ")")
  message(paste(rep("-", 50), collapse = ""))
  
  source("generate_DE_plots.R")
  
  for (comp in ALL_COMPARISONS) {
    label_out <- paste0(comp$label, pval_suffix)
    message("\n▶ MA Plot: ", label_out)
    
    tryCatch({
      run_deseq_analysis_with_ma_plot(
        samples = comp$samples,
        conditions = comp$conditions,
        base_dirs = BASE_DIRS,
        tx2gene_filename = TX2GENE_FILENAME,
        output_dir = file.path(OUTPUT_DIR, label_out, "MA"),
        comparison_label = label_out,
        lfc_cutoff = LOG2FC_THRESHOLD,
        padj_cutoff = PADJ_CUTOFF,
        use_pvalue = use_pvalue,
        pvalue_cutoff = PVALUE_CUTOFF
      )
    }, error = function(e) {
      message("❌ Error in MA plot for ", label_out, ": ", e$message)
    })
  }
}

# ============================================================================
# 4. GO ENRICHMENT
# ============================================================================
run_all_go_analysis <- function(use_pvalue) {
  pval_suffix <- if (use_pvalue) "_pvalue" else ""
  
  message("\n", paste(rep("-", 50), collapse = ""))
  message("STEP 4: Running GO Enrichment (", ifelse(use_pvalue, "p-value", "padj"), ")")
  message(paste(rep("-", 50), collapse = ""))
  
  source("go_analysis.R")
  
  for (comp in ALL_COMPARISONS) {
    label_out <- paste0(comp$label, pval_suffix)
    message("\n▶ GO Analysis: ", label_out)
    
    tryCatch({
      run_go_analysis(
        base_dirs = BASE_DIRS,
        samples_named_vector = comp$named_vector,
        output_dir = OUTPUT_DIR,
        label = label_out,
        tx2gene_filename = TX2GENE_FILENAME,
        log2fc_threshold = LOG2FC_THRESHOLD,
        pvalue_cutoff = if (use_pvalue) PVALUE_CUTOFF else PADJ_CUTOFF,
        use_pvalue = use_pvalue,
        top_per_direction = 3,
        ref_level = comp$ref_level
      )
    }, error = function(e) {
      message("❌ Error in GO analysis for ", label_out, ": ", e$message)
    })
  }
}

# ============================================================================
# 5. DEG ANALYSIS WITH HEATMAP
# ============================================================================
run_all_deg_analysis <- function(use_pvalue) {
  pval_suffix <- if (use_pvalue) "_pvalue" else ""
  
  message("\n", paste(rep("-", 50), collapse = ""))
  message("STEP 5: Running DEG Analysis (", ifelse(use_pvalue, "p-value", "padj"), ")")
  message(paste(rep("-", 50), collapse = ""))
  
  source("deg_analysis.R")
  
  for (comp in ALL_COMPARISONS) {
    label_out <- paste0(comp$label, pval_suffix)
    message("\n▶ DEG Analysis: ", label_out)
    
    tryCatch({
      run_deg_analysis(
        samples = comp$samples,
        conditions = comp$conditions,
        base_dirs = BASE_DIRS,
        tx2gene_filename = TX2GENE_FILENAME,
        output_dir = OUTPUT_DIR,
        label = label_out,
        use_pvalue = use_pvalue,
        pvalue_cutoff = PVALUE_CUTOFF,
        padj_cutoff = PADJ_CUTOFF,
        log2fc_threshold = LOG2FC_THRESHOLD,
        baseMean_min = 10,
        top_n_heatmap = 30
      )
    }, error = function(e) {
      message("❌ Error in DEG analysis for ", label_out, ": ", e$message)
    })
  }
}

# ============================================================================
# 6. RRHO ANALYSIS
# ============================================================================

# Step 6a: Generate RRHO vectors for each comparison
run_all_rrho_vectors <- function(use_pvalue) {
  pval_suffix <- if (use_pvalue) "_pvalue" else ""
  
  message("\n", paste(rep("-", 50), collapse = ""))
  message("STEP 6a: Generating RRHO Vectors (", ifelse(use_pvalue, "p-value", "padj"), ")")
  message(paste(rep("-", 50), collapse = ""))
  
  source("rrho_analysis.R")
  
  for (comp in ALL_COMPARISONS) {
    label_out <- paste0(comp$label, pval_suffix)
    message("\n▶ RRHO Vector: ", label_out)
    
    tryCatch({
      run_deseq_to_rrho_vector(
        sample_names = comp$samples,
        conditions = comp$conditions,
        base_dirs = BASE_DIRS,
        tx2gene_filename = TX2GENE_FILENAME,
        label = label_out,
        output_dir = OUTPUT_DIR,
        ref_level = comp$ref_level
      )
    }, error = function(e) {
      message("❌ Error generating RRHO vector for ", label_out, ": ", e$message)
    })
  }
}

# Step 6b: Compare RRHO vectors between experiments
run_all_rrho_comparisons <- function(use_pvalue) {
  pval_suffix <- if (use_pvalue) "_pvalue" else ""
  
  message("\n", paste(rep("-", 50), collapse = ""))
  message("STEP 6b: Running RRHO Comparisons (", ifelse(use_pvalue, "p-value", "padj"), ")")
  message(paste(rep("-", 50), collapse = ""))
  
  source("rrho_compute.R")
  
  # Define RRHO comparison pairs
  rrho_pairs <- list(
    list(
      comp1 = "vTRAP_BO_vTRAP_Trained",
      comp2 = "iPKRTRAP_BO_iPKRTRAP_Trained",
      label1 = "vTRAP",
      label2 = "iPKR"
    ),
    list(
      comp1 = "vTRAP_BO_vTRAP_Trained",
      comp2 = "vTRAP_BO_PT4h_vTRAP_Trained_PT4h",
      label1 = "vTRAP_PT1h",
      label2 = "vTRAP_PT4h"
    )
  )
  
  for (pair in rrho_pairs) {
    label_out <- paste0(pair$label1, "_vs_", pair$label2, pval_suffix)
    message("\n▶ RRHO Comparison: ", label_out)
    
    vec1_path <- file.path(OUTPUT_DIR, paste0(pair$comp1, pval_suffix), "RRHO", 
                           paste0(pair$comp1, pval_suffix, "_RRHO_vector.rds"))
    vec2_path <- file.path(OUTPUT_DIR, paste0(pair$comp2, pval_suffix), "RRHO", 
                           paste0(pair$comp2, pval_suffix, "_RRHO_vector.rds"))
    
    if (file.exists(vec1_path) && file.exists(vec2_path)) {
      tryCatch({
        run_rrho2_analysis(
          vec1_path = vec1_path,
          vec2_path = vec2_path,
          output_dir = file.path(OUTPUT_DIR, "RRHO", label_out),
          label1 = pair$label1,
          label2 = pair$label2
        )
      }, error = function(e) {
        message("❌ Error in RRHO comparison for ", label_out, ": ", e$message)
      })
    } else {
      message("⚠ RRHO vectors not found for ", label_out, " - skipping")
      message("  Expected: ", vec1_path)
      message("  Expected: ", vec2_path)
    }
  }
}

# Combined function to run both steps
run_all_rrho <- function(use_pvalue) {
  # First generate vectors for all comparisons
  run_all_rrho_vectors(use_pvalue)
  
  # Then compare the vectors
  run_all_rrho_comparisons(use_pvalue)
}

# ============================================================================
# MAIN PIPELINE
# ============================================================================
main <- function() {
  start_time <- Sys.time()
  
  # Step 1: Correlation (run once, independent of p-value choice)
  run_all_correlations()
  
  # Steps 2-6: Run for padj if requested
  if (RUN_PADJ) {
    message("\n", paste(rep("=", 70), collapse = ""))
    message("Running analyses with ADJUSTED P-VALUE (padj)")
    message(paste(rep("=", 70), collapse = ""))
    
    run_all_pca_volcano(use_pvalue = FALSE)
    run_all_ma_plots(use_pvalue = FALSE)
    run_all_go_analysis(use_pvalue = FALSE)
    run_all_deg_analysis(use_pvalue = FALSE)
    run_all_rrho(use_pvalue = FALSE) 
  }
  
  # Steps 2-6: Run for p-value if requested
  if (RUN_PVALUE) {
    message("\n", paste(rep("=", 70), collapse = ""))
    message("Running analyses with RAW P-VALUE")
    message(paste(rep("=", 70), collapse = ""))
    
    run_all_pca_volcano(use_pvalue = TRUE)
    run_all_ma_plots(use_pvalue = TRUE)
    run_all_go_analysis(use_pvalue = TRUE)
    run_all_deg_analysis(use_pvalue = TRUE)
    run_all_rrho(use_pvalue = TRUE)
  }
  
  # Summary
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  message("\n", paste(rep("=", 70), collapse = ""))
  message("PIPELINE COMPLETE")
  message(paste(rep("=", 70), collapse = ""))
  message("Output directory: ", OUTPUT_DIR)
  message("Total runtime: ", round(duration, 2), " minutes")
  message("Comparisons processed: ", length(ALL_COMPARISONS))
  message(paste(rep("=", 70), collapse = ""), "\n")
}

# Run the pipeline
main()
