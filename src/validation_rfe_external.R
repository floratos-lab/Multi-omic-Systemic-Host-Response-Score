# =============================================================================
# Script:  validation_rfe_external.R
# Purpose: Apply the MoSS gene signature to external validation cohorts.
#          Handles flexible input formats, cross-cohort gene ID harmonization,
#          and standardized output reporting.
# NOTE:    External validation cohort must have all 14-genes present to proceed.
#          MARCHF2, LINC00861, RALB, LRP10, IL17RA, GMPR, GLT1D1, RAB32, MCTP1,
#          NFIL3, TCF7, ALOX5, MXD1, FCMR.
# =============================================================================
# Functions defined:
#
#   load_validation_data(config)
#     Dispatches to a format-specific loader based on config$data_source
#     ("excel", "csv", "tsv", "geo", "rds"). Returns a named list with
#     elements: counts (genes × samples matrix), metadata (data frame with
#     sample_id and cohort), and gene_id_type.
#
#       load_from_excel(config)  — reads count matrix from an Excel file using
#         user-specified gene ID column and count start column.
#       load_from_csv(config)    — same structure for CSV input.
#       load_from_tsv(config)    — same structure for tab-delimited input.
#       load_from_geo(config)    — downloads from GEO via GEOquery, or uses
#         a pre-cached RDS count matrix and optional metadata CSV.
#       load_from_rds(config)    — accepts a bare matrix/data frame or a named
#         list with $counts and $metadata elements.
#
#   standardize_gene_ids(counts, id_type)
#     Strips Ensembl version suffixes (e.g., ".1") when id_type is "ensembl".
#     Resolves duplicate row IDs by retaining the row with the highest mean
#     expression across samples. Returns the cleaned count matrix.
#
#   map_model_genes_enhanced(model_coef_df, target_id_type, count_gene_ids)
#     Maps model gene symbols and Entrez IDs to the ID space used by the
#     validation cohort (target_id_type: "symbol", "entrez", or "ensembl").
#     Uses multi-stage fallback strategies:
#       symbol  → direct matching
#       entrez  → provided Entrez ID; fallback: symbol → Entrez via
#                 org.Hs.eg.db
#       ensembl → Entrez → Ensembl; fallback: symbol → Ensembl via
#                 org.Hs.eg.db
#     Returns a data frame recording each gene's mapped ID, mapping method,
#     and whether it was found in the count matrix.
#
#   create_validation_coefficients(mapping_results, model_coefficients)
#     Assembles the final coefficient table in the validation cohort's ID
#     space, preserving the intercept from the original model. Returns a data
#     frame with columns: feature, coefficient, gene_symbol.
#
#   preprocess_rnaseq_mofa_style(raw_mat, top_n, scale_views, keep_genes)
#     Applies count threshold filtering (≥2 counts in ≥2 samples), VST
#     normalization, zero-variance row removal, and retention of the top_n
#     highest-variance genes. If keep_genes is provided, those genes are
#     force-included regardless of variance rank, ensuring signature genes
#     survive the HVG filter. Row-scales if scale_views = TRUE. Returns a
#     scaled genes × samples matrix.
#
#   perform_qc(counts, model_genes, cohort_name)
#     Reports dataset dimensions, library size statistics (mean, median,
#     range), gene detection rates, and signature gene availability.
#     Returns a list with lib_sizes, genes_detected, genes_present, and
#     genes_missing.
#
#   plot_qc(counts, qc_results, cohort_name, output_file)
#     Saves a 4-panel PNG: library size bar chart, library size histogram,
#     gene detection histogram, and PCA of the top 500 variable genes
#     (log2-transformed counts).
#
#   predict_moss_scores(preprocessed_mat, model_coefs)
#     Computes per-sample MoSS scores as intercept + weighted sum of available
#     signature gene expression values. Reports any genes missing from the
#     preprocessed matrix. Returns a named numeric vector of MoSS scores.
#
#   plot_validation_results(moss_scores, metadata, cohort_name, output_file)
#     Prints summary statistics (N, mean, SD, median, IQR, range) and saves
#     a 4-panel PNG: histogram, density plot, boxplot with jitter, and Q-Q
#     plot. Returns a data frame of predictions joined to metadata.
#
#   run_validation(validation_config, model_coefficients, output_dir, 
#                  save_log)
#     Main orchestration function. Executes the full pipeline: load → ID
#     standardization → gene mapping → QC → preprocessing → prediction → output.
#
#     validation_config fields: cohort_name, data_source, file_path or geo_id,
#       gene_id_column, count_start_column, gene_id_type, metadata_file
#       (optional), counts_file (optional, for geo).
#     model_coefficients: data frame with columns gene_symbol, entrez_id,
#       coefficient — produced by rfe_gene_selection.R.
#
#     Outputs written to output_dir/{cohort_name}/:
#       {cohort}_gene_mapping.csv, {cohort}_qc.png,
#       {cohort}_moss_distribution.png, {cohort}_moss_predictions.csv,
#       {cohort}_final_coefficients.csv, {cohort}_summary.csv,
#       {cohort}_validation_log.txt (if save_log = TRUE).
#
# Dependencies: DESeq2, ggplot2, pheatmap, gridExtra, readxl, org.Hs.eg.db,
#               AnnotationDbi, GEOquery, tidyverse, methods.R
#   (methods.R provides: normalizeRNASEQwithVST(), removeNoVarRows())
# Note: dplyr::select() conflicts may arise with AnnotationDbi; use
#       dplyr::select() explicitly if needed.
# =============================================================================

# Required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(gridExtra)
  library(readxl)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GEOquery)
  library(tidyverse)  # Load last to ensure dplyr functions take precedence
})

# Note: If you encounter select() function conflicts, use dplyr::select() explicitly

# Source preprocessing functions
source("methods.R")

################################################################################
# DATA LOADING FUNCTIONS
################################################################################

#' Load validation cohort data from various sources
load_validation_data <- function(config) {

  cat("\n================================================================================\n")
  cat("LOADING VALIDATION DATA:", config$cohort_name, "\n")
  cat("================================================================================\n\n")

  data_source <- config$data_source

  if (data_source == "excel") {
    result <- load_from_excel(config)
  } else if (data_source == "csv") {
    result <- load_from_csv(config)
  } else if (data_source == "tsv") {
    result <- load_from_tsv(config)
  } else if (data_source == "geo") {
    result <- load_from_geo(config)
  } else if (data_source == "rds") {
    result <- load_from_rds(config)
  } else {
    stop("Unknown data source: ", data_source)
  }

  cat(sprintf("✓ Loaded data: %d genes × %d samples\n",
              nrow(result$counts), ncol(result$counts)))

  return(result)
}

#' Load data from Excel file
load_from_excel <- function(config) {
  cat("Loading from Excel file:", config$file_path, "\n")

  raw_data <- read_excel(config$file_path)
  gene_ids <- raw_data[[config$gene_id_column]]
  count_start_col <- config$count_start_column
  counts <- as.matrix(raw_data[, count_start_col:ncol(raw_data)])
  rownames(counts) <- gene_ids
  mode(counts) <- "numeric"

  valid_genes <- !is.na(gene_ids) & gene_ids != ""
  counts <- counts[valid_genes, ]

  metadata <- data.frame(
    sample_id = colnames(counts),
    cohort = config$cohort_name,
    stringsAsFactors = FALSE
  )

  return(list(
    counts = counts,
    metadata = metadata,
    gene_id_type = config$gene_id_type
  ))
}

#' Load data from CSV file
load_from_csv <- function(config) {
  cat("Loading from CSV file:", config$file_path, "\n")

  raw_data <- read.csv(config$file_path, stringsAsFactors = FALSE)
  gene_ids <- raw_data[[config$gene_id_column]]
  count_start_col <- config$count_start_column
  counts <- as.matrix(raw_data[, count_start_col:ncol(raw_data)])
  rownames(counts) <- gene_ids
  mode(counts) <- "numeric"

  valid_genes <- !is.na(gene_ids) & gene_ids != ""
  counts <- counts[valid_genes, ]

  metadata <- data.frame(
    sample_id = colnames(counts),
    cohort = config$cohort_name,
    stringsAsFactors = FALSE
  )

  return(list(
    counts = counts,
    metadata = metadata,
    gene_id_type = config$gene_id_type
  ))
}

#' Load data from TSV file
load_from_tsv <- function(config) {
  cat("Loading from TSV file:", config$file_path, "\n")

  raw_data <- read.delim(config$file_path, stringsAsFactors = FALSE)
  gene_ids <- raw_data[[config$gene_id_column]]
  count_start_col <- config$count_start_column
  counts <- as.matrix(raw_data[, count_start_col:ncol(raw_data)])
  rownames(counts) <- gene_ids
  mode(counts) <- "numeric"

  valid_genes <- !is.na(gene_ids) & gene_ids != ""
  counts <- counts[valid_genes, ]

  metadata <- data.frame(
    sample_id = colnames(counts),
    cohort = config$cohort_name,
    stringsAsFactors = FALSE
  )

  return(list(
    counts = counts,
    metadata = metadata,
    gene_id_type = config$gene_id_type
  ))
}

#' Load data from GEO
load_from_geo <- function(config) {
  cat("Loading from GEO:", config$geo_id, "\n")

  if (!is.null(config$counts_file) && file.exists(config$counts_file)) {
    cat("  Using pre-downloaded count matrix:", config$counts_file, "\n")
    counts <- readRDS(config$counts_file)
  } else {
    cat("  Downloading from GEO...\n")
    gse <- getGEO(config$geo_id, GSEMatrix = TRUE)
    expr_mat <- exprs(gse[[1]])
    counts <- expr_mat
  }

  if (!is.null(config$metadata_file) && file.exists(config$metadata_file)) {
    cat("  Using provided metadata:", config$metadata_file, "\n")
    metadata <- read.csv(config$metadata_file, stringsAsFactors = FALSE)
  } else {
    cat("  Extracting metadata from GEO...\n")
    gse <- getGEO(config$geo_id, GSEMatrix = TRUE)
    metadata <- pData(gse[[1]])
    metadata$sample_id <- rownames(metadata)
    metadata$cohort <- config$cohort_name
  }

  return(list(
    counts = counts,
    metadata = metadata,
    gene_id_type = config$gene_id_type
  ))
}

#' Load data from RDS file
load_from_rds <- function(config) {
  cat("Loading from RDS file:", config$file_path, "\n")

  data <- readRDS(config$file_path)

  if (is.matrix(data) || is.data.frame(data)) {
    counts <- as.matrix(data)
    metadata <- data.frame(
      sample_id = colnames(counts),
      cohort = config$cohort_name,
      stringsAsFactors = FALSE
    )
  } else if (is.list(data)) {
    counts <- data$counts
    metadata <- data$metadata
  } else {
    stop("Unknown RDS data structure")
  }

  return(list(
    counts = counts,
    metadata = metadata,
    gene_id_type = config$gene_id_type
  ))
}

################################################################################
# ENHANCED GENE ID HANDLING
################################################################################

#' Standardize gene IDs (remove version numbers, handle duplicates)
standardize_gene_ids <- function(counts, id_type) {

  cat("\nStandardizing gene IDs (type:", id_type, ")...\n")

  # Remove version numbers for Ensembl IDs
  if (id_type == "ensembl") {
    cat("  Removing Ensembl version numbers...\n")
    rownames(counts) <- gsub("\\..*", "", rownames(counts))
  }

  # Handle duplicate IDs
  if (any(duplicated(rownames(counts)))) {
    cat("  Resolving duplicate gene IDs...\n")
    dup_ids <- unique(rownames(counts)[duplicated(rownames(counts))])
    cat(sprintf("    Found %d duplicate IDs\n", length(dup_ids)))

    for (id in dup_ids) {
      dup_idx <- which(rownames(counts) == id)
      mean_expr <- rowMeans(counts[dup_idx, , drop = FALSE])
      keep_idx <- dup_idx[which.max(mean_expr)]
      remove_idx <- dup_idx[dup_idx != keep_idx]
      counts <- counts[-remove_idx, ]
    }
    cat("    ✓ Duplicates resolved\n")
  }

  return(counts)
}

#' Enhanced gene mapping using both symbols and Entrez IDs
#'
#' @param model_coef_df Data frame with gene_symbol, entrez_id, coefficient
#' @param target_id_type Target ID type in validation cohort
#' @param count_gene_ids Gene IDs present in count matrix
#' @return Data frame with successful mappings
map_model_genes_enhanced <- function(model_coef_df, target_id_type, count_gene_ids) {

  cat("\n================================================================================\n")
  cat("ENHANCED GENE MAPPING\n")
  cat("================================================================================\n\n")

  # Extract model genes (exclude intercept)
  model_genes <- model_coef_df %>% filter(gene_symbol != "(Intercept)")
  n_genes <- nrow(model_genes)

  cat(sprintf("Model genes to map: %d\n", n_genes))
  cat(sprintf("Target ID type: %s\n", target_id_type))
  cat(sprintf("Genes in count matrix: %d\n\n", length(count_gene_ids)))

  # Initialize mapping results
  mapping_results <- data.frame(
    gene_symbol = model_genes$gene_symbol,
    entrez_id = model_genes$entrez_id,
    coefficient = model_genes$coefficient,
    mapped_id = NA_character_,
    mapping_method = NA_character_,
    in_count_matrix = FALSE,
    stringsAsFactors = FALSE
  )

  if (target_id_type == "symbol") {
    # Direct symbol matching
    cat("Strategy: Direct symbol matching\n")
    mapping_results$mapped_id <- mapping_results$gene_symbol
    mapping_results$mapping_method <- "direct_symbol"
    mapping_results$in_count_matrix <- mapping_results$mapped_id %in% count_gene_ids

  } else if (target_id_type == "entrez") {
    # Use provided Entrez IDs
    cat("Strategy: Using provided Entrez IDs\n")
    mapping_results$mapped_id <- mapping_results$entrez_id
    mapping_results$mapping_method <- "provided_entrez"
    mapping_results$in_count_matrix <- mapping_results$mapped_id %in% count_gene_ids

    # For genes not found, try mapping symbol to entrez
    not_found <- !mapping_results$in_count_matrix
    if (any(not_found)) {
      cat(sprintf("  %d genes not found, trying symbol->entrez mapping...\n", sum(not_found)))
      alternative_entrez <- mapIds(
        org.Hs.eg.db,
        keys = mapping_results$gene_symbol[not_found],
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
      )

      # Update if alternative mapping found in count matrix
      for (i in which(not_found)) {
        alt_id <- alternative_entrez[mapping_results$gene_symbol[i]]
        if (!is.na(alt_id) && alt_id %in% count_gene_ids) {
          mapping_results$mapped_id[i] <- alt_id
          mapping_results$mapping_method[i] <- "symbol_to_entrez"
          mapping_results$in_count_matrix[i] <- TRUE
        }
      }
    }

  } else if (target_id_type == "ensembl") {
    # Try multiple strategies for Ensembl
    cat("Strategy: Multi-stage mapping to Ensembl\n")

    # Strategy 1: Map from provided Entrez IDs
    cat("  Stage 1: Entrez -> Ensembl\n")
    ensembl_from_entrez <- mapIds(
      org.Hs.eg.db,
      keys = mapping_results$entrez_id,
      column = "ENSEMBL",
      keytype = "ENTREZID",
      multiVals = "first"
    )

    mapping_results$mapped_id <- ensembl_from_entrez
    mapping_results$mapping_method <- "entrez_to_ensembl"
    mapping_results$in_count_matrix <- mapping_results$mapped_id %in% count_gene_ids

    # Strategy 2: For unmapped, try symbol -> ensembl
    not_found <- is.na(mapping_results$mapped_id) | !mapping_results$in_count_matrix
    if (any(not_found)) {
      cat(sprintf("  Stage 2: Symbol -> Ensembl for %d genes\n", sum(not_found)))
      ensembl_from_symbol <- mapIds(
        org.Hs.eg.db,
        keys = mapping_results$gene_symbol[not_found],
        column = "ENSEMBL",
        keytype = "SYMBOL",
        multiVals = "first"
      )

      for (i in which(not_found)) {
        alt_id <- ensembl_from_symbol[mapping_results$gene_symbol[i]]
        if (!is.na(alt_id) && alt_id %in% count_gene_ids) {
          mapping_results$mapped_id[i] <- alt_id
          mapping_results$mapping_method[i] <- "symbol_to_ensembl"
          mapping_results$in_count_matrix[i] <- TRUE
        }
      }
    }

  } else {
    stop("Unknown target ID type: ", target_id_type)
  }

  # Summary
  cat("\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("MAPPING SUMMARY\n")
  cat(paste(rep("-", 80), collapse = ""), "\n\n")

  n_mapped <- sum(!is.na(mapping_results$mapped_id))
  n_in_matrix <- sum(mapping_results$in_count_matrix)
  n_failed <- n_genes - n_in_matrix

  cat(sprintf("Total model genes: %d\n", n_genes))
  cat(sprintf("Successfully mapped to %s: %d (%.1f%%)\n",
              target_id_type, n_mapped, 100 * n_mapped / n_genes))
  cat(sprintf("Found in count matrix: %d (%.1f%%)\n",
              n_in_matrix, 100 * n_in_matrix / n_genes))
  cat(sprintf("Failed to map: %d (%.1f%%)\n\n",
              n_failed, 100 * n_failed / n_genes))

  if (n_failed > 0) {
    failed_genes <- mapping_results$gene_symbol[!mapping_results$in_count_matrix]
    cat("Failed genes:\n")
    cat(sprintf("  %s\n", paste(failed_genes, collapse = ", ")))
  }

  # Show mapping methods used
  method_table <- table(mapping_results$mapping_method[mapping_results$in_count_matrix])
  if (length(method_table) > 0) {
    cat("\nMapping methods used:\n")
    for (method in names(method_table)) {
      cat(sprintf("  %s: %d genes\n", method, method_table[method]))
    }
  }

  cat("\n")

  return(mapping_results)
}

#' Create model coefficients for validation cohort
create_validation_coefficients <- function(mapping_results, model_coefficients) {

  # Get intercept from original model
  intercept_coef <- model_coefficients$coefficient[model_coefficients$gene_symbol == "(Intercept)"]

  # Include intercept
  intercept_row <- data.frame(
    feature = "(Intercept)",
    coefficient = intercept_coef,
    gene_symbol = "(Intercept)",
    stringsAsFactors = FALSE
  )

  # Add successfully mapped genes
  mapped_genes <- mapping_results %>%
    filter(in_count_matrix) %>%
    dplyr::select(gene_symbol, mapped_id, coefficient) %>%
    rename(feature = mapped_id) %>%
    dplyr::select(feature, coefficient, gene_symbol)

  validation_coefs <- rbind(
    intercept_row,
    mapped_genes
  )

  return(validation_coefs)
}

################################################################################
# PREPROCESSING
################################################################################

preprocess_rnaseq_mofa_style <- function(raw_mat, top_n = 4000, scale_views = TRUE,
                                         keep_genes = NULL) {
  thresh <- 2
  keep_genes_thresh <- rowSums(raw_mat >= thresh) >= 2
  raw_mat <- raw_mat[keep_genes_thresh, ]

  vst <- normalizeRNASEQwithVST(raw_mat)
  vst <- removeNoVarRows(vst)

  gene_vars <- sort(apply(vst, 1, var), decreasing = TRUE)
  top_genes <- names(gene_vars)[1:min(top_n, length(gene_vars))]

  if (!is.null(keep_genes)) {
    keep_genes <- keep_genes[keep_genes %in% rownames(vst)]
    top_genes <- unique(c(top_genes, keep_genes))
  }

  vst <- vst[top_genes, ]

  if (scale_views) {
    vst <- t(scale(t(vst)))
  }

  return(vst)
}

################################################################################
# QC FUNCTIONS
################################################################################

perform_qc <- function(counts, model_genes, cohort_name) {

  cat("\n================================================================================\n")
  cat("DATA QUALITY CONTROL:", cohort_name, "\n")
  cat("================================================================================\n\n")

  cat("Dataset dimensions:\n")
  cat(sprintf("  Genes: %d\n", nrow(counts)))
  cat(sprintf("  Samples: %d\n", ncol(counts)))

  lib_sizes <- colSums(counts, na.rm = TRUE)
  cat(sprintf("\nLibrary sizes:\n"))
  cat(sprintf("  Mean: %.1f million\n", mean(lib_sizes, na.rm = TRUE) / 1e6))
  cat(sprintf("  Median: %.1f million\n", median(lib_sizes, na.rm = TRUE) / 1e6))
  cat(sprintf("  Range: %.1f - %.1f million\n",
              min(lib_sizes, na.rm = TRUE) / 1e6,
              max(lib_sizes, na.rm = TRUE) / 1e6))

  genes_detected <- rowSums(counts > 0)
  cat(sprintf("\nGene detection:\n"))
  cat(sprintf("  Genes with counts in all samples: %d\n",
              sum(genes_detected == ncol(counts))))
  cat(sprintf("  Genes with counts in >50%% samples: %d\n",
              sum(genes_detected > ncol(counts)/2)))

  cat(sprintf("\nSignature gene availability:\n"))
  genes_present <- model_genes[model_genes %in% rownames(counts)]
  genes_missing <- model_genes[!model_genes %in% rownames(counts)]

  cat(sprintf("  Present: %d/%d\n", length(genes_present), length(model_genes)))
  if (length(genes_missing) > 0) {
    cat(sprintf("  Missing: %s\n", paste(genes_missing, collapse = ", ")))
  }

  return(list(
    lib_sizes = lib_sizes,
    genes_detected = genes_detected,
    genes_present = genes_present,
    genes_missing = genes_missing
  ))
}

plot_qc <- function(counts, qc_results, cohort_name, output_file) {

  cat("\nGenerating QC plots...\n")

  lib_sizes <- qc_results$lib_sizes / 1e6
  lib_sizes <- lib_sizes[!is.na(lib_sizes)]

  sample_data <- data.frame(
    sample = names(lib_sizes),
    lib_size = lib_sizes,
    stringsAsFactors = FALSE
  )

  p1 <- ggplot(sample_data, aes(x = reorder(sample, lib_size), y = lib_size)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_hline(yintercept = median(lib_sizes, na.rm = TRUE),
               linetype = "dashed", color = "red") +
    labs(title = paste("Library Sizes -", cohort_name),
         x = "Sample", y = "Library Size (millions)") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p2 <- ggplot(sample_data, aes(x = lib_size)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    geom_vline(xintercept = median(lib_sizes, na.rm = TRUE),
               linetype = "dashed", color = "red", size = 1) +
    labs(title = "Library Size Distribution",
         x = "Library Size (millions)", y = "Number of Samples") +
    theme_bw()

  detection_data <- data.frame(samples_detected = qc_results$genes_detected)
  p3 <- ggplot(detection_data, aes(x = samples_detected)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    labs(title = "Gene Detection Distribution",
         x = "Number of Samples with Counts > 0", y = "Number of Genes") +
    theme_bw()

  counts_log <- log2(as.matrix(counts) + 1)
  gene_vars <- apply(counts_log, 1, var, na.rm = TRUE)
  counts_log <- counts_log[gene_vars > 0 & !is.na(gene_vars), ]
  top_var_genes <- names(sort(gene_vars[gene_vars > 0 & !is.na(gene_vars)],
                              decreasing = TRUE)[1:min(500, sum(gene_vars > 0))])
  pca_result <- prcomp(t(counts_log[top_var_genes, ]), scale. = TRUE, center = TRUE)

  pca_data <- data.frame(
    sample = rownames(pca_result$x),
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2]
  )

  var_explained <- summary(pca_result)$importance[2, 1:2] * 100

  p4 <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    labs(title = "PCA of Samples",
         x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
         y = sprintf("PC2 (%.1f%% variance)", var_explained[2])) +
    theme_bw()

  png(output_file, width = 3000, height = 3000, res = 300)
  grid.arrange(p1, p2, p3, p4, ncol = 2)
  dev.off()

  cat(sprintf("✓ QC plots saved: %s\n", output_file))
}

################################################################################
# PREDICTION AND PLOTTING
################################################################################

predict_moss_scores <- function(preprocessed_mat, model_coefs) {

  gene_coefs <- model_coefs %>% filter(feature != "(Intercept)")
  intercept <- model_coefs %>% filter(feature == "(Intercept)") %>% pull(coefficient)

  available_genes <- gene_coefs$feature[gene_coefs$feature %in% rownames(preprocessed_mat)]
  missing_genes <- setdiff(gene_coefs$feature, rownames(preprocessed_mat))

  cat(sprintf("\nPrediction gene availability: %d/%d\n",
              length(available_genes), nrow(gene_coefs)))
  if (length(missing_genes) > 0) {
    cat("Missing:", paste(missing_genes, collapse = ", "), "\n")
  }

  gene_expr <- preprocessed_mat[available_genes, , drop = FALSE]
  gene_weights <- gene_coefs %>%
    filter(feature %in% available_genes) %>%
    pull(coefficient)

  moss_scores <- as.vector(intercept + t(gene_weights) %*% gene_expr)
  names(moss_scores) <- colnames(preprocessed_mat)

  return(moss_scores)
}

plot_validation_results <- function(moss_scores, metadata, cohort_name, output_file) {

  results_df <- data.frame(
    sample = names(moss_scores),
    moss_score = moss_scores,
    stringsAsFactors = FALSE
  )

  if (!is.null(metadata) && "sample_id" %in% colnames(metadata)) {
    results_df <- results_df %>%
      left_join(metadata, by = c("sample" = "sample_id"))
  }

  cat("\n================================================================================\n")
  cat("VALIDATION RESULTS:", cohort_name, "\n")
  cat("================================================================================\n\n")
  cat(sprintf("  N samples: %d\n", length(moss_scores)))
  cat(sprintf("  Mean MoSS: %.4f (SD: %.4f)\n", mean(moss_scores), sd(moss_scores)))
  cat(sprintf("  Median: %.4f (IQR: %.4f)\n", median(moss_scores), IQR(moss_scores)))
  cat(sprintf("  Range: %.4f to %.4f\n\n", min(moss_scores), max(moss_scores)))

  # Plots
  p1 <- ggplot(results_df, aes(x = moss_score)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = mean(moss_scores), color = "red",
               linetype = "dashed", size = 1) +
    labs(title = paste("MoSS Distribution -", cohort_name),
         subtitle = sprintf("Mean: %.3f, SD: %.3f", mean(moss_scores), sd(moss_scores)),
         x = "MoSS Score", y = "Count") +
    theme_bw()

  p2 <- ggplot(results_df, aes(x = moss_score)) +
    geom_density(fill = "steelblue", alpha = 0.5, color = "black") +
    geom_vline(xintercept = mean(moss_scores), color = "red",
               linetype = "dashed", size = 1) +
    labs(title = "MoSS Density", x = "MoSS Score", y = "Density") +
    theme_bw()

  p3 <- ggplot(results_df, aes(y = moss_score, x = cohort_name)) +
    geom_boxplot(fill = "steelblue", alpha = 0.7) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
    labs(title = "MoSS Distribution", y = "MoSS Score", x = "") +
    theme_bw()

  p4 <- ggplot(results_df, aes(sample = moss_score)) +
    stat_qq(color = "steelblue", size = 2, alpha = 0.6) +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(title = "Q-Q Plot", x = "Theoretical", y = "Sample") +
    theme_bw()

  png(output_file, width = 3000, height = 3000, res = 300)
  grid.arrange(p1, p2, p3, p4, ncol = 2)
  dev.off()

  cat(sprintf("✓ Plots saved: %s\n", output_file))

  return(results_df)
}

################################################################################
# MAIN EXTERNAL VALIDATION FUNCTION
################################################################################

#' Run external validation of the MoSS gene signature
#'
#' @param validation_config  Named list describing the validation cohort data
#'   (cohort_name, data_source, file_path / geo_id, gene_id_column,
#'   count_start_column, gene_id_type, optionally metadata_file / counts_file).
#' @param model_coefficients Data frame with columns gene_symbol, entrez_id,
#'   coefficient.
#' @param output_dir         Directory for output files (created if needed).
#' @param save_log           Logical; if TRUE, tee console output to a log file.
run_validation <- function(validation_config,
                           model_coefficients,
                           output_dir = NULL,
                           save_log = TRUE) {

  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("MOSS EXTERNAL VALIDATION - ENHANCED GENE MAPPING\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  cohort_name <- validation_config$cohort_name

  if (is.null(output_dir)) {
    output_dir <- file.path("validation_results", cohort_name)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Set up logging if requested
  if (save_log) {
    log_file <- file.path(output_dir, paste0(cohort_name, "_validation_log.txt"))
    sink(log_file, split = TRUE)  # split = TRUE shows output in console AND saves to file
    cat("=================================================================================\n")
    cat("MOSS EXTERNAL VALIDATION LOG\n")
    cat("=================================================================================\n")
    cat(sprintf("Cohort: %s\n", cohort_name))
    cat(sprintf("Date: %s\n", Sys.time()))
    cat(sprintf("Log file: %s\n", log_file))
    cat("=================================================================================\n\n")
  }

  # Wrap everything in tryCatch to ensure sink is closed even if there's an error
  tryCatch({

    cat(sprintf("Cohort: %s\n", cohort_name))
    cat(sprintf("Output: %s\n\n", output_dir))

    # Track original model size
    original_n_genes <- sum(model_coefficients$gene_symbol != "(Intercept)")
    cat(sprintf("Original model: %d genes\n\n", original_n_genes))

    # Load validation data
    validation_data <- load_validation_data(validation_config)
    counts <- validation_data$counts
    metadata <- validation_data$metadata
    val_id_type <- validation_data$gene_id_type

    # Standardize gene IDs
    counts <- standardize_gene_ids(counts, val_id_type)

    # Enhanced gene mapping
    mapping_results <- map_model_genes_enhanced(
      model_coefficients,
      val_id_type,
      rownames(counts)
    )

    # Save mapping results
    write.csv(mapping_results,
              file.path(output_dir, paste0(cohort_name, "_gene_mapping.csv")),
              row.names = FALSE)

    # Create validation coefficients with mapped IDs
    validation_coefs <- create_validation_coefficients(mapping_results, model_coefficients)

    # Get mapped gene IDs for QC and preprocessing
    mapped_gene_ids <- mapping_results %>%
      filter(in_count_matrix) %>%
      pull(mapped_id)

    # QC
    qc_results <- perform_qc(counts, mapped_gene_ids, cohort_name)
    plot_qc(counts, qc_results, cohort_name,
            file.path(output_dir, paste0(cohort_name, "_qc.png")))

    # Preprocess
    cat("\nPreprocessing validation cohort...\n")
    cat("  - Threshold filtering (thresh=2, min 2 samples)\n")
    cat("  - VST normalization\n")
    cat("  - Top 4000 variable genes + signature genes\n")
    cat("  - Centering and scaling\n")

    preprocessed <- preprocess_rnaseq_mofa_style(
      counts,
      top_n = 4000,
      scale_views = TRUE,
      keep_genes = mapped_gene_ids
    )

    cat(sprintf("  ✓ Preprocessed: %d genes\n", nrow(preprocessed)))

    # Check shared genes after preprocessing
    shared_genes <- intersect(mapped_gene_ids, rownames(preprocessed))
    cat(sprintf("\nShared genes after preprocessing: %d/%d\n",
                length(shared_genes), original_n_genes))

    # Incomplete Input (Not all 14 genes are present)
    if (length(shared_genes) < original_n_genes) {
      missing_genes <- mapping_results$gene_symbol[!mapping_results$in_count_matrix]
      dropped_in_preprocess <- setdiff(
        mapping_results$gene_symbol[mapping_results$in_count_matrix],
        mapping_results$gene_symbol[mapping_results$mapped_id %in% shared_genes]
      )
      all_missing <- unique(c(missing_genes, dropped_in_preprocess))
    
      stop(
        "\n\n",
        paste(rep("=", 80), collapse = ""), "\n",
        "PARTIAL SIGNATURE — CANNOT PROCEED\n",
        paste(rep("=", 80), collapse = ""), "\n\n",
        sprintf("  %d of %d signature genes were found after preprocessing.\n",
                length(shared_genes), original_n_genes),
        sprintf("  Missing genes: %s\n\n", paste(all_missing, collapse = ", ")),
        "  Ensure the validation dataset has all 14 genes present:\n",
        "  MARCHF2, LINC00861, RALB, LRP10, IL17RA, GMPR, GLT1D1,\n",
        "  RAB32, MCTP1, NFIL3, TCF7, ALOX5, MXD1, FCMR.\n",
        "  If you are certain all genes are present in your dataset,\n",
        "  check your gene ID formats.\n",
        paste(rep("=", 80), collapse = ""), "\n",
        call. = FALSE
      )
    } else {
      cat("\n✓ Full signature available - using original coefficients\n")
      final_model_coefs <- validation_coefs
    }

    # Predict MoSS scores
    cat("\nPredicting MoSS scores...\n")
    moss_scores <- predict_moss_scores(preprocessed, final_model_coefs)

    # Plot and save results
    results_df <- plot_validation_results(
      moss_scores,
      metadata,
      cohort_name,
      file.path(output_dir, paste0(cohort_name, "_moss_distribution.png"))
    )

    # Add gene symbols to results
    results_df_detailed <- results_df
    if ("gene_symbol" %in% colnames(final_model_coefs)) {
      results_df_detailed$genes_used <- paste(
        final_model_coefs$gene_symbol[final_model_coefs$feature != "(Intercept)"],
        collapse = ", "
      )
    }

    # Save outputs
    write.csv(results_df_detailed,
              file.path(output_dir, paste0(cohort_name, "_moss_predictions.csv")),
              row.names = FALSE)

    write.csv(final_model_coefs,
              file.path(output_dir, paste0(cohort_name, "_final_coefficients.csv")),
              row.names = FALSE)

    # Summary
    validation_summary <- data.frame(
      metric = c("Cohort", "Original signature size", "Genes successfully mapped",
                 "Genes after preprocessing", "Final model size",
                 "Genes used", "Mean MoSS", "SD MoSS", "Median MoSS", "IQR MoSS"),
      value = c(cohort_name, original_n_genes,
                sum(mapping_results$in_count_matrix),
                length(shared_genes),
                sum(final_model_coefs$feature != "(Intercept)"),
                paste(final_model_coefs$gene_symbol[final_model_coefs$feature != "(Intercept)"],
                      collapse = ", "),
                mean(moss_scores), sd(moss_scores),
                median(moss_scores), IQR(moss_scores))
    )
    write.csv(validation_summary,
              file.path(output_dir, paste0(cohort_name, "_summary.csv")),
              row.names = FALSE)

    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    cat("VALIDATION COMPLETE!\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
    cat("Output files:\n")
    cat(sprintf("  1. %s\n", paste0(cohort_name, "_gene_mapping.csv")))
    cat(sprintf("  2. %s\n", paste0(cohort_name, "_qc.png")))
    cat(sprintf("  3. %s\n", paste0(cohort_name, "_moss_distribution.png")))
    cat(sprintf("  4. %s\n", paste0(cohort_name, "_moss_predictions.csv")))
    cat(sprintf("  5. %s\n", paste0(cohort_name, "_final_coefficients.csv")))
    cat(sprintf("  6. %s\n", paste0(cohort_name, "_summary.csv")))
    if (save_log) {
      cat(sprintf("  7. %s\n", paste0(cohort_name, "_validation_log.txt")))
    }
    cat("\n")

    if (save_log) {
      cat(sprintf("✓ Complete validation log saved to: %s\n\n", log_file))
    }

    result <- list(
      cohort = cohort_name,
      counts = counts,
      preprocessed = preprocessed,
      mapping_results = mapping_results,
      moss_scores = moss_scores,
      results_df = results_df_detailed,
      final_model_coefs = final_model_coefs,
      qc = qc_results,
      log_file = if(save_log) log_file else NULL
    )

    return(result)

  }, finally = {
    # Always close the sink connection, even if there's an error
    if (save_log && sink.number() > 0) {
      sink()
    }
  })
}
