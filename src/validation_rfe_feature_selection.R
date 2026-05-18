# =============================================================================
# Script:  validation_rfe_feature_selection.R
# Purpose: Recursive feature elimination (RFE) with linear models to identify
#          compact gene and protein signatures that predict MOFA Factor 1 in
#          the discovery cohort (RESERVE-U-2-TOR, Tororo) and validate them
#          externally in the validation cohort (RESERVE-U-1-EBB, Entebbe).
#          Merges the logic of the original rfe_gene_selection.R and
#          rfe_protein_selection.R into a single parameterized pipeline.
# =============================================================================
# Script structure:
#   1. Preprocessing functions and run_rfe_pipeline() are defined.
#   2. Discovery and validation cohort data are loaded and preprocessed
#      (RNA-seq and Olink, both cohorts; MOFA factor scores for both cohorts).
#   3. run_rfe_pipeline() is called once for the gene (RNA) view, then once
#      for the protein (Olink) view. Both calls follow the same nine-step
#      internal workflow (see below).
#
# Inputs:
#   - data/cohort1_tororo_all_274_raw_ensembl.csv   [controlled access via dbGaP]
#   - data/cohort2_entebbe_all_128_raw_ensembl.csv
#       Raw RNA-seq count matrices for discovery and validation cohorts.
#   - data/cohort1_tororo_all_287_olink_proteomics.csv
#   - data/cohort2_entebbe_all_278_olink_proteomics.csv
#       Olink proteomics matrices for both cohorts.
#   - results/MOFA_runs/Tororo_woGSVA_selected/model.hdf5
#       Trained MOFA model (discovery cohort); provides Factor 1 training labels.
#   - results/MOFA_runs/Entebbe_woGSVA_denovo/model.hdf5
#       De novo MOFA model (validation cohort); provides Factor 1 test labels.
#   - results/validation_runs/non0_rna_Factor1_lambdaMin_woGSEA.csv
#   - results/validation_runs/non0_protein_Factor1_lambdaMin_woGSEA.csv
#       LASSO-selected features at lambda.min (produced by
#       validation_lasso_postprocess.R).
#
# Functions defined:
#
#   preprocess_rnaseq_mofa_style(raw_mat, top_n, scale_views)
#     VST normalization, zero-variance removal, top-N gene selection, optional
#     z-scoring. Mirrors the preprocessing in MOFA_analysis.R.
#
#   preprocess_olinc_mofa_style(raw_mat, top_n, scale_views, olink_remove)
#     Removes flagged Olink proteins, zero-variance removal, optional z-scoring.
#
#   run_rfe_pipeline(X_train_sf, X_test_sf, Y_train, Y_test, lasso_file,
#                    feature_col, feature_type, output_dir, factor_name,
#                    optimal_size, signature_sizes)
#     General-purpose RFE pipeline parameterized by feature type ("gene" or
#     "protein"). Accepts data already in samples x features format.
#     NOTE: RNA matrices must be transposed (t()) before passing; protein
#     matrices are already in samples x features format.
#     Internal steps:
#       1. Filter features to the LASSO-selected set available in both cohorts.
#       2. Rank features by importance using RFE with lmFuncs (10-fold CV,
#          5 repeats, sizes 5–25).
#       3. Train a linear model for each signature size and record RMSE and R².
#       4. Assess feature concordance across sizes: Jaccard similarity matrix,
#          core features (sizes 5–10), per-size composition, and stability
#          frequency across all tested sizes.
#       5. Select signature size (post-hoc decision: size 15 by default).
#       6. Retrain final linear model on the full discovery cohort at the
#          selected size.
#       7. Apply final model to the validation cohort; compute Pearson r, R²,
#          RMSE, and MAE.
#       8. Save all results to output_dir: feature rankings, RMSE-by-size
#          table, final signature, model coefficients, predictions, validation
#          summary, and complete RDS object.
#       9. Generate and save RMSE-by-size curve and external validation scatter.
#     Returns complete_results list invisibly.
#
# Outputs (written to results/validation_RFE/gene/ and .../protein/):
#   - rfe_<type>_rankings_factor1.csv     Feature importance rankings from RFE
#   - performance_by_size_factor1.csv     RMSE and R² for each signature size
#   - final_signature_factor1.csv         Selected signature features and ranks
#   - final_model_coefficients_factor1.csv LM coefficients
#   - external_predictions_factor1.csv    True vs. predicted Factor 1 (Entebbe)
#   - validation_summary_factor1.csv      Summary of train/external metrics
#   - complete_results_factor1.rds        Full results list
#   - rmse_by_size_factor1.png            RMSE curve with optimal size marked
#   - external_validation_scatter_factor1.png Predicted vs. true scatter
#
# Additional outputs (gene pipeline only):
#   - results/validation_RFE/gene/tororo_preprocessed_rna.rds
#   - results/validation_RFE/gene/tororo_factor1_scores.rds
#
# Dependencies: tidyverse, MOFA2, caret, here
#               Sources src/methods.R
# =============================================================================

library(tidyverse)
library(MOFA2)
library(caret)
library(here)

source(here("src", "methods.R"))

# =============================================================================
# Preprocessing functions
# =============================================================================

preprocess_rnaseq_mofa_style <- function(raw_mat, top_n = 4000, scale_views = TRUE) {
  vst <- normalizeRNASEQwithVST(raw_mat)
  vst <- removeNoVarRows(vst)
  vst <- vst[names(sort(apply(vst, 1, var), decreasing = TRUE))[1:top_n], ]
  if (scale_views) vst <- t(scale(t(vst)))
  return(vst)
}

preprocess_olinc_mofa_style <- function(raw_mat, top_n = 4000, scale_views = TRUE, olink_remove = NULL) {
  if (!is.null(olink_remove))
    raw_mat <- raw_mat[setdiff(rownames(raw_mat), olink_remove), ]
  olinc_d <- removeNoVarRows(raw_mat)
  if (scale_views) olinc_d <- t(scale(t(olinc_d)))
  return(olinc_d)
}

# =============================================================================
# Helper
# =============================================================================

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

# =============================================================================
# RFE pipeline (parameterized for gene or protein)
# X_train_sf / X_test_sf must be in samples x features format before calling.
# For RNA, pass t(X_train_rna); for protein, pass X_train_protein directly.
# =============================================================================

run_rfe_pipeline <- function(X_train_sf,
                              X_test_sf,
                              Y_train,
                              Y_test,
                              lasso_file,
                              feature_col,
                              feature_type,
                              output_dir,
                              factor_name     = "Factor1",
                              optimal_size    = 15,
                              signature_sizes = seq(5, 25, by = 1)) {

  ft <- feature_type
  FT <- paste0(toupper(substr(ft, 1, 1)), substr(ft, 2, nchar(ft)))

  # ------------------------------------------------------------------
  # STEP 1: Filter to LASSO-selected features
  # ------------------------------------------------------------------
  cat(sprintf("\n╔══ STEP 1: Filtering to LASSO-Selected %ss ══╗\n\n", FT))

  lasso_df        <- read.csv(lasso_file)
  lasso_features  <- as.character(lasso_df[[feature_col]])
  cat(sprintf("LASSO selected %ss (lambda.min): %d\n\n", ft, length(lasso_features)))

  avail_train     <- intersect(lasso_features, colnames(X_train_sf))
  avail_test      <- intersect(lasso_features, colnames(X_test_sf))
  common_features <- intersect(avail_train, avail_test)
  cat(sprintf("%ss available in Tororo: %d\n", FT, length(avail_train)))
  cat(sprintf("%ss available in Entebbe: %d\n", FT, length(avail_test)))
  cat(sprintf("%ss available in both cohorts: %d\n\n", FT, length(common_features)))

  X_train_sub <- X_train_sf[, common_features, drop = FALSE]
  X_test_sub  <- X_test_sf[, common_features, drop = FALSE]

  common_tr   <- intersect(rownames(X_train_sub), rownames(Y_train))
  X_train_sub <- X_train_sub[common_tr, , drop = FALSE]
  y_train     <- Y_train[common_tr, factor_name]

  common_te   <- intersect(rownames(X_test_sub), rownames(Y_test))
  X_test_sub  <- X_test_sub[common_te, , drop = FALSE]
  y_test      <- Y_test[common_te, factor_name]

  cat(sprintf("Final dimensions:\n  Training: %d samples x %d %ss\n  Test: %d samples x %d %ss\n\n",
              nrow(X_train_sub), ncol(X_train_sub), ft,
              nrow(X_test_sub),  ncol(X_test_sub),  ft))

  colnames(X_train_sub) <- make.names(colnames(X_train_sub))
  colnames(X_test_sub)  <- make.names(colnames(X_test_sub))
  X_train_df <- as.data.frame(X_train_sub)
  X_test_df  <- as.data.frame(X_test_sub)

  # ------------------------------------------------------------------
  # STEP 2: RFE ranking
  # ------------------------------------------------------------------
  cat(sprintf("╔══ STEP 2: Running RFE to Rank %ss ══╗\n\n", FT))

  set.seed(123)
  cv_folds <- createMultiFolds(y_train, k = 10, times = 5)
  cat("Created", length(cv_folds), "folds (10-fold CV, 5 repeats)\n\n")

  rfe_ctrl <- rfeControl(functions = lmFuncs, method = "cv",
                         index = cv_folds, verbose = FALSE)

  cat(sprintf("Running RFE for %s ranking (this may take several minutes)...\n\n", ft))
  rfe_result <- rfe(x = X_train_df, y = y_train,
                    sizes = signature_sizes, rfeControl = rfe_ctrl, metric = "RMSE")

  feature_rankings <- rfe_result$variables %>%
    group_by(var) %>%
    summarise(mean_importance = mean(Overall), .groups = "drop") %>%
    arrange(desc(mean_importance))

  cat(sprintf("RFE ranking complete! Top 10 %ss:\n", ft))
  print(head(feature_rankings, 10))

  # ------------------------------------------------------------------
  # STEP 3: Train LM for each signature size
  # ------------------------------------------------------------------
  cat(sprintf("\n╔══ STEP 3: Training LM for Each Signature Size ══╗\n\n"))
  cat("Testing signature sizes:", paste(signature_sizes, collapse = ", "), "\n\n")

  results_by_size <- list()
  for (size in signature_sizes) {
    selected   <- feature_rankings$var[1:size]
    train_data <- data.frame(y = y_train, X_train_df[, selected, drop = FALSE])
    lm_model   <- lm(y ~ ., data = train_data)
    train_pred <- stats::predict(lm_model)
    train_rmse <- sqrt(mean((y_train - train_pred)^2))
    train_r2   <- cor(y_train, train_pred)^2
    results_by_size[[as.character(size)]] <- list(
      size = size, features = selected, model = lm_model,
      train_rmse = train_rmse, train_r2 = train_r2
    )
    cat(sprintf("Size %2d: RMSE = %.4f, R² = %.4f\n", size, train_rmse, train_r2))
  }

  # ------------------------------------------------------------------
  # STEP 4: Concordance analysis across signature sizes
  # ------------------------------------------------------------------
  cat(sprintf("\n╔══ STEP 4: Assessing %s Concordance Across Signature Sizes ══╗\n\n", FT))

  feature_lists <- lapply(signature_sizes, function(s)
    gsub("^X", "", results_by_size[[as.character(s)]]$features))
  names(feature_lists) <- paste0("size_", signature_sizes)

  sim_mat <- outer(seq_along(signature_sizes), seq_along(signature_sizes),
                   Vectorize(function(i, j) jaccard(feature_lists[[i]], feature_lists[[j]])))
  rownames(sim_mat) <- colnames(sim_mat) <- paste0("Size", signature_sizes)
  cat(sprintf("Jaccard Similarity Matrix (%s overlap between sizes):\n\n", ft))
  print(round(sim_mat, 3))

  cat(sprintf("\n\n%s composition by signature size:\n%s\n\n", FT, strrep("=", 60)))
  for (size in signature_sizes) {
    feats <- gsub("^X", "", results_by_size[[as.character(size)]]$features)
    cat(sprintf("Size %2d (RMSE=%.4f, R²=%.4f):\n  %ss: %s\n\n",
                size, results_by_size[[as.character(size)]]$train_rmse,
                results_by_size[[as.character(size)]]$train_r2,
                FT, paste(feats, collapse = ", ")))
  }

  core_features <- Reduce(intersect, feature_lists[1:6])
  cat(sprintf("\nCore %ss (present in sizes 5–10): %s\n  Count: %d\n\n",
              ft, paste(core_features, collapse = ", "), length(core_features)))

  cat(sprintf("\n%ss added at each size:\n", FT))
  for (i in 2:length(signature_sizes)) {
    new_f <- setdiff(feature_lists[[i]], feature_lists[[i - 1]])
    if (length(new_f) > 0)
      cat(sprintf("  Size %d → %d: %s\n",
                  signature_sizes[i - 1], signature_sizes[i], paste(new_f, collapse = ", ")))
  }

  feat_freq  <- sort(table(unlist(feature_lists)), decreasing = TRUE)
  universal  <- names(feat_freq[feat_freq == length(signature_sizes)])
  cat(sprintf("\nTop 10 most frequently selected %ss:\n", ft))
  print(head(feat_freq, 10))
  cat(sprintf("\n%ss appearing in ALL tested sizes (5–25):\n", FT))
  if (length(universal) > 0)
    cat(sprintf("  %s\n  Count: %d\n", paste(universal, collapse = ", "), length(universal)))
  else
    cat(sprintf("  None (no %ss appear in all sizes)\n", ft))

  # ------------------------------------------------------------------
  # STEP 5: Select signature size
  # ------------------------------------------------------------------
  cat(sprintf("\n╔══ STEP 5: Selecting Signature Size ══╗\n\n"))

  rmse_vals     <- sapply(signature_sizes, function(s) results_by_size[[as.character(s)]]$train_rmse)
  min_rmse_size <- signature_sizes[which.min(rmse_vals)]
  cat(sprintf("RMSE-optimal size (data-driven): %d %ss (RMSE = %.4f)\n",
              min_rmse_size, ft, min(rmse_vals)))
  cat(sprintf("Selected size (post-hoc decision): %d %ss\n", optimal_size, ft))
  cat(sprintf("Training RMSE at size %d: %.4f\n", optimal_size,
              results_by_size[[as.character(optimal_size)]]$train_rmse))
  cat(sprintf("Training R² at size %d: %.4f\n\n", optimal_size,
              results_by_size[[as.character(optimal_size)]]$train_r2))

  # ------------------------------------------------------------------
  # STEP 6: Retrain final model at selected size
  # ------------------------------------------------------------------
  cat(sprintf("╔══ STEP 6: Training Final Model at Size %d ══╗\n\n", optimal_size))

  final_features       <- results_by_size[[as.character(optimal_size)]]$features
  original_feature_ids <- gsub("^X", "", final_features)

  cat(sprintf("Final %s signature (%d %ss):\n", ft, optimal_size, ft))
  print(original_feature_ids)

  train_data_final   <- data.frame(y = y_train, X_train_df[, final_features, drop = FALSE])
  final_model        <- lm(y ~ ., data = train_data_final)
  final_coefficients <- coef(final_model)
  cat("\nFinal model coefficients:\n"); print(round(final_coefficients, 4))

  # ------------------------------------------------------------------
  # STEP 7: External validation on Entebbe cohort
  # ------------------------------------------------------------------
  cat(sprintf("\n╔══ STEP 7: External Validation on Entebbe Cohort ══╗\n\n"))

  X_test_final     <- X_test_df[, final_features, drop = FALSE]
  test_predictions <- stats::predict(final_model, newdata = as.data.frame(X_test_final))
  test_correlation <- cor(y_test, test_predictions)
  test_r2          <- test_correlation^2
  test_rmse        <- sqrt(mean((y_test - test_predictions)^2))
  test_mae         <- mean(abs(y_test - test_predictions))

  cat(sprintf("  Pearson correlation: %.4f\n  R²: %.4f\n  RMSE: %.4f\n  MAE: %.4f\n\n",
              test_correlation, test_r2, test_rmse, test_mae))

  # ------------------------------------------------------------------
  # STEP 8: Save results
  # ------------------------------------------------------------------
  cat("╔══ STEP 8: Saving Results ══╗\n\n")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  rankings_out <- feature_rankings
  rankings_out[[paste0("original_", ft, "_id")]] <- gsub("^X", "", rankings_out$var)
  write.csv(rankings_out,
            file.path(output_dir, sprintf("rfe_%s_rankings_factor1.csv", ft)),
            row.names = FALSE)

  perf_by_size <- data.frame(
    signature_size = signature_sizes,
    train_rmse = sapply(signature_sizes, function(s) results_by_size[[as.character(s)]]$train_rmse),
    train_r2   = sapply(signature_sizes, function(s) results_by_size[[as.character(s)]]$train_r2)
  )
  write.csv(perf_by_size,
            file.path(output_dir, "performance_by_size_factor1.csv"), row.names = FALSE)

  final_sig_df <- data.frame(r_name = final_features, original_id = original_feature_ids,
                              rank = seq_along(final_features))
  colnames(final_sig_df)[2] <- paste0("original_", ft, "_id")
  write.csv(final_sig_df,
            file.path(output_dir, "final_signature_factor1.csv"), row.names = FALSE)

  write.csv(data.frame(feature = names(final_coefficients),
                        coefficient = as.numeric(final_coefficients)),
            file.path(output_dir, "final_model_coefficients_factor1.csv"), row.names = FALSE)

  write.csv(data.frame(sample_id = rownames(X_test_final),
                        true_factor1 = y_test, predicted_factor1 = test_predictions),
            file.path(output_dir, "external_predictions_factor1.csv"), row.names = FALSE)

  write.csv(data.frame(
    metric = c("Signature Size", "Training RMSE", "Training R²",
               "External Correlation", "External R²", "External RMSE", "External MAE"),
    value  = c(optimal_size,
               results_by_size[[as.character(optimal_size)]]$train_rmse,
               results_by_size[[as.character(optimal_size)]]$train_r2,
               test_correlation, test_r2, test_rmse, test_mae)),
    file.path(output_dir, "validation_summary_factor1.csv"), row.names = FALSE)

  saveRDS(list(
    feature_rankings     = rankings_out,
    results_by_size      = results_by_size,
    optimal_size         = optimal_size,
    final_features       = final_features,
    original_feature_ids = original_feature_ids,
    final_model          = final_model,
    final_coefficients   = final_coefficients,
    external_validation  = list(predictions = test_predictions, true_values = y_test,
                                correlation = test_correlation, r2 = test_r2,
                                rmse = test_rmse, mae = test_mae)
  ), file.path(output_dir, "complete_results_factor1.rds"))

  # ------------------------------------------------------------------
  # STEP 9: Visualizations
  # ------------------------------------------------------------------
  cat("╔══ STEP 9: Creating Visualizations ══╗\n\n")

  p1 <- ggplot(perf_by_size, aes(x = signature_size, y = train_rmse)) +
    geom_line(color = "steelblue", size = 1.2) +
    geom_point(color = "steelblue", size = 3) +
    geom_vline(xintercept = optimal_size, linetype = "dashed", color = "red", size = 1) +
    annotate("text", x = optimal_size, y = max(perf_by_size$train_rmse),
             label = paste("Optimal:", optimal_size, sprintf("%ss", ft)),
             hjust = -0.1, color = "red") +
    labs(title = "Training RMSE by Signature Size (Factor 1)", subtitle = "Tororo cohort",
         x = sprintf("Signature Size (Number of %ss)", FT), y = "Training RMSE") +
    theme_minimal() + theme(plot.title = element_text(face = "bold", size = 14))
  ggsave(file.path(output_dir, "rmse_by_size_factor1.png"), plot = p1,
         width = 10, height = 6, dpi = 300)

  p2 <- ggplot(data.frame(true = y_test, predicted = test_predictions),
               aes(x = true, y = predicted)) +
    geom_point(alpha = 0.6, size = 3, color = "steelblue") +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    annotate("text", x = min(y_test), y = max(test_predictions),
             label = sprintf("r = %.3f\nR² = %.3f\nRMSE = %.3f", test_correlation, test_r2, test_rmse),
             hjust = 0, vjust = 1, size = 5) +
    labs(title = "External Validation: Factor 1",
         subtitle = sprintf("%s signature: %d %ss | Entebbe cohort", FT, optimal_size, ft),
         x = "True Factor 1 Score", y = "Predicted Factor 1 Score") +
    theme_minimal() + theme(plot.title = element_text(face = "bold", size = 14))
  ggsave(file.path(output_dir, "external_validation_scatter_factor1.png"), plot = p2,
         width = 8, height = 8, dpi = 300)

  cat(sprintf("\n═══ FINAL RESULTS: %s ═══\n\n", toupper(ft)))
  cat(sprintf("Optimal Signature: %d %ss\n  IDs: %s\n\n",
              optimal_size, ft, paste(original_feature_ids, collapse = ", ")))
  cat(sprintf("Training Performance (Tororo):\n  RMSE: %.4f\n  R²: %.4f\n\n",
              results_by_size[[as.character(optimal_size)]]$train_rmse,
              results_by_size[[as.character(optimal_size)]]$train_r2))
  cat(sprintf("External Validation (Entebbe):\n  Correlation: %.4f\n  R²: %.4f\n  RMSE: %.4f\n  MAE: %.4f\n\n",
              test_correlation, test_r2, test_rmse, test_mae))
  cat(sprintf("Files saved in: %s\n\n", output_dir))

  invisible(list(
    feature_rankings = rankings_out, results_by_size = results_by_size,
    final_features = final_features, final_model = final_model,
    y_train = y_train
  ))
}

# =============================================================================
# Load and preprocess data
# =============================================================================
cat("\n╔══ Loading and Preprocessing Data ══╗\n\n")

# RNA
cat("Loading Tororo RNA...\n")
X_train_rna <- prepRNAseq(here("data", "cohort1_tororo_all_274_raw_ensembl_only.csv"), verbose = TRUE)
X_train_rna <- preprocess_rnaseq_mofa_style(X_train_rna)

cat("Loading Entebbe RNA...\n")
X_test_rna  <- prepRNAseq(here("data", "cohort2_entebbe_all_128_raw_ensembl.csv"), verbose = TRUE)
X_test_rna  <- preprocess_rnaseq_mofa_style(X_test_rna)

# Protein
olinc_remove <- c("il1alpha", "il2", "il33", "il4", "il13", "prcp", "ltbp2",
                  "sod1", "itgam", "fap", "mfap5")

cat("Loading Tororo Olink proteomics...\n")
X_train_protein <- preprocess_olinc_mofa_style(
  read.csv(here("data", "cohort1_tororo_all_287_olink_proteomics.csv"), row.names = 1, check.names = FALSE),
  olink_remove = olinc_remove)

cat("Loading Entebbe Olink proteomics...\n")
X_test_protein  <- preprocess_olinc_mofa_style(
  read.csv(here("data", "cohort1_entebbe_all_278_olink_proteomics.csv"), row.names = 1, check.names = FALSE),
  olink_remove = olinc_remove)

# MOFA factor scores
cat("Loading MOFA factor scores...\n")
model_train <- load_model(here("results", "MOFA_runs", "Tororo_woGSVA_selected", "model.hdf5"))
Y_train     <- get_factors(model_train)[[1]]

model_test  <- load_model(here("results", "MOFA_runs", "Entebbe_woGSVA_denovo", "model.hdf5"))
Y_test      <- get_factors(model_test)[[1]]

cat("\nData loaded successfully!\n")
cat("  Tororo RNA:      ", nrow(X_train_rna),     "genes x",    ncol(X_train_rna),     "samples\n")
cat("  Entebbe RNA:     ", nrow(X_test_rna),      "genes x",    ncol(X_test_rna),      "samples\n")
cat("  Tororo Protein:  ", nrow(X_train_protein), "proteins x", ncol(X_train_protein), "samples\n")
cat("  Entebbe Protein: ", nrow(X_test_protein),  "proteins x", ncol(X_test_protein),  "samples\n\n")

# =============================================================================
# Run RFE pipeline: Gene (RNA)
# RNA is genes x samples; transpose to samples x genes before passing.
# =============================================================================
gene_results <- run_rfe_pipeline(
  X_train_sf   = t(X_train_rna),
  X_test_sf    = t(X_test_rna),
  Y_train      = Y_train,
  Y_test       = Y_test,
  lasso_file   = here("results", "validation_runs", "non0_rna_Factor1_lambdaMin_woGSEA.csv"),
  feature_col  = "gene",
  feature_type = "gene",
  output_dir   = here("results", "validation_RFE", "gene")
)

# Additional saves for gene pipeline: preprocessed RNA matrix and training factor scores
saveRDS(X_train_rna,
        here("results", "validation_RFE", "gene", "tororo_preprocessed_rna.rds"))
saveRDS(gene_results$y_train,
        here("results", "validation_RFE", "gene", "tororo_factor1_scores.rds"))

# =============================================================================
# Run RFE pipeline: Protein
# Protein matrix is already in samples x features format.
# =============================================================================
run_rfe_pipeline(
  X_train_sf   = X_train_protein,
  X_test_sf    = X_test_protein,
  Y_train      = Y_train,
  Y_test       = Y_test,
  lasso_file   = here("results", "validation_runs", "non0_protein_Factor1_lambdaMin_woGSEA.csv"),
  feature_col  = "protein",
  feature_type = "protein",
  output_dir   = here("results", "validation_RFE", "protein")
)
