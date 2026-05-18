# =============================================================================
# Script:  validation_lasso.R
# Purpose: Cross-cohort LASSO validation of MOFA factors. Trains LASSO models
#          in the discovery cohort (RESERVE-U-2-TOR, Tororo) and evaluates
#          transferability to the validation cohort (RESERVE-U-1-EBB, Entebbe)
#          using RNA-seq and Olink proteomics separately. Extracts non-zero
#          LASSO features at lambda.min and computes Pearson correlations
#          between predicted and de novo MOFA Factor 1 scores.
# =============================================================================
# Inputs (paths relative to repository root; update via here()):
#   - data/cohort1_tororo_all_274_raw_ensembl_only.csv   [controlled access via dbGaP]
#       Raw RNA-seq counts, discovery cohort.
#   - data/cohort2_entebbe_all_128_raw_ensembl.csv
#       Raw RNA-seq counts, validation cohort.
#   - data/cohort1_tororo_all_287_olink_proteomics.csv
#       Olink proteomics, discovery cohort.
#   - data/cohort2_entebbe_all_278_olink_proteomics.csv
#       Olink proteomics, validation cohort.
#   - results/MOFA_runs/<run_folder>/model.hdf5
#       Trained MOFA model from the discovery cohort (produced by MOFA_analysis.R).
#   - results/MOFA_runs/Entebbe_woGSVA_denovo/model.hdf5
#       De novo MOFA model trained on the validation cohort (for correlation).
#
# Functions defined:
#
#   preprocess_rnaseq_mofa_style(raw_mat, top_n, scale_views)
#     Applies VST normalization, removes zero-variance genes, selects the top
#     `top_n` (default 4,000) most variable genes, and optionally z-scores.
#     Mirrors the preprocessing used in MOFA_analysis.R.
#
#   preprocess_olinc_mofa_style(raw_mat, top_n, scale_views, olink_remove)
#     Removes flagged Olink proteins, drops zero-variance features, and
#     optionally z-scores the proteomics matrix.
#
#   run_internal_cv_lasso(X_full, Y_full, factor_name, outcome_type,
#                         n_repeats, train_frac, seed_range)
#     Runs `n_repeats` (default=20) repeated 70/30 splits of the discovery cohort and
#     evaluates LASSO prediction of a single MOFA factor. Returns mean, min,
#     max, and SD of R-squared (continuous) or accuracy (binary) across splits.
#     Seeds are logged to a .log file for reproducibility.
#
#   predict_on_external(X_train, Y_train, X_test, factor_name, outcome_type,
#                       n_repeats, seed_range)
#     Trains `n_repeats` (default=20) LASSO models on the full discovery cohort and
#     applies each to the validation cohort. Selects the best model by minimum
#     internal CV loss at lambda.min and returns predictions and model object.
#
#   extract_lambda_min_nonzero(model, type_label)
#     Extracts non-zero LASSO coefficients at lambda.min from a cv.glmnet
#     object. Returns a data frame with columns: name (feature), coefficient.
#
# Key steps:
#   1. Preprocess RNA-seq and Olink data for both cohorts (MOFA-style).
#   2. Extract MOFA factor scores from the discovery cohort model as training
#      labels (Y).
#   3. Run repeated internal CV LASSO (discovery cohort only) and external
#      LASSO prediction (discovery → validation) for Factors 1 and 2,
#      separately for RNA and protein views.
#   4. Extract non-zero features at lambda.min for each view × factor
#      combination and write to CSV.
#   5. Compute Pearson correlations between predicted and de novo MOFA Factor 1
#      scores in the validation cohort (RNA and protein views) and write to CSV.
#
# Outputs (written to results/validation_runs/):
#   - intVal_Factor<n>_<view>_woGSEA.rds         Internal CV results (per factor/view)
#   - extVal_Factor<n>_<view>_woGSEA.rds         External validation results
#   - lasso_cv_seeds_<factor>.log                Seeds used per internal CV run
#   - external_lasso_seeds_<factor>.log          Seeds used per external prediction run
#   - non0_rna_Factor<n>_lambdaMin_woGSEA.csv    Non-zero RNA features at lambda.min
#   - non0_protein_Factor<n>_lambdaMin_woGSEA.csv  Non-zero protein features at lambda.min
#   - correlation_entebbe.csv                    Pearson r: predicted vs. de novo Factor 1
#
# Dependencies: tidyverse, glmnet, MOFA2, here
#               Sources src/methods.R
#
# =============================================================================
# # Example usage:
#   library(here)
#   source(here("src", "methods.R"))
#
#   # --- Preprocess discovery RNA-seq ---
#   X_train_rna <- prepRNAseq(here("data", "cohort1_tororo_all_274_raw_ensembl_only.csv"),
#                              verbose = TRUE)
#   X_train_rna <- preprocess_rnaseq_mofa_style(X_train_rna)
#
#   # --- Preprocess validation RNA-seq ---
#   X_test_rna <- prepRNAseq(here("data", "cohort2_entebbe_all_128_raw_ensembl.csv"),
#                             verbose = TRUE)
#   X_test_rna <- preprocess_rnaseq_mofa_style(X_test_rna)
#
#   # --- Load MOFA factor scores from discovery model as training labels ---
#   model <- load_model(here("results", "MOFA_runs", "<run_folder>", "model.hdf5"))
#   Y_train_rna <- get_factors(model)[[1]]
#
#   # --- Internal CV and external prediction (repeat for Factor2 and protein) ---
#   f1_rna_int <- run_internal_cv_lasso(t(X_train_rna), Y_train_rna, "Factor1")
#   f1_rna_ext <- predict_on_external(t(X_train_rna), Y_train_rna,
#                                     t(X_test_rna),  "Factor1")
#   saveRDS(f1_rna_int,
#           here("results", "validation_runs", "intVal_Factor1_rna_woGSEA.rds"))
#   saveRDS(f1_rna_ext,
#           here("results", "validation_runs", "extVal_Factor1_rna_woGSEA.rds"))
#
# =============================================================================

library(tidyverse)
library(glmnet)
library(MOFA2)
library(here)
source(here("src", "methods.R"))

# Note: MOFA factor matrix should be of shape (samples x factors)
#       i.e., rownames(Y) must match rownames(X), and columns are Factor1, Factor2, etc.

# ------------------------------
# Function: preprocess_rnaseq_mofa_style
# ------------------------------
preprocess_rnaseq_mofa_style <- function(raw_mat, top_n = 4000, scale_views = TRUE) {
  vst <- normalizeRNASEQwithVST(raw_mat)
  vst <- removeNoVarRows(vst)
  vst <- vst[names(sort(apply(vst, 1, var), decreasing = TRUE))[1:top_n], ]
  if (scale_views) {
    vst <- t(scale(t(vst)))
  }
  return(vst)
}

preprocess_olinc_mofa_style <- function(raw_mat, top_n = 4000, scale_views = TRUE, olink_remove = NULL) {
  if (!is.null(olink_remove)) {
    raw_mat[setdiff(rownames(raw_mat), olinc_remove), ]
  }
  olinc_d <- removeNoVarRows(raw_mat)
  if (scale_views) {
    olinc_d <- t(scale(t(olinc_d)))
  }
  return(olinc_d)
}

# ------------------------------
# Function: run_internal_cv_lasso
# ------------------------------
run_internal_cv_lasso <- function(X_full, Y_full, factor_name = "Factor1", outcome_type = "continuous", n_repeats = 20, train_frac = 0.7, seed_range = c(100, 900)) {
  seeds <- sample(seed_range[1]:seed_range[2], n_repeats, replace = FALSE)
  writeLines(as.character(seeds), con = paste0("lasso_cv_seeds_", factor_name, ".log"))
  
  perf_vec <- numeric(n_repeats)
  family_type <- ifelse(outcome_type == "binary", "binomial", "gaussian")
  
  for (i in 1:n_repeats) {
    set.seed(seeds[i])
    sample_ids <- rownames(X_full)
    train_ids <- sample(sample_ids, size = floor(train_frac * length(sample_ids)))
    test_ids  <- setdiff(sample_ids, train_ids)
    
    X_train <- X_full[train_ids, , drop = FALSE]
    Y_train <- Y_full[train_ids, factor_name]
    X_test  <- X_full[test_ids, , drop = FALSE]
    Y_test  <- Y_full[test_ids, factor_name]
    
    nzv <- apply(X_train, 2, var, na.rm = TRUE) > 0
    X_train <- X_train[, nzv]
    X_test  <- X_test[, colnames(X_train)]
    
    cv_fit <- cv.glmnet(x = X_train, y = Y_train, alpha = 1, family = family_type)
    pred_test <- predict(cv_fit, newx = X_test, s = "lambda.min", type = ifelse(outcome_type == "binary", "response", "link"))
    
    if (outcome_type == "continuous") {
      ss_total <- sum((Y_test - mean(Y_test))^2)
      ss_resid <- sum((Y_test - pred_test)^2)
      perf_vec[i] <- 1 - (ss_resid / ss_total)
    } else {
      pred_bin <- as.numeric(pred_test > 0.5)
      perf_vec[i] <- mean(pred_bin == Y_test)
    }
  }
  
  cat("\nPerformance summary over", n_repeats, "runs (", outcome_type, "):")
  cat("\n  Mean:", round(mean(perf_vec), 4))
  cat("\n  Min:", round(min(perf_vec), 4))
  cat("\n  Max:", round(max(perf_vec), 4))
  cat("\n  SD:", round(sd(perf_vec), 4), "\n")
  
  return(list(mean = mean(perf_vec), min = min(perf_vec), max = max(perf_vec), sd = sd(perf_vec), all = perf_vec))
}

# ------------------------------
# Function: predict_on_external
# ------------------------------
predict_on_external <- function(X_train, Y_train, X_test, factor_name = "Factor1", outcome_type = "continuous", n_repeats = 20, seed_range = c(100, 900)) {
  seeds <- sample(seed_range[1]:seed_range[2], n_repeats, replace = FALSE)
  writeLines(as.character(seeds), con = paste0("external_lasso_seeds_", factor_name, ".log"))
  
  common_genes <- Reduce(intersect, list(colnames(X_train), colnames(X_test)))
  X_train <- X_train[, common_genes]
  X_test  <- X_test[, common_genes]
  
  # Ensure X and Y training datasets have matching sample_ids
  sample_ids <- rownames(X_train)
  Y_train <- Y_train[sample_ids, ]
  nzv <- apply(X_train, 2, var, na.rm = TRUE) > 0
  X_train <- X_train[, nzv]
  X_test  <- X_test[, colnames(X_train)]
  
  models <- list()
  losses <- numeric(n_repeats)
  predictions <- list()
  
  family_type <- ifelse(outcome_type == "binary", "binomial", "gaussian")
  for (i in 1:n_repeats) {
    set.seed(seeds[i])
    cv_fit <- cv.glmnet(x = X_train, y = Y_train[, factor_name], alpha = 1, family = family_type)
    pred <- predict(cv_fit, newx = X_test, s = "lambda.min", type = ifelse(outcome_type == "binary", "response", "link"))
    
    models[[i]] <- cv_fit
    
    # Get CV error at lambda.min
    lambda_min_idx <- which(cv_fit$lambda == cv_fit$lambda.min)
    losses[i] <- cv_fit$cvm[lambda_min_idx]
    
    predictions[[i]] <- pred
  }
  
  best_run <- which.min(losses)
  cat("\nExternal prediction (best internal CV loss at lambda.min):")
  cat("\nBest run:", best_run)
  cat("\nCV loss at lambda.min:", round(losses[best_run], 4), "\n")
  
  return(list(
    model = models[[best_run]],
    pred = predictions[[best_run]],
    loss = losses[best_run],
    all_loss = losses,
    all_preds = predictions
  ))
}

# ------------------------------
# Function: extract_lambda_min_nonzero
# ------------------------------
# Extracts non-zero coefficients at lambda.min from a cv.glmnet object.
# Returns a data frame with columns: name (feature), coefficient.
# ------------------------------
extract_lambda_min_nonzero <- function(model, type_label = "gene") {
  coeff_matrix <- as.matrix(model$glmnet.fit$beta)
  lambda_min <- model$lambda.min
  lambda_idx <- which.min(abs(model$glmnet.fit$lambda - lambda_min))
  lambda_name <- colnames(coeff_matrix)[lambda_idx]

  cat(sprintf("\n%s - Lambda min index: %s\n", type_label, lambda_name))

  coeffs <- coeff_matrix[, lambda_idx]
  nonzero <- coeffs[coeffs != 0]
  df <- data.frame(
    name = names(nonzero),
    coefficient = nonzero
  )
  return(df)
}

# ------------------------------
# Load external validation results
# ------------------------------
f1_rna_ext     <- readRDS(here("results", "validation_runs", "extVal_Factor1_rna_woGSEA.rds"))
f2_rna_ext     <- readRDS(here("results", "validation_runs", "extVal_Factor2_rna_woGSEA.rds"))
f1_protein_ext <- readRDS(here("results", "validation_runs", "extVal_Factor1_protein_woGSEA.rds"))
f2_protein_ext <- readRDS(here("results", "validation_runs", "extVal_Factor2_protein_woGSEA.rds"))

# ------------------------------
# Extract non-zero features at lambda.min
# ------------------------------

# Factor 1 - RNA
rna1_min_df <- extract_lambda_min_nonzero(f1_rna_ext$model, type_label = "Factor1 RNA")
colnames(rna1_min_df)[1] <- "gene"
write.csv(rna1_min_df,
          here("results", "validation_runs", "non0_rna_Factor1_lambdaMin_woGSEA.csv"),
          row.names = FALSE)

# Factor 1 - Protein
protein1_min_df <- extract_lambda_min_nonzero(f1_protein_ext$model, type_label = "Factor1 Protein")
colnames(protein1_min_df)[1] <- "protein"
write.csv(protein1_min_df,
          here("results", "validation_runs", "non0_protein_Factor1_lambdaMin_woGSEA.csv"),
          row.names = FALSE)

# Factor 2 - RNA
rna2_min_df <- extract_lambda_min_nonzero(f2_rna_ext$model, type_label = "Factor2 RNA")
colnames(rna2_min_df)[1] <- "gene"
write.csv(rna2_min_df,
          here("results", "validation_runs", "non0_rna_Factor2_lambdaMin_woGSEA.csv"),
          row.names = FALSE)

# Factor 2 - Protein
protein2_min_df <- extract_lambda_min_nonzero(f2_protein_ext$model, type_label = "Factor2 Protein")
colnames(protein2_min_df)[1] <- "protein"
write.csv(protein2_min_df,
          here("results", "validation_runs", "non0_protein_Factor2_lambdaMin_woGSEA.csv"),
          row.names = FALSE)

# ------------------------------
# Pearson correlation: predicted vs. de novo Factor 1 (Entebbe)
# ------------------------------
denovo_model <- load_model(here("results", "MOFA_runs", "Entebbe_woGSVA_denovo", "model.hdf5"))
fac_mat      <- get_factors(denovo_model)[[1]]

common_rna     <- intersect(rownames(fac_mat), rownames(f1_rna_ext$pred))
common_protein <- intersect(rownames(fac_mat), rownames(f1_protein_ext$pred))

cor_matrix <- data.frame(
  view    = c("RNA", "Protein"),
  pearson_r = c(
    cor(fac_mat[common_rna,     "Factor1"], f1_rna_ext$pred[common_rna, ]),
    cor(fac_mat[common_protein, "Factor1"], f1_protein_ext$pred[common_protein, ])
  )
)

write.csv(cor_matrix,
          here("results", "validation_runs", "correlation_entebbe.csv"),
          row.names = FALSE)
