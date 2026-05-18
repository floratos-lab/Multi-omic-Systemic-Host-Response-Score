# =============================================================================
# Script:  MOFA_analysis.R
# Purpose: Preprocess multi-omic input data and train a Multi-Omics Factor
#          Analysis (MOFA) model on the discovery cohort (RESERVE-U-2-TOR).
# =============================================================================
# IMPORTANT: cohort1_tororo_all_274_raw_ensembl.csv needs to be accessed via dbGaP (see README)
#
# Inputs:
#   - data/cohort1_tororo_all_298_clinical.csv
#       Clinical metadata for the discovery cohort.
#   - data/cohort1_tororo_all_274_raw_ensembl.csv  [controlled access via dbGaP]
#       Raw RNA-seq counts; variable rnaf_name is commented out by default.
#   - data/cohort1_tororo_all_274_immune_cells_relative_cibersortx.csv
#       CIBERSORTx relative immune cell fraction estimates.
#   - data/cohort1_tororo_all_287_olink_proteomics.csv
#       Olink proteomics panel for the discovery cohort.
#   - results/tororo/GSVA/cohort_tororo_GSVA_reactome.tsv
#   - results/tororo/GSVA/cohort_tororo_GSVA_msigdb_c2.tsv
#   - results/tororo/GSVA/cohort_tororo_GSVA_msigdb_c5.tsv
#       GSVA single-sample enrichment scores (produced by ssEnrichmentScores.R).
#
# Key steps:
#   1. Load and harmonize all data views; remove flagged Olink proteins and
#      low-variance features.
#   2. Apply variance-stabilizing transformation (VST) to RNA-seq counts;
#      select top 4,000 most variable genes.
#   3. Optionally z-score scale each view (scale_views = TRUE by default).
#   4. Pad missing samples across views to satisfy MOFA package requirements.
#   5. Run n_runs independent MOFA training iterations (default: 20), each
#      with a distinct random seed, saving each model as model.hdf5 in a
#      timestamped output folder.
#   6. Select the best model by highest ELBO; export factor score matrix and
#      ELBO comparison table.
#
# Outputs (written to results/MOFA_runs/):
#   - model.hdf5                        Trained MOFA model (per run)
#   - factor_scores.csv                 Sample-level factor loadings (best run)
#   - elbo_scores_<n>_trials.csv        ELBO values across all runs
#   - sessionInfo.txt                   R session information (per run)
#
# Dependencies: MOFA2, here, tidyverse, purrr, tibble
#               Sources src/methods.R
# =============================================================================

library(MOFA2)
library(here)

source(here("src", "methods.R"), chdir=TRUE)
output_dir <- here("results", "MOFA_runs")

# Location of data files -- UPDATE ACCORDINGLY
# NOTE: Raw RNA-sequence data requires controlled access authorization via dbGaP (see README)
# rnaf_name <- here("data", "cohort1_tororo_all_274_raw_ensembl.csv") 
clinf_name <- here("data", "cohort1_tororo_all_298_clinical.csv")
cibersortf_name <- here("data", "cohort1_tororo_all_274_immune_cells_relative_cibersortx.csv")
olincf_name <- here("data", "cohort1_tororo_all_287_olink_proteomics.csv")
# The following datafiles are created by `ssEnrishmentScores.R`.
reactome_fname <- here("results", "tororo/GSVA/cohort_tororo_GSVA_reactome.tsv")
msigc2_fname <- here("results", "tororo/GSVA/cohort_tororo_GSVA_msigdb_c2.tsv")
msigc5_fname <- here("results", "tororo/GSVA/cohort_tororo_GSVA_msigdb_c5.tsv")

# list of proteins to exclude from analysis, per Olink.
olinc_remove = c("il1alpha", "il2", "il33", "il4", "il13", "prcp", "ltbp2", 
		"sod1", "itgam", "fap", "mfap5")

# Indicates if the data should be scaled. According to the MOFA paper,
# scaling is not needed. Scaling the data, however, may be helpful
# for some of the plots; e.g., in heatmaps of gene expression data,
# clusters are more distinguishable visually when using z-scores 
scale_views <- TRUE

# Specify the number of MOFA runs to perform. Each run uses a different
# random seed.
#
# If n_runs > 1, multiple models are trained. The "best" model (i.e.,
# the one with the smallest ELBO value) can then be selected for
# downstream analysis.
n_runs <- 20

# Read in and prep the data matrices
cnt_d <- prepRNAseq(rnaf_name, verbose=TRUE)
clin_d <- read.csv(clinf_name, row.names=1, check.names=FALSE)
ciber_d <- data.frame(t(read.csv(cibersortf_name, row.names=1, check.names=FALSE)), check.names=FALSE)
colnames(ciber_d) <- gsub("_", "-", colnames(ciber_d))
olinc_d <- data.frame(t(read.csv(olincf_name, row.names=1, check.names=FALSE)), check.names=FALSE)
olinc_d <- olinc_d[setdiff(rownames(olinc_d), olinc_remove), , drop = FALSE]
reactome_d <- read.table(reactome_fname, header=TRUE, row.names = 1, check.names=FALSE, quote="", sep="\t") 
msigc2_d <- read.table(msigc2_fname, header=TRUE, row.names = 1, check.names=FALSE, quote="", sep="\t")
msigc5_d <- read.table(msigc5_fname, header=TRUE, row.names = 1, check.names=FALSE, quote="", sep="\t")

# Per MOFA paper/vignette recommendations, perform variance stabilization
# and remove features with low variance.
cnt_vst = normalizeRNASEQwithVST(cnt_d)
cnt_vst = removeNoVarRows(cnt_vst)
ciber_d = removeNoVarRows(ciber_d)
olinc_d = removeNoVarRows(olinc_d)
reactome_d = removeNoVarRows(reactome_d)
msigc2_d = removeNoVarRows(msigc2_d)
msigc5_d = removeNoVarRows(msigc5_d)
n_feat <- 4000
gexp_d <- cnt_vst[names(sort(apply(cnt_vst, 1, var), decreasing=TRUE))[1:n_feat], , drop = FALSE]

# Scale the data, if desired
if (scale_views){
	gexp_d <- t(scale(t(gexp_d)))
	ciber_d <- t(scale(t(ciber_d)))
	olinc_d <- t(scale(t(olinc_d)))
	reactome_d <- t(scale(t(reactome_d)))
	msigc2_d <- t(scale(t(msigc2_d)))
	msigc5_d <- t(scale(t(msigc5_d)))
}

# Prep data matrices by filling in missing columns. This is a requirement of
# the MOFA package.
sample_ids <- Reduce(union, lapply(list(gexp_d, olinc_d, ciber_d, reactome_d,
						msigc2_d, msigc5_d), colnames))
sample_ids <- sample_ids[sample_ids != ""] # Ensuring only valid samples are included (remove NA)
gexp_d <- fillEmptyCols(gexp_d, sample_ids)
ciber_d <- fillEmptyCols(ciber_d, sample_ids)
olinc_d <- fillEmptyCols(olinc_d, sample_ids)
reactome_d <- fillEmptyCols(reactome_d, sample_ids)
msigc2_d <- fillEmptyCols(msigc2_d, sample_ids)
msigc5_d <- fillEmptyCols(msigc5_d, sample_ids)

# row name fix due to format reading errors
# replace with a clean ASCII version
rn_ciber <- rownames(ciber_d)
rn_ciber[rn_ciber == "B cells na\xefve"] <- "B cells naive"
rownames(ciber_d) <- rn_ciber

# Create mofa object, comprising all data matrices (aka, "views")
mdd <- create_mofa(list(as.matrix(gexp_d), as.matrix(olinc_d), as.matrix(ciber_d),
				as.matrix(reactome_d), as.matrix(msigc2_d), as.matrix(msigc5_d)))
# Define the view names (for use later on, with the visualizations)
views_names(mdd) <- c("gexp", "olinc", "cibersort", "reactome", "msigc2", "msigc5")

# Set data options
data_opts <- get_default_data_options(mdd)
data_opts$use_float32 <- FALSE
# Set model options 
model_opts <- get_default_model_options(mdd)

run_folders <- c() # to keep track of the output folders for each run, in case of multiple runs
for (run_iter in 1:n_runs){
	# If there are to be more than one runs, generate a seeds vector
	if (n_runs > 1)
		seeds <- sample.int(10*n_runs, size=n_runs)
	# Set training options
	train_opts <- get_default_training_options(mdd)
	train_opts$drop_factor_threshold <- 0.01 # Drop factors that explain less than 1% variance
	train_opts$convergence_mode <- "slow" # To improve quality of fit
	train_opts$maxiter <- 4000 # increase from default value of 1000, just in case more iters are required for convergence
	if (run_iter > 1)
		train_opts$seed <- seeds[run_iter]
	
	# Prepare the MOFA object and add the metadata data frame, with one row per sample.
	# Per package requirements, include columns named "sample" and "group", indicating 
	# sample id and group membership, for each sample.
	mo <- prepare_mofa(object = mdd, data_options = data_opts, 
			model_options = model_opts,	training_options = train_opts)
	meta_d <- cbind(data.frame(group = rep("group1", nrow(clin_d)), sample = rownames(clin_d)), clin_d)
	rownames(meta_d) <- meta_d$sample
	meta_d <- meta_d[samples_names(mo)[[1]], , drop = FALSE]
	samples_metadata(mo) <- meta_d
	
	# Tag the name of the output folder with the current date and time, for logging
	# purposes. 
	res_folder <- paste(output_dir,
			ifelse(n_runs > 1, "multi_runs", ""), format(Sys.time(), "%Y_%b_%d_%H_%M_%S"), sep="/")
	run_folders <- c(run_folders, res_folder) # keep track of the output folders for each run, in case of multiple runs
	res_file <- paste(res_folder, "model.hdf5", sep="/")
	log_file <- paste(res_folder, "mofa_out.txt", sep="/")
	dir.create(res_folder, recursive=TRUE)
	
	# Train the MOFA model and save the output of the "run_mofa" method
	sink(file=log_file)
	mo.tr <- run_mofa(mo, res_file, use_basilisk=TRUE)
	sink()
	
	# Logging - write out the session info
	session_file <- paste(res_folder, "sessionInfo.txt", sep="/")
	writeLines(capture.output(sessionInfo()), session_file)	
}

# For multiple runs, perform model selection using the lowest ELBO.
# Extract the factor score matrix from the selected model and save the
# ELBO values from all n runs.
if (n_runs > 1) {
	library(tidyverse)
	library(purrr)
	library(tibble)
	
	# Function to load MOFA model and return final ELBO
	get_final_elbo <- function(model_path) {
	  model <- load_model(model_path)
	  elbo_vals <- model@training_stats$elbo
	  valid_elbo <- elbo_vals[!is.nan(elbo_vals)]
	  final_elbo <- tail(valid_elbo, 1)
	  return(final_elbo)
	}
	
	# Loop over runs and collect ELBOs
	elbo_results <- map_dfr(run_folders, function(folder) {
	  model_file <- file.path(folder, "model.hdf5")
	  if (file.exists(model_file)) {
	    tryCatch({
	      elbo <- get_final_elbo(model_file)
	      tibble(run_folder = folder, final_elbo = elbo)
	    }, error = function(e) {
	      message("Error loading model in folder: ", folder)
	      return(NULL)
	    })
	  } else {
	    message("No model.hdf5 found in folder: ", folder)
	    return(NULL)
	  }
	})
	
	# Identify the best run (highest ELBO)
	best_run <- elbo_results[[which.max(elbo_results$final_elbo), "run_folder"]]
	write.csv(elbo_results,
			  file = file.path(best_run, sprintf("elbo_scores_%d_trials.csv", n_runs)))
	
	# Extract factors & factor scores
	# Load the best model
	best_model <- load_model(file.path(best_run, "model.hdf5"))
	# Extract factor scores (Z matrix) from the model
	factors <- get_factors(best_model, factors = "all", groups = "all")$group1
	# Save the factor matrix
	write.csv(factors, file = file.path(best_run, "factor_scores.csv"))

	cat("Best MOFA run is:", best_run, "\n")
}
