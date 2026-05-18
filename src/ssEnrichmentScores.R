# =============================================================================
# Script:  ssEnrichmentScores.R
# Purpose: Compute single-sample gene set enrichment scores (GSVA or ssGSEA)
#          from VST-normalized RNA-seq data, using gene set collections
#          converted from the MOFAdata package. Output matrices are used as
#          input views in MOFA_analysis.R.
# =============================================================================
# IMPORTANT: cohort1_tororo_all_274_raw_ensembl.csv needs to be accessed via dbGaP (see README)
# 
# Inputs:
#   - data/cohort1_tororo_all_274_raw_ensembl.csv   [controlled access via dbGaP]
#       Raw RNA-seq count matrix for the discovery cohort (Ensembl IDs).
#   - annotations/reactomeGS.rda
#   - annotations/MSigDB_C2_human.rda
#   - annotations/MSigDB_C5_human.rda
#       NCBI-keyed gene set membership matrices produced by convertMOFAdata.R.
#
# Key steps:
#   1. Load and preprocess RNA-seq counts: Ensembl-to-NCBI ID conversion, low-
#      count gene removal, and VST normalization (via methods.R utilities).
#   2. For each gene set collection (Reactome, MSigDB C2, MSigDB C5), build a
#      gene set list and compute single-sample enrichment scores using GSVA
#      (default) or ssGSEA. Gene sets with fewer than 5 or more than 500 genes
#      are excluded.
#   3. Write each enrichment score matrix to a tab-separated file.
#
# Outputs (written to results/<cohort>/<method>/):
#   - cohort_<cohort>_<method>_reactome.tsv      Reactome enrichment scores
#   - cohort_<cohort>_<method>_msigdb_c2.tsv     MSigDB C2 enrichment scores
#   - cohort_<cohort>_<method>_msigdb_c5.tsv     MSigDB C5 enrichment scores
#   - sessionInfo.txt                             R session information
#
# User configuration:
#   - enrichment_method: set to "GSVA" (default) or "ssGSEA"
#   - cohort:            descriptive label used to name the results subfolder
#
# Dependencies: GSVA, here
#               Sources src/methods.R
# =============================================================================

library(here)

source(here("src", "methods.R"), chdir=TRUE)

# Location of the converted MOFAdata files specifying the gene sets.
# See src/convertMOFAdata.R for context.
gs_collection_names <- c("reactome", "msigdb_c2", "msigdb_c5")
mofad_files <- setNames(c(here("annotations", "reactomeGS.rda"), 
				here("annotations", "MSigDB_C2_human.rda"), 
				here("annotations", "MSigDB_C5_human.rda")), 
		gs_collection_names)

# Compute single-sample gene set enrichment scores using GSVA or ssGSEA.
#
# Arguments:
#   dd:        VST-normalized gene expression matrix
#   method:    Enrichment method to use ("GSVA" or "ssGSEA")
#   gset:      Gene set collection
#   save:      Logical; whether to write results to file
#   out_file:  Output file path (used if save = TRUE)
#   verbose:   Logical; print progress messages
#
# Returns:
#   A matrix with gene sets as rows and samples as columns. Each entry
#   is the enrichment score of a gene set in a given sample.
calcSSscore <- function(dd, method=c("GSVA", "ssGSEA"), 
		gset = gs_collection_names, save=FALSE, out_file = NULL,
		verbose = FALSE){
	
	# Load required library
	require(GSVA)
	
	# Check that converted MOFAdata files exist
	nothing <- sapply(mofad_files, function(fname){
				if (!file.exists(fname))
					stop("Converted MOFAdata files are missing")
			})
	
	# Prepare the gene sets lists, as required by GSVA
	gs_mat <- get(load(mofad_files[gset[1]]))
	gs_list <- lapply(1:nrow(gs_mat), function(idx){
				return(colnames(gs_mat)[gs_mat[idx, ] != 0])
			})
	names(gs_list) = rownames(gs_mat)
	
	# Set up the parameters object to use. Only consider gene sets with a 
	# minimum size of 5 and maximum size of 500.
	if (method[1] == "GSVA")
		params <- gsvaParam(dd, gs_list, minSize = 5, maxSize = 500)
	else if (method[1] == "ssGSEA")
		params <- ssgseaParam(dd, gs_list, minSize = 5, maxSize = 500)
	else
		stop("Invalid value in 'method' argument.")
	
	# Run the analysis
	gse_scores <- gsva(params, verbose=verbose)
	
	if (save)
		write.table(gse_scores, file=out_file, quote=FALSE, sep="\t", row.names=TRUE)
	
	return(gse_scores)
}

# Location of RNA-seq data file. 
rnaseq_file <- here("data", "cohort1_tororo_all_274_raw_ensembl.csv") # controlled access via dbGaP.
# Descriptive name of cohort, to be used when creating results folder
cohort <- "tororo"

# Open file, VST-normalize, convert IDs from Ensembl to NCBI, and remove genes 
# with no variation
cnt_d <- prepRNAseq(rnaseq_file, verbose=TRUE)
cnt_vst <- normalizeRNASEQwithVST(cnt_d)
cnt_vst <- removeNoVarRows(cnt_vst)

# Specify single-sample enrichment method to use (GSVA or ssGSEA)
enrichment_method <- "GSVA"

# Where to save the results
res_folder <- here("results", cohort, enrichment_method)
dir.create(res_folder, recursive=TRUE)

for (gs_set in gs_collection_names){
	f_name <- paste(paste("cohort", cohort, enrichment_method, gs_set, sep="_"),
			"tsv", sep=".")
	res_file <- paste(res_folder, f_name, sep="/")
	gs_scores <- calcSSscore(as.matrix(cnt_vst), method=enrichment_method, 
			gset = gs_set, save=TRUE, out_file = res_file)
}

# Logging - write out the session info
session_file <- paste(res_folder, "sessionInfo.txt", sep="/")
writeLines(capture.output(sessionInfo()), session_file)

