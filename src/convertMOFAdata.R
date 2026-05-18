
# =============================================================================
# Script:  convertMOFAdata.R
# Purpose: Convert gene set collection tables from the MOFAdata package
#          (https://github.com/bioFAM/MOFAdata/) from Ensembl gene IDs to
#          NCBI gene IDs, for compatibility with the analysis pipeline.
# =============================================================================
# Background:
#   The MOFAdata package provides gene set membership matrices intended for
#   use with MOFA's gene set enrichment functionality. Columns in these
#   matrices are keyed by Ensembl gene IDs. This script re-keys them to
#   NCBI gene IDs via the convertMOFAGeneSetData() function, which wraps
#   the ensemblToNCBI() utility defined in src/methods.R.
#
# Usage:
#   Call convertMOFAGeneSetData(in_file, out_file) with:
#     in_file  -- path to the original MOFAdata .rda file
#     out_file -- path where the NCBI-keyed .rda file will be saved
#
#   This script is intended to be run once. Precomputed versions of both
#   the original and converted tables are already available under the
#   data/annotations/ directory.
#
# Key steps (within convertMOFAGeneSetData):
#   1. Load the gene set matrix from the input .rda file.
#   2. Remove gene sets with missing (NA) row names (QC step).
#   3. Map Ensembl column IDs to NCBI gene IDs via ensemblToNCBI().
#   4. Drop columns that could not be mapped; remove gene sets that become
#      empty after remapping.
#   5. Save the converted matrix to out_file, preserving the original
#      variable name.
#
# Inputs:
#   - Any MOFAdata gene set .rda file (e.g., Reactome or MSigDB collections)
#
# Outputs:
#   - NCBI-keyed .rda file written to out_file
#
# Dependencies: here
#               Sources src/methods.R
# =============================================================================

library(here)

# Location of the methods.R file; this is included in the GitHub repository. 
source(here("src", "methods.R"))

convertMOFAGeneSetData <- function(in_file, out_file){
	if (!file.exists(in_file))
		stop(paste("Input file", in_file, "does not exist."))
	
	gs_mat_name = load(in_file)
	gs_mat = get(gs_mat_name)
	
	# A bit of QC, as at least one gene sets in the reactome collection that they
	# provide is not assigned a name
	gs_mat = gs_mat[!(is.na(rownames(gs_mat))),]
	
	# Replace Ensembl ids with NCBI gene ids
	gids = ensemblToNCBI(colnames(gs_mat))
	gs_mat = gs_mat[, names(gids)]
	colnames(gs_mat) = gids
	
	# Remove empty gene sets, in case any were introduced
	gs_mat = gs_mat[apply(gs_mat, 1, sum) > 0,]
	
	# Save the converted table using its original variable name
	assign(gs_mat_name, gs_mat)
	save(list=gs_mat_name, file=out_file)
}

