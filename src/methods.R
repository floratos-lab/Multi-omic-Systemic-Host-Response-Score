# =============================================================================
# Script:  methods.R
# Purpose: Shared utility functions sourced by other scripts in this pipeline.
#          Not intended to be run directly.
# =============================================================================
# Functions defined:
#
#   ensemblToNCBI(ensIds)
#     Maps a character vector of Ensembl gene IDs (human) to NCBI Entrez gene
#     IDs. Strips version suffixes, discards unmapped IDs, and resolves
#     one-to-many mappings by retaining the first occurrence. Returns a named
#     character vector (names = Ensembl IDs, values = Entrez IDs), or NULL if
#     no valid IDs are found.
#
#   prepRNAseq(f_name, verbose, thresh, out_file)
#     Loads a raw RNA-seq count matrix (CSV, genes x samples), converts row
#     Ensembl IDs to NCBI gene IDs via ensemblToNCBI(), resolves duplicate
#     mappings by retaining the Ensembl ID with the highest median count, and
#     drops genes whose counts do not exceed `thresh` in any sample. Optionally
#     saves the resulting matrix and prints progress messages. Returns the
#     filtered, NCBI-keyed count matrix.
#
#   removeNoVarRows(dmat, min_var)
#     Removes rows from a matrix whose variance is at or below `min_var`
#     (default: 0). Returns the filtered matrix.
#
#   fillEmptyCols(dmat, sids)
#     Adds NA-filled columns to a matrix for any sample IDs in `sids` that are
#     absent from the matrix, then reorders columns to match `sids` exactly.
#     Used to align data views to a common sample set before MOFA input.
#
#   normalizeRNASEQwithVST(readCount, filter_thresh)
#     Applies DESeq2 variance-stabilizing transformation (VST) to a raw count
#     matrix (genes x samples). Returns the VST-transformed data frame.
#
# Internal state:
#   utils.env  -- package-level environment (created once) caching org.Hs.eg.db
#                 key lookups for gene symbols, aliases, and Ensembl IDs, used
#                 by ensemblToNCBI().
#
# Dependencies: annotate, org.Hs.eg.db, AnnotationDbi, DESeq2, here
# =============================================================================

# Load required libraries
library(annotate)
library(org.Hs.eg.db)
library(here)

# Variables to be used by the methods that map back and forth betweer ensembl ids
# and gene ids/gene symbols
if (!exists("utils.env")){
	utils.env <- new.env()
	utils.env$gene_symbols <- keys(org.Hs.eg.db,  keytype="SYMBOL")
	utils.env$gene_aliases <- keys(org.Hs.eg.db,  keytype="ALIAS")
	utils.env$ensembl_ids <- keys(org.Hs.eg.db,  keytype="ENSEMBL")
}


# Takes as input a vector "ensIds" of character strings representing ENSEMBL ids 
# for human genes, e.g., ensIds = c("ENSG00000180346", "ENSG00000157404.15", ...), 
# and returns a vector "res"of character strings representing NCBI gene IDs, where:
# 	res[i] = NCBI ID correponding to ENSEMBL id ensIds[i]
#	names(res)[i] = ensIds[i]
# NOTE: if ensIds[i] does not map to a gene id, res[i] is removed. Accordingly, 
# the lenght of "res"can be shorter than the length of "ensIds"

ensemblToNCBI <- function(ensIds){
	if(sum(grepl("ENSG", ensIds)) != length(ensIds)){
		warning("Argument list contains non-ENSEMBL human gene id entries")
		return(NULL)
	}
		
	# Remove the "version" portion from the ENSEMBL ids, if there.
	ensIds <- sapply(strsplit(ensIds, ".", fixed=TRUE), function(t){return(t[1])})
	# Get rid of invalid ensembl ids
	ensIds <- ensIds[ensIds %in% utils.env$ensembl_ids]
	if (length(ensIds) == 0)
		return(NULL)
	map <- as.matrix(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=ensIds,
							keytype="ENSEMBL", columns=c("ENTREZID"))))
	# Remove duplicate mappings, if any
	t <- which(duplicated(map[,1]))
	if (length(t) > 0)
		map <- map[-t, , drop=FALSE]
	res <- map[,2]
	names(res) <- map[, 1]
	return(res)
}


# Load an RNA-seq data file and replace Ensembl IDs with NCBI gene IDs.
# Ensembl IDs that do not map uniquely to a gene ID are discarded.
#
# Rows with total read counts across all samples less than or equal to
# "thresh" are removed. This helps filter out low-quality features.
#
# If "out_file" is provided, the resulting gene expression matrix is
# written to that location. If "verbose" is TRUE, progress messages are
# printed.
#
# The goal is to remove pseudogenes and other low-quality features.
#
# Return the resulting matrix.
prepRNAseq <- function(f_name, verbose = FALSE, thresh = 2, out_file = NULL){
	dd <- read.csv(f_name, row.names=1, check.names=FALSE)
	
	# convert Ensembl ids to gene ids.
	eids <- ensemblToNCBI(rownames(dd))
	nomap <- setdiff(rownames(dd), names(eids))
	if (verbose)
		writeLines(paste("Number of Ensembl ids in input file= ", nrow(dd), 
						"\nNumber of Ensembl ids mapping to gene ids =", length(eids), 
						"\nNumber of unamapped Ensembl ids=", length(nomap), "\n"))
	
	# Check for duplicates, where two different Ensembl ids map to the same gene id. 
	# In such cases, keep the one with the highest median read count across samples.
	dups <- duplicated(eids)
	if (sum(dups) > 0){
		if (verbose)
			writeLines(paste("Among the mapped Ensembl ids,there are", sum(dups), 
							"duplicates, i.e., ids that\ndo not uniquely map to", 
							"a gene id. For each set of duplicates, we retain\nthe",
							"Ensembl id with the highest median read count across samples.\n"))
		dups <- unique(eids[duplicated(eids)])
		drop <- unlist(lapply(dups, function(eid){
							xx <- apply(dd[names(eids)[eids==eid],], 1, median)
							return(names(sort(xx, decreasing=TRUE))[2:length(xx)])
						}))
		eids <- eids[setdiff(names(eids), drop)]
	}
	
	# Create a new expression matrix, using only entries that map to gene ids
	ddn <- dd[names(eids), , drop = FALSE]
	rownames(ddn) <- eids

	# Drop genes with no more than "thresh" counts in each sample
	keep <- apply(ddn, 1, function(r){return(sum(r > thresh ) > 0)})
	ddn <- ddn[keep, , drop = FALSE]
	
	if (verbose)
		writeLines(paste("Gene ids with with a read count of", thresh, "or fewer across *ALL* samples are removed:", 
						"\n\tNumber of gene ids (rows) comprising the final read count matrix =", sum(keep), 
						"\n\tProportion of original reads retained in the final read count matrix= ", sum(ddn)/sum(dd)))
	
	if (!is.null(out_file))
		save(ddn, file=out_file)
	
	return(ddn)
}


# remove rows with variance below a given threshold
removeNoVarRows <- function(dmat, min_var = 0){
	drop <- (apply(dmat, 1, var) <= min_var)
	return(dmat[!drop, , drop = FALSE])
}

# Given a data matrix and a set of sample ids, generate NA colunms for 
# those samples not already in the matrix
fillEmptyCols <- function(dmat, sids){
	new_cols <- setdiff(sids, colnames(dmat))
	nc_mat <- matrix(NA, nrow(dmat), length(new_cols))
	rownames(nc_mat) <- rownames(dmat)
	colnames(nc_mat) <- new_cols
	dmat <- cbind(dmat, nc_mat)
	return(dmat[, sids, drop = FALSE])
}

# A function to transform RNA-Seq data with VST in the DESeq2 package
# readCount: RNA-Seq rawcounts in a matrix or in a data frame form
#            Rows are genes and columns are samples
# Returns the vst-transformed data matrix
normalizeRNASEQwithVST <- function(readCount, filter_thresh = 1) {
    
	# Load required library
	require(DESeq2)
    
	# make a design matrix for DESeq2 data
	condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
	# Data preparation for DESeq2 format
	deSeqData <- DESeqDataSetFromMatrix(countData = readCount, colData = condition, design = ~0)
    
	# VST
	vsd <- vst(deSeqData)
	transCnt <- data.frame(assay(vsd), check.names = FALSE)
    
	return(transCnt)
    
}
