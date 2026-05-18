# annotations/ — Gene Set Annotation Files

This folder contains gene set annotation files used by the enrichment scoring
and MOFA+ analysis steps in `src/`.

## Files in this folder

These files were converted from the originals in `annotations/MOFAdata/` using
`src/convertMOFAdata.R::convertMOFAGeneSetData()`, which re-keys gene set
membership from Ensembl IDs to NCBI (Entrez) IDs. The converted files are used
by `gsEnrichmentAnalysis()` in `src/methods.R`.

- `reactomeGS.rda` — Reactome gene set collection (Entrez IDs)
- `MSigDB_C2_human.rda` — MSigDB C2 curated gene sets (Entrez IDs)
- `MSigDB_C5_human.rda` — MSigDB C5 GO gene sets (Entrez IDs)

## MOFAdata/ subfolder

`annotations/MOFAdata/` contains the original source files downloaded from the
[MOFAdata GitHub repository](https://github.com/bioFAM/MOFAdata/tree/master/data).
These files use Ensembl gene IDs and are the inputs to `convertMOFAdata.R`.

- `MSigDB_v6.0_C2_human.RData`
- `MSigDB_v6.0_C5_human.RData`
- `reactomeGS.rda`
