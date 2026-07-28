#!/usr/bin/env Rscript
# Entrypoint for the PRIMARY bulk RNA-seq analysis (see ADR-0003): BayesPrism
# deconvolution of the 24-patient GSE176078-matched bulk cohort against the exact
# Phase 1 scRNA-seq signature artifact. This cohort has no built-in tumour/normal
# split (every sample is tumour tissue from a patient also in the Phase 1 scRNA
# data), so unlike the TCGA-BRCA secondary pipeline (run_bulk_pipeline_tcga.R) there
# is no DESeq2/edgeR/enrichment step here - its role in this project is providing
# same-patient deconvolution proportions for Phase 5's concordance benchmark, not a
# two-group DE comparison. See docs/bulk-pipeline.md.
#
# Usage: Rscript scripts/run_bulk_pipeline_matched.R
# (run from the repository root, with the omicsatlas conda environment active)

suppressPackageStartupMessages({
  library(OmicsAtlasBulk)
})

# Run from the repository root (matches the Makefile targets this script is
# invoked from, and Phase 1's `python -m omicsatlas.scrna.pipeline` convention).
repo_root <- getwd()

counts_path <- file.path(
  repo_root, "data", "raw", "bulk", "gse176078_matched",
  "GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt"
)
if (!file.exists(counts_path)) {
  stop(
    "Matched-bulk counts file not found at ", counts_path, ". Run `make bulk-data-fetch` ",
    "and `make bulk-data-extract` first."
  )
}

# Resolves the exact Phase 1 scRNA-seq signature artifact path via Python's
# signature_path() contract (src/omicsatlas/scrna/artifact.py) rather than
# hardcoding or reimplementing that path logic in R - see DATA_SCHEMA.md.
resolve_signature_path <- function() {
  out <- system2(
    "python",
    c("-c", "from omicsatlas.scrna.artifact import signature_path; print(signature_path())"),
    stdout = TRUE
  )
  trimws(out[length(out)])
}

message("Loading matched-bulk count matrix...")
bulk_se <- read_matched_bulk_counts(counts_path)
bulk_se <- filter_low_count_genes(bulk_se, min_total_count = 10)
message(sprintf("  %d genes x %d samples after low-count filtering.", nrow(bulk_se), ncol(bulk_se)))

signature_h5ad_path <- resolve_signature_path()
if (!file.exists(signature_h5ad_path)) {
  stop(
    "Phase 1 scRNA-seq signature artifact not found at ", signature_h5ad_path, ". ",
    "Run `make scrna-signature` first (see docs/scrna-pipeline.md)."
  )
}
message("Loading Phase 1 scRNA-seq signature artifact from ", signature_h5ad_path, "...")
reference_sce <- load_scrna_signature(signature_h5ad_path)

message("Running BayesPrism deconvolution (this takes several minutes)...")
proportions <- run_bayesprism_deconvolution(
  bulk_se = bulk_se,
  reference_sce = reference_sce
)

output_dir <- file.path(repo_root, "data", "processed", "bulk", "gse176078_matched")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "bayesprism_cell_type_proportions.csv")
utils::write.csv(proportions, output_path, row.names = TRUE)
message("Wrote deconvolution proportions to ", output_path)

print(round(proportions, 3))
