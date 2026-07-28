#!/usr/bin/env Rscript
# Entrypoint for the SECONDARY bulk RNA-seq analysis (see ADR-0003): a
# generalisability check on the larger, independent TCGA-BRCA cohort (30 Primary
# Tumor + 15 Solid Tissue Normal). Runs the full DESeq2 + edgeR + comparison +
# enrichment + report pipeline; its Primary Tumor vs. Solid Tissue Normal design is
# the natural two-group comparison the matched-bulk cohort lacks (see
# run_bulk_pipeline_matched.R). Not fed into Phase 5's concordance metric - see
# docs/bulk-pipeline.md.
#
# Usage: Rscript scripts/run_bulk_pipeline_tcga.R
# (run from the repository root, with the omicsatlas conda environment active;
# needs live network access for KEGG enrichment)

suppressPackageStartupMessages({
  library(OmicsAtlasBulk)
})

repo_root <- getwd()

manifest_path <- file.path(repo_root, "data", ".provenance", "tcga_brca_subset.json")
if (!file.exists(manifest_path)) {
  stop("TCGA-BRCA fetch manifest not found at ", manifest_path, ". Run `make bulk-data-fetch` first.")
}

output_dir <- file.path(repo_root, "data", "processed", "bulk", "tcga_brca")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading TCGA-BRCA STAR-Counts subset...")
bulk_se <- read_tcga_brca_counts(manifest_path)
bulk_se <- filter_low_count_genes(bulk_se, min_total_count = 10)
message(sprintf("  %d genes x %d samples after low-count filtering.", nrow(bulk_se), ncol(bulk_se)))
message(sprintf(
  "  %d Primary Tumor, %d Solid Tissue Normal.",
  sum(SummarizedExperiment::colData(bulk_se)$sample_type == "Primary Tumor"),
  sum(SummarizedExperiment::colData(bulk_se)$sample_type == "Solid Tissue Normal")
))

message("Running DESeq2...")
deseq2_res <- run_deseq2(bulk_se, condition_column = "sample_type")
message(sprintf("  %d significant genes.", sum(deseq2_res$is_significant)))

message("Running edgeR...")
edger_res <- run_edger(bulk_se, condition_column = "sample_type")
message(sprintf("  %d significant genes.", sum(edger_res$is_significant)))

message("Comparing DESeq2 vs. edgeR...")
comparison <- compare_de_methods(deseq2_res, edger_res)
message(sprintf(
  "  %d genes significant in both, %d DESeq2-only, %d edgeR-only.",
  length(comparison$overlap$both), length(comparison$overlap$deseq2_only),
  length(comparison$overlap$edger_only)
))

gene_symbols <- SummarizedExperiment::rowData(bulk_se)$gene_symbol

run_enrichment_for <- function(res, label) {
  sig_genes <- unique(gene_symbols[res$is_significant])
  ranked <- stats::setNames(res$log2FoldChange, gene_symbols)
  ranked <- ranked[!is.na(ranked) & !duplicated(names(ranked))]
  message(sprintf("Running enrichment on %s significant genes (%d genes, live KEGG)...", label, length(sig_genes)))
  run_enrichment(sig_genes, ranked, run_kegg = TRUE)
}
enrichment_deseq2 <- run_enrichment_for(deseq2_res, "DESeq2")
enrichment_edger <- run_enrichment_for(edger_res, "edgeR")

utils::write.csv(deseq2_res, file.path(output_dir, "deseq2_results.csv"), row.names = TRUE)
utils::write.csv(edger_res, file.path(output_dir, "edger_results.csv"), row.names = TRUE)
utils::write.csv(comparison$discrepancy, file.path(output_dir, "de_method_discrepancy.csv"), row.names = FALSE)

report_path <- file.path(output_dir, "tcga_brca_report.html")
message("Rendering report to ", report_path, "...")
render_bulk_report(
  deseq2_res = deseq2_res,
  edger_res = edger_res,
  comparison = comparison,
  enrichment_deseq2 = enrichment_deseq2,
  enrichment_edger = enrichment_edger,
  dataset_label = "TCGA-BRCA subset (30 tumor + 15 normal, secondary/generalisability analysis)",
  output_file = report_path
)

message("Done. Outputs written to ", output_dir)
