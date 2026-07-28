# Builds the small synthetic bulk RNA-seq SummarizedExperiment bundled with the
# package (`data(example_bulk_se)`), used throughout roxygen2 @examples blocks and
# testthat tests instead of the real ~500MB+ alignment inputs — deterministic (fixed
# seed) so every example/test run is reproducible.
#
# 20 genes x 8 samples, 2 conditions x 4 replicates, with 5 genes given a real mean
# shift between conditions so DESeq2/edgeR/enrichment examples have a genuine (if
# tiny) signal to find, not pure noise.

set.seed(1)

n_genes <- 20
n_samples_per_condition <- 4
gene_ids <- sprintf("ENSG%011dEX", seq_len(n_genes))
gene_symbols <- c(
  "TP53", "BRCA1", "BRCA2", "ERBB2", "ESR1", "PGR", "MKI67", "EGFR", "PTEN", "MYC",
  "CCND1", "CDKN2A", "RB1", "PIK3CA", "AKT1", "GATA3", "FOXA1", "KRT8", "KRT18", "VIM"
)

sample_ids <- c(sprintf("control_%d", 1:n_samples_per_condition), sprintf("treated_%d", 1:n_samples_per_condition))
condition <- factor(rep(c("control", "treated"), each = n_samples_per_condition), levels = c("control", "treated"))

baseline <- stats::rgamma(n_genes, shape = 4, rate = 0.5) + 5
shifted_genes <- 1:5 # first 5 genes carry a real condition effect
counts <- matrix(0L, nrow = n_genes, ncol = length(sample_ids), dimnames = list(gene_ids, sample_ids))
for (j in seq_along(sample_ids)) {
  gene_means <- baseline
  if (condition[j] == "treated") {
    gene_means[shifted_genes] <- gene_means[shifted_genes] * 4
  }
  counts[, j] <- stats::rpois(n_genes, lambda = gene_means)
}

col_data <- S4Vectors::DataFrame(condition = condition, row.names = sample_ids)
row_data <- S4Vectors::DataFrame(gene_symbol = gene_symbols, row.names = gene_ids)

example_bulk_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = counts),
  colData = col_data,
  rowData = row_data
)

usethis::use_data(example_bulk_se, overwrite = TRUE)
