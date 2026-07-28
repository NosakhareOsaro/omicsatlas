#' Run DESeq2 differential expression
#'
#' Wraps DESeq2's standard `DESeqDataSetFromMatrix()` -> `DESeq()` -> `results()`
#' workflow for a simple two-level comparison.
#'
#' @param se A `SummarizedExperiment` with a raw-count assay.
#' @param condition_column Name of a two-level `colData(se)` column defining the
#'   comparison groups (e.g. `"condition"`).
#' @param assay_name Which assay in `se` holds counts. Default `"counts"`.
#' @param alpha Significance threshold passed to `DESeq2::results()` for its
#'   independent-filtering step, and used for the `is_significant` column below.
#'
#' @return A data frame (one row per gene) with DESeq2's standard columns
#'   (`baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`) plus
#'   `is_significant` (`padj < alpha`, `FALSE` if `padj` is `NA`).
#'
#' Non-integer counts (e.g. RSEM/Salmon-style estimated counts, as in the
#' GSE176078-matched bulk dataset — see `ADR-0003`) are rounded before use; DESeq2
#' requires integer counts.
#'
#' @examples
#' data(example_bulk_se)
#' res <- run_deseq2(example_bulk_se, condition_column = "condition")
#' head(res[order(res$padj), ])
#' @export
run_deseq2 <- function(se, condition_column, assay_name = "counts", alpha = 0.05) {
  counts_mat <- round(as.matrix(SummarizedExperiment::assay(se, assay_name)))
  storage.mode(counts_mat) <- "integer"

  col_data <- as.data.frame(SummarizedExperiment::colData(se))
  col_data[[condition_column]] <- factor(col_data[[condition_column]])
  if (nlevels(col_data[[condition_column]]) != 2) {
    stop(
      "`condition_column` must have exactly 2 levels, got ",
      nlevels(col_data[[condition_column]]), ": ",
      paste(levels(col_data[[condition_column]]), collapse = ", ")
    )
  }

  design <- stats::as.formula(paste("~", condition_column))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_mat, colData = col_data, design = design)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- as.data.frame(DESeq2::results(dds, alpha = alpha))

  res$is_significant <- !is.na(res$padj) & res$padj < alpha
  res
}
