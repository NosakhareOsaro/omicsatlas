#' Run edgeR differential expression
#'
#' Wraps edgeR's standard `DGEList()` -> `calcNormFactors()` -> `estimateDisp()` ->
#' `glmQLFit()` -> `glmQLFTest()` quasi-likelihood workflow for a simple two-level
#' comparison — run on the *same* count matrix `run_deseq2()` uses, so the two
#' methods' results are directly comparable (see `compare_de_methods()`).
#'
#' @param se A `SummarizedExperiment` with a raw-count assay.
#' @param condition_column Name of a two-level `colData(se)` column defining the
#'   comparison groups (e.g. `"condition"`).
#' @param assay_name Which assay in `se` holds counts. Default `"counts"`.
#' @param alpha Significance threshold used for the `is_significant` column below
#'   (`FDR < alpha`).
#'
#' @return A data frame (one row per gene) with edgeR's standard `topTags` columns
#'   (`logFC`, `logCPM`, `F`, `PValue`, `FDR`) plus `is_significant`.
#'
#' Non-integer counts (e.g. RSEM/Salmon-style estimated counts — see `ADR-0003`) are
#' rounded before use, matching `run_deseq2()`.
#'
#' @examples
#' data(example_bulk_se)
#' res <- run_edger(example_bulk_se, condition_column = "condition")
#' head(res[order(res$FDR), ])
#' @export
run_edger <- function(se, condition_column, assay_name = "counts", alpha = 0.05) {
  counts_mat <- round(as.matrix(SummarizedExperiment::assay(se, assay_name)))

  col_data <- as.data.frame(SummarizedExperiment::colData(se))
  col_data[[condition_column]] <- factor(col_data[[condition_column]])
  if (nlevels(col_data[[condition_column]]) != 2) {
    stop(
      "`condition_column` must have exactly 2 levels, got ",
      nlevels(col_data[[condition_column]]), ": ",
      paste(levels(col_data[[condition_column]]), collapse = ", ")
    )
  }

  dge <- edgeR::DGEList(counts = counts_mat, group = col_data[[condition_column]])
  dge <- edgeR::calcNormFactors(dge)

  design_matrix <- stats::model.matrix(stats::as.formula(paste("~", condition_column)), data = col_data)
  dge <- edgeR::estimateDisp(dge, design_matrix)
  fit <- edgeR::glmQLFit(dge, design_matrix)
  qlf <- edgeR::glmQLFTest(fit, coef = ncol(design_matrix))

  res <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  res$is_significant <- !is.na(res$FDR) & res$FDR < alpha
  res
}
