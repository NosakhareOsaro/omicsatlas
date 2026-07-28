test_that("run_deseq2 returns one row per gene with expected columns", {
  data(example_bulk_se, package = "OmicsAtlasBulk", envir = environment())

  res <- run_deseq2(example_bulk_se, condition_column = "condition")

  expect_equal(nrow(res), nrow(example_bulk_se))
  for (col in c("baseMean", "log2FoldChange", "pvalue", "padj", "is_significant")) {
    expect_true(col %in% colnames(res))
  }
})

test_that("run_deseq2 recovers the genes with a real injected condition effect", {
  data(example_bulk_se, package = "OmicsAtlasBulk", envir = environment())

  res <- run_deseq2(example_bulk_se, condition_column = "condition")

  # The fixture's first 5 genes (see data-raw/example_bulk_se.R) have a real 4x
  # mean-expression shift between conditions; they should show a large positive
  # log2FoldChange, unlike the unshifted genes.
  shifted <- res[1:5, ]
  unshifted <- res[6:20, ]
  expect_true(mean(shifted$log2FoldChange) > mean(unshifted$log2FoldChange))
})

test_that("run_deseq2 rounds non-integer counts rather than erroring", {
  data(example_bulk_se, package = "OmicsAtlasBulk", envir = environment())
  se_fractional <- example_bulk_se
  SummarizedExperiment::assay(se_fractional, "counts") <-
    SummarizedExperiment::assay(se_fractional, "counts") + 0.37

  expect_no_error(run_deseq2(se_fractional, condition_column = "condition"))
})

test_that("run_deseq2 requires exactly 2 condition levels", {
  data(example_bulk_se, package = "OmicsAtlasBulk", envir = environment())
  se_three_levels <- example_bulk_se
  SummarizedExperiment::colData(se_three_levels)$condition <- factor(
    rep(c("a", "b", "c"), length.out = ncol(se_three_levels))
  )

  expect_error(run_deseq2(se_three_levels, condition_column = "condition"), "exactly 2 levels")
})
