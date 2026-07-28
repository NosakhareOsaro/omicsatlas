# render_bulk_report() shells out to the real Quarto CLI to render the bundled
# inst/quarto/bulk_report.qmd template - not mocked, so this test genuinely proves
# the template renders against real pipeline output. Skipped (rather than faked) if
# the Quarto CLI isn't available in the environment running the tests, matching the
# skip_on_ci()/skip_if_offline() pattern already used for the KEGG enrichment test.

skip_if_no_quarto <- function() {
  if (!requireNamespace("quarto", quietly = TRUE) || !quarto::quarto_available()) {
    testthat::skip("Quarto CLI not available")
  }
}

test_that("render_bulk_report() renders a real HTML report from real pipeline output", {
  skip_if_no_quarto()
  data(example_bulk_se, package = "OmicsAtlasBulk", envir = environment())

  deseq2_res <- run_deseq2(example_bulk_se, condition_column = "condition")
  edger_res <- run_edger(example_bulk_se, condition_column = "condition")
  comparison <- compare_de_methods(deseq2_res, edger_res)

  sig_genes <- SummarizedExperiment::rowData(example_bulk_se)$gene_symbol[deseq2_res$is_significant]
  ranked <- stats::setNames(
    deseq2_res$log2FoldChange,
    SummarizedExperiment::rowData(example_bulk_se)$gene_symbol
  )
  enrichment_deseq2 <- run_enrichment(sig_genes, ranked, run_kegg = FALSE)

  output_file <- tempfile(fileext = ".html")

  result <- render_bulk_report(
    deseq2_res = deseq2_res,
    edger_res = edger_res,
    comparison = comparison,
    enrichment_deseq2 = enrichment_deseq2,
    dataset_label = "example_bulk_se (test)",
    output_file = output_file
  )

  expect_identical(result, output_file)
  expect_true(file.exists(output_file))
  expect_gt(file.info(output_file)$size, 0)

  html <- paste(readLines(output_file, warn = FALSE), collapse = "\n")
  expect_true(grepl("example_bulk_se \\(test\\)", html))
  expect_true(grepl("DESeq2", html))
})

test_that("render_bulk_report() requires the Quarto CLI to be available", {
  testthat::local_mocked_bindings(
    quarto_available = function(...) FALSE,
    .package = "quarto"
  )

  expect_error(
    render_bulk_report(
      deseq2_res = data.frame(),
      edger_res = data.frame(),
      comparison = list(),
      dataset_label = "irrelevant",
      output_file = tempfile(fileext = ".html")
    ),
    "Quarto CLI is not available"
  )
})
