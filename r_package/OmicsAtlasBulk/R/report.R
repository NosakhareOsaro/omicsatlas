#' Render an auto-generated Quarto HTML report from bulk pipeline outputs
#'
#' Renders a full HTML report - a differential expression summary and top genes for
#' both methods, the DESeq2/edgeR overlap comparison (Venn diagram, discrepancy
#' table), functional enrichment results, and (optionally) BayesPrism deconvolution
#' proportions - entirely from the R objects the pipeline itself already produced
#' ([run_deseq2()], [run_edger()], [compare_de_methods()], [run_enrichment()],
#' [run_bayesprism_deconvolution()]). The bundled Quarto template
#' (`inst/quarto/bulk_report.qmd`) has no hardcoded numbers or narrative baked in -
#' every table, count, and plot it shows is computed live from whatever objects are
#' passed in here, so re-running this against fresh pipeline output can never leave
#' stale numbers behind.
#'
#' Requires the Quarto CLI to be installed and on `PATH` (checked via
#' `quarto::quarto_available()`); the `quarto` R package is a `Suggests` dependency
#' rather than `Imports` since it wraps an external CLI this package doesn't itself
#' install, following the standard R packaging convention for that case.
#'
#' @param deseq2_res Result of [run_deseq2()].
#' @param edger_res Result of [run_edger()], on the same data.
#' @param comparison Result of `compare_de_methods(deseq2_res, edger_res)`.
#' @param enrichment_deseq2 Result of [run_enrichment()] on DESeq2's significant
#'   genes (a list with `go`/`kegg`/`gsea`), or `NULL` to omit that section.
#' @param enrichment_edger Result of [run_enrichment()] on edgeR's significant
#'   genes, or `NULL` to omit that section.
#' @param deconvolution_props Result of [run_bayesprism_deconvolution()] (a
#'   samples-by-cell-type proportion matrix), or `NULL` to omit that section.
#' @param dataset_label Short human-readable label for the dataset this report
#'   covers (e.g. `"GSE176078-matched bulk (primary)"`), used in the report title.
#' @param output_file Path to write the rendered HTML report to.
#' @param quiet Passed to `quarto::quarto_render()`. Default `TRUE`.
#'
#' @return `output_file`, invisibly.
#'
#' @examples
#' \donttest{
#' # Requires the Quarto CLI to be installed; see tests/testthat/test-report.R for a
#' # complete, real, executed run on the bundled example_bulk_se fixture.
#' data(example_bulk_se)
#' deseq2_res <- run_deseq2(example_bulk_se, condition_column = "condition")
#' edger_res <- run_edger(example_bulk_se, condition_column = "condition")
#' comparison <- compare_de_methods(deseq2_res, edger_res)
#' render_bulk_report(
#'   deseq2_res, edger_res, comparison,
#'   dataset_label = "example_bulk_se",
#'   output_file = tempfile(fileext = ".html")
#' )
#' }
#' @export
render_bulk_report <- function(deseq2_res,
                                edger_res,
                                comparison,
                                enrichment_deseq2 = NULL,
                                enrichment_edger = NULL,
                                deconvolution_props = NULL,
                                dataset_label,
                                output_file,
                                quiet = TRUE) {
  if (!requireNamespace("quarto", quietly = TRUE)) {
    stop("The 'quarto' package is required to render this report. Install it with ",
      "install.packages(\"quarto\").",
      call. = FALSE
    )
  }
  if (!quarto::quarto_available()) {
    stop("The Quarto CLI is not available on PATH; install it to render this report.",
      call. = FALSE
    )
  }

  work_dir <- tempfile("bulk_report_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  data_path <- file.path(work_dir, "report_data.rds")
  saveRDS(
    list(
      deseq2_res = deseq2_res,
      edger_res = edger_res,
      comparison = comparison,
      enrichment_deseq2 = enrichment_deseq2,
      enrichment_edger = enrichment_edger,
      deconvolution_props = deconvolution_props,
      dataset_label = dataset_label
    ),
    data_path
  )

  template <- system.file("quarto", "bulk_report.qmd", package = "OmicsAtlasBulk")
  if (!nzchar(template)) {
    stop("Could not find the bundled bulk_report.qmd template.", call. = FALSE)
  }
  rendered_qmd <- file.path(work_dir, "bulk_report.qmd")
  file.copy(template, rendered_qmd)

  output_filename <- basename(output_file)
  quarto::quarto_render(
    input = rendered_qmd,
    execute_params = list(report_data_path = data_path),
    output_file = output_filename,
    quiet = quiet
  )

  file.copy(file.path(work_dir, output_filename), output_file, overwrite = TRUE)
  invisible(output_file)
}
