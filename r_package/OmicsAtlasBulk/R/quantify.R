#' Quantify gene-level counts with featureCounts
#'
#' Wraps `Rsubread::featureCounts()` (Bioconductor, in-process — no external binary,
#' unlike [run_star_align()]; see `ADR-0004`) to quantify aligned reads against a GTF
#' annotation, returning a `SummarizedExperiment` rather than featureCounts' native
#' list-of-matrices output, matching this package's `SummarizedExperiment`-based
#' public API convention.
#'
#' @param bam_files Character vector of BAM file paths.
#' @param annotation_gtf Path to a GTF annotation file.
#' @param sample_names Optional character vector naming each BAM file's sample (used
#'   as the result's column names); defaults to the BAM file basenames with `.bam`/
#'   `.BAM` stripped.
#' @param ... Additional arguments passed through to `Rsubread::featureCounts()`
#'   (e.g. `isPairedEnd`, `GTF.featureType`, `GTF.attrType`).
#'
#' @return A `SummarizedExperiment` with the gene-by-sample count matrix in
#'   `assays$counts` and featureCounts' per-gene annotation (`Length`, etc.) in
#'   `rowData`.
#'
#' @examples
#' # A real, small BAM + matching GTF are bundled with the package (built once via
#' # Rsubread::align() against a tiny synthetic reference — alignment itself has a
#' # ~3-minute fixed per-call overhead in Rsubread regardless of input size, so
#' # rebuilding it on every example/test run isn't practical; the *result* is
#' # real and bundled instead, and featureCounts() itself runs on it directly here,
#' # not stubbed).
#' bam <- system.file("extdata", "toy_aligned.bam", package = "OmicsAtlasBulk")
#' gtf <- system.file("extdata", "toy_annotation.gtf", package = "OmicsAtlasBulk")
#' se <- run_featurecounts(bam, gtf, sample_names = "toy_sample", isPairedEnd = FALSE)
#' SummarizedExperiment::assay(se, "counts")
#' @export
run_featurecounts <- function(bam_files, annotation_gtf, sample_names = NULL, ...) {
  fc <- Rsubread::featureCounts(
    files = bam_files,
    annot.ext = annotation_gtf,
    isGTFAnnotationFile = TRUE,
    ...
  )

  if (is.null(sample_names)) {
    sample_names <- sub("\\.(bam|BAM)$", "", basename(bam_files))
  }
  colnames(fc$counts) <- sample_names

  row_data <- S4Vectors::DataFrame(fc$annotation, row.names = fc$annotation$GeneID)

  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = fc$counts),
    rowData = row_data,
    colData = S4Vectors::DataFrame(sample = sample_names, row.names = sample_names)
  )
}
