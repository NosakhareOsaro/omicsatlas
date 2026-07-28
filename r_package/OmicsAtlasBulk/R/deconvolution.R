#' Run BayesPrism deconvolution against the OmicsAtlas scRNA-seq signature
#'
#' Deconvolves bulk RNA-seq samples into cell-type proportions using BayesPrism
#' (installed from GitHub, pinned to a specific commit — see `ADR-0004`), with the
#' reference built from the exact, versioned Phase 1 scRNA-seq signature artifact
#' (see the parent repository's `DATA_SCHEMA.md` and
#' `src/omicsatlas/scrna/artifact.py`'s `signature_path()` contract) — never
#' recomputed or approximated here.
#'
#' `reference_sce` takes an already-loaded `SingleCellExperiment` (typically via
#' `zellkonverter::readH5AD(signature_h5ad_path)`) rather than a file path, so this
#' function stays testable/portable — callers control exactly how/when the
#' (potentially large) real artifact is loaded, and tests can substitute a small
#' synthetic reference without needing a real `.h5ad` file at all. The pipeline
#' entrypoint (see the parent repository) is what resolves the real path via Python's
#' `signature_path()` contract and loads it before calling this.
#'
#' @param bulk_se `SummarizedExperiment` with a raw-count assay (bulk RNA-seq
#'   samples).
#' @param reference_sce A `SingleCellExperiment` loaded from the Phase 1 scRNA-seq
#'   signature artifact. Its default assay is used as the reference counts (the
#'   artifact's `.X`, i.e. QC-passed raw counts — *not* `scran_normalized`; BayesPrism
#'   wants raw counts).
#' @param assay_name Which assay in `bulk_se` holds raw counts. Default `"counts"`.
#' @param cell_type_col `colData` column in `reference_sce` giving cell-type labels.
#'   Default `"singler_label"` (the Phase 1 signature's SingleR annotation column —
#'   see `DATA_SCHEMA.md`).
#' @param cell_state_col `colData` column in `reference_sce` giving cell-state labels
#'   (a finer-grained tier within each type, as BayesPrism recommends). Default
#'   `"leiden"` (the Phase 1 signature's cluster column).
#' @param n_cores Passed to `BayesPrism::run.prism()`. Default `1`.
#' @param ... Additional arguments passed to `BayesPrism::new.prism()` (e.g.
#'   `outlier.cut`, `outlier.fraction`).
#'
#' @return A samples-by-cell-type matrix of estimated proportions
#'   (`BayesPrism::get.fraction(..., state.or.type = "type")`).
#'
#' No `key` (tumour-reference) is passed to `BayesPrism::new.prism()` — the Phase 1
#' signature's SingleR reference has no malignant/normal-epithelial split to key on
#' (a documented limitation carried from `ADR-0002`/`DATA_SCHEMA.md`), so this runs
#' BayesPrism's plain cell-type/cell-state mode rather than inventing a malignant key.
#'
#' @examples
#' \donttest{
#' # BayesPrism's Gibbs sampling has real runtime (tens of seconds even on tiny
#' # data), so this is marked \donttest rather than run on every R CMD check. See
#' # tests/testthat/test-deconvolution.R for a real, executed run against a small
#' # synthetic reference with known ground-truth mixture proportions.
#' data(example_bulk_se)
#' # reference_sce would normally come from zellkonverter::readH5AD() on the real
#' # Phase 1 artifact; omitted here since building one from scratch needs several
#' # lines - see the test file for a complete worked construction.
#' }
#' @export
run_bayesprism_deconvolution <- function(bulk_se,
                                          reference_sce,
                                          assay_name = "counts",
                                          cell_type_col = "singler_label",
                                          cell_state_col = "leiden",
                                          n_cores = 1,
                                          ...) {
  ref_counts <- t(as.matrix(SummarizedExperiment::assay(reference_sce)))
  cell_type_labels <- as.character(SummarizedExperiment::colData(reference_sce)[[cell_type_col]])
  cell_state_labels <- as.character(SummarizedExperiment::colData(reference_sce)[[cell_state_col]])

  bulk_counts <- t(round(as.matrix(SummarizedExperiment::assay(bulk_se, assay_name))))

  prism <- BayesPrism::new.prism(
    reference = ref_counts,
    input.type = "count.matrix",
    cell.type.labels = cell_type_labels,
    cell.state.labels = cell_state_labels,
    key = NULL,
    mixture = bulk_counts,
    ...
  )
  result <- BayesPrism::run.prism(prism = prism, n.cores = n_cores)
  BayesPrism::get.fraction(bp = result, which.theta = "final", state.or.type = "type")
}

#' Load the OmicsAtlas scRNA-seq signature artifact
#'
#' Thin wrapper around `zellkonverter::readH5AD()`, isolated in its own function so
#' it's the one place that knows how the Phase 1 `.h5ad` artifact gets into R.
#'
#' @param signature_h5ad_path Path to the Phase 1 scRNA-seq signature `.h5ad`
#'   artifact (see `src/omicsatlas/scrna/artifact.py`'s `signature_path()` in the
#'   parent repository for the canonical path/version contract — this function does
#'   not hardcode or guess that path, callers resolve and pass it explicitly).
#'
#' @return A `SingleCellExperiment`.
#' @export
load_scrna_signature <- function(signature_h5ad_path) {
  zellkonverter::readH5AD(signature_h5ad_path)
}
