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
#' @param cell_state_col `colData` column in `reference_sce` giving raw cell-state
#'   labels (a finer-grained tier within each type, as BayesPrism recommends).
#'   Default `"leiden"` (the Phase 1 signature's cluster column). The actual state
#'   label passed to BayesPrism is `paste(cell_type, raw_state)`, not this column's
#'   raw values directly - BayesPrism requires states to nest within types, and
#'   `leiden` clusters (computed independently of `singler_label`) don't nest on
#'   their own (verified on the real signature: a single Leiden cluster spans
#'   nearly every SingleR type). Combining the two guarantees valid nesting by
#'   construction, for any `cell_type_col`/`cell_state_col` pair, without reshaping
#'   the underlying artifact. Combinations with fewer than `min_state_cells` cells
#'   are further pooled into a per-type `"other"` state - see `min_state_cells`.
#' @param min_state_cells Minimum number of cells a `cell_type`/`cell_state_col`
#'   combination must have to be kept as its own BayesPrism cell state; smaller
#'   combinations are pooled into a single `"<cell_type>_other"` state per type.
#'   Default `30` - chosen from the real Phase 1 signature's state-size distribution
#'   (99 raw type x state combinations; 53% had fewer than 5 cells; right-skewed,
#'   with 1-3 dominant clusters carrying nearly all cells in every type). `N = 30`
#'   collapses that down to 42 states while pooling only ~2% of cells; the range
#'   20-50 all land within a few states of each other, so 30 isn't a finely-tuned
#'   edge case (see `ADR-0004`). This only changes `cell_state_labels` granularity -
#'   `cell_type_labels`, and therefore the per-type proportions BayesPrism reports,
#'   are unaffected.
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
                                          min_state_cells = 30,
                                          n_cores = 1,
                                          ...) {
  ref_counts <- t(as.matrix(SummarizedExperiment::assay(reference_sce)))
  cell_type_labels <- as.character(SummarizedExperiment::colData(reference_sce)[[cell_type_col]])
  cell_state_raw <- as.character(SummarizedExperiment::colData(reference_sce)[[cell_state_col]])
  # BayesPrism requires cell.state.labels to nest within cell.type.labels (each
  # state belongs to exactly one type) - errors with "one or more cell states
  # belong to multiple cell types" otherwise. The documented Phase 1 schema's
  # `leiden` clusters are computed independently of `singler_label` cell types and
  # don't nest - verified empirically on the real signature artifact (a single
  # Leiden cluster spans nearly every SingleR type). Combining type and raw state
  # into one label nests by construction while still carrying leiden's within-type
  # heterogeneity signal, without reshaping the underlying artifact at all. On the
  # real signature this produces many near-singleton states unsuited to stable
  # profile estimation (and blew Gibbs sampling's estimated runtime out to 13+
  # hours on a first real run) - collapse_rare_cell_states() pools those below
  # min_state_cells into a per-type "other" state (see ADR-0004).
  cell_state_labels <- collapse_rare_cell_states(cell_type_labels, cell_state_raw, min_state_cells)

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

#' Pool rare cell-type/cell-state combinations into a per-type "other" state
#'
#' `paste(cell_type, raw_state)` guarantees BayesPrism's required type/state nesting
#' (see `run_bayesprism_deconvolution()`) but, on the real Phase 1 signature,
#' produces many near-singleton combinations unsuited to stable expression-profile
#' estimation - see `ADR-0004` for the real state-size distribution this was measured
#' against and why `min_cells` defaults to `30` there. Combinations with fewer than
#' `min_cells` cells are replaced with a single `"<cell_type>_other"` label per type;
#' combinations at or above `min_cells` pass through as `"<cell_type>_<raw_state>"`
#' unchanged. Rarity is judged from each combination's total size across all cells
#' passed in, not recomputed after pooling.
#'
#' @param cell_type_labels Character vector, one per cell.
#' @param cell_state_raw Character vector, one per cell (the raw `cell_state_col`
#'   values, before combining with type).
#' @param min_cells Minimum cell count for a `cell_type`/`cell_state_raw` combination
#'   to be kept as its own state.
#'
#' @return Character vector of final cell-state labels, one per cell.
#' @keywords internal
collapse_rare_cell_states <- function(cell_type_labels, cell_state_raw, min_cells) {
  combined <- paste(cell_type_labels, cell_state_raw, sep = "_")
  combo_sizes <- table(combined)
  # unname()/as.integer(): table's `[` keeps dim/dimnames on the lookup result, and
  # ifelse() then copies that shape onto its output - unwrapped, the returned labels
  # would silently carry array dim/dimnames instead of being a plain character vector.
  sizes <- unname(as.integer(combo_sizes[combined]))
  ifelse(sizes < min_cells, paste0(cell_type_labels, "_other"), combined)
}

#' Load the OmicsAtlas scRNA-seq signature artifact
#'
#' Thin wrapper around `zellkonverter::readH5AD()`, isolated in its own function so
#' it's the one place that knows how the Phase 1 `.h5ad` artifact gets into R. Uses
#' `reader = "R"` (zellkonverter's native HDF5 reader) rather than the default
#' Python/basilisk-backed reader: `.h5ad` is a well-defined on-disk format readable
#' without going through Python at all, and the default reader tried to bootstrap an
#' entirely new, from-source-compiled Python via basilisk/reticulate the first time
#' this was run for real - unnecessary and far slower than just reading the file.
#'
#' @param signature_h5ad_path Path to the Phase 1 scRNA-seq signature `.h5ad`
#'   artifact (see `src/omicsatlas/scrna/artifact.py`'s `signature_path()` in the
#'   parent repository for the canonical path/version contract — this function does
#'   not hardcode or guess that path, callers resolve and pass it explicitly).
#'
#' @return A `SingleCellExperiment`.
#' @export
load_scrna_signature <- function(signature_h5ad_path) {
  zellkonverter::readH5AD(signature_h5ad_path, reader = "R")
}
