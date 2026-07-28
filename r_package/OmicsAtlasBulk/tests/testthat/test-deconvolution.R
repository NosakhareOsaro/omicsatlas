# run_bayesprism_deconvolution() wraps a real BayesPrism Gibbs-sampling run, which has
# genuine runtime (tens of seconds) even on tiny data. Rather than mock BayesPrism
# itself (this package's whole point in this module is orchestrating the real
# algorithm correctly), this builds a small but *realistic* synthetic reference and
# bulk mixture in-memory - independent, uncorrelated random counts make BayesPrism's
# outlier-filtering step reject every gene (verified empirically while prototyping
# this module), so cell types get distinct mean expression profiles, and bulk samples
# are Poisson draws around a *weighted mixture* of those same profiles with known
# ground-truth proportions. That known ground truth lets the test assert BayesPrism
# actually recovers the right answer, not just that it runs without erroring. The
# whole (real, executed) run is done once and its result checked from multiple
# angles, rather than repeated across several test_that() blocks.

make_synthetic_reference_and_bulk <- function(seed = 1) {
  withr::local_seed(seed)

  n_genes <- 40
  n_cells_per_type <- 15
  cell_types <- rep(c("TypeA", "TypeB", "TypeC"), each = n_cells_per_type)
  gene_ids <- paste0("gene", seq_len(n_genes))

  type_means <- list(
    TypeA = exp(stats::rnorm(n_genes, mean = 3, sd = 1)),
    TypeB = exp(stats::rnorm(n_genes, mean = 3, sd = 1)),
    TypeC = exp(stats::rnorm(n_genes, mean = 3, sd = 1))
  )

  ref_counts <- matrix(
    0, nrow = n_genes, ncol = length(cell_types),
    dimnames = list(gene_ids, paste0("cell", seq_along(cell_types)))
  )
  for (i in seq_along(cell_types)) {
    ref_counts[, i] <- stats::rpois(n_genes, lambda = type_means[[cell_types[i]]])
  }

  reference_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = ref_counts),
    colData = S4Vectors::DataFrame(
      singler_label = cell_types,
      leiden = cell_types,
      row.names = colnames(ref_counts)
    )
  )

  true_props <- rbind(
    sample1 = c(TypeA = 0.6, TypeB = 0.3, TypeC = 0.1),
    sample2 = c(TypeA = 0.2, TypeB = 0.6, TypeC = 0.2),
    sample3 = c(TypeA = 0.1, TypeB = 0.1, TypeC = 0.8)
  )
  bulk_counts <- matrix(
    0, nrow = n_genes, ncol = nrow(true_props),
    dimnames = list(gene_ids, rownames(true_props))
  )
  for (i in seq_len(nrow(true_props))) {
    mixed_mean <- true_props[i, "TypeA"] * type_means$TypeA +
      true_props[i, "TypeB"] * type_means$TypeB +
      true_props[i, "TypeC"] * type_means$TypeC
    bulk_counts[, i] <- stats::rpois(n_genes, lambda = mixed_mean * 50)
  }

  bulk_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = bulk_counts),
    colData = S4Vectors::DataFrame(sample = colnames(bulk_counts), row.names = colnames(bulk_counts))
  )

  list(reference_sce = reference_sce, bulk_se = bulk_se, true_props = true_props)
}

test_that("run_bayesprism_deconvolution recovers known ground-truth mixture proportions", {
  testthat::skip_on_cran()
  fixture <- make_synthetic_reference_and_bulk()

  estimated <- run_bayesprism_deconvolution(
    bulk_se = fixture$bulk_se,
    reference_sce = fixture$reference_sce,
    outlier.cut = 0.5,
    outlier.fraction = 0.5
  )

  expect_true(is.matrix(estimated))
  expect_setequal(rownames(estimated), rownames(fixture$true_props))
  expect_setequal(colnames(estimated), colnames(fixture$true_props))

  estimated <- estimated[rownames(fixture$true_props), colnames(fixture$true_props)]

  # Proportions in each row should sum to ~1.
  expect_equal(unname(rowSums(estimated)), rep(1, nrow(estimated)), tolerance = 1e-6)

  # BayesPrism should recover something close to the true composition per sample -
  # not exact (it's a probabilistic estimate), but well within a generous margin.
  max_abs_error <- max(abs(estimated - fixture$true_props))
  expect_lt(max_abs_error, 0.15)

  # Each sample's dominant true cell type should also be its dominant estimated one.
  true_dominant <- colnames(fixture$true_props)[apply(fixture$true_props, 1, which.max)]
  estimated_dominant <- colnames(estimated)[apply(estimated, 1, which.max)]
  expect_equal(estimated_dominant, true_dominant)
})

test_that("load_scrna_signature() calls zellkonverter::readH5AD() on the given path", {
  called_with <- NULL
  fake_sce <- SingleCellExperiment::SingleCellExperiment()
  testthat::local_mocked_bindings(
    readH5AD = function(path, ...) {
      called_with <<- path
      fake_sce
    },
    .package = "zellkonverter"
  )

  result <- load_scrna_signature("/some/path/brca_scrna_signature.h5ad")

  expect_identical(called_with, "/some/path/brca_scrna_signature.h5ad")
  expect_identical(result, fake_sce)
})
