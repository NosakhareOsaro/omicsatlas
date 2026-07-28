# Uses the package's bundled toy_aligned.bam / toy_annotation.gtf (see
# data-raw/build_quantify_fixtures.R) — a real BAM from a real Rsubread::align() run
# against a tiny synthetic reference, built once and bundled rather than
# re-aligned at test time (alignment has a ~3-minute fixed per-call overhead in
# Rsubread regardless of input size). featureCounts() itself runs directly on it
# here, genuinely, not stubbed.

bam <- system.file("extdata", "toy_aligned.bam", package = "OmicsAtlasBulk")
gtf <- system.file("extdata", "toy_annotation.gtf", package = "OmicsAtlasBulk")

test_that("bundled fixtures exist", {
  expect_true(file.exists(bam))
  expect_true(file.exists(gtf))
})

test_that("run_featurecounts returns a SummarizedExperiment with the expected counts", {
  se <- run_featurecounts(bam, gtf, sample_names = "toy_sample", isPairedEnd = FALSE)

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(colnames(se), "toy_sample")
  expect_equal(rownames(se), c("TOYGENE1", "TOYGENE2", "TOYGENE3"))
  # Each of the 3 synthetic reads was constructed to fall entirely inside exactly
  # one gene's region — a real, verifiable ground truth, not an arbitrary count.
  counts <- SummarizedExperiment::assay(se, "counts")
  expect_equal(as.integer(counts[, "toy_sample"]), c(1L, 1L, 1L))
})

test_that("run_featurecounts defaults sample_names to the BAM basename", {
  se <- run_featurecounts(bam, gtf, isPairedEnd = FALSE)

  expect_equal(colnames(se), "toy_aligned")
})

test_that("run_featurecounts carries featureCounts' gene annotation into rowData", {
  se <- run_featurecounts(bam, gtf, sample_names = "toy_sample", isPairedEnd = FALSE)

  row_data <- SummarizedExperiment::rowData(se)
  expect_true("Length" %in% colnames(row_data))
  expect_equal(row_data["TOYGENE1", "Length"], 100)
})
