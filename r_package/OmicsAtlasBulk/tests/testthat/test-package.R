test_that("package namespace loads", {
  expect_true(isNamespaceLoaded("OmicsAtlasBulk") || requireNamespace("OmicsAtlasBulk", quietly = TRUE))
})
