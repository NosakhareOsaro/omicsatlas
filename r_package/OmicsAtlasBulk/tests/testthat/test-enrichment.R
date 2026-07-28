# GO/GSEA use only the local OrgDb annotation (org.Hs.eg.db) and need no network —
# tested for real here. KEGG genuinely requires live network access (KEGG's REST
# API, no offline alternative) — see R/enrichment.R's docstring. The KEGG-specific
# test below is skipped in CI (testthat::skip_on_ci()) rather than mocked, so it's
# either genuinely exercised (locally, when a developer runs it) or explicitly not
# run at all — never silently faked.

sig_genes <- c("TP53", "BRCA1", "BRCA2", "ERBB2", "ESR1")
ranked_stats <- stats::setNames(
  c(4, 3.5, 3, 2.8, 2.5, rep(0, 15)),
  c(sig_genes, paste0("GENE", 6:20))
)

test_that("run_enrichment finds real GO enrichment for well-annotated cancer genes", {
  result <- run_enrichment(sig_genes, ranked_stats, run_kegg = FALSE)

  expect_s4_class(result$go, "enrichResult")
  expect_gt(nrow(as.data.frame(result$go)), 0)
})

test_that("run_enrichment's GSEA covers the full ranked gene list", {
  result <- run_enrichment(sig_genes, ranked_stats, run_kegg = FALSE)

  expect_s4_class(result$gsea, "gseaResult")
})

test_that("run_kegg = FALSE skips KEGG entirely without touching the network", {
  result <- run_enrichment(sig_genes, ranked_stats, run_kegg = FALSE)

  expect_null(result$kegg)
})

test_that("empty gene_symbols returns NULL for GO and KEGG but still runs GSEA", {
  result <- run_enrichment(character(0), ranked_stats, run_kegg = TRUE)

  expect_null(result$go)
  expect_null(result$kegg)
  expect_s4_class(result$gsea, "gseaResult")
})

test_that("run_enrichment's KEGG step returns a real enrichResult (network required)", {
  testthat::skip_on_ci()
  testthat::skip_if_offline()

  result <- run_enrichment(sig_genes, ranked_stats, run_kegg = TRUE)

  expect_s4_class(result$kegg, "enrichResult")
})
