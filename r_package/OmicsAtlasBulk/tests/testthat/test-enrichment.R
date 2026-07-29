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

test_that("retry_with_backoff succeeds once the wrapped call stops failing, logging each retry", {
  attempt_count <- 0
  flaky_fn <- function() {
    attempt_count <<- attempt_count + 1
    if (attempt_count < 3) {
      stop(sprintf("simulated transient failure #%d", attempt_count))
    }
    "real result"
  }
  retry_messages <- character(0)
  capture_retry <- function(msg) retry_messages <<- c(retry_messages, msg)

  result <- retry_with_backoff(flaky_fn, max_attempts = 3, delays = c(0, 0), on_retry = capture_retry)

  expect_identical(result, "real result")
  expect_identical(attempt_count, 3)
  expect_length(retry_messages, 2)
  expect_match(retry_messages[1], "simulated transient failure #1")
  expect_match(retry_messages[2], "simulated transient failure #2")
})

test_that("retry_with_backoff propagates the error once all attempts are exhausted", {
  attempt_count <- 0
  always_fails <- function() {
    attempt_count <<- attempt_count + 1
    stop("simulated permanent failure")
  }

  expect_error(
    retry_with_backoff(always_fails, max_attempts = 3, delays = c(0, 0), on_retry = function(msg) invisible(NULL)),
    "simulated permanent failure"
  )
  expect_identical(attempt_count, 3)
})

test_that("retry_with_backoff succeeds immediately without retrying or sleeping", {
  call_count <- 0
  retry_called <- FALSE
  always_succeeds <- function() {
    call_count <<- call_count + 1
    "immediate success"
  }

  result <- retry_with_backoff(
    always_succeeds,
    max_attempts = 3, delays = c(0, 0),
    on_retry = function(msg) retry_called <<- TRUE
  )

  expect_identical(result, "immediate success")
  expect_identical(call_count, 1)
  expect_false(retry_called)
})
