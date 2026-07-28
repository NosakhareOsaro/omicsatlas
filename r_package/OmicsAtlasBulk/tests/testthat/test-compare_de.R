# compare_de_methods()'s overlap/discrepancy logic is tested two ways: against real
# run_deseq2()/run_edger() output on the bundled fixture (an integration check — on
# this fixture the two methods happen to agree perfectly, which is itself a real,
# validly-tested outcome, not a skipped case), and against hand-constructed
# data.frames with deliberate disagreement (a unit check that the discrepancy table
# is actually correct when methods *do* disagree, independent of whether DESeq2/edgeR
# happen to agree on any particular real dataset).

test_that("compare_de_methods works end-to-end on real run_deseq2()/run_edger() output", {
  data(example_bulk_se, package = "OmicsAtlasBulk", envir = environment())
  deseq2_res <- run_deseq2(example_bulk_se, condition_column = "condition")
  edger_res <- run_edger(example_bulk_se, condition_column = "condition")

  comparison <- compare_de_methods(deseq2_res, edger_res)

  expect_type(comparison$overlap, "list")
  expect_s3_class(comparison$venn_plot, "ggplot")
  expect_equal(ncol(comparison$upset_matrix), 2)
  expect_s3_class(comparison$discrepancy, "data.frame")
  # The 5 genes with a real injected effect (see data-raw/example_bulk_se.R) should
  # be found by both methods.
  expect_true(all(paste0("ENSG0000000000", 1:5, "EX") %in% comparison$overlap$both))
})

test_that("compare_de_methods requires matching gene sets", {
  data(example_bulk_se, package = "OmicsAtlasBulk", envir = environment())
  deseq2_res <- run_deseq2(example_bulk_se, condition_column = "condition")
  mismatched <- deseq2_res[-1, ]

  expect_error(compare_de_methods(deseq2_res, mismatched), "identical")
})

test_that("compare_de_methods correctly separates both/deseq2_only/edger_only", {
  genes <- paste0("GENE", 1:6)
  deseq2_res <- data.frame(
    log2FoldChange = c(2, 2, 2, 0, 0, -2),
    padj = c(0.01, 0.01, 0.01, 0.5, 0.5, 0.01),
    is_significant = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE),
    row.names = genes
  )
  edger_res <- data.frame(
    logFC = c(2, 2, 0, 0, -2, -2),
    FDR = c(0.01, 0.01, 0.5, 0.5, 0.01, 0.5),
    is_significant = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE),
    row.names = genes
  )
  # GENE1, GENE2: significant in both. GENE3: DESeq2 only. GENE5: edgeR only.
  # GENE4: neither. GENE6: DESeq2 only.

  comparison <- compare_de_methods(deseq2_res, edger_res)

  expect_setequal(comparison$overlap$both, c("GENE1", "GENE2"))
  expect_setequal(comparison$overlap$deseq2_only, c("GENE3", "GENE6"))
  expect_setequal(comparison$overlap$edger_only, "GENE5")

  expect_equal(nrow(comparison$discrepancy), 3)
  expect_setequal(comparison$discrepancy$gene_id, c("GENE3", "GENE5", "GENE6"))
  gene3_row <- comparison$discrepancy[comparison$discrepancy$gene_id == "GENE3", ]
  expect_equal(gene3_row$significant_in, "DESeq2 only")
  expect_equal(gene3_row$deseq2_log2FoldChange, 2)
  expect_equal(gene3_row$edger_logFC, 0)
})

test_that("compare_de_methods handles zero discrepant genes cleanly", {
  genes <- paste0("GENE", 1:3)
  same_sig <- data.frame(is_significant = c(TRUE, FALSE, TRUE), row.names = genes)
  deseq2_res <- cbind(same_sig, log2FoldChange = c(2, 0, -2), padj = c(0.01, 0.5, 0.01))
  edger_res <- cbind(same_sig, logFC = c(2, 0, -2), FDR = c(0.01, 0.5, 0.01))

  comparison <- compare_de_methods(deseq2_res, edger_res)

  expect_equal(nrow(comparison$discrepancy), 0)
  expect_setequal(comparison$overlap$both, c("GENE1", "GENE3"))
})
