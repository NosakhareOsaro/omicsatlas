#' Run GO, KEGG, and GSEA functional enrichment
#'
#' Runs over-representation analysis (GO, KEGG) on a significant-gene set and gene
#' set enrichment analysis (GSEA, on GO terms) on a full ranked gene list, via
#' `clusterProfiler`. Intended to run once per DE method's significant-gene set (see
#' [run_deseq2()]/[run_edger()]/[compare_de_methods()]) so results from both methods
#' can be reported side by side rather than only enriching one method's output.
#'
#' GO enrichment (`enrichGO`/`gseGO`) uses the local `OrgDb` annotation and needs no
#' network access. **KEGG enrichment (`enrichKEGG`) calls KEGG's REST API and
#' genuinely requires live network access** — clusterProfiler has no offline KEGG
#' data source. This package's own test suite exercises the GO/GSEA paths for real
#' but skips the KEGG path in CI (`testthat::skip_on_ci()`) rather than either
#' silently mocking it or letting a CI run depend on an external service; KEGG
#' enrichment is validated for real only in the actual pipeline run against the
#' project's real bulk RNA-seq data (see `docs/bulk-pipeline.md`).
#'
#' @param gene_symbols Character vector of significant gene symbols (used for GO/KEGG
#'   over-representation).
#' @param ranked_stats Named numeric vector, gene symbol -> ranking statistic
#'   (typically `log2FoldChange`/`logFC`), covering every *tested* gene (not just the
#'   significant ones) — used for GSEA.
#' @param ont GO ontology: `"BP"`, `"MF"`, `"CC"`, or `"ALL"`. Default `"BP"`.
#' @param organism_db An `OrgDb` annotation object. Default
#'   `org.Hs.eg.db::org.Hs.eg.db`.
#' @param pvalue_cutoff p-value cutoff for the GO/KEGG over-representation tests.
#'   Default `0.05`. GSEA is run with no cutoff (`1`) so callers can inspect/filter
#'   the full ranked result themselves.
#' @param run_kegg If `FALSE`, skip the network-dependent KEGG step entirely and
#'   return `kegg = NULL`. Default `TRUE`.
#'
#' @return A list with `go`, `kegg`, and `gsea` elements (each a clusterProfiler
#'   `enrichResult`/`gseaResult`, or `NULL` if that step had no input or was
#'   skipped).
#'
#' @examples
#' data(example_bulk_se)
#' res <- run_deseq2(example_bulk_se, condition_column = "condition")
#' sig_genes <- SummarizedExperiment::rowData(example_bulk_se)$gene_symbol[res$is_significant]
#' ranked <- stats::setNames(
#'   res$log2FoldChange,
#'   SummarizedExperiment::rowData(example_bulk_se)$gene_symbol
#' )
#' # KEGG needs live network access; skipped here so this example runs offline.
#' enrichment <- run_enrichment(sig_genes, ranked, run_kegg = FALSE)
#' enrichment$go
#' @export
run_enrichment <- function(gene_symbols,
                            ranked_stats,
                            ont = "BP",
                            organism_db = org.Hs.eg.db::org.Hs.eg.db,
                            pvalue_cutoff = 0.05,
                            run_kegg = TRUE) {
  go_result <- if (length(gene_symbols) > 0) {
    clusterProfiler::enrichGO(
      gene = gene_symbols,
      OrgDb = organism_db,
      keyType = "SYMBOL",
      ont = ont,
      pvalueCutoff = pvalue_cutoff
    )
  } else {
    NULL
  }

  kegg_result <- NULL
  if (run_kegg && length(gene_symbols) > 0) {
    entrez_map <- clusterProfiler::bitr(
      gene_symbols,
      fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism_db
    )
    if (nrow(entrez_map) > 0) {
      kegg_result <- clusterProfiler::enrichKEGG(
        gene = entrez_map$ENTREZID, organism = "hsa", pvalueCutoff = pvalue_cutoff
      )
    }
  }

  ranked_sorted <- sort(ranked_stats, decreasing = TRUE)
  gsea_result <- clusterProfiler::gseGO(
    geneList = ranked_sorted,
    OrgDb = organism_db,
    keyType = "SYMBOL",
    ont = ont,
    pvalueCutoff = 1,
    verbose = FALSE
  )

  list(go = go_result, kegg = kegg_result, gsea = gsea_result)
}
