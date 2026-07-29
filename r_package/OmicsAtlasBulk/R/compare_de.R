#' Compare DESeq2 and edgeR differential expression results
#'
#' Compares the significant-gene sets from [run_deseq2()] and [run_edger()] run on
#' the same data, producing an overlap summary, Venn and UpSet-style
#' visualisations, and an explicit discrepancy table — genes exactly one method
#' calls significant — rather than silently picking one method's result as "the"
#' answer. `deseq2_res` and `edger_res` must come from the same
#' `SummarizedExperiment` (same genes, same row order).
#'
#' @param deseq2_res Result of [run_deseq2()].
#' @param edger_res Result of [run_edger()], on the same data.
#'
#' @return A list with:
#' \describe{
#'   \item{overlap}{List of gene-ID character vectors: `both`, `deseq2_only`,
#'     `edger_only`.}
#'   \item{venn_plot}{A `ggplot` Venn diagram of the two significant-gene sets.}
#'   \item{upset_matrix}{A binary membership matrix (`UpSetR::fromList()` output),
#'     renderable via `UpSetR::upset()`.}
#'   \item{discrepancy}{Data frame, one row per gene significant in exactly one
#'     method, with each method's `log2FoldChange`/`logFC` and significance flag
#'     side by side.}
#' }
#'
#' @examples
#' data(example_bulk_se)
#' deseq2_res <- run_deseq2(example_bulk_se, condition_column = "condition")
#' edger_res <- run_edger(example_bulk_se, condition_column = "condition")
#' comparison <- compare_de_methods(deseq2_res, edger_res)
#' comparison$overlap
#' @export
compare_de_methods <- function(deseq2_res, edger_res) {
  if (!identical(rownames(deseq2_res), rownames(edger_res))) {
    stop(
      "`deseq2_res` and `edger_res` must have identical, identically-ordered gene ",
      "IDs (i.e. come from run_deseq2()/run_edger() on the same SummarizedExperiment)."
    )
  }

  gene_ids <- rownames(deseq2_res)
  deseq2_sig <- gene_ids[deseq2_res$is_significant]
  edger_sig <- gene_ids[edger_res$is_significant]

  both <- intersect(deseq2_sig, edger_sig)
  deseq2_only <- setdiff(deseq2_sig, edger_sig)
  edger_only <- setdiff(edger_sig, deseq2_sig)

  venn_plot <- ggVennDiagram::ggVennDiagram(
    list(DESeq2 = deseq2_sig, edgeR = edger_sig)
  ) + ggplot2::scale_fill_gradient(low = "white", high = "steelblue")

  upset_matrix <- UpSetR::fromList(list(DESeq2 = deseq2_sig, edgeR = edger_sig))

  discrepant_genes <- union(deseq2_only, edger_only)
  discrepancy <- data.frame(
    gene_id = discrepant_genes,
    significant_in = ifelse(discrepant_genes %in% deseq2_only, "DESeq2 only", "edgeR only"),
    deseq2_log2FoldChange = deseq2_res[discrepant_genes, "log2FoldChange"],
    deseq2_padj = deseq2_res[discrepant_genes, "padj"],
    edger_logFC = edger_res[discrepant_genes, "logFC"],
    edger_FDR = edger_res[discrepant_genes, "FDR"],
    row.names = NULL
  )

  list(
    overlap = list(both = both, deseq2_only = deseq2_only, edger_only = edger_only),
    venn_plot = venn_plot,
    upset_matrix = upset_matrix,
    discrepancy = discrepancy
  )
}

#' Extract the fold-change column from a DESeq2 or edgeR result
#'
#' [run_deseq2()] results have a `log2FoldChange` column; [run_edger()] results have
#' `logFC` instead. Callers that need "the fold-change column, whichever DE method
#' this came from" (e.g. building a ranked gene list for GSEA) should use this
#' rather than hardcoding one column name — a real incident (see `ADR-0006`)
#' silently produced `NULL` for the un-hardcoded method by doing exactly that.
#' Errors explicitly if neither column is present, rather than returning `NULL` for
#' some future third DE method with yet another naming convention.
#'
#' @param res A [run_deseq2()] or [run_edger()] result (or any data frame/DataFrame
#'   with a `log2FoldChange` or `logFC` column).
#'
#' @return The fold-change column, as given (whatever numeric type/class `res`
#'   itself uses for its columns).
#'
#' @examples
#' data(example_bulk_se)
#' deseq2_res <- run_deseq2(example_bulk_se, condition_column = "condition")
#' edger_res <- run_edger(example_bulk_se, condition_column = "condition")
#' extract_fold_change(deseq2_res)
#' extract_fold_change(edger_res)
#' @export
extract_fold_change <- function(res) {
  if (!is.null(res$log2FoldChange)) {
    return(res$log2FoldChange)
  }
  if (!is.null(res$logFC)) {
    return(res$logFC)
  }
  stop(
    "`res` has neither a `log2FoldChange` (DESeq2) nor `logFC` (edgeR) column - ",
    "cannot determine which fold-change column to use."
  )
}
