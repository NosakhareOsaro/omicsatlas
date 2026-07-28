#' Load the GSE176078-matched bulk RNA-seq count matrix
#'
#' Reads the raw, patient-matched bulk RNA-seq count matrix fetched by
#' `scripts/fetch_data.py` (`gse176078_bulk_matched` - see `ADR-0003`) into a
#' `SummarizedExperiment`. Values are RSEM-style non-integer estimated counts (not
#' literal raw counts, despite the source filename) - [run_deseq2()]/[run_edger()]
#' already round them.
#'
#' @param counts_path Path to the raw
#'   `GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt` file: tab-separated, a
#'   `Genes` column of gene symbols, then one column per patient ID (e.g.
#'   `CID3586`) - the same patient IDs as the Phase 1 scRNA-seq signature's
#'   `obs['orig.ident']`.
#'
#' @return A `SummarizedExperiment` with a `counts` assay (genes x samples) and
#'   `rowData$gene_symbol`. No `colData` beyond sample names is set here - this
#'   cohort has no single, sample-level tumour/normal split to build a DESeq2/edgeR
#'   design from (every sample is tumour tissue); its role in this project is
#'   feeding [run_bayesprism_deconvolution()], not a two-group DE comparison. See
#'   `docs/bulk-pipeline.md`.
#'
#' @examples
#' \donttest{
#' # Requires the real fetched file (scripts/fetch_data.py gse176078_bulk_matched).
#' se <- read_matched_bulk_counts("data/raw/bulk/gse176078_matched/GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt")
#' }
#' @export
read_matched_bulk_counts <- function(counts_path) {
  raw <- utils::read.delim(counts_path, check.names = FALSE, row.names = 1)
  counts_mat <- as.matrix(raw)
  storage.mode(counts_mat) <- "double"

  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts_mat),
    rowData = S4Vectors::DataFrame(gene_symbol = rownames(counts_mat), row.names = rownames(counts_mat))
  )
}

#' Load the TCGA-BRCA STAR-Counts subset into a SummarizedExperiment
#'
#' Reads the per-sample GDC STAR-Counts TSVs referenced by the fetch manifest
#' `scripts/fetch_data.py` wrote (`tcga_brca_subset` - see `ADR-0003`), using each
#' file's `unstranded` column as the raw count (GDC's recommended default when
#' library strandedness isn't otherwise pinned down). Gene rows are keyed by
#' Ensembl gene ID with its version suffix stripped (`gene_name`/symbol is not
#' unique across all ~60k GENCODE features in this file - 110 duplicated symbols in
#' the fetched subset, checked directly); `rowData$gene_symbol` carries `gene_name`
#' for [run_enrichment()], duplicates and all.
#'
#' @param manifest_path Path to the fetch manifest JSON, e.g.
#'   `data/.provenance/tcga_brca_subset.json`.
#' @param repo_root Root directory the manifest's `dest` paths (e.g.
#'   `data/raw/bulk/tcga_brca/...`) are relative to. Default: two directories up
#'   from `manifest_path` (right for the manifest's real
#'   `data/.provenance/tcga_brca_subset.json` location - the repository root -
#'   but overridable, e.g. for tests using a fixture manifest elsewhere).
#'
#' @return A `SummarizedExperiment` with a `counts` assay (genes x samples),
#'   `rowData$gene_symbol`, and `colData` columns `case_submitter_id` and
#'   `sample_type` (`"Primary Tumor"`/`"Solid Tissue Normal"` - GDC's own values,
#'   suitable directly as a [run_deseq2()]/[run_edger()] `condition_column`).
#'
#' @examples
#' \donttest{
#' # Requires the real fetched files (scripts/fetch_data.py tcga_brca_subset).
#' se <- read_tcga_brca_counts("data/.provenance/tcga_brca_subset.json")
#' }
#' @export
read_tcga_brca_counts <- function(manifest_path, repo_root = dirname(dirname(manifest_path))) {
  manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  files_info <- manifest$files

  read_one <- function(file_info) {
    path <- file.path(repo_root, file_info$dest)
    df <- utils::read.delim(path, skip = 1, check.names = FALSE)
    df <- df[!grepl("^N_", df$gene_id), ]
    list(
      gene_id = sub("\\..*$", "", df$gene_id),
      gene_name = df$gene_name,
      counts = df$unstranded
    )
  }

  first <- read_one(files_info[[1]])
  gene_ids <- first$gene_id

  counts_mat <- matrix(
    0, nrow = length(gene_ids), ncol = length(files_info),
    dimnames = list(gene_ids, vapply(files_info, `[[`, character(1), "file_id"))
  )
  counts_mat[, 1] <- first$counts
  for (i in seq_along(files_info)[-1]) {
    one <- read_one(files_info[[i]])
    if (!identical(one$gene_id, gene_ids)) {
      stop(
        "Gene ID order differs between TCGA-BRCA files (file ", i, "); cannot ",
        "assemble a shared count matrix."
      )
    }
    counts_mat[, i] <- one$counts
  }

  col_data <- S4Vectors::DataFrame(
    case_submitter_id = vapply(files_info, `[[`, character(1), "case_submitter_id"),
    sample_type = vapply(files_info, `[[`, character(1), "sample_type"),
    row.names = colnames(counts_mat)
  )
  row_data <- S4Vectors::DataFrame(gene_symbol = first$gene_name, row.names = gene_ids)

  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts_mat), rowData = row_data, colData = col_data
  )
}

#' Filter out genes with very low total counts
#'
#' A standard pre-filtering step recommended before DESeq2/edgeR (both packages'
#' own vignettes suggest it): drops genes whose counts summed across all samples
#' fall below `min_total_count`, before either method's own internal
#' independent-filtering runs. Mainly useful for the TCGA-BRCA subset, whose
#' ~60,660-feature GENCODE annotation includes many pseudogenes/lncRNAs with
#' near-zero counts across the whole cohort.
#'
#' @param se A `SummarizedExperiment` with a raw-count assay.
#' @param min_total_count Minimum summed count (across all samples) for a gene to
#'   be kept. Default `10`.
#' @param assay_name Which assay in `se` holds counts. Default `"counts"`.
#'
#' @return `se`, subset to the genes passing the filter.
#'
#' @examples
#' data(example_bulk_se)
#' filter_low_count_genes(example_bulk_se, min_total_count = 20)
#' @export
filter_low_count_genes <- function(se, min_total_count = 10, assay_name = "counts") {
  totals <- rowSums(SummarizedExperiment::assay(se, assay_name))
  se[totals >= min_total_count, ]
}
