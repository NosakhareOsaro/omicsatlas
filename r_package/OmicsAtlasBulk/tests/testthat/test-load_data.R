test_that("read_matched_bulk_counts reads a GSE176078-style matrix", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  writeLines(
    c(
      "Genes\tCID0001\tCID0002",
      "DDX11L1\t0\t32.39",
      "WASH7P\t780.93\t32.36"
    ),
    tmp
  )

  se <- read_matched_bulk_counts(tmp)

  expect_equal(dim(se), c(2, 2))
  expect_identical(colnames(se), c("CID0001", "CID0002"))
  expect_identical(rownames(se), c("DDX11L1", "WASH7P"))
  expect_identical(SummarizedExperiment::rowData(se)$gene_symbol, c("DDX11L1", "WASH7P"))
  expect_equal(SummarizedExperiment::assay(se, "counts")["WASH7P", "CID0001"], 780.93)
})

test_that("read_tcga_brca_counts reads the bundled fixture manifest and TSVs", {
  fixture_dir <- system.file("extdata", "tcga_fixture", package = "OmicsAtlasBulk")
  manifest_path <- file.path(fixture_dir, "manifest.json")

  se <- read_tcga_brca_counts(manifest_path, repo_root = fixture_dir)

  # N_unmapped/N_multimapping/N_noFeature/N_ambiguous summary rows are dropped.
  expect_equal(dim(se), c(5, 2))
  expect_identical(colnames(se), c("fixture-file-1", "fixture-file-2"))
  # Ensembl IDs are stripped of their version suffix.
  expect_identical(
    rownames(se),
    c(
      "ENSG00000000001", "ENSG00000000002", "ENSG00000000003",
      "ENSG00000000004", "ENSG00000000005"
    )
  )
  # gene_symbol carries gene_name, duplicates included (GENEB appears twice).
  expect_identical(
    SummarizedExperiment::rowData(se)$gene_symbol,
    c("GENEA", "GENEB", "GENEB", "GENEC", "GENED")
  )
  expect_equal(SummarizedExperiment::assay(se, "counts")["ENSG00000000005", "fixture-file-1"], 100)
  expect_equal(SummarizedExperiment::assay(se, "counts")["ENSG00000000004", "fixture-file-2"], 1)

  col_data <- SummarizedExperiment::colData(se)
  expect_identical(as.character(col_data$case_submitter_id), c("TCGA-AA-0001", "TCGA-AA-0002"))
  expect_identical(as.character(col_data$sample_type), c("Primary Tumor", "Solid Tissue Normal"))
})

test_that("read_tcga_brca_counts errors if gene order differs between files", {
  fixture_dir <- withr::local_tempdir()
  writeLines(
    c(
      "# gene-model: GENCODE v36",
      "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded",
      "ENSG00000000001.1\tGENEA\tprotein_coding\t10\t5\t5\t1.0\t1.0\t1.0",
      "ENSG00000000002.1\tGENEB\tprotein_coding\t20\t10\t10\t2.0\t2.0\t2.0"
    ),
    file.path(fixture_dir, "a.tsv")
  )
  writeLines(
    c(
      "# gene-model: GENCODE v36",
      "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded",
      "ENSG00000000002.1\tGENEB\tprotein_coding\t20\t10\t10\t2.0\t2.0\t2.0",
      "ENSG00000000001.1\tGENEA\tprotein_coding\t10\t5\t5\t1.0\t1.0\t1.0"
    ),
    file.path(fixture_dir, "b.tsv")
  )
  manifest_path <- file.path(fixture_dir, "manifest.json")
  jsonlite::write_json(
    list(files = list(
      list(file_id = "a", case_submitter_id = "X", sample_type = "Primary Tumor", dest = "a.tsv"),
      list(file_id = "b", case_submitter_id = "Y", sample_type = "Primary Tumor", dest = "b.tsv")
    )),
    manifest_path,
    auto_unbox = TRUE
  )

  expect_error(
    read_tcga_brca_counts(manifest_path, repo_root = fixture_dir),
    "Gene ID order differs"
  )
})

test_that("filter_low_count_genes drops genes below the total-count threshold", {
  data(example_bulk_se, package = "OmicsAtlasBulk", envir = environment())

  totals <- rowSums(SummarizedExperiment::assay(example_bulk_se, "counts"))
  threshold <- stats::median(totals)

  filtered <- filter_low_count_genes(example_bulk_se, min_total_count = threshold)

  expect_true(nrow(filtered) < nrow(example_bulk_se))
  expect_true(all(rowSums(SummarizedExperiment::assay(filtered, "counts")) >= threshold))
})
