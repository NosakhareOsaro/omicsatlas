# Never assumes STAR is actually installed — see R/align.R and ADR-0004. The
# "real execution" tests below use a tiny stub shell script standing in for the STAR
# binary, generated at test time (not bundled, to avoid any git-checkout
# executable-bit concerns), so the system2() call, exit-status handling, and output
# path contract are still genuinely exercised, not just the dry-run arg construction.

r1 <- system.file("extdata", "toy_reads_R1.fastq.gz", package = "OmicsAtlasBulk")
r2 <- system.file("extdata", "toy_reads_R2.fastq.gz", package = "OmicsAtlasBulk")

make_stub_star <- function(dir, exit_code = 0L, touch_output = TRUE) {
  stub_path <- file.path(dir, "STAR_stub.sh")
  lines <- c(
    "#!/bin/sh",
    "# Minimal stand-in for STAR: parses --outFileNamePrefix and, on success,",
    "# touches the BAM STAR would have produced there.",
    "prefix=\"\"",
    "while [ \"$#\" -gt 0 ]; do",
    "  if [ \"$1\" = \"--outFileNamePrefix\" ]; then prefix=\"$2\"; fi",
    "  shift",
    "done",
    if (touch_output) "touch \"${prefix}Aligned.sortedByCoord.out.bam\"" else "",
    paste0("exit ", exit_code)
  )
  writeLines(lines, stub_path)
  Sys.chmod(stub_path, mode = "0755")
  stub_path
}

test_that("dry_run constructs single-end command without running anything", {
  args <- run_star_align(
    fastq_files = r1,
    genome_dir = "/fake/genome",
    output_prefix = "/fake/out_",
    dry_run = TRUE
  )

  expect_true("--genomeDir" %in% args)
  expect_true("/fake/genome" %in% args)
  expect_true(r1 %in% args)
  expect_true("--readFilesCommand" %in% args) # .gz input
  expect_true("zcat" %in% args)
})

test_that("dry_run constructs paired-end command with both FASTQ files", {
  args <- run_star_align(
    fastq_files = c(r1, r2),
    genome_dir = "/fake/genome",
    output_prefix = "/fake/out_",
    dry_run = TRUE
  )

  expect_true(r1 %in% args)
  expect_true(r2 %in% args)
})

test_that("invalid fastq_files length errors clearly", {
  expect_error(
    run_star_align(character(0), "/fake/genome", "/fake/out_", dry_run = TRUE),
    "length 1"
  )
  expect_error(
    run_star_align(c("a", "b", "c"), "/fake/genome", "/fake/out_", dry_run = TRUE),
    "length 1"
  )
})

test_that("a successful stub STAR run returns the expected BAM path and it exists", {
  tmp <- withr::local_tempdir()
  stub <- make_stub_star(tmp, exit_code = 0L)
  output_prefix <- file.path(tmp, "sample1_")

  bam_path <- run_star_align(
    fastq_files = r1,
    genome_dir = tmp,
    output_prefix = output_prefix,
    star_path = stub
  )

  expect_equal(bam_path, paste0(output_prefix, "Aligned.sortedByCoord.out.bam"))
  expect_true(file.exists(bam_path))
})

test_that("a non-zero STAR exit status raises an informative error", {
  tmp <- withr::local_tempdir()
  stub <- make_stub_star(tmp, exit_code = 1L, touch_output = FALSE)

  expect_error(
    run_star_align(
      fastq_files = r1,
      genome_dir = tmp,
      output_prefix = file.path(tmp, "sample1_"),
      star_path = stub
    ),
    "non-zero status"
  )
})
