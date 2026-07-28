#' Run STAR alignment
#'
#' Wraps a real STAR command-line invocation for RNA-seq read alignment against a
#' pre-built genome index. This is a genuine wrapper around the external `STAR`
#' binary, not a reimplementation — see `ADR-0004` in the parent repository for why
#' STAR stays external-binary-wrapped (rather than, say, `Rsubread::align()`), and
#' for why this function is only tested against a tiny synthetic FASTQ + reference
#' fixture: neither of this project's two approved bulk RNA-seq datasets
#' (`ADR-0003`) has practically available raw FASTQs to align for real.
#'
#' @param fastq_files Character vector of 1 (single-end) or 2 (paired-end) FASTQ
#'   file paths.
#' @param genome_dir Path to a STAR genome index directory (built with
#'   `STAR --runMode genomeGenerate`).
#' @param output_prefix Path prefix for STAR's output files (STAR appends its own
#'   suffixes, e.g. `Aligned.sortedByCoord.out.bam`).
#' @param star_path Path to the STAR executable. Defaults to `"STAR"`, resolved via
#'   `PATH`.
#' @param extra_args Character vector of additional STAR command-line arguments,
#'   appended after the standard ones.
#' @param dry_run If `TRUE`, construct and return the command-line arguments without
#'   running STAR — used by tests and by the example below, which don't assume STAR
#'   is installed.
#'
#' @return If `dry_run = TRUE`, the constructed character vector of STAR arguments
#'   (invisibly). Otherwise, the path to the aligned, coordinate-sorted BAM file STAR
#'   produces (`<output_prefix>Aligned.sortedByCoord.out.bam`), invisibly; errors if
#'   STAR exits with a non-zero status.
#'
#' @examples
#' # Construct (but do not run) a STAR command line, using the package's bundled toy
#' # paired-end FASTQ fixture:
#' r1 <- system.file("extdata", "toy_reads_R1.fastq.gz", package = "OmicsAtlasBulk")
#' r2 <- system.file("extdata", "toy_reads_R2.fastq.gz", package = "OmicsAtlasBulk")
#' run_star_align(
#'   fastq_files = c(r1, r2),
#'   genome_dir = tempdir(),
#'   output_prefix = file.path(tempdir(), "toy_"),
#'   dry_run = TRUE
#' )
#' @export
run_star_align <- function(fastq_files,
                            genome_dir,
                            output_prefix,
                            star_path = "STAR",
                            extra_args = character(),
                            dry_run = FALSE) {
  if (!length(fastq_files) %in% c(1L, 2L)) {
    stop("`fastq_files` must have length 1 (single-end) or 2 (paired-end), got ", length(fastq_files))
  }

  fastq_arg <- paste(fastq_files, collapse = ",")
  read_files_command <- if (any(grepl("\\.gz$", fastq_files))) c("--readFilesCommand", "zcat") else character()

  args <- c(
    "--runMode", "alignReads",
    "--genomeDir", genome_dir,
    "--readFilesIn", fastq_files,
    read_files_command,
    "--outFileNamePrefix", output_prefix,
    "--outSAMtype", "BAM", "SortedByCoordinate",
    extra_args
  )

  if (dry_run) {
    return(invisible(args))
  }

  status <- system2(star_path, args = args)
  if (status != 0) {
    stop("STAR exited with non-zero status ", status, " (command: ", star_path, " ", paste(args, collapse = " "), ")")
  }

  invisible(paste0(output_prefix, "Aligned.sortedByCoord.out.bam"))
}
