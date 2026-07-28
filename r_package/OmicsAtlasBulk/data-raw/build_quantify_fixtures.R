# Builds inst/extdata/toy_aligned.bam + inst/extdata/toy_annotation.gtf: a real,
# small BAM produced by actually running Rsubread::align() against a tiny synthetic
# 600bp reference, with a matching 3-gene GTF. Run once, not at test/example time —
# Rsubread::align() has a ~3-minute fixed per-call overhead regardless of input
# size, which isn't practical to pay on every test run, but the *result* is a real
# alignment, not a stub. See R/quantify.R.
library(Rsubread)

set.seed(42)
work_dir <- tempfile("quantify_fixtures_")
dir.create(work_dir)

bases <- c("A", "C", "G", "T")
ref_seq <- paste(sample(bases, 600, replace = TRUE), collapse = "")
ref_path <- file.path(work_dir, "toy_ref.fa")
writeLines(c(">chr1", ref_seq), ref_path)

idx_prefix <- file.path(work_dir, "toy_idx")
buildindex(basename = idx_prefix, reference = ref_path, indexSplit = FALSE)

read_starts <- c(10, 200, 400)
read_len <- 50
r1_lines <- c()
for (i in seq_along(read_starts)) {
  s <- read_starts[i]
  seq <- substr(ref_seq, s, s + read_len - 1)
  r1_lines <- c(r1_lines, paste0("@read", i), seq, "+", strrep("I", read_len))
}
fq_path <- file.path(work_dir, "toy_reads.fastq")
writeLines(r1_lines, fq_path)

bam_prefix <- file.path(work_dir, "toy_aligned")
align(
  index = idx_prefix, readfile1 = fq_path,
  output_file = paste0(bam_prefix, ".BAM"), phredOffset = 64, nthreads = 1
)

gene_regions <- list(c(1, 100), c(150, 260), c(350, 460))
gtf_lines <- c()
for (i in seq_along(gene_regions)) {
  g <- gene_regions[[i]]
  gtf_lines <- c(gtf_lines, paste(
    "chr1", "toy", "exon", g[1], g[2], ".", "+", ".",
    sprintf("gene_id \"TOYGENE%d\"; transcript_id \"TOYTX%d\";", i, i),
    sep = "\t"
  ))
}
gtf_path <- file.path(work_dir, "toy_annotation.gtf")
writeLines(gtf_lines, gtf_path)

extdata_dir <- "inst/extdata"
dir.create(extdata_dir, showWarnings = FALSE, recursive = TRUE)
file.copy(paste0(bam_prefix, ".BAM"), file.path(extdata_dir, "toy_aligned.bam"), overwrite = TRUE)
file.copy(gtf_path, file.path(extdata_dir, "toy_annotation.gtf"), overwrite = TRUE)

cat("Wrote", file.path(extdata_dir, "toy_aligned.bam"), "and toy_annotation.gtf\n")
