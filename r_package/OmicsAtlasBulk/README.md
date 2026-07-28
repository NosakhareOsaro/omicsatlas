# OmicsAtlasBulk

The bulk RNA-seq arm of the [OmicsAtlas](https://github.com/NosakhareOsaro/omicsatlas)
project: alignment orchestration (STAR), quantification
(`Rsubread::featureCounts()`), differential expression (DESeq2 and edgeR run in
parallel, with an explicit method-comparison step rather than silently picking one),
functional enrichment (clusterProfiler GO/KEGG/GSEA), and BayesPrism deconvolution
against the project's scRNA-seq signature.

Runs against two datasets — see
[`adr/ADR-0003-bulk-data-choice.md`](../../adr/ADR-0003-bulk-data-choice.md) in the
parent repository for the full reasoning: a GSE176078-matched bulk cohort (24
patients exact-ID-matched to the project's scRNA-seq data — the **primary** dataset
feeding the project's CrossOmicsConcordance benchmark) and a TCGA-BRCA subset (a
**secondary**, generalisability-only analysis).

Architecture and Bioconductor-readiness decisions are recorded in
[`adr/ADR-0004-r-package-architecture.md`](../../adr/ADR-0004-r-package-architecture.md).

## Installation

This package lives inside the OmicsAtlas monorepo, not as a standalone repository:

```r
# from the repo root
devtools::install("r_package/OmicsAtlasBulk")
```

## Example

```r
library(OmicsAtlasBulk)

# A small synthetic dataset (20 genes x 8 samples, 2 conditions) is bundled for
# examples and tests, rather than requiring the real ~500MB+ alignment inputs.
data(example_bulk_se)
example_bulk_se
```
