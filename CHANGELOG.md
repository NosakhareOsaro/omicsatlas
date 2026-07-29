# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Repository scaffold for the multi-modal pipeline (scrna, bulk, spatial, atac,
  concordance, browser) and supporting directories.
- ADR-0001: license (MIT) and overall architecture decision.
- Python packaging (`pyproject.toml`, hatchling, Python 3.11), pre-commit hooks
  (ruff, black, mypy), conda environment, Docker image (pinned by digest), and a
  lint/test GitHub Actions CI workflow.
- `CONTRIBUTING.md`, issue templates, and PR template.
- scRNA-seq pipeline (`src/omicsatlas/scrna/`): checksum-verified GSE176078 fetcher,
  MAD-based per-sample QC, per-sample Scrublet doublet detection, scran normalisation
  via an rpy2/R bridge, HVG/PCA/UMAP preprocessing, Leiden clustering with an
  empirical resolution sweep, SingleR annotation against
  `celldex::HumanPrimaryCellAtlasData()`, a scVelo module (fixture-tested; not
  scientifically applicable to GSE176078), Milo differential neighbourhood abundance
  via `pertpy`/`pydeseq2`, and a pipeline entrypoint producing a versioned signature
  artifact (`make scrna-signature` / `make scrna-signature-fixture`).
- A basic Streamlit UMAP browser (`src/omicsatlas/browser/app.py`) for the scRNA
  signature artifact.
- ADR-0002: scRNA-seq pipeline design (dataset choice, QC/doublet methodology,
  SingleR reference and its limitations, scVelo/Milo scope decisions, the signature
  artifact contract).
- `DATA_SCHEMA.md`: the scRNA-seq signature artifact's schema contract for Phases 3
  and 5.
- `docs/scrna-pipeline.md`: rationale for every biological/statistical choice in the
  scRNA-seq pipeline, plus empirical results from a real end-to-end run against a
  6-patient GSE176078 subset.
- CI now builds a full conda/mamba environment (R/Bioconductor included) instead of a
  plain Python install, so it can exercise the R bridge.
- `OmicsAtlasBulk` R package (`r_package/OmicsAtlasBulk/`): STAR/featureCounts
  wrappers, DESeq2 + edgeR differential expression, `compare_de_methods()`,
  clusterProfiler-based enrichment, BayesPrism deconvolution against the Phase 1
  scRNA-seq signature, an auto-generated Quarto report, and two pipeline
  entrypoints - `scripts/run_bulk_pipeline_matched.R` (primary: the GSE176078-matched
  bulk cohort) and `scripts/run_bulk_pipeline_tcga.R` (secondary: TCGA-BRCA). See
  ADR-0003 (dataset choice) and ADR-0004 (package architecture).
- ADR-0003: bulk RNA-seq data choice - GSE176078-matched bulk (24 of the 26 Phase 1
  scRNA-seq patients) as the primary dataset feeding Phase 5's CrossOmicsConcordance
  benchmark, TCGA-BRCA as a secondary generalisability check.
- ADR-0004: R package architecture and Bioconductor readiness, later extended with
  the BayesPrism cell-state collapse and `n_cores` rationale below.

### Fixed

- Several rpy2 sparse-matrix and R S4-object handling bugs found while running the
  scRNA-seq pipeline against real (large, sparse) data rather than only small dense
  test fixtures — see individual commits for detail. Most notably: `compute_size_factors`,
  `scran_normalize`, `run_singler`, and `run_scrublet_per_sample` all used to
  silently mishandle or densify sparse input, which would have needed tens of GB of
  RAM (or crashed outright) on the real ~100k-cell GSE176078 matrix.
- Three bugs found on the first real (non-fixture) run of the matched-bulk pipeline:
  `read_tcga_brca_counts()`'s default `repo_root` resolved one directory too shallow
  (`dirname(dirname(manifest_path))` instead of three levels up); `load_scrna_signature()`
  used zellkonverter's default Python/basilisk-backed `.h5ad` reader, which tried to
  bootstrap an entirely new from-source Python on first real use, now
  `reader = "R"`; and `run_bulk_pipeline_matched.R`'s Python `signature_path()`
  resolution passed an unquoted shell argument to `system2()`, now `shQuote()`-wrapped.
- `run_bayesprism_deconvolution()`'s `cell_state_labels` didn't nest within
  `cell_type_labels` as BayesPrism requires (`leiden` clusters are computed
  independently of `singler_label` types) - fixed by combining them into one label.
  That fix alone then produced 99 BayesPrism cell states on the real 21,585-cell
  Phase 1 signature, 53% with fewer than 5 cells, driving Gibbs sampling's estimated
  runtime to 13+ hours; fixed by collapsing combinations below `min_state_cells`
  (default 30) into a per-type `"other"` state (99 → 42 states; see ADR-0004).
- Even after that collapse, the real run's Gibbs sampling estimate was still 7hrs
  15mins at the package's default `n_cores = 1`; `run_bulk_pipeline_matched.R` now
  explicitly passes `n_cores = 8` (BayesPrism parallelizes per bulk sample, so this
  scales close to linearly), bringing the real run down to ~113 minutes wall-clock.
  Mid-run, an unrelated leftover process - a stray, pre-fix, single-core run of the
  same script, started before this session began and almost certainly a holdover
  from the VS Code crash that interrupted the previous session - was found still
  executing and competing for the same CPU cores, and was killed; it accounted for
  part of the run's elapsed time and could otherwise have silently overwritten this
  run's output hours later.
