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

### Fixed

- Several rpy2 sparse-matrix and R S4-object handling bugs found while running the
  scRNA-seq pipeline against real (large, sparse) data rather than only small dense
  test fixtures — see individual commits for detail. Most notably: `compute_size_factors`,
  `scran_normalize`, `run_singler`, and `run_scrublet_per_sample` all used to
  silently mishandle or densify sparse input, which would have needed tens of GB of
  RAM (or crashed outright) on the real ~100k-cell GSE176078 matrix.
