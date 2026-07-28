# OmicsAtlas — Phase Status

Tracks progress against the phase plan. Updated at the start/end of each phase.

- [x] **Phase 0 — Foundations**: repo scaffold, ADR-0001 (license + architecture),
      pyproject.toml, pre-commit hooks, conda env, Dockerfile, CI skeleton,
      CONTRIBUTING, issue/PR templates.
- [x] **Phase 1 — scRNA-seq pipeline**: QC, Scrublet, scran normalisation, HVG/PCA/UMAP,
      Leiden clustering, SingleR annotation, scVelo, Milo, Streamlit UMAP browser.
- [ ] **Phase 2 — Bulk RNA-seq pipeline (R package)**: STAR/featureCounts, DESeq2 +
      edgeR, clusterProfiler, BayesPrism deconvolution, pkgdown site.
- [ ] **Phase 3 — Spatial transcriptomics arm**: Visium ingestion, SCTransform,
      Moran's I, RCTD deconvolution, Squidpy, LIANA/NicheNet.
- [ ] **Phase 4 — ATAC-seq pipeline (Snakemake)**: Bowtie2, MACS3, DiffBind, chromVAR,
      HOMER, deepTools, Seurat v5 bridge integration.
- [ ] **Phase 5 — CrossOmicsConcordance (novel contribution)**: metric ADR, independent
      installable package, benchmark study, PyPI + Bioconda publication.
- [ ] **Phase 6 — Unified browser & integration**: cellxgene deployment, unified
      Streamlit/Dash app, Nextflow DSL2 master workflow.
- [ ] **Phase 7 — Documentation, manuscript, release**: Quarto site, methods paper,
      Zenodo DOI, CITATION.cff finalised, README rewrite, v1.0.0 tag.

## Current status

Phase 1 complete (2026-07-28): full scRNA-seq pipeline (fetcher with checksum
provenance, MAD-based QC, per-sample Scrublet, scran normalisation, HVG/PCA/UMAP,
Leiden with an empirical resolution sweep, SingleR annotation, scVelo, Milo, a basic
Streamlit UMAP browser, and the versioned signature artifact + schema contract Phases
3/5 will depend on) built, tested (fixture-only in CI, matching the constraint from
Phase 0 review), and run end-to-end against real GSE176078 data (6-patient subset,
21,585 cells — see `docs/scrna-pipeline.md` for results; full 26-patient cohort is
HPC-scale, documented as a follow-up). ADR-0002 records the design decisions. CI green
on GitHub Actions at
[NosakhareOsaro/omicsatlas](https://github.com/NosakhareOsaro/omicsatlas). Pausing for
review before starting Phase 2 (bulk RNA-seq pipeline, R package).

Phase 0 complete (2026-07-28): repo scaffolded, ADR-0001 accepted (MIT license,
monorepo architecture, Python 3.11), packaging/lint/CI/Docker all verified locally
and green on GitHub Actions.
