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

## Known limitations / follow-ups

- Phase 1's generic SingleR reference produces some biologically implausible
  cell-type labels (`Chondrocytes`, `Smooth_muscle_cells` at ~9% each in Phase 2
  matched-bulk deconvolution) in tissue types the reference wasn't designed for.
  Tracked as a candidate Phase 5 follow-up (breast-specific reference or CNV-based
  malignant/normal refinement) since it affects labels feeding both RCTD and
  BayesPrism deconvolution.

## Current status

Phase 2 in progress (2026-07-29): matched-bulk (primary) arm complete and reviewed.
`OmicsAtlasBulk`'s full pipeline code (align/quantify/DESeq2/edgeR/enrichment/
BayesPrism deconvolution, both entrypoint scripts, pkgdown site) was already written
and fixture-tested; this session ran `scripts/run_bulk_pipeline_matched.R` against
the real GSE176078-matched bulk data (24 patients) and the real Phase 1 signature
(21,585 cells) for the first time. That first real run surfaced three bugs (a
`repo_root` path off-by-one in `read_tcga_brca_counts()`, a BayesPrism cell-state
nesting requirement `leiden` clusters don't satisfy on their own, and unquoted
shell args in the TCGA path-resolution `system2()` call), all fixed, then a second
issue: the type/state nesting fix produced 99 near-singleton BayesPrism cell states
on the real signature, driving Gibbs sampling's estimated runtime to 13+ hours.
Fixed by collapsing sub-`min_state_cells` (default 30) combinations into a per-type
`"other"` state (99 → 42 states; see ADR-0004), and running with `n_cores = 8`
(BayesPrism parallelizes per bulk sample) instead of the package default of 1,
bringing the real run down to **~113 minutes wall-clock** (Gibbs sampling estimates
of 45 + 26 minutes, plus data loading/alignment overhead and some CPU contention -
see below). Mid-run, an unrelated leftover process (a stray, pre-fix `n_cores = 1`
run of the same script, started before this session began, almost certainly a
holdover from the VS Code crash that interrupted the previous session) was found
still executing on the same machine and killed; it explained part of the observed
runtime and could otherwise have silently overwritten this run's output hours
later. Result: 24 samples, cell-type proportions summing to 1 in every sample,
`Epithelial_cells` dominant in 22/24 (biologically expected for tumour bulk tissue);
non-trivial `Chondrocytes`/`Smooth_muscle_cells` proportions trace to Phase 1's
generic pan-tissue SingleR reference (`celldex::HumanPrimaryCellAtlasData()`), not a
Phase 2 defect. TCGA-BRCA secondary analysis (DESeq2/edgeR/clusterProfiler/Quarto
report) still to come before Phase 2 as a whole is complete.

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
