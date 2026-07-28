# ADR-0002: scRNA-seq Pipeline Design and Dataset Choice

## Status

Accepted

## Context

Phase 1 builds the first of four modality pipelines: scRNA-seq QC, doublet detection,
normalisation, clustering, cell-type annotation, RNA velocity, and differential
abundance testing, producing the annotated object that Phase 3 (RCTD) and Phase 5
(BayesPrism, CrossOmicsConcordance) will import as a reference cell-type signature. I
need to settle which public dataset to use, how QC/doublet thresholds are set, which
SingleR reference to annotate against, and how the two brief-scoped tools that don't
trivially fit this dataset (scVelo, Milo) are scoped, before writing any pipeline code —
these decisions constrain everything downstream, including Phases 3 and 5.

## Decision

### Dataset: GEO GSE176078 (Wu et al. 2021, *Nature Genetics*)

"A single-cell and spatially resolved atlas of human breast cancers" (Wu, Al-Eryani,
Roden et al., Swarbrick lab) — 26 primary breast tumours (11 ER+, 5 HER2+, 10 TNBC),
~100,000 cells, processed count matrix and published cell-type metadata available
directly from GEO (`GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz`, ~533MB). Only the raw
FASTQs are EGA-controlled-access; the processed matrix I need has no such restriction.

I chose this over the brief's other option (10x PBMC) because the project's central
premise is a *same-tissue* cross-modality comparison — CrossOmicsConcordance (Phase 5)
scores agreement between cell-type composition vectors derived from scRNA-seq, spatial,
and bulk deconvolution *of the same tissue*. PBMC has no epithelial, stromal, or tumour
compartment; pairing it with breast-cancer spatial (Phase 3) and bulk (Phase 2) data
would make the three modalities biologically incomparable before any concordance
analysis could even start. GSE176078 is also the field's de facto standard breast-cancer
scRNA reference (three clinical subtypes, richly annotated), which gives Phase 5 an
external published-label baseline to sanity-check this pipeline's own Leiden/SingleR
calls against.

**Noted, not acted on:** GSE176078 also ships a bulk RNA-seq raw-count matrix from the
same 26 patients (`GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz`). Matched-patient
bulk data would be a stronger concordance comparator than TCGA-BRCA's unrelated cohort.
The brief specifies TCGA-BRCA for Phase 2, so I'm not deciding this now — I've recorded
it in `data/README.md` as an option for Phase 2 review.

### QC: per-sample MAD thresholds, not a single pooled threshold

Outlier detection runs **per sample** (per patient/10x lane), not on the pooled
100k-cell object — the 26 samples vary enough in sequencing depth that a single global
threshold would systematically penalise lower-depth patients rather than flagging
genuinely bad cells within each sample.

- `log1p(total_counts)` or `log1p(n_genes_by_counts)` more than 5 MADs from the
  per-sample median → outlier. Log-transformed because both distributions are heavy
  right-tailed; 5 MADs is the standard sc-best-practices cutoff, permissive enough to
  keep the biological heterogeneity expected across three clinical subtypes while still
  catching genuine extremes.
- `pct_counts_mt` more than 3 MADs above the per-sample median, **or** above an
  absolute 20% ceiling → outlier. 20% rather than the PBMC-typical 5–10%: this is
  enzymatically dissociated solid tumour tissue, where dissociation stress elevates
  mitochondrial fraction even in healthy cells, and is consistent with the threshold Wu
  et al. used for this same tissue type.
- `n_genes_by_counts < 200` → dropped regardless of MAD, as an absolute floor against
  empty droplets/debris even in an already-low-quality sample.
- Genes expressed in fewer than 3 cells are dropped (`sc.pp.filter_genes(min_cells=3)`).

### Doublet detection: Scrublet, run per-sample with scaled expected rate

Scrublet's simulated doublets are only meaningful within a single 10x GEM well, so it
runs separately per sample, never on the pooled object. `expected_doublet_rate` is
scaled per sample by loaded cell count (~0.8% per 1,000 cells, approximating 10x's
published loading/doublet-rate table, capped at 20%) rather than applying Scrublet's
flat 6% default uniformly regardless of sample size. The primary doublet call uses
Scrublet's automatic bimodal threshold detection; where that fails (unimodal score
histogram, which happens in low-doublet-rate samples), a fixed `doublet_score > 0.25`
fallback is used, and every fallback is logged per-sample rather than applied silently.

### Cell-type annotation: SingleR against `celldex::HumanPrimaryCellAtlasData()`

The standard, most widely used pan-tissue SingleR reference, covering the lineages
present in a breast tumour microenvironment (epithelial, T/B/myeloid immune subsets,
fibroblasts, endothelial cells). No standard, citable breast-cancer-specific SingleR
reference exists, so a well-established pan-tissue reference is the defensible choice
for broad-lineage annotation.

**Documented limitation:** this reference has no "malignant epithelial" category —
SingleR calls tumour and normal epithelial cells both "epithelial." Separating malignant
from normal epithelium needs copy-number inference (e.g. inferCNV), which is out of
scope for Phase 1. This doesn't block Phase 3/5 concordance work, since the
deconvolution methods used there typically operate at the same lineage-level
resolution, but it is a real limitation, stated here rather than left implicit.

### Clustering: Leiden with an empirical resolution sweep, silhouette-selected

Resolution is not fixed in advance. The pipeline sweeps a range of Leiden resolutions,
scores each with silhouette on the PCA embedding, and selects by silhouette while
sanity-checking cluster count against Wu et al.'s published major-lineage count. Exact
pairwise silhouette is O(n²) and infeasible at ~100k cells, so silhouette is computed on
a fixed random subsample (`sample_size=5000, random_state=0`), a standard practice at
this data scale — the sweep table and chosen resolution are recorded empirically in
`docs/scrna-pipeline.md` once run, not asserted here.

### scVelo: implemented, but not scientifically applicable to this dataset

RNA velocity needs spliced/unspliced count layers (from velocyto/STARsolo against raw
FASTQs). GSE176078's processed download is a standard count matrix with no such layers,
and the raw FASTQs are EGA-controlled-access. The `velocity.py` module is implemented
and tested against a synthetic fixture with injected spliced/unspliced layers, for
pipeline completeness and because the brief scopes scVelo into Phase 1, but it is not
run against GSE176078 and its output is not part of the signature artifact —
`docs/scrna-pipeline.md` states this plainly.

### Milo: meaningful here, implemented via `pertpy`

Unlike scVelo, Milo has a real grouping variable to test: the three clinical subtypes
(ER+/HER2+/TNBC) recorded in GSE176078's published per-patient metadata. Milo tests
differential neighbourhood abundance across those subtypes, a genuine, documentable
result. Implemented via `pertpy` (scverse's Python-native Milo port) rather than an
rpy2 wrapper around R's `miloR` — scran and SingleR stay on the R bridge because they
have no solid Python equivalent, but Milo does, so adding another R round-trip for it
would be unjustified complexity for no benefit.

### The signature artifact is a stable, versioned, explicitly-imported build product

The final annotated object (QC-filtered, doublet-removed, normalised, clustered,
SingleR-labelled) is saved as a versioned `.h5ad` under a path returned by
`signature_path(version=...)` in `src/omicsatlas/scrna/artifact.py`. Its schema (what's
in `.obs`, `.var`, `.obsm`, `.uns`) is documented in `DATA_SCHEMA.md`. Phases 3 and 5
import this exact versioned path explicitly — they do not regenerate the signature
themselves with potentially different parameters. The `.h5ad` file itself is a build
product, not a committed file (`data/processed/` is gitignored, per Phase 0); what's
committed is the schema contract, the version identifier, and `make scrna-signature`
(real data) / `make scrna-signature-fixture` (synthetic fixture, for CI-safe contract
testing) as the reproducible commands that produce it.

### CI runs a real conda/mamba environment, not just plain Python

Phase 0's CI (`actions/setup-python` only) cannot run the scran/SingleR rpy2 bridge —
there's no R interpreter. Rather than mock or skip the R-dependent code paths in CI
(which would leave that integration completely untested in the one place that gates
merges), CI now builds the same `environment/env.yml` via a conda/mamba action, so it
genuinely exercises the R bridge against the synthetic fixtures. This still has zero
dependency on the real dataset — only on installing the R/Bioconductor packages
themselves, which is a build-environment concern, not a data-provenance one.

`bioconductor-singler` has no `osx-arm64` conda build as of this writing (only
`linux-64`/`osx-64`/`linux-aarch64`). Local development on this Apple Silicon Mac uses
`CONDA_SUBDIR=osx-64` (Rosetta 2) to get a working build; CI runs on `ubuntu-latest`,
where native `linux-64` builds exist for the full R stack, so this is a local-dev-only
workaround, not a CI concern.

## Consequences

- Every downstream phase (3, 5) has a single, explicit, versioned signature to import —
  no ambiguity about which QC/clustering parameters produced the cell-type calls they're
  benchmarking against.
- CI is now meaningfully slower (conda/mamba + Bioconductor resolution) than Phase 0's
  plain pip install, in exchange for actually testing the R bridge rather than assuming
  it works.
- scVelo's inclusion is partly for brief-completeness rather than scientific payoff on
  this specific dataset; the manuscript/benchmark work in Phase 5 won't lean on velocity
  results from this pipeline.
- The malignant/normal epithelial ambiguity in SingleR's output is an accepted,
  documented limitation carried forward into Phase 3/5, not resolved here.

## Alternatives considered

- **10x PBMC** as the scRNA-seq dataset: rejected — biologically incompatible with the
  breast-cancer spatial/bulk/ATAC data the rest of the project uses, which would break
  the concordance benchmark's core premise.
- **A breast-cancer-specific curated SingleR reference** (e.g. built from Wu et al.'s
  own published annotations): rejected for Phase 1 — using the paper's own labels as
  the reference to annotate its own cells is circular for validating this pipeline's
  independent annotation; the published labels are instead kept as an external
  sanity-check baseline, not the reference itself.
- **rpy2-wrapped `miloR`** for differential abundance: rejected in favour of `pertpy`'s
  native Python Milo implementation, since it avoids an unnecessary additional R
  round-trip when a solid native alternative exists.
- **Pooled (non-per-sample) QC and doublet detection**: rejected — both MAD thresholds
  and Scrublet's doublet simulation are only statistically valid within a single
  sample/GEM well; pooling across 26 samples of varying depth and composition would
  bias both.
