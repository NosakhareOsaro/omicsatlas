# scRNA-seq pipeline: rationale

This documents *why* each biological/statistical choice in the scRNA-seq pipeline
(`src/omicsatlas/scrna/`) was made, not how to run it — see `README.md` and the
module docstrings for usage. Full design context and alternatives considered are in
`adr/ADR-0002-scrna-pipeline-design.md`; this page is the narrower, results-focused
companion to it.

## Dataset

GEO [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) — Wu et
al. 2021 (*Nature Genetics*), "A single-cell and spatially resolved atlas of human
breast cancers." 26 primary tumours (11 ER+, 5 HER2+, 10 TNBC), ~100,064 cells,
29,733 genes. I chose this over the alternative in the brief (10x PBMC) because the
whole point of this project is a same-tissue cross-modality comparison — see
ADR-0002 for the full reasoning.

## QC: MAD thresholds

Computed **per sample** (per patient/10x lane), not pooled across the cohort — the 26
samples vary enough in sequencing depth that a pooled threshold would systematically
penalise lower-depth patients rather than flagging genuinely bad cells within each
sample.

- **Library size / complexity**: `log1p(total_counts)` or `log1p(n_genes_by_counts)`
  more than 5 MADs from the per-sample median. 5 MADs is the sc-best-practices
  convention (Theis lab QC tutorial) — permissive enough to keep the real biological
  heterogeneity expected across three clinical subtypes, while still catching genuine
  extremes. I use two-sided MAD here (not just a lower bound) because both
  unusually-low and unusually-high library size/complexity are suspicious in
  practice (empty droplets vs. multiplets/debris aggregates).
- **Mitochondrial fraction**: more than 3 MADs above the per-sample median, *or*
  above an absolute 20% ceiling. I used 20% rather than the PBMC-typical 5–10%
  because this is enzymatically dissociated solid tumour tissue — dissociation
  stress elevates mitochondrial fraction even in healthy cells, and 20% is consistent
  with the threshold Wu et al. themselves used for this same tissue type in their own
  QC (see their Methods).
- **Absolute floor**: `n_genes_by_counts < 200` dropped regardless of MAD, to catch
  empty droplets/debris even in an already-low-quality sample where the MAD baseline
  itself might be depressed.

Per-sample QC pass/fail counts from the real run (6-patient subset — see
"Empirical run" below) are in the table there.

## Doublet detection: Scrublet, per sample

Scrublet's simulated doublets are only statistically meaningful within a single 10x
GEM well, so it's run separately per sample, never on the pooled object — pooling
would let Scrublet "detect" cross-sample composition differences as spurious
doublets.

`expected_doublet_rate` is scaled per sample by loaded cell count (~0.8% per 1,000
cells, approximating 10x's published loading/doublet-rate table, capped at 20%)
rather than Scrublet's flat 6% default applied uniformly regardless of sample size —
a 600-cell sample and a 6,000-cell sample from the same 10x chemistry do not have the
same expected doublet rate.

I validated the per-sample Scrublet wrapper on synthetic data before trusting it on
real data: at the project's standard ~40-cell/sample test-fixture scale, Scrublet's
KNN/PCA classifier has essentially no power (injected-doublet vs. normal-cell scores
were statistically indistinguishable) — it only produces a meaningful separation once
a sample has a few hundred cells, which is why the fixture used for detection-quality
testing (`tests/test_doublets.py`) is deliberately larger than the project's other
fixtures. This matters for reading real-sample results in the table below: smaller
patients in the subset should be expected to have less reliable doublet calls than
larger ones, for the same underlying statistical reason.

## Normalisation: scran deconvolution

`scran::computeSumFactors` (deconvolution-based size factor estimation) rather than
scanpy's simple total-count normalisation, run via an rpy2 bridge to R/Bioconductor —
see ADR-0002 for why this stays on the R bridge instead of a Python port. Kept sparse
throughout (both the R-side matrix and the Python-side log-normalised layer) —
GSE176078 is large and sparse enough (~100k cells, ~6% density) that a dense
float64 copy would need tens of GB of RAM; this was a real bug I hit and fixed while
preparing the empirical run below, not a preemptive optimisation.

## Clustering: Leiden with an empirical resolution sweep

Resolution is not fixed in code. `sweep_leiden_resolutions` runs Leiden at each
candidate resolution and scores it with silhouette on the PCA embedding (subsampled
for tractability — exact pairwise silhouette is O(n²), infeasible at real-dataset
scale), and `choose_best_resolution` picks the highest-scoring one. The actual sweep
and chosen resolution from the real run are in "Empirical run" below.

## Cell-type annotation: SingleR against `celldex::HumanPrimaryCellAtlasData()`

The standard, most widely used pan-tissue SingleR reference, covering the lineages
present in a breast tumour microenvironment (epithelial, T/B/myeloid immune subsets,
fibroblasts, endothelial cells). No standard, citable breast-cancer-specific SingleR
reference exists, so a well-established pan-tissue reference is the defensible
choice for broad-lineage annotation.

**Known limitation, not solved here**: this reference has no "malignant epithelial"
category — SingleR calls tumour and normal epithelial cells both "epithelial." It
cannot distinguish them; that needs copy-number inference (e.g. inferCNV) from
expression, out of scope for this pipeline. This doesn't block the Phase 3/5
concordance work, since the deconvolution methods used there typically operate at
the same lineage-level resolution, but it's a real limitation worth stating plainly
rather than leaving implicit.

## scVelo and Milo: one is meaningful here, one isn't

**scVelo is implemented but not run meaningfully against this dataset.** RNA velocity
needs spliced/unspliced count layers (from velocyto/STARsolo against raw FASTQs);
GSE176078's processed download is a standard count matrix with neither layer, and the
raw FASTQs are EGA-controlled-access. The module is tested against a synthetic
fixture with injected spliced/unspliced layers, for pipeline completeness, but its
output isn't part of the signature artifact and shouldn't be treated as a real
biological result on this data.

**Milo is meaningful here.** The cohort carries a real grouping variable — clinical
subtype (ER+/HER2+/TNBC) — so Milo tests differential neighbourhood abundance across
subtypes, a genuine, documentable question. Implemented via `pertpy`'s Python-native
Milo port with the `pydeseq2` solver, keeping it off the R bridge entirely (scran and
SingleR stay on the R bridge because they have no solid Python equivalent; Milo does,
so an extra R round-trip for it would be unjustified).

## Empirical run

Run against a **6-patient subset** of GSE176078 (2 patients per clinical subtype:
ER+ `CID4067`/`CID4398`, HER2+ `CID3921`/`CID4066`, TNBC `CID3963`/`CID4515`), not the
full 26-patient/~100,064-cell cohort. This is a deliberate scope decision, not a
shortcut: the full cohort is HPC-scale (matches the project's stated reproducibility
policy — full runs are documented as needing HPC/cloud, small/laptop-scale runs
validate the pipeline itself). This subset is still real GSE176078 data, still
multi-patient, still spans all three clinical subtypes — it validates the pipeline
end-to-end on real data, which is what this section is for; it is not a substitute
for a full-cohort scientific analysis.

Command:

```python
from omicsatlas.scrna.pipeline import run_pipeline, GSE176078_EXTRACTED_DIR
from omicsatlas.scrna.artifact import REPO_ROOT

raw_dir = REPO_ROOT / "data" / "raw" / "scrna" / "gse176078" / GSE176078_EXTRACTED_DIR
result = run_pipeline(
    raw_dir,
    patient_subset=["CID4067", "CID4398", "CID3921", "CID4066", "CID3963", "CID4515"],
)
```

Ran end-to-end in **5m 10s** on a single laptop (no HPC).

### QC / doublet retention per sample

| Sample | Subtype | Raw cells (GEO metadata) | Cells in signature (post-QC, post-doublet-removal) | Retained |
|---|---|---|---|---|
| CID3921 | HER2+ | 3,024 | 2,586 | 85.5% |
| CID3963 | TNBC  | 3,527 | 2,924 | 82.9% |
| CID4066 | HER2+ | 5,309 | 4,736 | 89.2% |
| CID4067 | ER+   | 3,764 | 3,705 | 98.4% |
| CID4398 | ER+   | 4,451 | 3,770 | 84.7% |
| CID4515 | TNBC  | 4,149 | 3,864 | 93.1% |
| **Total** | | **24,224** | **21,585** | **89.1%** |

None of the 6 samples needed Scrublet's fixed-threshold fallback — automatic bimodal
threshold detection succeeded for every sample at this scale (0/21,585 cells flagged
`doublet_threshold_fallback_used`).

One real quirk worth recording honestly: `scran::computeSumFactors` logged 8
`"encountered non-positive size factor estimates"` warnings during its internal
pooling/deconvolution computation on this run. These are warnings from scran's
intermediate calculations, not the final output — `scran_normalize`'s own validation
(which raises if any *final* returned size factor is non-positive or non-finite) did
not trigger, and the saved artifact's `scran_size_factor` column confirms every final
value is strictly positive (min `7.3e-7`, mean `1.0` as expected of scran's
normalisation, max `21.0`). I'm not suppressing or explaining away the warning here —
just noting it didn't propagate to an invalid result on this run, and it's worth
watching for on future/larger runs.

### Leiden resolution sweep

Silhouette computed on a 5,000-cell subsample of the PCA embedding (see
`cluster.py`/ADR-0002 for why exact silhouette isn't tractable at this scale):

| Resolution | Clusters | Silhouette |
|---|---|---|
| 0.2 | 8  | **0.117** (chosen) |
| 0.4 | 10 | 0.100 |
| 0.6 | 11 | 0.097 |
| 0.8 | 13 | 0.087 |
| 1.0 | 16 | 0.085 |
| 1.2 | 16 | 0.078 |
| 1.4 | 19 | 0.076 |

Silhouette decreases monotonically as resolution (and cluster count) increases on
this subset — unsurprising, since finer partitions of a continuous expression
manifold tend to produce less-separated clusters by this metric. `0.2` (8 clusters)
was selected as the highest-scoring resolution.

### SingleR labels

| Label | Cells |
|---|---|
| T_cells | 8,733 |
| Epithelial_cells | 4,996 |
| NK_cell | 1,835 |
| Tissue_stem_cells | 1,162 |
| Endothelial_cells | 980 |
| B_cell | 958 |
| Macrophage | 915 |
| Fibroblasts | 650 |
| Monocyte | 611 |
| Chondrocytes | 393 |
| DC | 125 |
| Smooth_muscle_cells | 118 |
| (11 further labels, each <50 cells) | 109 |

A sensible mix for a breast tumour microenvironment — dominated by T cells and
epithelial (tumour + normal, per the documented limitation above) cells, with the
expected supporting immune (NK, B, macrophage, monocyte, DC) and stromal
(endothelial, fibroblast, smooth muscle) populations alongside them. This is a
reasonable sanity check that the pipeline behaves sensibly on real data — not a
validated ground-truth comparison against Wu et al.'s own published labels, which
would need a proper agreement analysis (exactly the kind of cross-method comparison
Phase 5's CrossOmicsConcordance metric is for).

### Leiden clusters are subtype-skewed, consistent with tumour-intrinsic epithelial programs

Cross-tabulating Leiden cluster against clinical subtype shows several clusters
strongly dominated by one subtype (cluster 7: 3,683/3,714 ER+; cluster 5: 2,298/2,360
ER+; cluster 1: 1,978/1,982 TNBC), while others are more evenly mixed across subtypes
(cluster 2, the largest, spans all three: 729 ER+ / 3,481 HER2+ / 2,881 TNBC). Cross-
referencing against the SingleR labels, the subtype-dominated clusters are
epithelial-enriched and the mixed clusters are immune-enriched (e.g. cluster 2 is
majority T_cells) — consistent with the expected biology: tumour epithelial cells
carry strong patient/subtype-specific transcriptional programs, while immune cell
states are shared more broadly across the tumour microenvironment regardless of
subtype. This is a coarse, single-subset observation, not a rigorously tested claim,
but it's the kind of sanity check that gives some confidence the pipeline is
producing biologically sensible structure rather than noise.

### Scope note

This 6-patient/21,585-cell run validates the pipeline end-to-end against real data.
The full 26-patient/~100,064-cell cohort was not run for this phase — see "Empirical
run" intro above for why. Re-running against the full cohort (`run_pipeline(raw_dir)`
with no `patient_subset`) on suitable HPC/cloud compute would be a natural follow-up
before treating any of the numbers above as a finished scientific result rather than
a pipeline validation.
