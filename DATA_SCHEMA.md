# scRNA-seq Signature Artifact Schema

This documents the schema of the versioned scRNA-seq signature artifact produced by
`src/omicsatlas/scrna/pipeline.py` and returned by `signature_path()` in
`src/omicsatlas/scrna/artifact.py`. Phases 3 (RCTD) and 5 (BayesPrism,
CrossOmicsConcordance) import this exact versioned path explicitly — if you change
this schema, bump `CURRENT_SIGNATURE_VERSION` in `artifact.py` rather than silently
changing what a given version identifier means.

The artifact itself (`.h5ad`) is a build product, not committed to git
(`data/processed/` is gitignored). Regenerate it with `make scrna-signature` (real
GSE176078 data) or `make scrna-signature-fixture` (synthetic fixture, for schema
contract testing without the real download).

The columns below apply to both. The fixture artifact additionally carries a few
underscore-prefixed test-only columns (`_qc_case`, `_true_blob`, `_synthetic_doublet`)
from the synthetic fixture generator — these aren't part of the contract and won't
appear in the real artifact; downstream consumers should ignore them.

## `.X`

Not used directly downstream — raw QC-passed, doublet-removed counts. Consumers should
use `.layers['scran_normalized']` for expression values.

## `.layers`

| Key | Contents |
|---|---|
| `scran_normalized` | scran size-factor normalised, log1p-transformed expression (see `normalize.py`). This is the layer used for clustering and SingleR annotation. |

## `.obs`

| Column | Type | Contents |
|---|---|---|
| `orig.ident` | str | Sample/patient ID (10x GEM well) — GSE176078's own column name, kept as-is. |
| `subtype` | str | Clinical subtype: `ER+`, `HER2+`, or `TNBC`. |
| `total_counts`, `n_genes_by_counts`, `pct_counts_mt` | float | QC metrics from `scanpy.pp.calculate_qc_metrics` (see `qc.py`). |
| `qc_outlier`, `qc_fail_reason` | bool, str | Per-sample MAD-based QC outlier flag and reason(s); all cells in the saved artifact have `qc_outlier == False` (outliers are filtered before saving). |
| `doublet_score`, `predicted_doublet` | float, bool | Per-sample Scrublet output; all cells in the saved artifact have `predicted_doublet == False`. |
| `scran_size_factor` | float | Per-cell scran size factor used to produce `layers['scran_normalized']`. |
| `leiden` | str (categorical) | Leiden cluster assignment at the empirically chosen resolution (see `.uns['leiden_chosen_resolution']`). |
| `singler_label` | str | SingleR predicted cell type. For `make scrna-signature` (real data), from `celldex::HumanPrimaryCellAtlasData()` — see ADR-0002 for the documented malignant/normal-epithelial limitation of this reference. For `make scrna-signature-fixture`, from a data-driven placeholder reference (`annotate.build_placeholder_reference`) with meaningless `placeholder_type_N` labels — the fixture artifact exists to prove this schema, not to produce real cell-type calls. |

## `.var`

| Column | Type | Contents |
|---|---|---|
| `mt` | bool | Whether the gene is a mitochondrial gene (`MT-` prefix), used for QC. |
| `highly_variable` | bool | Whether the gene was selected as a highly variable gene for PCA/clustering. |

## `.obsm`

| Key | Contents |
|---|---|
| `X_pca` | PCA embedding used for neighbours/UMAP/Leiden/silhouette scoring. |
| `X_umap` | 2D UMAP embedding, for visualisation (the browser reads this directly). |

## `.uns`

| Key | Contents |
|---|---|
| `leiden_resolution_sweep` | DataFrame with `resolution`, `n_clusters`, `silhouette` columns — one row per resolution tried in the empirical sweep (see `docs/scrna-pipeline.md` for the actual numbers from the real run). |
| `leiden_chosen_resolution` | The resolution selected by `choose_best_resolution()` (highest silhouette). |
| `hvg` | scanpy's standard HVG metadata. |
| `neighbors` | scanpy's standard neighbour-graph metadata. |

## Known limitations (carried forward, not resolved in this artifact)

- `singler_label` cannot distinguish malignant from normal epithelial cells (no such
  category in `HumanPrimaryCellAtlasData()`). See ADR-0002.
- No RNA velocity results are included — scVelo is not scientifically applicable to
  GSE176078's processed matrix (no spliced/unspliced layers). See ADR-0002.
- Milo differential-abundance results (`src/omicsatlas/scrna/milo.py`) are a separate
  analysis output, not part of this artifact — they test a hypothesis about the cohort,
  they don't annotate individual cells.
