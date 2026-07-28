"""Synthetic scRNA-seq AnnData fixtures used across the scrna test suite.

None of the QC/doublet/clustering/annotation tests touch the real GSE176078 download
— CI must never depend on network access to the real dataset. These builders produce
small, structured, deterministic AnnData objects instead: multiple synthetic "samples"
(mirroring the real dataset's per-sample `orig.ident`/`subtype` columns), a handful of
Gaussian expression blobs so clustering has real structure to find, and explicitly
labelled QC-outlier and doublet cells so tests have ground truth to assert against.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

MT_GENE_NAMES = [
    "MT-ND1",
    "MT-ND2",
    "MT-CO1",
    "MT-CO2",
    "MT-ATP8",
    "MT-ATP6",
]

SUBTYPES = ["ER+", "HER2+", "TNBC"]


def _sample_metadata(n_samples: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    sample_ids = [f"SYN{i:03d}" for i in range(n_samples)]
    # Skew subtype assignment so Milo has a real compositional signal to find later,
    # rather than a uniform random assignment with no structure.
    subtypes = [SUBTYPES[i % len(SUBTYPES)] for i in range(n_samples)]
    return pd.DataFrame(
        {"orig.ident": sample_ids, "subtype": subtypes},
        index=sample_ids,
    ).sample(
        frac=1, random_state=rng
    )  # shuffle so blob assignment isn't order-biased


def build_synthetic_scrna_adata(
    *,
    seed: int = 0,
    n_samples: int = 4,
    cells_per_sample: int = 50,
    n_genes: int = 500,
    n_blobs: int = 3,
    n_doublets_per_sample: int = 2,
    n_qc_outliers_per_sample: int = 2,
    include_velocity_layers: bool = False,
) -> ad.AnnData:
    """Build a small, structured synthetic scRNA-seq AnnData object.

    Schema mirrors the real GSE176078 pipeline's expectations: ``obs['orig.ident']``
    (sample/10x-lane id) and ``obs['subtype']`` (clinical subtype) are the same column
    names the real metadata.csv uses, so QC/doublet/Milo code written against this
    fixture works unchanged against the real data. ``obs['_qc_case']`` and
    ``obs['_synthetic_doublet']`` record ground truth for tests only — production code
    never reads them.

    Returns an object with ``(n_samples * (cells_per_sample + n_doublets_per_sample +
    n_qc_outliers_per_sample))`` cells and ``n_genes`` genes (first
    ``len(MT_GENE_NAMES)`` of which are named like real human mitochondrial genes).
    """
    rng = np.random.default_rng(seed)
    sample_meta = _sample_metadata(n_samples, seed)

    var_names = MT_GENE_NAMES + [f"GENE{i:04d}" for i in range(n_genes - len(MT_GENE_NAMES))]
    n_genes_total = len(var_names)

    # A handful of Gaussian expression blobs give clustering/silhouette real structure
    # to find, rather than pure noise. Blob log-means are offset on a random subset of
    # genes each, so blobs are separable but not trivially identical to one axis.
    blob_log_means = rng.uniform(0.5, 3.0, size=(n_blobs, n_genes_total))
    blob_signal_genes = [
        rng.choice(n_genes_total, size=n_genes_total // n_blobs, replace=False)
        for _ in range(n_blobs)
    ]
    baseline_log_mean = 0.3

    rows: list[np.ndarray] = []
    obs_records: list[dict] = []

    for sample_id, sample_row in sample_meta.iterrows():
        subtype = sample_row["subtype"]
        # Subtype-skewed blob weighting: gives Milo a real per-subtype compositional
        # signal instead of uniform random cluster membership.
        subtype_idx = SUBTYPES.index(subtype)
        blob_weights = np.ones(n_blobs)
        blob_weights[subtype_idx % n_blobs] *= 3.0
        blob_weights /= blob_weights.sum()

        normal_cells = []
        normal_blobs = []
        for _ in range(cells_per_sample):
            blob = int(rng.choice(n_blobs, p=blob_weights))
            log_means = np.full(n_genes_total, baseline_log_mean)
            log_means[blob_signal_genes[blob]] = blob_log_means[blob, blob_signal_genes[blob]]
            counts = rng.poisson(np.exp(log_means)).astype(np.float32)
            normal_cells.append(counts)
            normal_blobs.append(blob)
            obs_records.append(
                {
                    "orig.ident": sample_id,
                    "subtype": subtype,
                    "_qc_case": "normal",
                    "_true_blob": blob,
                }
            )
        rows.extend(normal_cells)

        for i in range(n_qc_outliers_per_sample):
            counts = np.array(normal_cells[0], copy=True)
            if i % 2 == 0:
                # Empty-droplet-like: collapse to a handful of counts total.
                counts[:] = 0
                counts[: min(5, n_genes_total)] = rng.poisson(1, size=min(5, n_genes_total))
                qc_case = "low_count_outlier"
            else:
                # Dissociation-stress-like: mitochondrial fraction dominates.
                counts[:] = rng.poisson(0.5, size=n_genes_total)
                counts[: len(MT_GENE_NAMES)] = rng.poisson(200, size=len(MT_GENE_NAMES))
                qc_case = "high_mito_outlier"
            rows.append(counts.astype(np.float32))
            # -1: an outlier's expression is overwritten, so it has no meaningful blob
            # origin — not the blob normal_cells[0] happened to come from.
            obs_records.append(
                {
                    "orig.ident": sample_id,
                    "subtype": subtype,
                    "_qc_case": qc_case,
                    "_true_blob": -1,
                }
            )

        for _ in range(n_doublets_per_sample):
            a, b = rng.choice(len(normal_cells), size=2, replace=False)
            doublet_counts = normal_cells[a] + normal_cells[b]
            rows.append(doublet_counts.astype(np.float32))
            # -1: a doublet is a mixture of two (possibly different) blobs, no single
            # true origin.
            obs_records.append(
                {
                    "orig.ident": sample_id,
                    "subtype": subtype,
                    "_qc_case": "doublet",
                    "_true_blob": -1,
                }
            )

    x = np.vstack(rows)
    obs = pd.DataFrame(obs_records)
    obs["_synthetic_doublet"] = obs["_qc_case"] == "doublet"
    obs.index = [f"{row['orig.ident']}_cell{i:05d}" for i, row in obs.iterrows()]

    var = pd.DataFrame(index=pd.Index(var_names, name="gene_symbol"))

    adata = ad.AnnData(X=x, obs=obs, var=var)
    # Generative blob parameters, attached for tests that need to build a matching
    # synthetic reference (e.g. SingleR annotation recovery) — not read by production
    # code.
    adata.uns["_blob_signal_genes"] = blob_signal_genes
    adata.uns["_blob_log_means"] = blob_log_means
    adata.uns["_baseline_log_mean"] = baseline_log_mean

    if include_velocity_layers:
        # Synthetic spliced/unspliced layers so the scVelo module has something to run
        # against structurally. Not biologically meaningful — GSE176078's processed
        # matrix has no such layers (see ADR-0002); this fixture exists purely to
        # exercise the code path.
        spliced_fraction = rng.uniform(0.6, 0.9, size=x.shape)
        adata.layers["spliced"] = (x * spliced_fraction).astype(np.float32)
        adata.layers["unspliced"] = (x * (1 - spliced_fraction)).astype(np.float32)

    return adata
