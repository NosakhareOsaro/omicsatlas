"""HVG selection, PCA, UMAP, and Leiden clustering with an empirical resolution sweep.

Resolution is not fixed in advance — ``sweep_leiden_resolutions`` scores each
candidate with silhouette on the PCA embedding (subsampled at scale, since exact
pairwise silhouette is O(n^2)) and the final chosen resolution is recorded empirically
in docs/scrna-pipeline.md once run against the real data. See ADR-0002.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_score

DEFAULT_RESOLUTIONS = (0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)
DEFAULT_SILHOUETTE_SAMPLE_SIZE = 5000


def preprocess_for_clustering(
    adata: ad.AnnData,
    *,
    layer_key: str = "scran_normalized",
    n_top_genes: int = 2000,
    n_pcs: int = 50,
    n_neighbors: int = 15,
    random_state: int = 0,
) -> ad.AnnData:
    """HVG selection -> PCA -> neighbours -> UMAP, in place. Uses ``layers[layer_key]``
    (log-normalised expression) as the working matrix rather than raw counts. Returns
    ``adata`` for chaining."""
    if layer_key not in adata.layers:
        raise ValueError(f"Missing layer {layer_key!r}; run scran_normalize() first.")

    working = adata.copy()
    working.X = working.layers[layer_key]

    n_top_genes = min(n_top_genes, working.n_vars)
    sc.pp.highly_variable_genes(working, n_top_genes=n_top_genes, flavor="seurat")
    working = working[:, working.var["highly_variable"]].copy()

    n_pcs = min(n_pcs, working.n_obs - 1, working.n_vars - 1)
    sc.pp.pca(working, n_comps=n_pcs, random_state=random_state)
    sc.pp.neighbors(
        working, n_neighbors=min(n_neighbors, working.n_obs - 1), random_state=random_state
    )
    sc.tl.umap(working, random_state=random_state)

    adata.obsm["X_pca"] = working.obsm["X_pca"]
    adata.obsm["X_umap"] = working.obsm["X_umap"]
    adata.obsp["connectivities"] = working.obsp["connectivities"]
    adata.obsp["distances"] = working.obsp["distances"]
    adata.uns["neighbors"] = working.uns["neighbors"]
    adata.uns["hvg"] = working.uns["hvg"]
    adata.var["highly_variable"] = adata.var_names.isin(working.var_names)
    return adata


def _silhouette(
    embedding: np.ndarray, labels: np.ndarray, *, sample_size: int, random_state: int
) -> float:
    n = embedding.shape[0]
    if len(set(labels)) < 2:
        return float("nan")
    effective_sample = min(sample_size, n)
    return float(
        silhouette_score(embedding, labels, sample_size=effective_sample, random_state=random_state)
    )


def sweep_leiden_resolutions(
    adata: ad.AnnData,
    *,
    resolutions: tuple[float, ...] = DEFAULT_RESOLUTIONS,
    silhouette_sample_size: int = DEFAULT_SILHOUETTE_SAMPLE_SIZE,
    random_state: int = 0,
) -> pd.DataFrame:
    """Run Leiden at each candidate resolution and score with silhouette on
    ``obsm['X_pca']`` (subsampled for tractability at scale — see module docstring).
    Does not mutate ``adata``. Returns a DataFrame with one row per resolution:
    ``resolution``, ``n_clusters``, ``silhouette``."""
    if "X_pca" not in adata.obsm:
        raise ValueError("Missing adata.obsm['X_pca']; run preprocess_for_clustering() first.")

    rows = []
    for resolution in resolutions:
        working = adata.copy()
        sc.tl.leiden(
            working,
            resolution=resolution,
            random_state=random_state,
            key_added="_leiden_sweep",
            flavor="igraph",
            n_iterations=2,
        )
        labels = working.obs["_leiden_sweep"].to_numpy()
        score = _silhouette(
            adata.obsm["X_pca"],
            labels,
            sample_size=silhouette_sample_size,
            random_state=random_state,
        )
        rows.append(
            {
                "resolution": resolution,
                "n_clusters": len(set(labels)),
                "silhouette": score,
            }
        )
    return pd.DataFrame(rows)


def run_leiden(
    adata: ad.AnnData, *, resolution: float, random_state: int = 0, obs_key: str = "leiden"
) -> ad.AnnData:
    """Run Leiden at a single, already-chosen resolution, in place."""
    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=random_state,
        key_added=obs_key,
        flavor="igraph",
        n_iterations=2,
    )
    return adata


def choose_best_resolution(sweep_results: pd.DataFrame) -> float:
    """Pick the resolution with the highest silhouette score (NaN-safe)."""
    valid = sweep_results.dropna(subset=["silhouette"])
    if valid.empty:
        raise ValueError("No resolution in the sweep produced a valid silhouette score.")
    best_position = int(np.argmax(valid["silhouette"].to_numpy()))
    resolutions = valid["resolution"].to_numpy()
    return float(resolutions[best_position])
