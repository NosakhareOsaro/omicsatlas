"""Tests for src/omicsatlas/scrna/cluster.py against the synthetic fixture."""

from __future__ import annotations

from omicsatlas.scrna.cluster import (
    choose_best_resolution,
    preprocess_for_clustering,
    run_leiden,
    sweep_leiden_resolutions,
)
from omicsatlas.scrna.normalize import scran_normalize
from tests.fixtures.synthetic_scrna import build_synthetic_scrna_adata


def _normalized_fixture():
    adata = build_synthetic_scrna_adata(
        n_samples=3, cells_per_sample=60, n_qc_outliers_per_sample=0, n_doublets_per_sample=0
    )
    scran_normalize(adata, pool_sizes=(10, 20, 30))
    return adata


def test_preprocess_for_clustering_adds_embeddings() -> None:
    adata = _normalized_fixture()

    preprocess_for_clustering(adata, n_top_genes=100, n_pcs=10)

    assert "X_pca" in adata.obsm
    assert "X_umap" in adata.obsm
    assert adata.obsm["X_pca"].shape[0] == adata.n_obs
    assert adata.obsm["X_umap"].shape == (adata.n_obs, 2)
    assert "connectivities" in adata.obsp
    assert adata.var["highly_variable"].sum() > 0


def test_preprocess_for_clustering_requires_normalized_layer() -> None:
    adata = build_synthetic_scrna_adata(n_samples=1, cells_per_sample=30)
    try:
        preprocess_for_clustering(adata)
        raised = False
    except ValueError:
        raised = True
    assert raised


def test_sweep_leiden_resolutions_returns_one_row_per_resolution() -> None:
    adata = _normalized_fixture()
    preprocess_for_clustering(adata, n_top_genes=100, n_pcs=10)

    sweep = sweep_leiden_resolutions(adata, resolutions=(0.2, 0.5, 0.8), silhouette_sample_size=100)

    assert list(sweep["resolution"]) == [0.2, 0.5, 0.8]
    assert (sweep["n_clusters"] >= 1).all()
    assert sweep["silhouette"].notna().any()


def test_choose_best_resolution_picks_max_silhouette() -> None:
    adata = _normalized_fixture()
    preprocess_for_clustering(adata, n_top_genes=100, n_pcs=10)
    sweep = sweep_leiden_resolutions(adata, resolutions=(0.2, 0.5, 0.8), silhouette_sample_size=100)

    best = choose_best_resolution(sweep)

    assert best == sweep.loc[sweep["silhouette"].idxmax(), "resolution"]


def test_run_leiden_assigns_cluster_labels_to_every_cell() -> None:
    adata = _normalized_fixture()
    preprocess_for_clustering(adata, n_top_genes=100, n_pcs=10)

    run_leiden(adata, resolution=0.5)

    assert "leiden" in adata.obs.columns
    assert adata.obs["leiden"].notna().all()
    assert adata.obs["leiden"].nunique() >= 1


def test_clustering_recovers_blob_structure_reasonably_well() -> None:
    """The fixture's 3 expression blobs should be at least partially recoverable —
    not a strict ARI assertion (Leiden cluster count won't necessarily match blob
    count exactly), but resolution 0.5 shouldn't collapse everything into one cluster
    given how separable the fixture's blobs are."""
    adata = _normalized_fixture()
    preprocess_for_clustering(adata, n_top_genes=100, n_pcs=10)
    run_leiden(adata, resolution=0.5)

    assert adata.obs["leiden"].nunique() > 1
