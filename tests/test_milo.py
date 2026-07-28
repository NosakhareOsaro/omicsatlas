"""Tests for src/omicsatlas/scrna/milo.py (Milo differential abundance via pertpy).

Runs the real pertpy/pydeseq2 pipeline against the synthetic fixture — no mocking.
Milo has a real grouping variable to test here (clinical subtype), unlike scVelo; see
ADR-0002.
"""

from __future__ import annotations

from omicsatlas.scrna.cluster import preprocess_for_clustering
from omicsatlas.scrna.milo import run_milo_differential_abundance
from omicsatlas.scrna.normalize import scran_normalize
from tests.fixtures.synthetic_scrna import build_synthetic_scrna_adata


def _clustered_fixture(seed: int = 5):
    adata = build_synthetic_scrna_adata(
        n_samples=6,
        cells_per_sample=80,
        n_genes=300,
        n_qc_outliers_per_sample=0,
        n_doublets_per_sample=0,
        seed=seed,
    )
    scran_normalize(adata, pool_sizes=(10, 20, 30))
    preprocess_for_clustering(adata, n_top_genes=100, n_pcs=10)
    return adata


def test_run_milo_differential_abundance_requires_pca() -> None:
    adata = build_synthetic_scrna_adata(n_samples=2, cells_per_sample=30)

    try:
        run_milo_differential_abundance(adata)
        raised = False
    except ValueError:
        raised = True
    assert raised


def test_run_milo_differential_abundance_returns_expected_structure() -> None:
    adata = _clustered_fixture()

    mdata = run_milo_differential_abundance(adata, n_neighbors=10, nhood_prop=0.2)

    assert "milo" in mdata.mod
    da_results = mdata["milo"].var
    for column in ("logFC", "PValue", "FDR", "SpatialFDR"):
        assert column in da_results.columns
    assert len(da_results) > 0
    assert da_results["PValue"].between(0, 1).all()
    assert da_results["SpatialFDR"].between(0, 1).all()
