"""Tests for src/omicsatlas/scrna/velocity.py.

Runs against the synthetic fixture's injected spliced/unspliced layers only — scVelo
is not scientifically applicable to GSE176078 (see ADR-0002 and the module docstring).
These tests exist to verify the code path is correct and tested, not that the results
are biologically meaningful.
"""

from __future__ import annotations

from omicsatlas.scrna.velocity import run_scvelo
from tests.fixtures.synthetic_scrna import build_synthetic_scrna_adata


def test_run_scvelo_requires_spliced_unspliced_layers() -> None:
    adata = build_synthetic_scrna_adata(n_samples=1, cells_per_sample=30)

    try:
        run_scvelo(adata)
        raised = False
    except ValueError:
        raised = True
    assert raised


def test_run_scvelo_adds_velocity_layer_and_graph() -> None:
    adata = build_synthetic_scrna_adata(
        n_samples=2,
        cells_per_sample=60,
        n_genes=300,
        n_qc_outliers_per_sample=0,
        n_doublets_per_sample=0,
        include_velocity_layers=True,
    )

    run_scvelo(adata, n_pcs=10, n_neighbors=15)

    assert "velocity" in adata.layers
    assert adata.layers["velocity"].shape[0] == adata.n_obs
    assert "velocity_graph" in adata.uns or "velocity_graph" in adata.obsp
