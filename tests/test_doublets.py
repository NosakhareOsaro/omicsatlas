"""Tests for src/omicsatlas/scrna/doublets.py.

Scrublet is a KNN/PCA-based classifier that needs a few hundred cells per sample
before it has any real statistical power — at the project's standard ~40-cell/sample
fixture size, its doublet scores are close to noise (verified empirically while
building this module). So this file splits its assertions: fast structural/schema
checks run against the standard small fixture (matching the rest of the suite's
scale), and a single dedicated test verifies actual detection quality against a
larger, still-synthetic, still-network-free fixture sized so Scrublet has enough
signal to work with.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from omicsatlas.scrna.doublets import expected_doublet_rate, run_scrublet_per_sample
from tests.fixtures.synthetic_scrna import build_synthetic_scrna_adata


def test_expected_doublet_rate_scales_with_cell_count() -> None:
    assert expected_doublet_rate(1000) == pytest.approx(0.008)
    assert expected_doublet_rate(2000) == pytest.approx(0.016)
    assert expected_doublet_rate(500) == pytest.approx(0.004)


def test_expected_doublet_rate_is_capped_and_floored() -> None:
    assert expected_doublet_rate(100_000) == 0.20
    assert expected_doublet_rate(0) >= 1e-4


def test_run_scrublet_per_sample_adds_expected_columns() -> None:
    adata = build_synthetic_scrna_adata(
        n_samples=2, cells_per_sample=40, n_doublets_per_sample=3, n_qc_outliers_per_sample=0
    )
    run_scrublet_per_sample(adata)

    assert "doublet_score" in adata.obs.columns
    assert "predicted_doublet" in adata.obs.columns
    assert "doublet_threshold_fallback_used" in adata.obs.columns
    assert adata.obs["doublet_score"].notna().all()
    assert adata.obs["predicted_doublet"].dtype == bool
    assert np.isfinite(adata.obs["doublet_score"]).all()
    assert (adata.obs["doublet_score"] >= 0).all()


def test_run_scrublet_per_sample_preserves_cell_order_and_count() -> None:
    adata = build_synthetic_scrna_adata(
        n_samples=2, cells_per_sample=30, n_doublets_per_sample=2, n_qc_outliers_per_sample=0
    )
    original_index = adata.obs.index.copy()
    original_n_obs = adata.n_obs

    run_scrublet_per_sample(adata)

    assert adata.n_obs == original_n_obs
    assert list(adata.obs.index) == list(original_index)


def test_run_scrublet_per_sample_works_on_sparse_input() -> None:
    """Regression check: a bare np.asarray() on adata.X[positions] used to silently
    mishandle scipy sparse input (wraps it in a 0-d object array rather than
    densifying, so counts_matrix.shape[0] raised IndexError) — found running against
    the real, sparse GSE176078 data. Scrublet accepts sparse input natively; the fix
    was to stop forcing it through np.asarray at all."""
    adata = build_synthetic_scrna_adata(
        n_samples=2, cells_per_sample=40, n_doublets_per_sample=3, n_qc_outliers_per_sample=0
    )
    adata.X = sp.csr_matrix(adata.X)

    run_scrublet_per_sample(adata)

    assert adata.obs["doublet_score"].notna().all()
    assert adata.obs["predicted_doublet"].dtype == bool


def test_run_scrublet_detects_injected_doublets_at_sufficient_scale() -> None:
    """At ~500 cells (well within Scrublet's intended operating range), injected
    doublets should score and get flagged at a materially higher rate than normal
    cells. This is the module's real correctness check; see module docstring for why
    it doesn't run at the project's standard tiny fixture size."""
    adata = build_synthetic_scrna_adata(
        n_samples=1,
        cells_per_sample=550,
        n_genes=1500,
        n_doublets_per_sample=30,
        n_qc_outliers_per_sample=0,
        seed=1,
    )

    run_scrublet_per_sample(adata)

    is_doublet = adata.obs["_qc_case"] == "doublet"
    is_normal = adata.obs["_qc_case"] == "normal"
    doublet_scores = adata.obs.loc[is_doublet, "doublet_score"]
    normal_scores = adata.obs.loc[is_normal, "doublet_score"]
    doublet_call_rate = adata.obs.loc[is_doublet, "predicted_doublet"].mean()
    normal_call_rate = adata.obs.loc[is_normal, "predicted_doublet"].mean()

    assert doublet_scores.mean() > normal_scores.mean()
    assert doublet_call_rate > normal_call_rate
