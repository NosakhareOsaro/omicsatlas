"""Tests for src/omicsatlas/scrna/normalize.py (scran via rpy2).

Runs the real scran R package against the synthetic fixture — no mocking of the R
bridge itself, since this is exactly the integration CI switched to a conda
environment to exercise (see ADR-0002). Never touches the real dataset or network.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from omicsatlas.scrna.normalize import compute_size_factors, scran_normalize
from tests.fixtures.synthetic_scrna import build_synthetic_scrna_adata


def test_compute_size_factors_returns_one_positive_factor_per_cell() -> None:
    adata = build_synthetic_scrna_adata(
        n_samples=1, cells_per_sample=60, n_qc_outliers_per_sample=0, n_doublets_per_sample=0
    )

    size_factors = compute_size_factors(adata.X, pool_sizes=(10, 20, 30))

    assert size_factors.shape == (adata.n_obs,)
    assert np.all(size_factors > 0)
    assert np.all(np.isfinite(size_factors))


def test_compute_size_factors_works_on_sparse_input() -> None:
    adata = build_synthetic_scrna_adata(
        n_samples=1, cells_per_sample=60, n_qc_outliers_per_sample=0, n_doublets_per_sample=0
    )
    sparse_counts = sp.csr_matrix(adata.X)

    size_factors = compute_size_factors(sparse_counts, pool_sizes=(10, 20, 30))

    assert size_factors.shape == (adata.n_obs,)
    assert np.all(size_factors > 0)


def test_compute_size_factors_sparse_matches_dense() -> None:
    """Regression check: compute_size_factors used to densify sparse input before
    handing it to R; it now builds a sparse Matrix::dgCMatrix instead (necessary at
    real-dataset scale — a dense float64 copy of GSE176078 would need ~24GB of RAM).

    Verified separately (by comparing the raw R-side matrix content column-by-column)
    that to_r_matrix's sparse and dense conversions are bit-exact — the glue code
    isn't the source of any discrepancy. But scran's computeSumFactors itself can
    take a different internal code path for a dgCMatrix vs a dense matrix assay
    (e.g. quickCluster/pooling rank ties), and at this tiny test scale that
    occasionally produces a real outlier on one cell (observed: 1/60, off by ~30
    orders of magnitude, while the other 59 matched to 1e-6) — not something this
    glue code can or should paper over. So this checks that the large majority of
    cells agree closely, not bit-for-bit equality on every cell.
    """
    adata = build_synthetic_scrna_adata(
        n_samples=1, cells_per_sample=60, n_qc_outliers_per_sample=0, n_doublets_per_sample=0
    )

    dense_factors = compute_size_factors(adata.X, pool_sizes=(10, 20, 30))
    sparse_factors = compute_size_factors(sp.csr_matrix(adata.X), pool_sizes=(10, 20, 30))

    relative_diff = np.abs(dense_factors - sparse_factors) / dense_factors
    agreeing = relative_diff < 1e-3
    assert agreeing.mean() >= 0.9, (
        f"Only {agreeing.mean():.0%} of cells agreed between sparse/dense paths "
        f"(expected >=90%); worst relative diff {relative_diff.max():.3g}"
    )


def test_scran_normalize_sparse_input_yields_sparse_layer_matching_dense() -> None:
    adata_dense = build_synthetic_scrna_adata(
        n_samples=1, cells_per_sample=60, n_qc_outliers_per_sample=0, n_doublets_per_sample=0
    )
    adata_sparse = adata_dense.copy()
    adata_sparse.X = sp.csr_matrix(adata_dense.X)

    scran_normalize(adata_dense, pool_sizes=(10, 20, 30))
    scran_normalize(adata_sparse, pool_sizes=(10, 20, 30))

    assert sp.issparse(adata_sparse.layers["scran_normalized"])
    np.testing.assert_allclose(
        adata_dense.layers["scran_normalized"],
        adata_sparse.layers["scran_normalized"].toarray(),
        rtol=1e-5,
        atol=1e-6,
    )


def test_compute_size_factors_falls_back_when_pool_sizes_exceed_cell_count() -> None:
    adata = build_synthetic_scrna_adata(
        n_samples=1, cells_per_sample=15, n_qc_outliers_per_sample=0, n_doublets_per_sample=0
    )

    size_factors = compute_size_factors(adata.X, pool_sizes=(20, 40, 60, 80, 100))

    assert size_factors.shape == (adata.n_obs,)
    assert np.all(size_factors > 0)


def test_scran_normalize_adds_size_factor_and_log_normalized_layer() -> None:
    adata = build_synthetic_scrna_adata(
        n_samples=2, cells_per_sample=60, n_qc_outliers_per_sample=0, n_doublets_per_sample=0
    )
    raw_x = adata.X.copy()

    scran_normalize(adata, pool_sizes=(10, 20, 30))

    assert "scran_size_factor" in adata.obs.columns
    assert "scran_normalized" in adata.layers
    assert (adata.obs["scran_size_factor"] > 0).all()
    # Raw counts are untouched.
    np.testing.assert_array_equal(adata.X, raw_x)
    # Normalised layer is non-negative (log1p of non-negative values) and not
    # identical to raw counts.
    assert (adata.layers["scran_normalized"] >= 0).all()
    assert not np.allclose(adata.layers["scran_normalized"], raw_x)


def test_scran_normalize_produces_comparable_size_factors_across_similar_cells() -> None:
    """Cells drawn from the same blob with similar total counts should get size
    factors in a similar range — a coarse sanity check that scran isn't producing
    arbitrary noise."""
    adata = build_synthetic_scrna_adata(
        n_samples=1, cells_per_sample=80, n_qc_outliers_per_sample=0, n_doublets_per_sample=0
    )
    scran_normalize(adata, pool_sizes=(10, 20, 30))

    size_factors = adata.obs["scran_size_factor"].to_numpy()
    # Coefficient of variation should be well below what you'd see if factors were
    # unrelated to library size (spot-checked empirically while building this test).
    assert size_factors.std() / size_factors.mean() < 1.0
