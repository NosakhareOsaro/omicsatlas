"""Sanity checks on the synthetic scRNA-seq fixture builder itself.

These aren't testing pipeline logic — they're testing that the fixture used by every
other scrna test actually has the structure its docstring promises (right shape,
doublets genuinely sum their parents, QC-outlier cases are genuinely outliers).
"""

from __future__ import annotations

import numpy as np

from tests.fixtures.synthetic_scrna import MT_GENE_NAMES, build_synthetic_scrna_adata


def test_shape_matches_requested_composition() -> None:
    adata = build_synthetic_scrna_adata(
        n_samples=4,
        cells_per_sample=50,
        n_genes=500,
        n_doublets_per_sample=2,
        n_qc_outliers_per_sample=2,
    )

    expected_cells = 4 * (50 + 2 + 2)
    assert adata.shape == (expected_cells, 500)


def test_mitochondrial_genes_present_and_prefixed() -> None:
    adata = build_synthetic_scrna_adata()

    for gene in MT_GENE_NAMES:
        assert gene in adata.var_names
        assert gene.startswith("MT-")


def test_deterministic_given_seed() -> None:
    a = build_synthetic_scrna_adata(seed=42)
    b = build_synthetic_scrna_adata(seed=42)

    np.testing.assert_array_equal(a.X, b.X)
    assert list(a.obs["_qc_case"]) == list(b.obs["_qc_case"])


def test_qc_outlier_cases_are_genuinely_extreme() -> None:
    adata = build_synthetic_scrna_adata(n_samples=2, cells_per_sample=30)
    total_counts = np.asarray(adata.X.sum(axis=1)).ravel()

    normal_mask = (adata.obs["_qc_case"] == "normal").to_numpy()
    low_count_mask = (adata.obs["_qc_case"] == "low_count_outlier").to_numpy()
    high_mito_mask = (adata.obs["_qc_case"] == "high_mito_outlier").to_numpy()

    assert total_counts[low_count_mask].max() < total_counts[normal_mask].min()

    mt_idx = [adata.var_names.get_loc(g) for g in MT_GENE_NAMES]
    mt_fraction = np.asarray(adata.X[:, mt_idx].sum(axis=1)).ravel() / total_counts
    assert mt_fraction[high_mito_mask].min() > mt_fraction[normal_mask].max()


def test_doublets_are_pairwise_sums_of_two_normal_cells_in_sample() -> None:
    adata = build_synthetic_scrna_adata(
        n_samples=1, cells_per_sample=20, n_doublets_per_sample=3, n_qc_outliers_per_sample=0
    )
    normal = adata[adata.obs["_qc_case"] == "normal"].X
    doublets = adata[adata.obs["_qc_case"] == "doublet"].X

    # Every doublet's total count must exceed any single normal cell's total count —
    # it's a sum of two real cells, not just noisier.
    normal_totals = np.asarray(normal.sum(axis=1)).ravel()
    doublet_totals = np.asarray(doublets.sum(axis=1)).ravel()
    assert doublet_totals.min() > normal_totals.max()


def test_velocity_layers_present_only_when_requested() -> None:
    without = build_synthetic_scrna_adata(include_velocity_layers=False)
    assert "spliced" not in without.layers
    assert "unspliced" not in without.layers

    with_layers = build_synthetic_scrna_adata(include_velocity_layers=True)
    assert "spliced" in with_layers.layers
    assert "unspliced" in with_layers.layers
    np.testing.assert_allclose(
        with_layers.layers["spliced"] + with_layers.layers["unspliced"],
        with_layers.X,
        rtol=1e-5,
    )


def test_subtype_skews_blob_composition_for_milo_signal() -> None:
    """Different subtypes should have measurably different expression-blob composition
    (via the dominant-blob signal genes), giving Milo something real to detect later."""
    adata = build_synthetic_scrna_adata(n_samples=6, cells_per_sample=40, seed=7)

    er_mean = adata[adata.obs["subtype"] == "ER+"].X.mean(axis=0)
    tnbc_mean = adata[adata.obs["subtype"] == "TNBC"].X.mean(axis=0)

    assert not np.allclose(er_mean, tnbc_mean, rtol=0.2)
