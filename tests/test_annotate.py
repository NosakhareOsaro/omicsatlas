"""Tests for src/omicsatlas/scrna/annotate.py (SingleR via rpy2).

Never downloads the real celldex reference — builds a tiny synthetic R
SummarizedExperiment reference instead (see
tests/fixtures/synthetic_r_reference.py), matching the AnnData fixture's own blob
generative structure so the test can check that SingleR actually recovers the right
label per cell, not just that the code path runs. See ADR-0002.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from omicsatlas.scrna.annotate import annotate_with_singler, run_singler
from omicsatlas.scrna.normalize import scran_normalize
from tests.fixtures.synthetic_r_reference import build_synthetic_singler_reference
from tests.fixtures.synthetic_scrna import build_synthetic_scrna_adata


def _fixture_with_reference(seed: int = 3):
    adata = build_synthetic_scrna_adata(
        n_samples=1,
        cells_per_sample=90,
        n_blobs=3,
        n_qc_outliers_per_sample=0,
        n_doublets_per_sample=0,
        seed=seed,
    )
    scran_normalize(adata, pool_sizes=(10, 20, 30))

    reference, labels = build_synthetic_singler_reference(
        gene_names=list(adata.var_names),
        n_blobs=3,
        blob_signal_genes=adata.uns["_blob_signal_genes"],
        blob_log_means=adata.uns["_blob_log_means"],
        baseline_log_mean=adata.uns["_baseline_log_mean"],
        seed=seed,
    )
    return adata, reference, labels


def test_run_singler_returns_one_label_per_cell() -> None:
    adata, reference, _ = _fixture_with_reference()

    labels = run_singler(adata.layers["scran_normalized"], list(adata.var_names), reference)

    assert len(labels) == adata.n_obs


def test_run_singler_recovers_the_true_blob_label_for_most_cells() -> None:
    adata, reference, labels = _fixture_with_reference()

    predicted = run_singler(adata.layers["scran_normalized"], list(adata.var_names), reference)

    true_blob = adata.obs["_true_blob"].to_numpy()
    predicted_label_names = [labels[b] for b in true_blob]

    accuracy = float((predicted == predicted_label_names).mean())
    # Not asserting near-100%: SingleR's correlation-based scoring on a small
    # synthetic reference won't be perfect, but it should be well above the ~33%
    # chance rate for 3 blobs.
    assert accuracy > 0.7


def test_run_singler_sparse_query_matches_dense() -> None:
    """Regression check: run_singler used to densify the query with a plain
    np.asarray call, which silently mishandles scipy sparse input (it doesn't
    densify — it wraps it in a 0-d object array). Fixed to route through
    normalize.to_r_matrix, which stays sparse. Both must give the same labels."""
    adata, reference, _ = _fixture_with_reference()
    dense_layer = adata.layers["scran_normalized"]
    sparse_layer = sp.csr_matrix(dense_layer)

    dense_labels = run_singler(dense_layer, list(adata.var_names), reference)
    sparse_labels = run_singler(sparse_layer, list(adata.var_names), reference)

    np.testing.assert_array_equal(dense_labels, sparse_labels)


def test_annotate_with_singler_writes_expected_obs_column() -> None:
    adata, reference, _ = _fixture_with_reference()

    annotate_with_singler(adata, reference)

    assert "singler_label" in adata.obs.columns
    assert adata.obs["singler_label"].notna().all()


def test_annotate_with_singler_requires_normalized_layer() -> None:
    adata = build_synthetic_scrna_adata(n_samples=1, cells_per_sample=30)
    _, reference, _ = _fixture_with_reference()

    try:
        annotate_with_singler(adata, reference)
        raised = False
    except ValueError:
        raised = True
    assert raised
