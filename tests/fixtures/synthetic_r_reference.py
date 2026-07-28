"""Builds a tiny synthetic R SummarizedExperiment for SingleR tests.

Never downloads the real celldex reference — CI must stay network-free (see
ADR-0002). This builds a minimal, in-memory reference with a handful of synthetic
"cell type" profiles distinguished by the same signal-gene structure the AnnData
fixture's expression blobs use, so tests can check that SingleR actually recovers the
right label for cells generated from a given blob, not just that the code path runs.
"""

from __future__ import annotations

import numpy as np

from omicsatlas.scrna.annotate import build_summarized_experiment_reference


def build_synthetic_singler_reference(
    *,
    gene_names: list[str],
    n_blobs: int,
    blob_signal_genes: list[np.ndarray],
    blob_log_means: np.ndarray,
    baseline_log_mean: float,
    labels: list[str] | None = None,
    profiles_per_blob: int = 5,
    seed: int = 0,
):
    """Return an R SummarizedExperiment with a ``logcounts`` assay and
    ``colData$label.main``, built from ``profiles_per_blob`` synthetic pseudo-bulk
    profiles per blob (log1p of the same per-blob signal-gene means the AnnData
    fixture uses, with a little noise so it isn't a degenerate single point). Uses
    ``annotate.build_summarized_experiment_reference`` for the actual R object
    construction — the same low-level constructor the production placeholder
    reference uses, so this test double stays representative.
    """
    labels = labels or [f"synthetic_type_{i}" for i in range(n_blobs)]
    rng = np.random.default_rng(seed)

    n_genes = len(gene_names)
    profiles = []
    profile_labels = []
    for blob in range(n_blobs):
        log_means = np.full(n_genes, baseline_log_mean)
        log_means[blob_signal_genes[blob]] = blob_log_means[blob, blob_signal_genes[blob]]
        for _ in range(profiles_per_blob):
            noisy = log_means + rng.normal(0, 0.05, size=n_genes)
            profiles.append(np.log1p(np.exp(noisy)))
            profile_labels.append(labels[blob])

    ref_matrix_gene_by_sample = np.asarray(profiles).T  # genes x ref-samples

    reference = build_summarized_experiment_reference(
        ref_matrix_gene_by_sample, gene_names, profile_labels
    )

    return reference, labels
