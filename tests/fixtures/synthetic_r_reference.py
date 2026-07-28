"""Builds a tiny synthetic R SummarizedExperiment for SingleR tests.

Never downloads the real celldex reference — CI must stay network-free (see
ADR-0002). This builds a minimal, in-memory reference with a handful of synthetic
"cell type" profiles distinguished by the same signal-gene structure the AnnData
fixture's expression blobs use, so tests can check that SingleR actually recovers the
right label for cells generated from a given blob, not just that the code path runs.
"""

from __future__ import annotations

import numpy as np


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
    fixture uses, with a little noise so it isn't a degenerate single point).
    """
    import rpy2.robjects as robjects
    from rpy2.robjects import numpy2ri
    from rpy2.robjects.packages import importr

    summarized_experiment = importr("SummarizedExperiment")
    s4vectors = importr("S4Vectors")

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

    # See annotate.py: R calls made *inside* the numpy2ri context also get their
    # return values silently converted back to numpy, stripping attributes like
    # rownames. Only the initial py2rpy conversion needs the context.
    with (robjects.default_converter + numpy2ri.converter).context():
        r_matrix_raw = robjects.conversion.get_conversion().py2rpy(ref_matrix_gene_by_sample)

    r_matrix = robjects.r["rownames<-"](r_matrix_raw, robjects.StrVector(gene_names))

    col_data = s4vectors.DataFrame(**{"label.main": robjects.StrVector(profile_labels)})
    reference = summarized_experiment.SummarizedExperiment(
        assays=robjects.ListVector({"logcounts": r_matrix}),
        colData=col_data,
    )

    return reference, labels
