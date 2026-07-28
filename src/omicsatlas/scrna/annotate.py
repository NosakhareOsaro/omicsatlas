"""SingleR cell-type annotation via an rpy2 bridge to R/Bioconductor.

Reference is ``celldex::HumanPrimaryCellAtlasData()`` in production (see ADR-0002 for
why, and for the documented malignant/normal-epithelial limitation this reference
carries). The reference is deliberately an injectable argument, not hardcoded inside
``annotate_with_singler`` — tests pass a tiny synthetic reference instead of
downloading the real celldex panel, keeping CI free of network dependencies (see
constraint on fixture-only tests in ADR-0002 / docs/scrna-pipeline.md). Any R
``SummarizedExperiment`` with a ``logcounts`` assay and a label column in ``colData``
works, real or synthetic.
"""

from __future__ import annotations

from typing import Any

import anndata as ad
import numpy as np

DEFAULT_LABEL_COLUMN = "label.main"


def load_human_primary_cell_atlas() -> Any:
    """Fetch celldex's HumanPrimaryCellAtlasData reference (network access to
    Bioconductor's ExperimentHub required — not used in CI/unit tests, only in the
    real end-to-end pipeline run). Returns an R SummarizedExperiment (untyped: rpy2
    objects have no useful static type)."""
    from rpy2.robjects.packages import importr

    celldex = importr("celldex")
    return celldex.HumanPrimaryCellAtlasData()


def run_singler(
    query_log_norm: np.ndarray,
    query_genes: list[str],
    reference: Any,
    *,
    label_column: str = DEFAULT_LABEL_COLUMN,
) -> np.ndarray:
    """Run SingleR, returning one predicted label per query cell (row of
    ``query_log_norm``, a cells-by-genes matrix). ``reference`` is an R
    SummarizedExperiment with a ``logcounts`` assay and ``colData[[label_column]]``."""
    import rpy2.robjects as robjects
    from rpy2.robjects import numpy2ri
    from rpy2.robjects.packages import importr

    singler = importr("SingleR")
    summarized_experiment = importr("SummarizedExperiment")

    # Only the numpy -> R conversion itself needs the numpy2ri context. Calling R
    # functions *inside* that context also converts their return values back to
    # numpy, silently stripping R attributes like rownames (found empirically: a
    # `rownames<-` call made inside the context returned a plain numpy.ndarray, and
    # SingleR then failed with "'test' must have row names" despite the assignment
    # having "worked"). So every R-side call below runs outside it.
    query_gene_by_cell = np.asarray(query_log_norm.T, dtype=np.float64)
    with (robjects.default_converter + numpy2ri.converter).context():
        r_query_raw = robjects.conversion.get_conversion().py2rpy(query_gene_by_cell)

    r_query = robjects.r["rownames<-"](r_query_raw, robjects.StrVector(query_genes))

    # colData() returns an S4 DFrame, not a plain rpy2 ListVector, so .rx2() isn't
    # available on it (found empirically) — use R's generic `[[` instead, which
    # dispatches correctly on DFrame.
    r_double_bracket = robjects.r["[["]
    ref_labels = r_double_bracket(summarized_experiment.colData(reference), label_column)

    result = singler.SingleR(test=r_query, ref=reference, labels=ref_labels)
    # SingleR's result is also an S4 DFrame (see the colData note above) — same `[[`
    # fix applies.
    labels = np.asarray(r_double_bracket(result, "labels"))

    return labels


def annotate_with_singler(
    adata: ad.AnnData,
    reference: Any,
    *,
    layer_key: str = "scran_normalized",
    label_column: str = DEFAULT_LABEL_COLUMN,
    obs_key: str = "singler_label",
) -> ad.AnnData:
    """Annotate ``adata`` in place with SingleR labels, read from ``layers[layer_key]``
    (expects log-normalised expression, e.g. from ``normalize.scran_normalize``).
    Returns ``adata`` for chaining."""
    if layer_key not in adata.layers:
        raise ValueError(f"Missing layer {layer_key!r}; run scran_normalize() first.")

    labels = run_singler(
        adata.layers[layer_key],
        list(adata.var_names),
        reference,
        label_column=label_column,
    )
    adata.obs[obs_key] = labels
    return adata
