"""scran-based normalisation via an rpy2 bridge to R/Bioconductor.

Uses ``scran::computeSumFactors`` (deconvolution-based size factor estimation) rather
than scanpy's simple total-count normalisation — see ADR-0002 for why this stays on
the R bridge instead of a Python port. Only the size-factor estimation happens in R;
the actual log-normalisation is done in Python for simplicity and to keep the R call
surface minimal.
"""

from __future__ import annotations

from typing import TypeAlias

import anndata as ad
import numpy as np
import scipy.sparse as sp

DEFAULT_POOL_SIZES = (20, 40, 60, 80, 100)

CountsMatrix: TypeAlias = np.ndarray | sp.spmatrix


def compute_size_factors(
    counts: CountsMatrix,
    *,
    pool_sizes: tuple[int, ...] = DEFAULT_POOL_SIZES,
    min_mean: float = 0.1,
) -> np.ndarray:
    """Compute per-cell size factors via scran's deconvolution method.

    ``counts`` is a cells-by-genes matrix (AnnData convention); it is transposed
    before being handed to R, since Bioconductor's convention is genes-by-cells.
    ``pool_sizes`` must not include any pool larger than the number of cells — callers
    working with small (e.g. fixture-scale) data should pass smaller pool sizes.
    """
    import rpy2.robjects as robjects
    from rpy2.robjects import numpy2ri
    from rpy2.robjects.packages import importr

    scran = importr("scran")
    single_cell_experiment = importr("SingleCellExperiment")
    biocgenerics = importr("BiocGenerics")

    n_cells = counts.shape[0]
    usable_pool_sizes = tuple(s for s in pool_sizes if s <= n_cells)
    if not usable_pool_sizes:
        usable_pool_sizes = (min(n_cells, max(pool_sizes[0], 2)),)

    dense_gene_by_cell = np.asarray(
        counts.T.toarray() if sp.issparse(counts) else counts.T,  # type: ignore[union-attr]
        dtype=np.float64,
    )

    with (robjects.default_converter + numpy2ri.converter).context():
        r_matrix = robjects.conversion.get_conversion().py2rpy(dense_gene_by_cell)
        # computeSumFactors dispatches on assay(x, "counts"), which requires a real
        # SingleCellExperiment — a bare matrix isn't accepted in current scran (found
        # empirically: "unable to find an inherited method for function 'assay'").
        sce = single_cell_experiment.SingleCellExperiment(
            assays=robjects.ListVector({"counts": r_matrix})
        )
        sce = scran.computeSumFactors(
            sce,
            sizes=robjects.IntVector(usable_pool_sizes),
            **{"min.mean": min_mean},
        )
        size_factors = np.asarray(biocgenerics.sizeFactors(sce))

    return size_factors


def scran_normalize(
    adata: ad.AnnData,
    *,
    pool_sizes: tuple[int, ...] = DEFAULT_POOL_SIZES,
    min_mean: float = 0.1,
    layer_key: str = "scran_normalized",
) -> ad.AnnData:
    """Add scran size factors (``obs['scran_size_factor']``) and a log-normalised
    layer (``layers[layer_key]``) to ``adata``, in place. Returns ``adata`` for
    chaining. Raw counts in ``adata.X`` are left untouched."""
    counts = adata.X
    size_factors = compute_size_factors(counts, pool_sizes=pool_sizes, min_mean=min_mean)

    if np.any(size_factors <= 0) or np.any(~np.isfinite(size_factors)):
        raise ValueError(
            "scran returned non-positive or non-finite size factors for at least one "
            "cell — usually a sign the input has cells with near-zero total counts "
            "that should have been removed by QC first."
        )

    dense_counts = counts.toarray() if sp.issparse(counts) else np.asarray(counts)
    normalized = dense_counts / size_factors[:, None]
    log_normalized = np.log1p(normalized)

    adata.obs["scran_size_factor"] = size_factors
    adata.layers[layer_key] = log_normalized.astype(np.float32)
    return adata
