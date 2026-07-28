"""scran-based normalisation via an rpy2 bridge to R/Bioconductor.

Uses ``scran::computeSumFactors`` (deconvolution-based size factor estimation) rather
than scanpy's simple total-count normalisation — see ADR-0002 for why this stays on
the R bridge instead of a Python port. Only the size-factor estimation happens in R;
the actual log-normalisation is done in Python for simplicity and to keep the R call
surface minimal.

Sparse input is kept sparse throughout, on both the R and Python sides. GSE176078 is
~100k cells x ~30k genes at ~6% density; densifying that (as an earlier version of
this module did, validated only against small dense test fixtures) would need tens of
GB of RAM for a float64 copy — found this the hard way running against the real data.
"""

from __future__ import annotations

from typing import Any, TypeAlias

import anndata as ad
import numpy as np
import scipy.sparse as sp

DEFAULT_POOL_SIZES = (20, 40, 60, 80, 100)

CountsMatrix: TypeAlias = np.ndarray | sp.spmatrix


def to_r_matrix(matrix: CountsMatrix) -> Any:
    """Convert a genes-by-cells matrix to an R matrix, keeping it sparse
    (``Matrix::dgCMatrix``) when the input is sparse rather than densifying."""
    import rpy2.robjects as robjects
    from rpy2.robjects import numpy2ri
    from rpy2.robjects.packages import importr

    if sp.issparse(matrix):
        matrix_pkg = importr("Matrix")
        csc = matrix.tocsc()  # type: ignore[union-attr]
        # robjects.IntVector/FloatVector convert numpy arrays directly (fast — no
        # Python-level element loop) and produce a plain R atomic vector. Routing
        # this through the numpy2ri conversion context instead (as an earlier version
        # did) produced an R object with a `dim` attribute, giving it class "array"
        # rather than "numeric" — Matrix::sparseMatrix's internal slot assignment
        # rejects that for the `x` slot ("is(value, 'numeric') is not TRUE"), found
        # empirically.
        i_r = robjects.IntVector(csc.indices + 1)  # R is 1-indexed
        p_r = robjects.IntVector(csc.indptr)
        x_r = robjects.FloatVector(csc.data)
        return matrix_pkg.sparseMatrix(i=i_r, p=p_r, x=x_r, dims=robjects.IntVector(csc.shape))

    dense = np.asarray(matrix, dtype=np.float64)
    with (robjects.default_converter + numpy2ri.converter).context():
        return robjects.conversion.get_conversion().py2rpy(dense)


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
    from rpy2.robjects.packages import importr

    scran = importr("scran")
    single_cell_experiment = importr("SingleCellExperiment")
    biocgenerics = importr("BiocGenerics")

    n_cells = counts.shape[0]
    usable_pool_sizes = tuple(s for s in pool_sizes if s <= n_cells)
    if not usable_pool_sizes:
        usable_pool_sizes = (min(n_cells, max(pool_sizes[0], 2)),)

    r_matrix = to_r_matrix(counts.T)
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
    return np.asarray(biocgenerics.sizeFactors(sce))


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

    if sp.issparse(counts):
        counts_csr = counts.tocsr().astype(np.float64)
        normalized = sp.diags(1.0 / size_factors) @ counts_csr
        # log1p(0) == 0, so applying it only to the stored (nonzero) values is exact
        # and keeps the layer sparse — avoids materialising a dense ~100k x ~30k array.
        normalized.data = np.log1p(normalized.data)
        log_normalized: np.ndarray | sp.spmatrix = normalized.astype(np.float32)
    else:
        dense_counts = np.asarray(counts, dtype=np.float64)
        log_normalized = np.log1p(dense_counts / size_factors[:, None]).astype(np.float32)

    adata.obs["scran_size_factor"] = size_factors
    adata.layers[layer_key] = log_normalized
    return adata
