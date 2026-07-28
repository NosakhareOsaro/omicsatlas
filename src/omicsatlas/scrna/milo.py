"""Milo differential neighbourhood abundance, via pertpy's Python-native port.

Tests differential abundance across the three clinical subtypes (ER+/HER2+/TNBC)
recorded in GSE176078's published per-patient metadata — a real, documentable
grouping variable, unlike scVelo (see ADR-0002 for why scVelo isn't meaningful on this
dataset but Milo is). Implemented via ``pertpy.tools.Milo`` rather than an rpy2 wrapper
around R's ``miloR``: scran and SingleR stay on the R bridge because they have no solid
Python equivalent, but Milo does, so an extra R round-trip here would be unjustified.
``da_nhoods`` uses pertpy's ``pydeseq2`` solver specifically to keep this module fully
Python-native — no rpy2 dependency at all.
"""

from __future__ import annotations

from typing import Literal

import anndata as ad
import mudata as mu
import scanpy as sc

DEFAULT_SAMPLE_COL = "orig.ident"
DEFAULT_NHOOD_PROP = 0.1


def run_milo_differential_abundance(
    adata: ad.AnnData,
    *,
    sample_col: str = DEFAULT_SAMPLE_COL,
    design: str = "~subtype",
    model_contrasts: str | None = None,
    use_rep: str = "X_pca",
    n_neighbors: int = 30,
    nhood_prop: float = DEFAULT_NHOOD_PROP,
    solver: Literal["edger", "pydeseq2"] = "pydeseq2",
) -> mu.MuData:
    """Run Milo differential neighbourhood abundance.

    Requires ``adata.obsm[use_rep]`` to already exist (e.g. from ``cluster.py``'s PCA
    step). Returns the MuData object; differential abundance results are in
    ``mdata['milo'].var`` (columns include ``logFC``, ``PValue``, ``SpatialFDR``).
    """
    import pertpy as pt

    if use_rep not in adata.obsm:
        raise ValueError(f"Missing adata.obsm[{use_rep!r}]; run PCA first (see cluster.py).")

    milo = pt.tl.Milo()
    mdata = milo.load(adata)
    sc.pp.neighbors(mdata["rna"], use_rep=use_rep, n_neighbors=n_neighbors)
    milo.make_nhoods(mdata["rna"], prop=nhood_prop)
    mdata = milo.count_nhoods(mdata, sample_col=sample_col)
    milo.da_nhoods(
        mdata,
        design=design,
        model_contrasts=model_contrasts,
        solver=solver,
    )
    return mdata
