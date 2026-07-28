"""scVelo RNA velocity — implemented for pipeline completeness, not scientifically
applicable to GSE176078.

RNA velocity needs spliced/unspliced count layers (from velocyto/STARsolo against raw
FASTQs). GSE176078's processed download is a standard count matrix with no such
layers, and the raw FASTQs are EGA-controlled-access. This module is tested against a
synthetic fixture with injected spliced/unspliced layers only — see ADR-0002. It is
not run against the real dataset and its output is not part of the signature artifact.
"""

from __future__ import annotations

import anndata as ad


def run_scvelo(
    adata: ad.AnnData,
    *,
    n_pcs: int = 10,
    n_neighbors: int = 15,
    mode: str = "deterministic",
) -> ad.AnnData:
    """Run scVelo's standard preprocessing + velocity estimation pipeline in place.

    Requires ``layers['spliced']`` and ``layers['unspliced']``. Adds
    ``layers['velocity']`` and ``obsp['velocity_graph']``. Returns ``adata`` for
    chaining.
    """
    import scvelo as scv

    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise ValueError(
            "Missing 'spliced'/'unspliced' layers — scVelo cannot run without them. "
            "GSE176078's processed matrix has no such layers; see ADR-0002."
        )

    scv.pp.filter_and_normalize(adata, min_shared_counts=1)
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
    scv.tl.velocity(adata, mode=mode)
    scv.tl.velocity_graph(adata)
    return adata
