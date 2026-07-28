"""Basic Streamlit UMAP browser for the scRNA-seq signature artifact.

Phase 1 scope only — will be merged into the unified multi-modal browser (UMAP +
spatial viewer + ATAC tracks + concordance panel) in Phase 6. Structured with pure,
independently testable functions (``load_adata``, ``build_umap_figure``,
``available_color_options``) kept separate from the thin Streamlit wiring in
``main()``, since a Streamlit app's own run loop can't be meaningfully unit-tested.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from omicsatlas.scrna.artifact import signature_path

DEFAULT_COLOR_OPTIONS = ["leiden", "singler_label", "subtype", "orig.ident"]


def load_adata(path: str | Path) -> ad.AnnData:
    return ad.read_h5ad(path)


def available_color_options(adata: ad.AnnData) -> list[str]:
    """Which of the standard colour-by columns are actually present in ``adata.obs``,
    in display order."""
    return [c for c in DEFAULT_COLOR_OPTIONS if c in adata.obs.columns]


def build_umap_figure(adata: ad.AnnData, color_by: str) -> Figure:
    """Build a UMAP scatter plot coloured by a categorical ``obs`` column."""
    if "X_umap" not in adata.obsm:
        raise ValueError("Missing adata.obsm['X_umap']; run the clustering pipeline first.")
    if color_by not in adata.obs.columns:
        raise ValueError(f"Missing obs column {color_by!r}.")

    coords = adata.obsm["X_umap"]
    categories = adata.obs[color_by].astype("category")
    n_categories = len(categories.cat.categories)
    cmap = plt.get_cmap("tab20", max(n_categories, 1))

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(coords[:, 0], coords[:, 1], c=categories.cat.codes, cmap=cmap, s=5)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(f"Coloured by {color_by}")

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="", color=cmap(i), label=str(category))
        for i, category in enumerate(categories.cat.categories)
    ]
    ax.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc="upper left", fontsize="small")
    fig.tight_layout()
    return fig


def main() -> None:
    import streamlit as st

    st.title("OmicsAtlas — scRNA-seq UMAP browser")
    st.caption(
        "Basic Phase 1 browser; will be merged into the unified multi-modal browser " "in Phase 6."
    )

    use_fixture = st.sidebar.checkbox("Use fixture artifact", value=False)
    path = signature_path(fixture=use_fixture)
    if not path.exists():
        st.error(
            f"No signature artifact found at {path}. Run `make scrna-signature` "
            "(or `make scrna-signature-fixture`) first."
        )
        return

    adata = load_adata(path)
    st.write(f"{adata.n_obs} cells × {adata.n_vars} genes")

    options = available_color_options(adata)
    if not options:
        st.error("No known colour-by columns found in this artifact's obs.")
        return

    color_by = st.sidebar.selectbox("Colour by", options)
    fig = build_umap_figure(adata, color_by)
    st.pyplot(fig)


if __name__ == "__main__":
    main()
