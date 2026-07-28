"""Tests for src/omicsatlas/browser/app.py's pure helper functions.

Streamlit's own run loop (main()) isn't unit-tested — see module docstring. These
tests cover load_adata/available_color_options/build_umap_figure directly against a
small, directly-constructed AnnData.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure

from omicsatlas.browser.app import available_color_options, build_umap_figure, load_adata


def _small_adata() -> ad.AnnData:
    n_cells = 30
    rng = np.random.default_rng(0)
    obs = pd.DataFrame(
        {
            "leiden": pd.Categorical([str(i % 3) for i in range(n_cells)]),
            "subtype": pd.Categorical(["ER+", "TNBC"] * (n_cells // 2)),
        }
    )
    adata = ad.AnnData(X=rng.random((n_cells, 5)), obs=obs)
    adata.obsm["X_umap"] = rng.random((n_cells, 2))
    return adata


def test_load_adata_roundtrips_through_h5ad(tmp_path: Path) -> None:
    adata = _small_adata()
    path = tmp_path / "test.h5ad"
    adata.write_h5ad(path)

    loaded = load_adata(path)

    assert loaded.shape == adata.shape
    assert "X_umap" in loaded.obsm


def test_available_color_options_filters_to_present_columns() -> None:
    adata = _small_adata()

    options = available_color_options(adata)

    assert options == ["leiden", "subtype"]


def test_available_color_options_empty_when_none_present() -> None:
    adata = ad.AnnData(X=np.zeros((5, 3)))

    assert available_color_options(adata) == []


def test_build_umap_figure_returns_figure_with_legend_per_category() -> None:
    adata = _small_adata()

    fig = build_umap_figure(adata, "leiden")

    assert isinstance(fig, Figure)
    legend = fig.axes[0].get_legend()
    assert legend is not None
    assert len(legend.get_texts()) == 3  # 3 leiden categories


def test_build_umap_figure_requires_umap_embedding() -> None:
    adata = _small_adata()
    del adata.obsm["X_umap"]

    with pytest.raises(ValueError):
        build_umap_figure(adata, "leiden")


def test_build_umap_figure_requires_known_column() -> None:
    adata = _small_adata()

    with pytest.raises(ValueError):
        build_umap_figure(adata, "not_a_real_column")
