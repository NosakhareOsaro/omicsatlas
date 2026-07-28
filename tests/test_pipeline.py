"""Tests for src/omicsatlas/scrna/pipeline.py.

``run_pipeline_on_fixture`` is run for real at a small scale — this is exactly the
"schema contract, no real data" test the fixture path exists for (see ADR-0002 /
DATA_SCHEMA.md). Never touches GSE176078 or the network (the placeholder reference
path, not celldex, is used).
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy.io import mmwrite
from scipy.sparse import random as sparse_random

import omicsatlas.scrna.pipeline as pipeline_module
from omicsatlas.scrna.pipeline import load_gse176078, run_pipeline_on_fixture


def _write_fake_geo_layout(raw_dir: Path, n_genes: int = 50, n_cells: int = 20) -> None:
    raw_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    matrix = sparse_random(n_genes, n_cells, density=0.3, random_state=rng, dtype=np.float64)
    matrix.data = np.round(matrix.data * 10) + 1  # positive integer-like counts
    mmwrite(raw_dir / "count_matrix_sparse.mtx", matrix)

    genes = [f"GENE{i:03d}" for i in range(n_genes)]
    (raw_dir / "count_matrix_genes.tsv").write_text("\n".join(genes) + "\n")

    samples = ["CID001", "CID002"]
    barcodes = [f"{samples[i % 2]}_cell{i:03d}" for i in range(n_cells)]
    (raw_dir / "count_matrix_barcodes.tsv").write_text("\n".join(barcodes) + "\n")

    metadata_lines = [",orig.ident,subtype,celltype_major"]
    for i, barcode in enumerate(barcodes):
        sample = samples[i % 2]
        subtype = "ER+" if sample == "CID001" else "TNBC"
        metadata_lines.append(f"{barcode},{sample},{subtype},Epithelial")
    (raw_dir / "metadata.csv").write_text("\n".join(metadata_lines) + "\n")


def test_load_gse176078_parses_geo_layout(tmp_path: Path) -> None:
    raw_dir = tmp_path / "Wu_etal_2021_BRCA_scRNASeq"
    _write_fake_geo_layout(raw_dir, n_genes=50, n_cells=20)

    adata = load_gse176078(raw_dir)

    assert adata.shape == (20, 50)
    assert "orig.ident" in adata.obs.columns
    assert "subtype" in adata.obs.columns
    assert set(adata.obs["orig.ident"].unique()) == {"CID001", "CID002"}


def test_run_pipeline_on_fixture_produces_valid_schema(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture_out = tmp_path / "fixture_signature.h5ad"
    monkeypatch.setattr(pipeline_module, "signature_path", lambda *a, **kw: fixture_out)

    # n_genes must comfortably clear qc.py's 200-genes/cell floor (calibrated for
    # real gene panels, not this test's scale) or every cell gets QC-flagged —
    # found empirically.
    result = run_pipeline_on_fixture(n_samples=2, cells_per_sample=60, n_genes=600)

    assert fixture_out.exists()
    reloaded = ad.read_h5ad(fixture_out)

    for col in ("leiden", "singler_label", "scran_size_factor", "orig.ident", "subtype"):
        assert col in reloaded.obs.columns
    assert "X_pca" in reloaded.obsm
    assert "X_umap" in reloaded.obsm
    assert "scran_normalized" in reloaded.layers
    assert "leiden_chosen_resolution" in reloaded.uns
    assert reloaded.n_obs == result.n_obs
