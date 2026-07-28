"""End-to-end scRNA-seq pipeline entrypoint: fetch -> QC -> doublets -> normalise ->
cluster -> annotate -> save the versioned signature artifact.

Two entrypoints: :func:`run_pipeline` against real GSE176078 data (``make
scrna-signature``) and :func:`run_pipeline_on_fixture` against the synthetic test
fixture (``make scrna-signature-fixture``), which exercises the identical assembly
code to prove the artifact schema contract end-to-end without the real 100k-cell
download. See ADR-0002 and DATA_SCHEMA.md.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd
import scipy.io as sio

from omicsatlas.scrna.annotate import (
    annotate_with_singler,
    build_placeholder_reference,
    load_human_primary_cell_atlas,
)
from omicsatlas.scrna.artifact import REPO_ROOT, signature_path
from omicsatlas.scrna.cluster import (
    choose_best_resolution,
    preprocess_for_clustering,
    run_leiden,
    sweep_leiden_resolutions,
)
from omicsatlas.scrna.doublets import run_scrublet_per_sample
from omicsatlas.scrna.normalize import scran_normalize
from omicsatlas.scrna.qc import run_qc

GSE176078_EXTRACTED_DIR = "Wu_etal_2021_BRCA_scRNASeq"


def load_gse176078(raw_dir: str | Path) -> ad.AnnData:
    """Load the extracted GSE176078 processed matrix into an AnnData object.

    Expects ``raw_dir`` to contain ``count_matrix_sparse.mtx`` (genes x cells),
    ``count_matrix_genes.tsv``, ``count_matrix_barcodes.tsv``, and ``metadata.csv``
    (the exact layout GEO's supplementary tarball extracts to).
    """
    raw_dir = Path(raw_dir)
    matrix = sio.mmread(raw_dir / "count_matrix_sparse.mtx").T.tocsr()  # -> cells x genes
    genes = pd.read_csv(raw_dir / "count_matrix_genes.tsv", header=None)[0].to_numpy()
    barcodes = pd.read_csv(raw_dir / "count_matrix_barcodes.tsv", header=None)[0].to_numpy()
    metadata = pd.read_csv(raw_dir / "metadata.csv", index_col=0)

    adata = ad.AnnData(
        X=matrix,
        obs=metadata.reindex(barcodes),
        var=pd.DataFrame(index=pd.Index(genes, name="gene_symbol")),
    )
    adata.var_names_make_unique()
    return adata


def _assemble_signature(
    adata: ad.AnnData,
    *,
    resolutions: tuple[float, ...],
    reference_builder: Callable[[ad.AnnData], Any],
) -> ad.AnnData:
    adata = run_qc(adata)
    run_scrublet_per_sample(adata)
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    scran_normalize(adata)
    preprocess_for_clustering(adata)

    sweep = sweep_leiden_resolutions(adata, resolutions=resolutions)
    best_resolution = choose_best_resolution(sweep)
    run_leiden(adata, resolution=best_resolution)
    # anndata's h5ad writer supports a DataFrame in uns directly, but not a list of
    # dicts (found empirically: TypeError writing uns/leiden_resolution_sweep).
    adata.uns["leiden_resolution_sweep"] = sweep
    adata.uns["leiden_chosen_resolution"] = best_resolution

    # reference_builder runs after normalisation so it can see layers['scran_normalized']
    # — the placeholder reference used for the fixture path is built from the data
    # itself, not a fixed external panel.
    reference = reference_builder(adata)
    annotate_with_singler(adata, reference)
    return adata


def run_pipeline(
    raw_dir: str | Path,
    *,
    version: str | None = None,
    resolutions: tuple[float, ...] = (0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
) -> ad.AnnData:
    """Run the full pipeline against the real, extracted GSE176078 data and write the
    versioned signature artifact. Returns the resulting AnnData."""
    adata = load_gse176078(raw_dir)
    adata = _assemble_signature(
        adata,
        resolutions=resolutions,
        reference_builder=lambda _adata: load_human_primary_cell_atlas(),
    )

    out_path = signature_path(version) if version else signature_path()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)
    return adata


def _fixture_reference_builder(adata: ad.AnnData) -> Any:
    return build_placeholder_reference(adata.layers["scran_normalized"], list(adata.var_names))


def run_pipeline_on_fixture(
    *,
    version: str | None = None,
    resolutions: tuple[float, ...] = (0.2, 0.5, 0.8),
    n_samples: int = 4,
    cells_per_sample: int = 120,
    n_genes: int = 1200,
) -> ad.AnnData:
    """Run the identical assembly code against the synthetic test fixture and write
    the fixture-scale signature artifact — proves the schema contract in CI without
    the real 100k-cell dataset. Annotates against a data-driven placeholder reference
    (see ``annotate.build_placeholder_reference``) rather than the real celldex panel,
    so this stays network-free. ``n_samples``/``cells_per_sample``/``n_genes`` default
    to ``make scrna-signature-fixture``'s size; tests pass smaller values to stay fast.
    """
    from tests.fixtures.synthetic_scrna import build_synthetic_scrna_adata

    adata = build_synthetic_scrna_adata(
        n_samples=n_samples, cells_per_sample=cells_per_sample, n_genes=n_genes
    )
    adata = _assemble_signature(
        adata, resolutions=resolutions, reference_builder=_fixture_reference_builder
    )

    out_path = signature_path(version, fixture=True) if version else signature_path(fixture=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)
    return adata


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--fixture":
        result = run_pipeline_on_fixture()
        print(f"Wrote fixture signature: {signature_path(fixture=True)} ({result.shape})")
    else:
        raw_dir = REPO_ROOT / "data" / "raw" / "scrna" / "gse176078" / GSE176078_EXTRACTED_DIR
        result = run_pipeline(raw_dir)
        print(f"Wrote signature: {signature_path()} ({result.shape})")
