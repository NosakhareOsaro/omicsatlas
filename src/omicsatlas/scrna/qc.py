"""Per-sample, MAD-based QC metrics and outlier filtering.

Thresholds and rationale are documented in ``adr/ADR-0002-scrna-pipeline-design.md``
and ``docs/scrna-pipeline.md``: 5 MADs (log-scale) on library size/complexity, 3 MADs
(or an absolute 20% ceiling) on mitochondrial fraction, computed per sample rather than
pooled across the whole cohort, plus an absolute floor of 200 genes/cell.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

DEFAULT_SAMPLE_COL = "orig.ident"
DEFAULT_MT_PREFIX = "MT-"
DEFAULT_N_MADS_COUNTS = 5.0
DEFAULT_N_MADS_MITO = 3.0
DEFAULT_MITO_CEILING_PCT = 20.0
DEFAULT_MIN_GENES = 200
DEFAULT_MIN_CELLS_PER_GENE = 3


def compute_qc_metrics(adata: ad.AnnData, *, mt_prefix: str = DEFAULT_MT_PREFIX) -> ad.AnnData:
    """Annotate ``adata.obs`` with total_counts / n_genes_by_counts / pct_counts_mt
    (and their log1p variants), in place. Returns ``adata`` for chaining."""
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)
    return adata


def _mad_outlier_two_sided(values: np.ndarray, n_mads: float) -> np.ndarray:
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    if mad == 0:
        return np.zeros(len(values), dtype=bool)
    return (values < median - n_mads * mad) | (values > median + n_mads * mad)


def _mad_outlier_upper(values: np.ndarray, n_mads: float) -> np.ndarray:
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    if mad == 0:
        return np.zeros(len(values), dtype=bool)
    return values > median + n_mads * mad


def flag_qc_outliers(
    adata: ad.AnnData,
    *,
    sample_col: str = DEFAULT_SAMPLE_COL,
    n_mads_counts: float = DEFAULT_N_MADS_COUNTS,
    n_mads_mito: float = DEFAULT_N_MADS_MITO,
    mito_ceiling_pct: float = DEFAULT_MITO_CEILING_PCT,
    min_genes: int = DEFAULT_MIN_GENES,
) -> ad.AnnData:
    """Flag QC outliers per-sample (not pooled) and record why, in place.

    Adds ``obs['qc_outlier']`` (bool) and ``obs['qc_fail_reason']`` (comma-joined
    strings, empty for cells that pass). Requires ``compute_qc_metrics`` to have
    already been run. Returns ``adata`` for chaining.
    """
    required = {
        "log1p_total_counts",
        "log1p_n_genes_by_counts",
        "pct_counts_mt",
        "n_genes_by_counts",
    }
    missing = required - set(adata.obs.columns)
    if missing:
        raise ValueError(f"Missing QC columns {missing}; call compute_qc_metrics() first.")

    n_cells = adata.n_obs
    outlier = np.zeros(n_cells, dtype=bool)
    reasons: list[list[str]] = [[] for _ in range(n_cells)]

    for _, idx in adata.obs.groupby(sample_col, observed=True).groups.items():
        positions = adata.obs.index.get_indexer(idx)

        counts_outlier = _mad_outlier_two_sided(
            adata.obs["log1p_total_counts"].to_numpy()[positions], n_mads_counts
        )
        genes_outlier = _mad_outlier_two_sided(
            adata.obs["log1p_n_genes_by_counts"].to_numpy()[positions], n_mads_counts
        )
        mito_values = adata.obs["pct_counts_mt"].to_numpy()[positions]
        mito_mad_outlier = _mad_outlier_upper(mito_values, n_mads_mito)
        mito_ceiling_outlier = mito_values > mito_ceiling_pct
        low_complexity = adata.obs["n_genes_by_counts"].to_numpy()[positions] < min_genes

        sample_outlier = (
            counts_outlier
            | genes_outlier
            | mito_mad_outlier
            | mito_ceiling_outlier
            | low_complexity
        )

        for local_i, global_pos in enumerate(positions):
            outlier[global_pos] = outlier[global_pos] or sample_outlier[local_i]
            if counts_outlier[local_i] or genes_outlier[local_i]:
                reasons[global_pos].append("library_size_mad")
            if mito_mad_outlier[local_i]:
                reasons[global_pos].append("mito_mad")
            if mito_ceiling_outlier[local_i]:
                reasons[global_pos].append("mito_ceiling")
            if low_complexity[local_i]:
                reasons[global_pos].append("low_complexity_floor")

    adata.obs["qc_outlier"] = outlier
    adata.obs["qc_fail_reason"] = [",".join(r) for r in reasons]
    return adata


def run_qc(
    adata: ad.AnnData,
    *,
    sample_col: str = DEFAULT_SAMPLE_COL,
    mt_prefix: str = DEFAULT_MT_PREFIX,
    n_mads_counts: float = DEFAULT_N_MADS_COUNTS,
    n_mads_mito: float = DEFAULT_N_MADS_MITO,
    mito_ceiling_pct: float = DEFAULT_MITO_CEILING_PCT,
    min_genes: int = DEFAULT_MIN_GENES,
    min_cells_per_gene: int = DEFAULT_MIN_CELLS_PER_GENE,
) -> ad.AnnData:
    """Compute QC metrics, flag per-sample outliers, and return a filtered copy.

    Does not mutate ``adata``; returns a new, QC-passed object with genes expressed
    in fewer than ``min_cells_per_gene`` cells also removed.
    """
    working = adata.copy()
    compute_qc_metrics(working, mt_prefix=mt_prefix)
    flag_qc_outliers(
        working,
        sample_col=sample_col,
        n_mads_counts=n_mads_counts,
        n_mads_mito=n_mads_mito,
        mito_ceiling_pct=mito_ceiling_pct,
        min_genes=min_genes,
    )
    filtered = working[~working.obs["qc_outlier"]].copy()
    sc.pp.filter_genes(filtered, min_cells=min_cells_per_gene)
    return filtered


def qc_summary(adata: ad.AnnData, *, sample_col: str = DEFAULT_SAMPLE_COL) -> pd.DataFrame:
    """Per-sample counts of cells kept/flagged, for the QC section of the pipeline doc."""
    if "qc_outlier" not in adata.obs.columns:
        raise ValueError("Missing 'qc_outlier'; call flag_qc_outliers() first.")
    return (
        adata.obs.groupby(sample_col, observed=True)["qc_outlier"]
        .agg(n_cells="size", n_outliers="sum")
        .assign(n_kept=lambda df: df["n_cells"] - df["n_outliers"])
    )
