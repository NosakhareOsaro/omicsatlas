"""Tests for src/omicsatlas/scrna/qc.py against the synthetic fixture.

No test here touches the real dataset — see tests/fixtures/synthetic_scrna.py.
"""

from __future__ import annotations

from omicsatlas.scrna.qc import compute_qc_metrics, flag_qc_outliers, qc_summary, run_qc
from tests.fixtures.synthetic_scrna import build_synthetic_scrna_adata


def test_compute_qc_metrics_adds_expected_columns() -> None:
    adata = build_synthetic_scrna_adata()
    compute_qc_metrics(adata)

    for col in ["total_counts", "n_genes_by_counts", "pct_counts_mt", "log1p_total_counts"]:
        assert col in adata.obs.columns


def test_flag_qc_outliers_catches_injected_low_count_and_high_mito_cases() -> None:
    adata = build_synthetic_scrna_adata(n_samples=3, cells_per_sample=40)
    compute_qc_metrics(adata)
    flag_qc_outliers(adata)

    low_count = adata.obs["_qc_case"] == "low_count_outlier"
    high_mito = adata.obs["_qc_case"] == "high_mito_outlier"
    normal = adata.obs["_qc_case"] == "normal"

    assert adata.obs.loc[low_count, "qc_outlier"].all()
    assert adata.obs.loc[high_mito, "qc_outlier"].all()
    # Normal cells should overwhelmingly pass; allow a small tolerance since MAD-based
    # thresholds will always flag a few genuine tail cells even in "normal" data.
    normal_flag_rate = adata.obs.loc[normal, "qc_outlier"].mean()
    assert normal_flag_rate < 0.1

    assert (adata.obs.loc[low_count, "qc_fail_reason"] != "").all()
    assert (adata.obs.loc[high_mito, "qc_fail_reason"] != "").all()


def test_flag_qc_outliers_requires_qc_metrics_first() -> None:
    adata = build_synthetic_scrna_adata()
    try:
        flag_qc_outliers(adata)
        raised = False
    except ValueError:
        raised = True
    assert raised


def test_run_qc_removes_flagged_cells_and_rare_genes() -> None:
    adata = build_synthetic_scrna_adata(n_samples=3, cells_per_sample=40)
    filtered = run_qc(adata)

    assert filtered.n_obs < adata.n_obs
    assert not filtered.obs["qc_outlier"].any()
    # Original object is untouched.
    assert "qc_outlier" not in adata.obs.columns


def test_qc_summary_reports_per_sample_counts() -> None:
    adata = build_synthetic_scrna_adata(n_samples=3, cells_per_sample=40)
    compute_qc_metrics(adata)
    flag_qc_outliers(adata)

    summary = qc_summary(adata)

    assert set(summary.index) == set(adata.obs["orig.ident"].unique())
    assert (summary["n_kept"] + summary["n_outliers"] == summary["n_cells"]).all()
