"""Per-sample Scrublet doublet detection.

Run separately per sample (10x GEM well) — Scrublet's simulated doublets are only
statistically meaningful within a single lane, never across pooled samples. Each
sample's ``expected_doublet_rate`` is scaled by its own cell count rather than using
Scrublet's flat default, approximating 10x's published loading/doublet-rate table.
Rationale is recorded in ``adr/ADR-0002-scrna-pipeline-design.md``.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import scrublet as scr

DEFAULT_SAMPLE_COL = "orig.ident"
DEFAULT_FALLBACK_THRESHOLD = 0.25
DOUBLET_RATE_PER_1000_CELLS = 0.008
MAX_EXPECTED_DOUBLET_RATE = 0.20
MIN_EXPECTED_DOUBLET_RATE = 1e-4


def expected_doublet_rate(n_cells: int) -> float:
    """Approximate 10x's published loading/doublet-rate table: ~0.8% doublets per
    1,000 cells loaded, capped at 20%."""
    rate = DOUBLET_RATE_PER_1000_CELLS * (n_cells / 1000)
    return float(np.clip(rate, MIN_EXPECTED_DOUBLET_RATE, MAX_EXPECTED_DOUBLET_RATE))


def _run_scrublet_one_sample(
    counts_matrix: np.ndarray,
    *,
    random_state: int,
    fallback_threshold: float,
    scrublet_kwargs: dict,
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Returns (doublet_scores, predicted_doublets, used_fallback_threshold)."""
    n_cells = counts_matrix.shape[0]
    scrub = scr.Scrublet(
        counts_matrix,
        expected_doublet_rate=expected_doublet_rate(n_cells),
        random_state=random_state,
    )
    scores, calls = scrub.scrub_doublets(verbose=False, **scrublet_kwargs)
    if calls is None:
        # Automatic bimodal threshold detection failed (common in low-doublet-rate or
        # very small samples) — fall back to a fixed threshold, logged per-sample
        # rather than applied silently. See ADR-0002.
        calls = scrub.call_doublets(threshold=fallback_threshold, verbose=False)
        return np.asarray(scores), np.asarray(calls), True
    return np.asarray(scores), np.asarray(calls), False


def run_scrublet_per_sample(
    adata: ad.AnnData,
    *,
    sample_col: str = DEFAULT_SAMPLE_COL,
    fallback_threshold: float = DEFAULT_FALLBACK_THRESHOLD,
    random_state: int = 0,
    scrublet_kwargs: dict | None = None,
) -> ad.AnnData:
    """Run Scrublet independently per sample, in place. Adds
    ``obs['doublet_score']``, ``obs['predicted_doublet']``, and
    ``obs['doublet_threshold_fallback_used']``. Returns ``adata`` for chaining."""
    scrublet_kwargs = scrublet_kwargs or {}
    n_cells = adata.n_obs
    doublet_scores = np.full(n_cells, np.nan)
    predicted = np.zeros(n_cells, dtype=bool)
    used_fallback = np.zeros(n_cells, dtype=bool)

    for _, idx in adata.obs.groupby(sample_col, observed=True).groups.items():
        positions = adata.obs.index.get_indexer(idx)
        counts_matrix = np.asarray(adata.X[positions])
        scores, calls, fallback_used = _run_scrublet_one_sample(
            counts_matrix,
            random_state=random_state,
            fallback_threshold=fallback_threshold,
            scrublet_kwargs=scrublet_kwargs,
        )
        doublet_scores[positions] = scores
        predicted[positions] = calls
        used_fallback[positions] = fallback_used

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted
    adata.obs["doublet_threshold_fallback_used"] = used_fallback
    return adata
