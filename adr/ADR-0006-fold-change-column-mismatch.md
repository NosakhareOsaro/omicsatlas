# ADR-0006: `logFC` vs `log2FoldChange` Column Mismatch in the TCGA-BRCA Enrichment Script

## Status

Accepted

## Context

The real TCGA-BRCA pipeline run (`scripts/run_bulk_pipeline_tcga.R`) failed twice
with `Error in names(object) <- nm : attempt to set an attribute on NULL`. ADR-0005
records the initial (incorrect) diagnosis: KEGG REST API flakiness. That diagnosis
was reached honestly but was wrong, and it's worth recording exactly how so the
next person reading history doesn't have to reconstruct the investigation from
commit messages alone.

**What actually happened:** the script's local `run_enrichment_for(res, label)`
helper hardcoded `res$log2FoldChange` for both DE methods:

```r
ranked <- stats::setNames(res$log2FoldChange, gene_symbols)
```

[`run_deseq2()`](../r_package/OmicsAtlasBulk/R/deseq2.R) results have a
`log2FoldChange` column, but [`run_edger()`](../r_package/OmicsAtlasBulk/R/edger.R)
results use edgeR's own `logFC` column name instead — a real, pre-existing,
documented difference between the two methods' output conventions (also handled
explicitly, per-column, in `compare_de_methods()`'s discrepancy table). So
`edger_res$log2FoldChange` was always `NULL`, and `stats::setNames(NULL, ...)`
throws exactly this error internally (`setNames()` calls `names(object) <- nm`,
and `object` is `NULL`).

This is a **deterministic, 100%-reproducible bug**, not a random or
network-dependent one. It was misdiagnosed as KEGG flakiness because of how the
failure presented: the crash happens on the very first line of
`run_enrichment_for()`, before that function's own `message()` call reporting
which method it's processing — so the log showed DESeq2's enrichment succeeding
completely (including its own live KEGG calls and GSEA run), then silently jumped
to the crash with no visible marker that a *second*, different call
(`run_enrichment_for(edger_res, "edgeR")`) had even started. Every isolated
reproduction attempt built from the visible log output naturally reconstructed
only the DESeq2 path (which was never broken), so it kept succeeding — a full
`options(error = ...)` traceback, captured only after reproducing the failure via
the actual script file rather than an inline `Rscript -e` snippet, was what
finally revealed `run_enrichment_for(edger_res, "edgeR")` as the true call site.

## Decision

Add `extract_fold_change(res)` to `R/compare_de.R` (alongside
`compare_de_methods()`, which already has to know about this same DESeq2/edgeR
column-naming difference): checks for `log2FoldChange` first, then `logFC`, and
**errors explicitly** if neither is present, rather than silently returning `NULL`
again for some future third DE method with yet another naming convention.
`run_enrichment_for()` now calls this instead of hardcoding a column name.

The retry-with-backoff logic from ADR-0005 (`retry_with_backoff()` wrapping
`enrichKEGG()`) is **unaffected by this finding and is being kept**. It's correctly
implemented and tested against exactly what it claims to handle (a live external
API call that can genuinely fail transiently) — it simply wasn't the fix for
*this* incident, because this incident's failure never actually reached the
KEGG-wrapped call for the edgeR path at all.

## Consequences

- `extract_fold_change()` is now the one place that knows about DESeq2/edgeR's
  differing fold-change column names for *ranking* purposes, alongside
  `compare_de_methods()`'s existing per-column handling for its discrepancy table
  (not consolidated into `extract_fold_change()` itself, since that table
  deliberately reports both methods' values side by side rather than picking one).
- This bug shipped despite `test-enrichment.R`'s fixture tests passing, because
  every fixture test only ever exercised `run_enrichment()` directly with
  manually-constructed `ranked_stats`, never through the script's
  `run_enrichment_for()` wrapper - a gap between package-level and
  script-level testing this incident surfaced.

## Alternatives considered

- **Rename `run_edger()`'s output column to `log2FoldChange` for consistency**:
  rejected - `logFC` is edgeR's own standard convention, and changing the
  package's public return schema to paper over a script-level bug would be a much
  larger, riskier change than fixing the actual bug at its actual location.
- **Fix inline in the script with a bare `if (!is.null(...)) ... else ...`**:
  rejected per review feedback - a named, exported, unit-tested helper function
  makes the DESeq2/edgeR column-naming knowledge reusable and gives the "neither
  column present" case an explicit, testable error instead of silently
  reproducing the same class of bug for a future third DE method.
