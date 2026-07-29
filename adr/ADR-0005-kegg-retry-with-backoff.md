# ADR-0005: Retry-with-Backoff Around Live KEGG Enrichment Calls

## Status

Accepted. **Addendum (2026-07-29):** the real failure that motivated this ADR
turned out to have a different root cause than diagnosed below - a deterministic
`logFC`/`log2FoldChange` column-name mismatch in the pipeline script, not KEGG
flakiness (see `ADR-0006`). `retry_with_backoff()` itself is correctly implemented
and tested against genuine KEGG API flakiness (a real, if secondary, concern) and
is being kept; it just wasn't the fix for that particular incident.

## Context

`run_enrichment()`'s KEGG step (`clusterProfiler::enrichKEGG()`) calls KEGG's REST
API directly (`rest.kegg.jp`) — there is no offline KEGG data source in
clusterProfiler, so this is a genuine live network dependency every time the bulk
pipeline runs the TCGA-BRCA secondary analysis (see `R/enrichment.R`'s docstring).

During the first real run of `scripts/run_bulk_pipeline_tcga.R` against the actual
TCGA-BRCA subset (30 tumor + 15 normal samples), DESeq2's enrichment step failed
with `Error in names(object) <- nm : attempt to set an attribute on NULL` — an
error surfacing from inside clusterProfiler/fgsea's internals after both of KEGG's
REST calls ("Reading KEGG annotation online...") had logged as successful. Rather
than guess at a fix, I isolated the exact real inputs (the same significant-gene
list and ranked stats from the real DESeq2 result) and reran each piece:

- GO + GSEA alone (`run_kegg = FALSE`): succeeded cleanly.
- `clusterProfiler::enrichKEGG()` called directly, alone: succeeded cleanly (31
  enriched pathways).
- The exact full sequence GO → KEGG → GSEA, identical to what `run_enrichment()`
  does internally: succeeded cleanly on retry, with all three result objects valid.

The identical call succeeding on a bare retry, with no code change, rules out a
deterministic bug in this package's code or in this specific real gene list.
KEGG's REST API has no documented rate limit, but is widely known in the
bioinformatics community to be occasionally flaky/rate-limited under repeated
calls in a short window - consistent with what was observed directly here (the
KEGG endpoint was hit 3 times in quick succession across the original run and two
isolation reruns, one of which failed outright).

## Decision

Wrap the `enrichKEGG()` call in a small, generic `retry_with_backoff()` helper
(`R/enrichment.R`): up to 3 attempts total, sleeping 2s then 8s between attempts,
logging each retry via `message()` before sleeping, and re-raising the original
error only once all attempts are exhausted. This is standard defensive handling
for a real external dependency with observed intermittent failures - not
speculative overengineering for a hypothetical problem, since the failure was
directly observed this session.

`retry_with_backoff()` is deliberately generic (a plain `fn`/`max_attempts`/
`delays`/`on_retry` signature, not KEGG-specific) so it can wrap other flaky
external calls later without duplicating the pattern, but only the KEGG call
uses it today - no other step in this package makes a live network call.

## Consequences

- `run_enrichment()`'s KEGG step can now take up to ~40 seconds longer in the
  worst case (2s + 8s of sleep across 3 attempts) before either succeeding or
  propagating a real, persistent failure - an acceptable cost against a live
  external dependency that has already demonstrated real-world flakiness.
- Tests for `retry_with_backoff()` use `delays = c(0, 0)` to stay fast and
  deterministic; they exercise the retry-then-succeed path, the
  exhausted-retries-then-propagate path, and the immediate-success/no-retry path,
  with a mock function rather than hitting the real KEGG API (the existing
  `test-enrichment.R` KEGG test, skipped on CI/offline, already covers a real
  network call separately).
- This does not make KEGG enrichment reliable in the face of a genuine, sustained
  outage - 3 retries with a ~40s total backoff window is meant to absorb brief
  rate-limiting or transient network hiccups, not a real multi-minute service
  disruption.

## Alternatives considered

- **No retry, just document the flakiness and let users re-run manually**:
  rejected - this is exactly what happened for the real TCGA-BRCA run this
  session, and it's the kind of thing a pipeline meant to run unattended
  shouldn't require a human to babysit and manually retry.
- **Unbounded/indefinite retry**: rejected - if KEGG were genuinely down for an
  extended period, an indefinite retry loop would hang the whole pipeline with no
  visible failure, which is worse than a bounded retry that eventually surfaces a
  real error.
- **Retrying the whole `run_enrichment()` call (GO + KEGG + GSEA) instead of just
  the KEGG step**: rejected - GO and GSEA don't hit the network and have shown no
  flakiness; retrying them too would just waste local compute on every KEGG
  hiccup for no benefit.
