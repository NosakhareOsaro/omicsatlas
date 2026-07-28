# ADR-0003: Bulk RNA-seq Data Choice — TCGA-BRCA vs. GSE176078-Matched Bulk

## Status

Accepted

## Context

Phase 2 needs a bulk RNA-seq dataset. The original project brief scoped TCGA-BRCA.
While building Phase 1's scRNA-seq fetcher, I found that GSE176078 (the same GEO
series as the Phase 1 scRNA-seq data) also ships a bulk RNA-seq count matrix
(`GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz`) from the same patient cohort.
This matters because Phase 5's CrossOmicsConcordance metric scores agreement between
cell-type composition vectors derived independently from scRNA-seq, spatial, and bulk
deconvolution *of the same tissue*. Bulk data from an unrelated cohort (TCGA-BRCA)
lets Phase 5 compare modalities on the same tissue *type*, but not the same tissue —
a materially weaker claim than bulk data from the literal same patients as the Phase
1 scRNA-seq signature.

I verified the GSE176078-matched bulk file directly before writing this ADR, rather
than deciding from the filename/description alone:

- It's a gzipped **tar** archive containing one file
  (`GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt`), 5.3MB — needs untarring, not
  just gunzipping, in the fetcher.
- Header row: `Genes` + 24 sample columns named `CID3586`, `CID3921`, `CID3941`, ...
  — **exactly the same patient IDs** used in Phase 1's `obs['orig.ident']`. Of the 26
  patients in the scRNA-seq cohort, **24 have a matched bulk sample**
  (`CID44991` and `CID45171` are missing from the bulk file).
- Values are non-integer (e.g. `32.3900253069881`) — RSEM/Salmon-style estimated
  counts, not raw integer read counts. DESeq2/edgeR both handle this via standard
  rounding (the same accommodation `tximport`-based workflows make), but it means
  the "raw_counts" in the filename isn't literally true and the pipeline needs to
  round before feeding either tool.
- Gene IDs are gene symbols (`DDX11L1`, `WASH7P`, ...), the same convention as the
  Phase 1 scRNA `var_names` — no ID-mapping needed between the two.

TCGA-BRCA, for comparison: RNA-seq gene-level counts from GDC's STAR-Counts
workflow, 60,660 genes across 1,098 samples, hg38, fully open-access (no controlled-
access application needed for the gene-level quantification data, only for
raw/germline-sensitive files this project doesn't use). A large, independent,
well-characterised cohort, but not connected to any patient in the Phase 1 scRNA-seq
data.

## Decision

Use **both**, with different roles:

- **GSE176078-matched bulk (24 patients) is the primary dataset feeding Phase 5's
  CrossOmicsConcordance benchmark.** The verified 24/26 exact-ID overlap with the
  Phase 1 scRNA-seq cohort means BayesPrism deconvolution (Phase 2) and the
  concordance comparison (Phase 5) can run on bulk samples from the literal same
  patients as the scRNA-seq signature — and, when Phase 3 lands, the same patients
  the spatial data (where available) comes from too. This is the strongest form of
  the "same tissue" claim the whole project's novel contribution depends on.
- **TCGA-BRCA is a secondary, generalisability analysis** — larger-cohort GO/KEGG/GSEA
  and DE results reported separately in the bulk pipeline's output, explicitly *not*
  fed into the CrossOmicsConcordance metric. It answers a different, still valuable
  question ("does this hold up in a much larger, independent cohort?") without
  silently diluting what the concordance metric can claim about same-tissue
  agreement.

## Consequences

- Phase 2's pipeline needs two entrypoints (matched-bulk primary, TCGA-BRCA
  secondary), not one — slightly more code than the brief's original single-dataset
  scope, but directly serves Phase 5's scientific soundness.
- The matched-bulk cohort (24 patients) is much smaller than TCGA-BRCA (1,098
  samples) — Phase 2's DE/enrichment results on it will have less statistical power
  than a TCGA-BRCA-only analysis would. This is expected and is exactly why TCGA-BRCA
  stays in the pipeline as the generalisability check, not as a replacement.
- The bulk fetcher must untar (not just gunzip) the matched-bulk file, and the
  DESeq2/edgeR input step must round the non-integer estimated counts — both are
  concrete implementation details this verification surfaced ahead of time rather
  than mid-implementation.
- `CID44991` and `CID45171` have scRNA-seq data but no matched bulk sample; they
  simply won't appear in the matched-bulk arm of Phase 2/5's analysis.

## Alternatives considered

- **TCGA-BRCA only, as originally scoped**: rejected as the sole dataset — would
  leave Phase 5's concordance benchmark comparing modalities across unrelated
  patients, which is a materially weaker claim than the project's own premise
  requires, once a same-patient alternative is known to exist and was verified
  available.
- **GSE176078-matched bulk only, dropping TCGA-BRCA**: rejected — would lose the
  generalisability check entirely and leave the bulk RNA-seq arm underpowered (24
  samples) with no larger-cohort validation, and would also depart from the original
  brief's scope more than necessary when both datasets are freely obtainable.
