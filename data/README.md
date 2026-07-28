# Data provenance

All raw data is fetched via scripted, version-pinned fetchers (`scripts/fetch_data.py`,
added alongside each pipeline phase) — never placed manually with no provenance record.

For each dataset used in this project, this file will record: source, accession number,
download date, fetcher script/commit, and MD5 checksum.

## Datasets (populated per phase)

| Modality | Source | Accession | Fetched | MD5 |
|---|---|---|---|---|
| scRNA-seq | 10x Genomics / Human Breast Cancer Atlas | TBD (Phase 1) | - | - |
| Bulk RNA-seq | TCGA-BRCA (GDC / recount3, hg38) | TBD (Phase 2) | - | - |
| Spatial (Visium) | 10x Genomics Datasets portal | TBD (Phase 3) | - | - |
| ATAC-seq | ENCODE breast tissue/cell line | TBD (Phase 4) | - | - |

## Licensing

Each dataset's original license/usage terms are documented here as it is added; this
project's own code is licensed separately (see root `LICENSE`).
