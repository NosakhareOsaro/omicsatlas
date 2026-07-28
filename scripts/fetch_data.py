"""Scripted, checksum-verified fetcher for OmicsAtlas source datasets.

Every dataset this project uses is downloaded here, never placed manually. Each
fetcher writes a small provenance manifest (accession, URL, download date, SHA256) to
``data/.provenance/<dataset>.json`` and regenerates the corresponding table in
``data/README.md`` from those manifests, so the provenance record can never drift out
of sync with what was actually downloaded. Re-running a fetcher is idempotent: if the
target file exists and its SHA256 matches the recorded manifest, it is not
re-downloaded.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PROVENANCE_DIR = REPO_ROOT / "data" / ".provenance"
DATA_README = REPO_ROOT / "data" / "README.md"

CHUNK_SIZE = 1024 * 1024  # 1 MiB


@dataclass(frozen=True)
class Dataset:
    """A single fetchable file and the metadata recorded about it."""

    key: str
    modality: str
    description: str
    accession: str
    url: str
    dest: Path
    license_note: str
    extra: dict[str, str] = field(default_factory=dict)


GSE176078_SCRNA = Dataset(
    key="gse176078_scrna",
    modality="scRNA-seq",
    description=(
        "Wu et al. 2021 (Nat Genet) human breast cancer atlas, processed scRNA-seq "
        "count matrix and published cell-type metadata (26 primary tumours: 11 ER+, "
        "5 HER2+, 10 TNBC)"
    ),
    accession="GSE176078",
    url=(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/"
        "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
    ),
    dest=REPO_ROOT
    / "data"
    / "raw"
    / "scrna"
    / "gse176078"
    / "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz",
    license_note=(
        "Public GEO deposit, no access restriction on the processed matrix (raw FASTQs "
        "are EGA-controlled-access under EGAS00001005173 and are not used by this "
        "project). See the publication for reuse/citation expectations."
    ),
)

GSE176078_BULK_MATCHED = Dataset(
    key="gse176078_bulk_matched",
    modality="Bulk RNA-seq (GSE176078-matched)",
    description=(
        "Wu et al. 2021 (Nat Genet) bulk RNA-seq (RSEM-style estimated counts) from "
        "the same patient cohort as the scRNA-seq data — 24 of the 26 scRNA-seq "
        "patients have a matched bulk sample (verified by exact orig.ident overlap; "
        "see ADR-0003). Primary dataset for the Phase 5 CrossOmicsConcordance "
        "benchmark. Downloaded file is a gzipped tar containing one .txt file — "
        "extract with `make bulk-data-extract` before use."
    ),
    accession="GSE176078",
    url=(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/"
        "GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz"
    ),
    dest=REPO_ROOT
    / "data"
    / "raw"
    / "bulk"
    / "gse176078_matched"
    / "GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz",
    license_note=(
        "Public GEO deposit, no access restriction (same series as the scRNA-seq "
        "data; see its license_note for the EGA-controlled-access raw FASTQ caveat, "
        "which doesn't apply to this processed count matrix)."
    ),
)

DATASETS: dict[str, Dataset] = {
    GSE176078_SCRNA.key: GSE176078_SCRNA,
    GSE176078_BULK_MATCHED.key: GSE176078_BULK_MATCHED,
}

# --- TCGA-BRCA: fetched via the GDC API (multi-file, not a single static URL), so it
# doesn't fit the Dataset/fetch() pattern above. See ADR-0003: secondary,
# generalisability-only dataset, not fed into the Phase 5 concordance benchmark. A
# documented, reproducible subset (not the full 1,098-sample cohort — HPC-scale,
# same pattern as Phase 1's real-run scope decision), split across both sample types
# so the secondary DE analysis has a genuine tumour-vs-normal comparison.
TCGA_BRCA_KEY = "tcga_brca_subset"
TCGA_BRCA_MODALITY = "Bulk RNA-seq (TCGA-BRCA)"
TCGA_BRCA_GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
TCGA_BRCA_GDC_DATA_URL = "https://api.gdc.cancer.gov/data"
TCGA_BRCA_N_TUMOR = 30
TCGA_BRCA_N_NORMAL = 15
TCGA_BRCA_DEST_DIR = REPO_ROOT / "data" / "raw" / "bulk" / "tcga_brca"
TCGA_BRCA_FIELDS = (
    "file_id,file_name,md5sum,file_size,cases.case_id,cases.submitter_id,"
    "cases.samples.sample_type"
)


def _gdc_filters(sample_type: str) -> dict:
    return {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-BRCA"]}},
            {
                "op": "in",
                "content": {"field": "data_type", "value": ["Gene Expression Quantification"]},
            },
            {
                "op": "in",
                "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]},
            },
            {"op": "in", "content": {"field": "experimental_strategy", "value": ["RNA-Seq"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": [sample_type]}},
        ],
    }


def query_gdc_files(sample_type: str, size: int) -> list[dict]:
    """Query GDC for TCGA-BRCA STAR-Counts files of a given sample type, sorted by
    file_id for reproducibility (re-running with the same size returns the same
    files). Returns the raw GDC 'hits' list."""
    import urllib.request as _urllib_request

    body = json.dumps(
        {
            "filters": _gdc_filters(sample_type),
            "fields": TCGA_BRCA_FIELDS,
            "format": "JSON",
            "size": str(size),
            "sort": "file_id:asc",
        }
    ).encode()
    request = _urllib_request.Request(
        TCGA_BRCA_GDC_FILES_URL, data=body, headers={"Content-Type": "application/json"}
    )
    with _urllib_request.urlopen(request) as response:
        payload = json.loads(response.read())
    return payload["data"]["hits"]


def md5_of(path: Path) -> str:
    # MD5 here only to match GDC's own published per-file checksum for verification —
    # not used for any security purpose.
    digest = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_tcga_brca_subset(*, now: str, force: bool = False) -> dict:
    """Fetch a reproducible TCGA-BRCA subset (30 Primary Tumor + 15 Solid Tissue
    Normal, sorted by file_id) from GDC's open-access STAR-Counts gene expression
    files. Idempotent per-file against GDC's own published MD5. Writes one raw file
    per sample to ``TCGA_BRCA_DEST_DIR`` and a manifest shaped like the Dataset-based
    manifests (so render_data_readme doesn't need to special-case it)."""
    hits = query_gdc_files("Primary Tumor", TCGA_BRCA_N_TUMOR) + query_gdc_files(
        "Solid Tissue Normal", TCGA_BRCA_N_NORMAL
    )

    TCGA_BRCA_DEST_DIR.mkdir(parents=True, exist_ok=True)
    file_records = []
    for hit in hits:
        dest = TCGA_BRCA_DEST_DIR / f"{hit['file_id']}_{hit['file_name']}"
        if not force and dest.exists() and md5_of(dest) == hit["md5sum"]:
            pass  # already correct, skip download
        else:
            download(f"{TCGA_BRCA_GDC_DATA_URL}/{hit['file_id']}", dest)
            actual_md5 = md5_of(dest)
            if actual_md5 != hit["md5sum"]:
                raise ValueError(
                    f"MD5 mismatch for {hit['file_name']}: GDC says {hit['md5sum']}, "
                    f"got {actual_md5}"
                )
        sample_type = hit["cases"][0]["samples"][0]["sample_type"]
        case_id = hit["cases"][0]["submitter_id"]
        file_records.append(
            {
                "file_id": hit["file_id"],
                "file_name": hit["file_name"],
                "md5sum": hit["md5sum"],
                "case_submitter_id": case_id,
                "sample_type": sample_type,
                "dest": os.path.relpath(dest, REPO_ROOT),
            }
        )

    # A single reproducible checksum over the sorted (file_id, md5) pairs — GDC
    # already guarantees each individual file's integrity via its own md5sum; this
    # exists so the one-row-per-dataset provenance table has a single value
    # representing "exactly this set of files".
    sorted_records = sorted(file_records, key=lambda r: r["file_id"])
    combined_input = "".join(f"{r['file_id']}:{r['md5sum']}" for r in sorted_records)
    combined = hashlib.sha256(combined_input.encode()).hexdigest()

    manifest = {
        "key": TCGA_BRCA_KEY,
        "modality": TCGA_BRCA_MODALITY,
        "description": (
            f"TCGA-BRCA STAR-Counts gene-level RNA-seq counts (GDC), reproducible "
            f"{TCGA_BRCA_N_TUMOR}-tumour + {TCGA_BRCA_N_NORMAL}-normal subset "
            f"sorted by file_id (open-access; full cohort is 1,098 samples, "
            f"HPC-scale — see ADR-0003). Secondary, generalisability-only dataset, "
            f"not fed into the Phase 5 concordance benchmark."
        ),
        "accession": "TCGA-BRCA",
        "url": TCGA_BRCA_GDC_FILES_URL,
        "dest": os.path.relpath(TCGA_BRCA_DEST_DIR, REPO_ROOT),
        "sha256": combined,
        "downloaded_at": now,
        "license_note": (
            "Open-access GDC data (gene-level quantification only; this project "
            "does not use any dbGaP-controlled-access TCGA files)."
        ),
        "files": file_records,
    }
    PROVENANCE_DIR.mkdir(parents=True, exist_ok=True)
    (PROVENANCE_DIR / f"{TCGA_BRCA_KEY}.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp_dest = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as response, tmp_dest.open("wb") as out_file:
        while chunk := response.read(CHUNK_SIZE):
            out_file.write(chunk)
    tmp_dest.replace(dest)


def load_manifest(dataset: Dataset) -> dict | None:
    manifest_path = PROVENANCE_DIR / f"{dataset.key}.json"
    if not manifest_path.exists():
        return None
    return json.loads(manifest_path.read_text())


def write_manifest(dataset: Dataset, checksum: str, downloaded_at: str) -> dict:
    manifest = {
        "key": dataset.key,
        "modality": dataset.modality,
        "description": dataset.description,
        "accession": dataset.accession,
        "url": dataset.url,
        "dest": os.path.relpath(dataset.dest, REPO_ROOT),
        "sha256": checksum,
        "downloaded_at": downloaded_at,
        "license_note": dataset.license_note,
    }
    PROVENANCE_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = PROVENANCE_DIR / f"{dataset.key}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def fetch(dataset: Dataset, *, now: str, force: bool = False) -> dict:
    """Fetch ``dataset``, verifying checksums and skipping a redundant download.

    Idempotency contract: if ``dataset.dest`` already exists and its SHA256 matches
    the checksum recorded in the existing manifest, no network request is made. If the
    file exists but the checksum doesn't match (or no manifest exists yet), it is
    (re-)downloaded and a fresh manifest is written.
    """
    existing_manifest = load_manifest(dataset)
    if not force and dataset.dest.exists() and existing_manifest is not None:
        actual = sha256_of(dataset.dest)
        if actual == existing_manifest.get("sha256"):
            return existing_manifest

    download(dataset.url, dataset.dest)
    checksum = sha256_of(dataset.dest)
    return write_manifest(dataset, checksum, downloaded_at=now)


def render_data_readme(manifests: list[dict]) -> str:
    """Render data/README.md from the current set of provenance manifests.

    This function is pure (manifests in, markdown out) so it can be unit tested
    without touching the filesystem or network.
    """
    header = (
        "# Data provenance\n\n"
        "All raw data is fetched via scripted, version-pinned fetchers "
        "(`scripts/fetch_data.py`) — never placed manually with no provenance record. "
        "This file is regenerated by that script from the manifests in "
        "`data/.provenance/`; do not hand-edit the table below.\n\n"
        "## Datasets\n\n"
        "| Modality | Description | Accession | Downloaded | SHA256 | License |\n"
        "|---|---|---|---|---|---|\n"
    )
    by_modality = {m["modality"]: m for m in manifests}
    rows = []
    for modality in (
        "scRNA-seq",
        "Bulk RNA-seq (GSE176078-matched)",
        "Bulk RNA-seq (TCGA-BRCA)",
        "Spatial (Visium)",
        "ATAC-seq",
    ):
        m = by_modality.get(modality)
        if m is None:
            rows.append(f"| {modality} | TBD | TBD | - | - | - |")
            continue
        sha_short = m["sha256"][:12] + "…"
        rows.append(
            f"| {modality} | {m['description']} | {m['accession']} | "
            f"{m['downloaded_at']} | `{sha_short}` | {m['license_note']} |"
        )
    footer = (
        "\n## Notes\n\n"
        "- The GSE176078-matched bulk dataset was chosen as the **primary** bulk "
        "RNA-seq input for the Phase 5 CrossOmicsConcordance benchmark over "
        "TCGA-BRCA, because 24 of its 26 scRNA-seq patients are exact-ID-matched "
        "to the same patients as the Phase 1 scRNA-seq signature — see "
        "`adr/ADR-0003-bulk-data-choice.md` for the full comparison. TCGA-BRCA is "
        "fetched too, as a secondary, generalisability-only analysis.\n"
    )
    return header + "\n".join(rows) + "\n" + footer


def regenerate_data_readme() -> None:
    manifests = []
    for dataset in DATASETS.values():
        manifest = load_manifest(dataset)
        if manifest is not None:
            manifests.append(manifest)
    tcga_manifest_path = PROVENANCE_DIR / f"{TCGA_BRCA_KEY}.json"
    if tcga_manifest_path.exists():
        manifests.append(json.loads(tcga_manifest_path.read_text()))
    DATA_README.write_text(render_data_readme(manifests))


def _now_iso() -> str:
    import datetime as _dt

    return _dt.datetime.now(_dt.UTC).strftime("%Y-%m-%d")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "dataset",
        nargs="?",
        default="all",
        choices=["all", *DATASETS.keys(), TCGA_BRCA_KEY],
        help="Which dataset to fetch (default: all).",
    )
    parser.add_argument(
        "--force", action="store_true", help="Re-download even if checksum matches."
    )
    args = parser.parse_args(argv)

    now = _now_iso()
    if args.dataset == "all":
        dataset_keys = list(DATASETS.keys())
    elif args.dataset == TCGA_BRCA_KEY:
        dataset_keys = []
    else:
        dataset_keys = [args.dataset]

    for key in dataset_keys:
        manifest = fetch(DATASETS[key], now=now, force=args.force)
        print(f"{key}: {manifest['dest']} sha256={manifest['sha256'][:12]}…")

    if args.dataset in (TCGA_BRCA_KEY, "all"):
        manifest = fetch_tcga_brca_subset(now=now, force=args.force)
        print(
            f"{TCGA_BRCA_KEY}: {manifest['dest']} ({len(manifest['files'])} files) "
            f"sha256={manifest['sha256'][:12]}…"
        )

    regenerate_data_readme()
    print(f"Regenerated {DATA_README.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
