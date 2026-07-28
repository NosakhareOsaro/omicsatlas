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

DATASETS: dict[str, Dataset] = {GSE176078_SCRNA.key: GSE176078_SCRNA}


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
    for modality in ("scRNA-seq", "Bulk RNA-seq", "Spatial (Visium)", "ATAC-seq"):
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
        "- GSE176078 (scRNA-seq) also includes a bulk RNA-seq raw-count matrix from "
        "the same 26-patient cohort "
        "(`GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz`). This would be a "
        "more tightly matched concordance comparator than TCGA-BRCA's unrelated "
        "cohort — flagged here as an option for Phase 2 review, not yet used.\n"
    )
    return header + "\n".join(rows) + "\n" + footer


def regenerate_data_readme() -> None:
    manifests = []
    for dataset in DATASETS.values():
        manifest = load_manifest(dataset)
        if manifest is not None:
            manifests.append(manifest)
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
        choices=["all", *DATASETS.keys()],
        help="Which dataset to fetch (default: all).",
    )
    parser.add_argument(
        "--force", action="store_true", help="Re-download even if checksum matches."
    )
    args = parser.parse_args(argv)

    targets = list(DATASETS.values()) if args.dataset == "all" else [DATASETS[args.dataset]]
    for dataset in targets:
        manifest = fetch(dataset, now=_now_iso(), force=args.force)
        print(f"{dataset.key}: {manifest['dest']} sha256={manifest['sha256'][:12]}…")

    regenerate_data_readme()
    print(f"Regenerated {DATA_README.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
