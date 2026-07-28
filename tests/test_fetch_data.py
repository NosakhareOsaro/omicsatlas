"""Tests for scripts/fetch_data.py.

Network access is always mocked here — these tests must never depend on reaching the
real GEO servers, matching the project's CI constraint that tests never require
network access to real datasets.
"""

from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path

import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"
_spec = importlib.util.spec_from_file_location("fetch_data", _SCRIPTS_DIR / "fetch_data.py")
fetch_data = importlib.util.module_from_spec(_spec)
sys.modules["fetch_data"] = fetch_data
_spec.loader.exec_module(fetch_data)


@pytest.fixture
def dataset(tmp_path: Path) -> fetch_data.Dataset:
    return fetch_data.Dataset(
        key="test_dataset",
        modality="scRNA-seq",
        description="Test dataset",
        accession="GSE000000",
        url="https://example.invalid/test.tar.gz",
        dest=tmp_path / "raw" / "test.tar.gz",
        license_note="Test license note",
    )


@pytest.fixture(autouse=True)
def _isolated_provenance_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fetch_data, "PROVENANCE_DIR", tmp_path / "provenance")


def test_sha256_of_matches_hashlib(tmp_path: Path) -> None:
    path = tmp_path / "hello.txt"
    content = b"hello world" * 100_000  # larger than CHUNK_SIZE-independent for good measure
    path.write_bytes(content)

    assert fetch_data.sha256_of(path) == hashlib.sha256(content).hexdigest()


def test_fetch_downloads_when_no_existing_manifest(
    dataset: fetch_data.Dataset, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = []

    def fake_download(url: str, dest: Path) -> None:
        calls.append(url)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"fake dataset contents")

    monkeypatch.setattr(fetch_data, "download", fake_download)

    manifest = fetch_data.fetch(dataset, now="2026-07-28")

    assert calls == [dataset.url]
    assert dataset.dest.exists()
    assert manifest["sha256"] == fetch_data.sha256_of(dataset.dest)
    assert manifest["accession"] == "GSE000000"
    assert manifest["downloaded_at"] == "2026-07-28"


def test_fetch_is_idempotent_when_checksum_matches(
    dataset: fetch_data.Dataset, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset.dest.parent.mkdir(parents=True, exist_ok=True)
    dataset.dest.write_bytes(b"already downloaded")
    fetch_data.write_manifest(
        dataset, checksum=fetch_data.sha256_of(dataset.dest), downloaded_at="2026-01-01"
    )

    def fail_if_called(url: str, dest: Path) -> None:
        raise AssertionError("download() should not be called when checksum matches")

    monkeypatch.setattr(fetch_data, "download", fail_if_called)

    manifest = fetch_data.fetch(dataset, now="2026-07-28")

    assert manifest["downloaded_at"] == "2026-01-01"  # untouched, no re-download


def test_fetch_redownloads_when_checksum_mismatched(
    dataset: fetch_data.Dataset, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset.dest.parent.mkdir(parents=True, exist_ok=True)
    dataset.dest.write_bytes(b"stale corrupted content")
    fetch_data.write_manifest(dataset, checksum="0" * 64, downloaded_at="2026-01-01")

    calls = []

    def fake_download(url: str, dest: Path) -> None:
        calls.append(url)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"freshly downloaded content")

    monkeypatch.setattr(fetch_data, "download", fake_download)

    manifest = fetch_data.fetch(dataset, now="2026-07-28")

    assert calls == [dataset.url]
    assert manifest["downloaded_at"] == "2026-07-28"
    assert manifest["sha256"] == fetch_data.sha256_of(dataset.dest)


def test_render_data_readme_includes_known_modality_and_placeholders() -> None:
    manifests = [
        {
            "modality": "scRNA-seq",
            "description": "Test breast cancer atlas",
            "accession": "GSE176078",
            "downloaded_at": "2026-07-28",
            "sha256": "a" * 64,
            "license_note": "Public GEO deposit",
        }
    ]

    readme = fetch_data.render_data_readme(manifests)

    assert "GSE176078" in readme
    assert "2026-07-28" in readme
    assert "aaaaaaaaaaaa…" in readme
    assert "| Bulk RNA-seq (GSE176078-matched) | TBD | TBD | - | - | - |" in readme
    assert "| Bulk RNA-seq (TCGA-BRCA) | TBD | TBD | - | - | - |" in readme
    assert "do not hand-edit the table below" in readme


def test_md5_of_matches_hashlib(tmp_path: Path) -> None:
    path = tmp_path / "hello.txt"
    content = b"hello tcga"
    path.write_bytes(content)

    assert fetch_data.md5_of(path) == hashlib.md5(content).hexdigest()  # noqa: S324


def _fake_gdc_hit(file_id: str, content: bytes, sample_type: str, submitter_id: str) -> dict:
    return {
        "file_id": file_id,
        "file_name": f"{file_id}.tsv",
        "md5sum": hashlib.md5(content).hexdigest(),  # noqa: S324
        "cases": [{"submitter_id": submitter_id, "samples": [{"sample_type": sample_type}]}],
    }


def _query_returning_only_tumor(tumor_hit: dict):
    def _query(sample_type: str, size: int) -> list[dict]:
        return [tumor_hit] if sample_type == "Primary Tumor" else []

    return _query


def test_fetch_tcga_brca_subset_downloads_and_verifies_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(fetch_data, "TCGA_BRCA_DEST_DIR", tmp_path / "tcga_brca")
    tumor_hit = _fake_gdc_hit("file-tumor-1", b"tumor content", "Primary Tumor", "TCGA-AA-0001")
    normal_hit = _fake_gdc_hit(
        "file-normal-1", b"normal content", "Solid Tissue Normal", "TCGA-AA-0002"
    )

    fake_content = {"file-tumor-1": b"tumor content", "file-normal-1": b"normal content"}

    def fake_query_gdc_files(sample_type: str, size: int) -> list[dict]:
        return [tumor_hit] if sample_type == "Primary Tumor" else [normal_hit]

    def fake_download(url: str, dest: Path) -> None:
        file_id = url.rsplit("/", 1)[-1]
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(fake_content[file_id])

    monkeypatch.setattr(fetch_data, "query_gdc_files", fake_query_gdc_files)
    monkeypatch.setattr(fetch_data, "download", fake_download)

    manifest = fetch_data.fetch_tcga_brca_subset(now="2026-07-28")

    assert manifest["modality"] == fetch_data.TCGA_BRCA_MODALITY
    assert manifest["accession"] == "TCGA-BRCA"
    assert len(manifest["files"]) == 2
    sample_types = {f["sample_type"] for f in manifest["files"]}
    assert sample_types == {"Primary Tumor", "Solid Tissue Normal"}
    assert len(manifest["sha256"]) == 64


def test_fetch_tcga_brca_subset_raises_on_md5_mismatch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(fetch_data, "TCGA_BRCA_DEST_DIR", tmp_path / "tcga_brca")
    tumor_hit = _fake_gdc_hit("file-tumor-1", b"expected content", "Primary Tumor", "TCGA-AA-0001")

    def fake_download(url: str, dest: Path) -> None:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"WRONG content")

    monkeypatch.setattr(fetch_data, "query_gdc_files", _query_returning_only_tumor(tumor_hit))
    monkeypatch.setattr(fetch_data, "download", fake_download)

    with pytest.raises(ValueError, match="MD5 mismatch"):
        fetch_data.fetch_tcga_brca_subset(now="2026-07-28")


def test_fetch_tcga_brca_subset_is_idempotent_per_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dest_dir = tmp_path / "tcga_brca"
    monkeypatch.setattr(fetch_data, "TCGA_BRCA_DEST_DIR", dest_dir)
    content = b"already correct"
    tumor_hit = _fake_gdc_hit("file-tumor-1", content, "Primary Tumor", "TCGA-AA-0001")
    dest_dir.mkdir(parents=True, exist_ok=True)
    (dest_dir / f"{tumor_hit['file_id']}_{tumor_hit['file_name']}").write_bytes(content)

    monkeypatch.setattr(fetch_data, "query_gdc_files", _query_returning_only_tumor(tumor_hit))

    def fail_if_called(url: str, dest: Path) -> None:
        raise AssertionError("download() should not be called when the file is already correct")

    monkeypatch.setattr(fetch_data, "download", fail_if_called)

    manifest = fetch_data.fetch_tcga_brca_subset(now="2026-07-28")

    assert len(manifest["files"]) == 1


def test_render_data_readme_includes_both_bulk_rows_when_present() -> None:
    manifests = [
        {
            "modality": "Bulk RNA-seq (GSE176078-matched)",
            "description": "Matched bulk",
            "accession": "GSE176078",
            "downloaded_at": "2026-07-28",
            "sha256": "b" * 64,
            "license_note": "Public GEO deposit",
        },
        {
            "modality": "Bulk RNA-seq (TCGA-BRCA)",
            "description": "TCGA-BRCA subset",
            "accession": "TCGA-BRCA",
            "downloaded_at": "2026-07-28",
            "sha256": "c" * 64,
            "license_note": "Open-access GDC data",
        },
    ]

    readme = fetch_data.render_data_readme(manifests)

    assert "Matched bulk" in readme
    assert "TCGA-BRCA subset" in readme
    assert "ADR-0003" in readme
