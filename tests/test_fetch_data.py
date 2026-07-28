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
    assert "| Bulk RNA-seq | TBD | TBD | - | - | - |" in readme
    assert "do not hand-edit the table below" in readme
