"""Tests for src/omicsatlas/scrna/artifact.py."""

from __future__ import annotations

from omicsatlas.scrna.artifact import CURRENT_SIGNATURE_VERSION, REPO_ROOT, signature_path


def test_signature_path_uses_current_version_by_default() -> None:
    path = signature_path()
    assert CURRENT_SIGNATURE_VERSION in path.parts
    assert path.name == "brca_scrna_signature.h5ad"


def test_signature_path_real_vs_fixture_are_distinct() -> None:
    real_path = signature_path()
    fixture_path = signature_path(fixture=True)

    assert real_path != fixture_path
    assert "real" in real_path.parts
    assert "fixture" in fixture_path.parts


def test_signature_path_respects_explicit_version() -> None:
    path = signature_path("v2")
    assert "v2" in path.parts


def test_signature_path_is_under_repo_root_data_processed() -> None:
    path = signature_path()
    assert path.is_relative_to(REPO_ROOT / "data" / "processed")
