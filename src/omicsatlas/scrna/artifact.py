"""The scRNA-seq signature artifact: a stable, versioned build product.

Phases 3 (RCTD) and 5 (BayesPrism, CrossOmicsConcordance) import this exact versioned
path explicitly — they do not regenerate the signature themselves with potentially
different parameters. The ``.h5ad`` file itself is a build product, not committed to
git (``data/processed/`` is gitignored); what's committed is this module (the version
identifier and path contract) and ``DATA_SCHEMA.md`` (the schema it must satisfy). See
ADR-0002.
"""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
CURRENT_SIGNATURE_VERSION = "v1"


def signature_path(version: str = CURRENT_SIGNATURE_VERSION, *, fixture: bool = False) -> Path:
    """Path to the versioned scRNA-seq signature artifact.

    ``fixture=True`` returns the path used by ``make scrna-signature-fixture`` (built
    from the synthetic test fixture, for CI-safe contract testing) instead of the real
    ``make scrna-signature`` output.
    """
    subdir = "fixture" if fixture else "real"
    return (
        REPO_ROOT / "data" / "processed" / "scrna" / subdir / version / "brca_scrna_signature.h5ad"
    )
