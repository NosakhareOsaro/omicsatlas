"""Smoke test confirming the package imports and the test harness is wired up."""

import omicsatlas


def test_package_importable():
    assert omicsatlas.__version__
