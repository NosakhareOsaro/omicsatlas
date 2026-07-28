# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Repository scaffold for the multi-modal pipeline (scrna, bulk, spatial, atac,
  concordance, browser) and supporting directories.
- ADR-0001: license (MIT) and overall architecture decision.
- Python packaging (`pyproject.toml`, hatchling, Python 3.11), pre-commit hooks
  (ruff, black, mypy), conda environment, Docker image (pinned by digest), and a
  lint/test GitHub Actions CI workflow.
- `CONTRIBUTING.md`, issue templates, and PR template.
