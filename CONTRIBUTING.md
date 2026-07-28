# Contributing

OmicsAtlas is currently maintained by a single author. This guide documents the
conventions I hold the project to, both for my own consistency and for anyone who wants
to propose a change.

## Development setup

```
mamba env create -f environment/env.yml   # Apple Silicon: prefix with CONDA_SUBDIR=osx-64
conda activate omicsatlas
pip install -e ".[dev]" rpy2==3.6.7 pertpy==1.0.3 scikit-learn==1.9.0
pre-commit install
```

Or via Docker: `docker build -f environment/Dockerfile -t omicsatlas .`

## Before opening a PR

- `pre-commit run --all-files` must pass (ruff, black, mypy, plus basic hygiene hooks).
- `pytest` must pass with the coverage gate met (currently 85%, enforced in
  `pyproject.toml`).
- New pipeline stages must land with tests in the same commit (or the immediately
  following one), and a doc page under `docs/` explaining *why* that tool/parameter
  choice was made, not just how to run it.
- Non-trivial design decisions get an ADR under `adr/`, following the template used in
  `adr/ADR-0001-license-and-architecture.md` (Context / Decision / Consequences /
  Alternatives considered), written *before* the implementation commit.

## Commit style

[Conventional Commits](https://www.conventionalcommits.org/): `feat:`, `fix:`, `docs:`,
`test:`, `refactor:`, `chore:`, `perf:`, `ci:`. One logical change per commit — a
reviewer should be able to read any single commit in under two minutes. Pipeline-stage
commits add their test in the same commit.

## Data

Raw data is never committed or manually placed. It is fetched via scripted,
version-pinned fetchers with recorded MD5 checksums — see `data/README.md` for the
provenance record and accession numbers.
