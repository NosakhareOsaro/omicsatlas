# ADR-0001: Project License and Overall Architecture

## Status

Accepted

## Context

OmicsAtlas is a multi-modal genomics platform combining four independent pipelines
(scRNA-seq, bulk RNA-seq, spatial transcriptomics, ATAC-seq) plus a novel cross-modality
concordance metric, `CrossOmicsConcordance`. I need to settle two foundational questions
before writing any pipeline code: what license governs the project, and how the codebase
is structured so that four fairly different bioinformatics toolchains (heavy Python
scientific stack, R/Bioconductor, Snakemake, Nextflow) and one independently publishable
sub-package can coexist without becoming an unmaintainable monolith.

I am the sole author and maintainer. The project is intended to be reproducible,
inspectable, and reusable by others — both as a working analysis platform and as the
source of a standalone PyPI/Bioconda package.

## Decision

### License: MIT

The project is licensed under the MIT License (see `LICENSE`).

Rationale:
- The core scientific Python ecosystem this project builds directly on top of
  (scanpy, squidpy, scvelo) is released under BSD-3-Clause, a license with essentially
  the same permissiveness and intent as MIT. Matching that convention keeps the project
  legally consistent with its dependencies and familiar to anyone coming from that
  ecosystem.
- MIT is the shortest, least ambiguous permissive license, which matters for a
  standalone package (`CrossOmicsConcordance`) I intend to publish to PyPI and Bioconda
  and want others to be able to adopt with zero friction.
- As sole author, there are no other contributors whose IP needs the more formal
  contribution/patent language Apache-2.0 provides. That protection is not buying
  anything here, only added length.

**Alternative considered: Apache-2.0.** Rejected. Apache-2.0's explicit patent grant is
valuable for projects worried about submarine patent claims from a large contributor
base — not a real risk for a single-author academic/portfolio tool. The extra length and
required "NOTICE" file handling isn't worth it for the reuse-friction reduction MIT gives
instead.

### Overall architecture

- **Monorepo, not multi-repo.** All four modality pipelines, the R package, the
  Nextflow/Snakemake orchestration, and the concordance package live in one repository
  under `src/omicsatlas/<modality>/`. A single clone reproduces the entire platform,
  which matters for the project's reproducibility goals. The cost (larger repo, mixed
  toolchains) is acceptable at this scale and for a single maintainer.
- **`concordance/` is structured as an independent installable sub-package** (its own
  `pyproject.toml`, tests, README, semantic version) even though it lives inside the
  monorepo. This is the novel-contribution artifact and needs to be usable by others
  without pulling in the rest of OmicsAtlas's heavy dependency stack (Scanpy, Seurat
  bridges, etc.). Structuring it this way from the start avoids a painful later
  extraction.
- **Language boundary is explicit, not hidden.** Python owns scRNA-seq, spatial, ATAC
  orchestration, and the concordance metric. R owns the bulk RNA-seq statistical stack
  (DESeq2/edgeR/clusterProfiler) as a proper R package under `r_package/`, because those
  tools are best-in-class in R and reimplementing them in Python would both be wasted
  effort and reduce trust in the results. Where Python needs to call R (e.g. BayesPrism
  deconvolution feeding the concordance benchmark), that boundary is crossed explicitly
  via subprocess/rpy2 wrappers in `src/omicsatlas/bulk/`, not silently.
- **Two workflow engines, scoped by what they're good at.** Nextflow DSL2
  (`pipelines/nextflow/`) orchestrates the end-to-end multi-modal workflow (resume
  support, HPC/cloud resource profiles). Snakemake (`pipelines/snakemake/`) is scoped to
  the ATAC-seq pipeline specifically, matching existing ENCODE-adjacent tooling
  conventions in that space. I am not forcing one engine to do both jobs.
- **`notebooks/` is explicitly exploratory-only** and never imported by production code —
  called out in its own README so this constraint doesn't erode over time.
- **Reproducibility is layered:** `environment/env.yml` (conda/mamba) pins the full
  scientific + R-bridge stack; `pyproject.toml` layers standard Python packaging on top
  for the installable `omicsatlas` and `crossomicsconcordance` packages;
  `environment/Dockerfile` pins a container image for full end-to-end reproduction
  without relying on the host having conda at all.
- **Python version: 3.11**, pinned across `pyproject.toml` and `environment/env.yml`.
  Chosen over 3.12 for the safest current compatibility with `rpy2` and the Bioconda
  package set the R-bridge and deconvolution tooling depend on.

## Consequences

- Reusing `CrossOmicsConcordance` elsewhere means `pip install`-ing a path/subdirectory
  or, once published, the PyPI package — not the whole monorepo. This is the entire
  point of Phase 5's packaging work.
- CI (added in this phase) needs to run at least two toolchains (Python lint/test, R CMD
  check) rather than one; slightly more CI complexity in exchange for one coherent repo
  history, which matters for this project's role as a portfolio evidence trail.
- Anyone building from source needs both conda/mamba and R available (or the Docker
  image) — documented in the root README's quick start.

## Alternatives considered

- **Multi-repo split** (separate repos per modality + concordance package): rejected —
  fragments the reproducibility story and the commit history that's part of this
  project's evidence value, for a benefit (independent CI/versioning) that only really
  matters at multi-contributor scale.
- **GPL-family license**: rejected — copyleft terms are unnecessary friction for a
  library-style artifact (`CrossOmicsConcordance`) I want maximally reusable, and would
  be inconsistent with the permissive licenses of the ecosystem it plugs into.
