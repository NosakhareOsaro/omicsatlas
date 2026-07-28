.PHONY: setup lint test run clean

setup:
	mamba env create -f environment/env.yml
	@echo "Now: conda activate omicsatlas"
	@echo "Then: pip install -e '.[dev]' rpy2==3.6.7 pertpy==1.0.3 scikit-learn==1.9.0 pydeseq2==0.5.4"
	@echo "Then: pre-commit install"
	@echo "Apple Silicon: prefix the first line with CONDA_SUBDIR=osx-64 (see environment/env.yml)."

lint:
	ruff check .
	black --check .
	mypy

test:
	pytest

run:
	@echo "No end-to-end entry point yet; added with the Nextflow master workflow in Phase 6."

clean:
	rm -rf .pytest_cache .mypy_cache .ruff_cache htmlcov .coverage
	find . -name "__pycache__" -type d -prune -exec rm -rf {} +
