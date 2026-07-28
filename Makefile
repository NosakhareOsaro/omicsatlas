.PHONY: setup lint test run clean

setup:
	mamba env create -f environment/env.yml
	@echo "Run 'conda activate omicsatlas && pre-commit install' next."

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
