# python:3.11-slim-bookworm, pinned by digest for reproducibility.
FROM python@sha256:b18992999dbe963a45a8a4da40ac2b1975be1a776d939d098c647482bcad5cba

WORKDIR /opt/omicsatlas

COPY pyproject.toml ./
COPY src ./src
COPY README.md ./

RUN pip install --no-cache-dir -e .

COPY . .

CMD ["python", "-c", "import omicsatlas; print(omicsatlas.__version__)"]
