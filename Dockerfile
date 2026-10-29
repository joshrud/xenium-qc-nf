# Container for every xenium-qc-nf step: the xenium_qc package (xqc CLI) plus MultiQC with the Xenium Extra plugin.
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 PIP_NO_CACHE_DIR=1 MPLBACKEND=Agg
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/xenium-qc
COPY pyproject.toml README.md ./
COPY xenium_qc ./xenium_qc
RUN pip install ".[multiqc]"

# Nextflow needs a writable HOME for matplotlib/numba caches when running as a non-root user
ENV HOME=/tmp NUMBA_CACHE_DIR=/tmp MPLCONFIGDIR=/tmp
CMD ["xqc", "--help"]
