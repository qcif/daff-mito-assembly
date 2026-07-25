# wf5 custom-logic image.
# Contains the Python interpreter + pinned deps from requirements.txt.
# Python source is NOT baked in — scripts are staged at runtime via
# Nextflow's bin/ auto-staging (see plan.md §2.2 container pattern).
#
# Rebuild only when requirements.txt or this Dockerfile changes (CI images.yml).

FROM python:3.12-slim

LABEL org.opencontainers.image.source="https://github.com/qcif/daff-mito-assembly"
LABEL org.opencontainers.image.description="wf5 custom logic runtime"

WORKDIR /opt/wf5

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Smoke-test: ensure biopython imported cleanly
RUN python -c "from Bio import SeqIO; print('biopython ok')"
