FROM python:3.12-slim

# Nextflow shells out to `ps` for task metrics; slim images do not ship procps.
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && rm /tmp/requirements.txt

# numba caches compiled code next to the package; the Batch task user may not own it.
ENV NUMBA_CACHE_DIR=/tmp/numba
