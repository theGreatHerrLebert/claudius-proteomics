# San José Worker Container
# Self-contained unit for processing a single PRIDE dataset
#
# Usage:
#   docker build -t san-jose-worker .
#   docker run -v /path/to/fragpipe:/opt/fragpipe \
#              -v /path/to/diann:/opt/diann \
#              -v /path/to/output:/output \
#              san-jose-worker PXD019086

FROM python:3.12-slim

LABEL maintainer="San José Project"
LABEL description="Distributed worker for peptide property extraction"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    unzip \
    openjdk-17-jre-headless \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies
COPY requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Create directory structure
RUN mkdir -p /app /data /output /fastas /config
WORKDIR /app

# Copy worker scripts and config
COPY scripts/ /app/scripts/
COPY config/ /app/config/
COPY entrypoint.sh /app/

RUN chmod +x /app/entrypoint.sh /app/scripts/*.py

# Environment variables (can be overridden at runtime)
ENV FRAGPIPE_PATH=/opt/fragpipe
ENV DIANN_PATH=/opt/diann/diann-linux
ENV FASTA_DIR=/fastas
ENV OUTPUT_DIR=/output
ENV DATA_DIR=/data
ENV JOB_QUEUE_URL=""
ENV RESULTS_UPLOAD_URL=""
ENV WORKER_ID=""

# Entrypoint
ENTRYPOINT ["/app/entrypoint.sh"]
