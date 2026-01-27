# San José Worker

Self-contained container for processing PRIDE datasets autonomously.

## Overview

The San José Worker is a distributable unit that:
1. Downloads raw data from PRIDE
2. Fetches appropriate FASTA (based on organism)
3. Runs FragPipe (spectrum-centric search)
4. Runs DIA-NN (peptide-centric search)
5. Packages and uploads results

Workers can run standalone or poll a job queue for continuous processing.

## Quick Start

### Docker

```bash
# Build
docker build -t san-jose-worker .

# Process single dataset
docker run \
  -v /path/to/fragpipe:/opt/fragpipe \
  -v /path/to/diann:/opt/diann \
  -v /path/to/output:/output \
  san-jose-worker PXD019086
```

### Singularity (HPC)

```bash
# Build
singularity build san-jose-worker.sif Singularity.def

# Process single dataset
singularity run \
  --bind /path/to/fragpipe:/opt/fragpipe \
  --bind /path/to/diann:/opt/diann \
  --bind /path/to/output:/output \
  san-jose-worker.sif PXD019086
```

## Prerequisites

You must provide your own installations of:

1. **FragPipe** (v24.0+)
   - Download from: https://fragpipe.nesvilab.org/
   - Requires academic license
   - Mount at `/opt/fragpipe`

2. **DIA-NN** (v2.3.2+)
   - Download from: https://github.com/vdemichev/DiaNN
   - Mount at `/opt/diann`

## Modes of Operation

### 1. Single Dataset

Process one PRIDE accession:

```bash
./entrypoint.sh PXD019086
```

### 2. Job Manifest

Process from a pre-defined job file:

```bash
./entrypoint.sh --job /path/to/job.yaml
```

### 3. Daemon Mode

Poll a job queue continuously:

```bash
export JOB_QUEUE_URL="s3://san-jose-jobs/jobs/"
export RESULTS_UPLOAD_URL="s3://san-jose-results/"
export WORKER_ID="worker-001"

./entrypoint.sh --daemon
```

## Job Queue Backends

### S3 Bucket

```
s3://san-jose-jobs/
├── jobs/
│   ├── pending/        # Jobs waiting to be claimed
│   │   └── PXD019086.yaml
│   ├── processing/     # Jobs being processed
│   │   └── worker-001/
│   │       └── PXD019086.yaml
│   └── completed/      # Finished jobs
│       └── PXD019086.yaml
└── results/            # Output data
    └── PXD019086/
        ├── fragpipe_psm.tsv
        └── diann_report.parquet
```

### REST API

Implement endpoints:
- `POST /claim` - Claim next available job
- `POST /upload` - Upload results
- `POST /complete` - Mark job complete

### Local Directory

For testing:

```bash
export JOB_QUEUE_URL="file:///path/to/jobs"
./entrypoint.sh --daemon
```

## Job Manifest Schema

```yaml
job_id: "PXD019086_20240127"
accession: "PXD019086"
organism: "human"

fasta:
  organism: "human"
  proteome_id: "UP000005640"
  url: "https://..."
  include_contaminants: true
  generate_decoys: true

raw_files:
  - name: "sample1.d.zip"
    url: "ftp://..."
    size_mb: 2500

search_config:
  fragpipe:
    workflow: "LFQ-MBR"
    threads: 16
  diann:
    threads: 16
    qvalue: 0.01

output:
  format: "parquet"
```

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `FRAGPIPE_PATH` | FragPipe installation | `/opt/fragpipe` |
| `DIANN_PATH` | DIA-NN binary | `/opt/diann/diann-linux` |
| `FASTA_DIR` | FASTA cache directory | `/fastas` |
| `DATA_DIR` | Raw data directory | `/data` |
| `OUTPUT_DIR` | Results output | `/output` |
| `JOB_QUEUE_URL` | Job queue URL (daemon mode) | - |
| `RESULTS_UPLOAD_URL` | Results destination | - |
| `WORKER_ID` | Unique worker identifier | hostname |

## Output Structure

```
/output/PXD019086/
├── fragpipe_psm.tsv      # FragPipe PSM results
├── diann_report.parquet  # DIA-NN results
├── job.yaml              # Job manifest copy
└── COMPLETED             # Completion marker
```

## Resource Requirements

| Stage | CPU | RAM | Disk |
|-------|-----|-----|------|
| Download | 1 | 1GB | varies |
| FragPipe | 16 | 64GB | 50GB |
| DIA-NN | 16 | 32GB | 20GB |

Typical runtime: 2-8 hours per dataset (depends on size).

## Distributing to Collaborators

1. Share this worker package
2. Collaborator installs FragPipe + DIA-NN (licensing)
3. Collaborator builds container
4. Provide job queue credentials or job manifest files
5. Collaborator runs worker, results upload automatically

## Troubleshooting

**FragPipe not found:**
```
Mount your installation: -v /path/to/fragpipe-24.0:/opt/fragpipe
```

**DIA-NN not found:**
```
Mount your installation: -v /path/to/diann-2.3:/opt/diann
```

**Out of memory:**
```
Reduce threads in job manifest or increase container memory limit
```

**FASTA download fails:**
```
Check internet connectivity or provide pre-downloaded FASTA
```
