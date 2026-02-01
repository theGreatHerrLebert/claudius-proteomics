# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**San José** is a reproducible, metadata-rich, bias-aware reference layer for timsTOF data on PRIDE.

> We are not collecting peptides. We are collecting **peptide observations in experimental context**.

The pipeline systematically reprocesses public timsTOF data from PRIDE through triple orthogonal validation (FragPipe + DIA-NN + Sage), extracts raw signal features via rustims/imspy, and builds versioned snapshots for model training.

**Two phases:**
- **Phase 1 (San José Database)**: Download → Triple-engine search → Raw extraction → Bias-aware database
- **Phase 2 (Modeling)**: Train prediction models (CCS, RT, MS2) on versioned database snapshots

Primary target: CCS prediction from timsTOF data, benchmarking against Meier et al. 2021 (R > 0.99).

## Common Commands

```bash
# Show all available targets
snakemake help

# Phase 1: Database building
snakemake add_dataset --config accession=PXD019086 --cores 16
snakemake add_all_datasets --cores 16
snakemake create_snapshot --config version=v1.0

# Search methods (triple orthogonal validation)
snakemake process_fragpipe --cores 16  # FragPipe (spectrum-centric)
snakemake process_diann --cores 16     # DIA-NN (peptide-centric)
snakemake process_sage --cores 16      # Sage (fast, independent Rust engine)
snakemake process_all --cores 16       # All three (recommended)

# Phase 2: Model training
snakemake train_model --config snapshot=v1.0 model=ccs_v1 --cores 16
snakemake evaluate_model --config snapshot=v1.0 model=ccs_v1 --cores 16

# HPC execution (Mogon2 SLURM)
snakemake add_all_datasets --profile profiles/mogon2

# Utilities
snakemake db_stats
snakemake list_snapshots
snakemake prepare_fasta --config accession=PXD019086
snakemake clean_temp
```

## Local Development Setup

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install snakemake
./scripts/setup_imspy.sh  # Build imspy from source (requires Rust)
```

Test mode is enabled by default in config (`test_mode.max_files: 3`).

## Architecture

### Technology Stack
- **Orchestration**: Snakemake v8+
- **Raw processing**: rustims/imspy (Rust backend + Python bindings)
- **Search engines**: FragPipe (spectrum-centric) + DIA-NN (peptide-centric) + Sage (fast Rust engine)
- **Database**: DuckDB/Parquet (columnar, versioned snapshots)
- **ML**: TensorFlow/keras via imspy-predictors
- **Containers**: Singularity (HPC) / Docker (local)
- **HPC**: SLURM (Mogon2 @ JGU Mainz)

### Map-Reduce Pattern

**Map (Runner)**: N runners process N datasets in parallel, no communication
```
PRIDE accession + software config + FASTA → Runner → standardized output
```

**Reduce (Aggregator)**: Combines runner outputs into unified database + snapshots

### Runner Output Schema

See `docs/RUNNER_OUTPUT_SCHEMA.md` for full specification.

```
{accession}/
├── manifest.json      # Metadata, versions, QC metrics
├── fragpipe.parquet   # FragPipe PSMs
├── diann.parquet      # DIA-NN results
├── sage.parquet       # Sage results
├── consensus.parquet  # Engine agreement (union, 2/3, 3/3)
└── spectra/
    ├── index.parquet  # Lightweight index
    └── blobs.bin      # PrecursorWithFragments blobs
```

### Modified Sequence Format (UniMod)

Canonical format across all outputs:
```
[N-term]-SEQUENCE[UNIMOD:ID]-[C-term]

Examples:
[]-PEPTIDEK-[]                      # Unmodified
[]-AAC[UNIMOD:4]PEPTIDEK-[]         # Carbamidomethyl (C)
[]-M[UNIMOD:35]PEPTIDEK-[]          # Oxidation (M)
[UNIMOD:1]-PEPTIDEK-[]              # N-term Acetyl
```

Join key: `(modified_sequence, charge)` → unique precursor identity

### Raw Data Storage (Parquet)

See `docs/RAW_DATA_ARCHITECTURE.md` for full specification.

Precursor data is stored in Parquet with list columns for variable-length arrays:
- Full 4D raw MS1 signal (`raw_rt`, `raw_mz`, `raw_mobility`, `raw_intensity`)
- 1D projections (XIC, mobilogram, isotope envelope)
- Merged fragment spectrum
- Search engine identifications

Benefits: columnar reads, predicate pushdown, batch iteration for training.

### Dashboard

Interactive precursor browser at `dashboard/`:
- **Backend**: FastAPI serving Parquet data (`dashboard/backend/`)
- **Frontend**: React + TypeScript + Tailwind + Deck.gl (`dashboard/frontend/`)

```bash
# Run backend (use venv python directly - no activate script)
.venv/bin/python dashboard/backend/main.py --store data/processed/PXD019086/precursor_index_v3.parquet --port 8000

# Run frontend - Option 1: Serve pre-built files (avoids file watcher limit)
cd dashboard/frontend/dist && python3 -m http.server 5173

# Run frontend - Option 2: Dev server (requires Node 20+ via nvm)
cd dashboard/frontend && nvm use 22 && npm run dev  # http://localhost:5173

# If dev server fails with ENOSPC (file watcher limit), increase the limit:
sudo sysctl fs.inotify.max_user_watches=524288
```

Open http://localhost:5173 in browser to access the precursor browser.

### Rule Modules (`rules/`)
| Module | Purpose |
|--------|---------|
| `fasta.smk` | FASTA resource management, organism mapping, decoy generation |
| `download.smk` | PRIDE download or local data linking |
| `fragpipe.smk` | FragPipe spectrum-centric search |
| `diann.smk` | DIA-NN peptide-centric search |
| `sage.smk` | Sage fast Rust search engine |
| `extract.smk` | imspy raw feature extraction |
| `database.smk` | Database insertion & aggregation |
| `snapshot.smk` | Versioned snapshot creation |
| `train.smk` | Model training & evaluation |

### Key Design Decisions
- **Triple orthogonal validation**: FragPipe + DIA-NN + Sage for consensus and disagreement analysis
- **Bias-awareness**: Track lab_id, dataset_id, organism, gradient_length, column_type, acquisition_mode
- **Dual extraction**: PSM identifications + raw signal features (full 4D data)
- **Database partitioning**: One dataset = one PRIDE accession
- **Containerization**: FragPipe/DIA-NN binaries mounted at runtime (licensing)
- **Human checkpoints**: Dataset selection, QC review, consensus rules, snapshot approval, release gates

## Configuration

Edit `config/config.yaml`:
- `datasets`: List of PRIDE accessions to process
- `dataset_metadata`: Maps accession → organism
- `organisms`: FASTA sources (local path or UniProt proteome ID)
- `fragpipe.path`: User-provided FragPipe installation (required)
- `diann.path`: User-provided DIA-NN installation (required)
- `sage.path`: Sage binary path (required)
- `san_jose.report_all`: Report ALL PSMs without FDR filtering (default: true)
- `test_mode.enabled`: Quick validation with limited files

## External Dependencies

**FragPipe**, **DIA-NN**, and **Sage** require separate setup:
- FragPipe: https://fragpipe.nesvilab.org/ (academic license)
- DIA-NN: https://github.com/vdemichev/DiaNN (academic license)
- Sage: https://github.com/lazear/sage (open source, Rust)

Set paths in `config/config.yaml` before running.

**Important**: All search engines configured to report ALL results (no FDR or PEP filtering). San José stores full union with per-engine scores; thresholds applied at query time for downstream analysis.

## Key Documentation

- `docs/SAN_JOSE_PITCH.md` - Vision, architecture, and roadmap for San José
- `docs/RAW_DATA_ARCHITECTURE.md` - Raw 4D data storage, Parquet schema, dashboard API
- `docs/RUNNER_OUTPUT_SCHEMA.md` - Parquet schemas, manifest format
- `docs/PRECURSOR_BLOB_DESIGN.md` - 4D precursor + fragment blob structure (legacy)
- `docs/DATASET_DEFINITION.md` - Dataset semantics and edge cases
- `CLAUDIUS-PROTEOMICS.md` - Full project plan and technical details
