# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CLAUDIUS-PROTEOMICS is a two-phase peptide property prediction pipeline:
- **Phase 1 (Database)**: Downloads public proteomics data from PRIDE, processes through FragPipe and DIA-NN, extracts raw signal features via imspy, and builds an extensible peptide property database
- **Phase 2 (Modeling)**: Trains prediction models (primarily CCS) on versioned database snapshots using imspy-predictors

Primary target: Collisional Cross Section (CCS) prediction from timsTOF data, benchmarking against Meier et al. 2021 (R > 0.99).

## Common Commands

```bash
# Show all available targets
snakemake help

# Phase 1: Database building
snakemake add_dataset --config accession=PXD019086 --cores 16
snakemake add_all_datasets --cores 16
snakemake create_snapshot --config version=v1.0

# Search methods (orthogonal validation)
snakemake process_only --cores 16    # FragPipe (spectrum-centric)
snakemake process_diann --cores 16   # DIA-NN (peptide-centric)
snakemake process_both --cores 16    # Both (recommended)

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
- **Search engines**: FragPipe/MSFragger (spectrum-centric) + DIA-NN v2.3+ (peptide-centric)
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

### Precursor Blob Design

See `docs/PRECURSOR_BLOB_DESIGN.md` for full specification.

Each precursor stored as compressed blob containing:
- Full 4D precursor signal (XIC, mobilogram, isotope envelope)
- Merged fragment spectrum (TimsFrame)
- Metadata (charge, m/z, collision energies)

### Rule Modules (`rules/`)
| Module | Purpose |
|--------|---------|
| `fasta.smk` | FASTA resource management, organism mapping, decoy generation |
| `download.smk` | PRIDE download or local data linking |
| `fragpipe.smk` | FragPipe spectrum-centric search |
| `diann.smk` | DIA-NN peptide-centric search |
| `extract.smk` | imspy raw feature extraction |
| `database.smk` | Database insertion & aggregation |
| `snapshot.smk` | Versioned snapshot creation |
| `train.smk` | Model training & evaluation |

### Key Design Decisions
- **Orthogonal validation**: Both FragPipe and DIA-NN run for high-confidence consensus
- **Dual extraction**: PSM identifications + raw signal features (full 4D data)
- **Database partitioning**: One dataset = one PRIDE accession
- **Containerization**: FragPipe/DIA-NN binaries mounted at runtime (licensing)

## Configuration

Edit `config/config.yaml`:
- `datasets`: List of PRIDE accessions to process
- `dataset_metadata`: Maps accession → organism
- `organisms`: FASTA sources (local path or UniProt proteome ID)
- `fragpipe.path`: User-provided FragPipe installation (required)
- `diann.path`: User-provided DIA-NN installation (required)
- `test_mode.enabled`: Quick validation with limited files

## External Dependencies

**FragPipe** and **DIA-NN** require separate download (academic licensing):
- FragPipe: https://fragpipe.nesvilab.org/
- DIA-NN: https://github.com/vdemichev/DiaNN

Set paths in `config/config.yaml` before running.

## Key Documentation

- `docs/RUNNER_OUTPUT_SCHEMA.md` - Parquet schemas, manifest format
- `docs/PRECURSOR_BLOB_DESIGN.md` - 4D precursor + fragment blob structure
- `docs/DATASET_DEFINITION.md` - Dataset semantics and edge cases
- `CLAUDIUS-PROTEOMICS.md` - Full project plan and vision
