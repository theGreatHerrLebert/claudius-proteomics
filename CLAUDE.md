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
# Runner Pipeline (6-step processing)
.venv/bin/python runner/run_dataset.py PXD019086 --config config/config.yaml

# Test mode (limited files for quick validation)
.venv/bin/python runner/run_dataset.py PXD019086 --test-mode --max-files 1

# Use local data (skip download)
.venv/bin/python runner/run_dataset.py PXD019086 --local-data data/raw/PXD019086

# Resume from checkpoint
.venv/bin/python runner/run_dataset.py PXD019086 --resume

# Create distributable archive (step 6)
.venv/bin/python runner/run_dataset.py PXD019086 --package --package-version 1.0

# Collection management
python scripts/add_to_collection.py --archive data/packages/PXD019086_v1.0.zip \
    --collection /path/to/collection --study meier_2021_ccs

# Snakemake commands (legacy/alternative)
snakemake process_fragpipe --cores 16  # FragPipe (spectrum-centric)
snakemake process_diann --cores 16     # DIA-NN (peptide-centric)
snakemake process_sage --cores 16      # Sage (fast, independent Rust engine)

# Phase 2: Model training
snakemake train_model --config snapshot=v1.0 model=ccs_v1 --cores 16

# HPC execution (Mogon2 SLURM)
snakemake add_all_datasets --profile profiles/mogon2
```

## Local Development Setup

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install snakemake
./scripts/setup_imspy.sh  # Build imspy from source (requires Rust)
```

Test mode is enabled by default in config (`test_mode.max_files: 1`).

## 6-Step Runner Pipeline

The runner (`runner/run_dataset.py`) processes a PRIDE dataset through 6 steps with checkpointing:

| Step | Name | Input | Output |
|------|------|-------|--------|
| 1 | **Download** | PRIDE accession | `data/raw/{acc}/*.d` |
| 2 | **Search** | .d files + FASTA | `data/processed/{acc}/{engine}/` |
| 3 | **Stratify** | Engine results | `precursor_index.parquet`, `consensus/` |
| 4 | **Extract** | .d files + index | `data/extracted/{acc}/`, raw 4D features |
| 5 | **Merge** | Index + features | `data/merged/{acc}/precursor_store.parquet` |
| 6 | **Package** | All outputs | `data/packages/{acc}_v{ver}.zip` (optional, `--package`) |

### Folder Structure

```
data/
├── raw/{accession}/                    # Step 1: Bruker .d folders
├── processed/{accession}/              # Steps 2-3
│   ├── fragpipe_output/                # FragPipe PSMs
│   ├── diann/                          # DIA-NN results
│   ├── sage/                           # Sage results
│   ├── consensus/stratified/           # Precursors by engine agreement
│   └── precursor_index.parquet         # Unified index
├── extracted/{accession}/              # Step 4: Raw 4D signal
│   └── {raw_file}.d/blobs.bin          # Raw 4D data blobs
├── merged/{accession}/                 # Step 5: Final output
│   ├── precursor_store.parquet         # IDs + raw features joined
│   └── manifest.json                   # QC metrics
├── packages/                           # Step 6: Distributable archives
│   └── {accession}_v{version}.zip
└── checkpoints/{accession}/            # Pipeline state
    └── state.json
```

### Archive Format (Step 6)

Self-contained zip archive for distribution/upload:
```
{accession}_v{version}.zip
└── {accession}/
    ├── manifest.json              # Metadata, QC, versions
    ├── precursor_store.parquet    # Main data (IDs + features)
    ├── precursor_index.parquet    # Engine-level details
    ├── consensus/                 # Stratified results
    ├── extracted/{raw}.d/         # Raw 4D blobs + index
    ├── engines/                   # Search engine outputs
    │   ├── fragpipe/              # Configs + combined + per-file results
    │   ├── diann/                 # report.parquet + log + stats
    │   └── sage/                  # results + fragments + settings
    └── summaries/                 # Step summaries
```

## Architecture

### Technology Stack
- **Orchestration**: Snakemake v8+ / Python runner
- **Raw processing**: rustims/imspy (Rust backend + Python bindings)
- **Search engines**: FragPipe (spectrum-centric) + DIA-NN (peptide-centric) + Sage (fast Rust engine)
- **Database**: DuckDB/Parquet (columnar, versioned snapshots)
- **ML**: TensorFlow/keras via imspy-predictors
- **Containers**: Singularity (HPC) / Docker (local)
- **HPC**: SLURM (Mogon2 @ JGU Mainz)

### Search Engine Configuration

All three engines process Bruker `.d` folders directly (no mzML conversion needed):

| Engine | Approach | Key Flags | Notes |
|--------|----------|-----------|-------|
| **FragPipe** | Spectrum-centric | `--no-fdr-filter` | FDR=1.0 to report all PSMs |
| **DIA-NN** | Peptide-centric | `--dda --fasta-search` | DDA mode + library-free search |
| **Sage** | Spectrum-centric | `report_psms: 1` | In config JSON |

**DIA-NN note**: Library generation is slow (~40 min) but only needed once per FASTA. The library is reused for subsequent searches.

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
- Merged fragment spectrum with intensity normalization
- Fragment ion annotations (b/y ion matching)
- Search engine identifications from all three engines

Benefits: columnar reads, predicate pushdown, batch iteration for training.

### Training Spectra

The `training_spectra.parquet` extends base precursor data with:
- **Normalized intensities**: base_peak, tic, sqrt, or log normalization
- **Fragment annotations**: b/y ion type, position, charge, theoretical m/z, error in ppm
- **Coverage metrics**: intensity_explained, sequence_coverage_b, sequence_coverage_y

Build with:
```bash
python scripts/build_training_spectra.py \
    --input data/processed/PXD019086/precursors.parquet \
    --output data/processed/PXD019086/training_spectra.parquet \
    --normalization base_peak --mz-tolerance 20.0
```

### Dashboard

Interactive precursor browser at `dashboard/`:
- **Backend**: FastAPI serving Parquet data (`dashboard/backend/`)
- **Frontend**: React + TypeScript + Tailwind + Deck.gl (`dashboard/frontend/`)

```bash
# Run backend - Single dataset mode
.venv/bin/python dashboard/backend/main.py --store data/merged/PXD019086/precursor_store.parquet --port 8000

# Run backend - Collection mode (browse multiple datasets)
.venv/bin/python dashboard/backend/main.py --collection /path/to/san_jose_collection/ --port 8000

# Run frontend - Option 1: Serve pre-built files (avoids file watcher limit)
cd dashboard/frontend/dist && python3 -m http.server 5173

# Run frontend - Option 2: Dev server (requires Node 20+ via nvm)
cd dashboard/frontend && nvm use 22 && npm run dev  # http://localhost:5173

# If dev server fails with ENOSPC (file watcher limit), increase the limit:
sudo sysctl fs.inotify.max_user_watches=524288
```

**Collection mode endpoints:**
- `GET /studies` - List all studies
- `GET /studies/{id}/datasets` - List datasets in a study
- `POST /datasets/{accession}/load` - Load dataset into memory
- `GET /collection` - Collection metadata

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

### Key Scripts (`scripts/`)
| Script | Purpose |
|--------|---------|
| `run_fragpipe.py` | FragPipe wrapper with San José mode (FDR override) |
| `build_precursor_index.py` | Bridges raw timsTOF data to search engine results |
| `add_to_collection.py` | Add archives to collection, rebuild manifest |
| `sequence_utils.py` | UNIMOD normalization across engines (FragPipe/DIA-NN/Sage) |
| `fragment_matching.py` | Theoretical b/y ion matching for spectrum annotation |
| `build_training_spectra.py` | Adds normalized intensities and fragment annotations |
| `precursor_store.py` | Parquet store utilities for precursor data |
| `precursor_store_parquet.py` | Full 4D extraction to Parquet with list columns |
| `analyze_overlap.py` | Engine agreement/disagreement analysis |

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
