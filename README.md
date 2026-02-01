# San José

> A reproducible, metadata-rich, bias-aware reference layer for timsTOF data on PRIDE.

## Vision

**We are not collecting peptides. We are collecting peptide observations in experimental context.**

San José systematically reprocesses every timsTOF dataset on PRIDE to build a versioned, citable dataset for CCS, RT, and MS2 model training — grounded in real experimental physics across thousands of labs.

```
PRIDE Archive                    San José Database              Models
┌─────────────┐                 ┌─────────────────┐           ┌─────────────┐
│ PXD019086   │ ──┐             │ precursors/     │           │ ccs_v1/     │
│ PXD010012   │   │  download   │   PXD019086     │ snapshot  │   model.h5  │
│ PXD...      │   ├──────────►  │   PXD010012     │ ────────► │   metrics   │
│             │   │  search     │ snapshots/      │  train    │             │
└─────────────┘ ──┘  extract    │   v1.0/         │           └─────────────┘
                                └─────────────────┘
```

See [CLAUDIUS-PROTEOMICS.md](CLAUDIUS-PROTEOMICS.md) for the full project plan and [docs/SAN_JOSE_PITCH.md](docs/SAN_JOSE_PITCH.md) for the vision.

## Key Design Principles

### Triple Orthogonal Validation

Each dataset is processed with three independent search engines:

| Engine | Approach | Strength |
|--------|----------|----------|
| **FragPipe** | Spectrum-centric | Gold standard for DDA |
| **DIA-NN** | Peptide-centric | Best for DIA, deep learning rescoring |
| **Sage** | Fast Rust engine | Independent validation, cheap at scale |

We store: full union, 2/3 agreement, 3/3 agreement, and engine-specific unique IDs. Disagreements are scientifically valuable.

### Bias-Awareness

PRIDE datasets cluster heavily by lab, sample prep, column chemistry, and gradient length. San José tracks:

```
lab_id, dataset_id, organism, gradient_length, column_type, acquisition_mode
```

This enables stratified sampling, bias analysis, and cross-lab validation.

### San José Mode (Report ALL Results)

All search engines are configured to report ALL PSMs without FDR filtering. San José stores the full union with per-engine scores; thresholds are applied at query time for downstream analysis.

### Human Checkpoints

| Checkpoint | Purpose |
|------------|---------|
| Dataset Selection | Curated PRIDE accessions, blacklist for problematic datasets |
| QC Dashboard | Auto-generated reports, anomaly detection, human review when flagged |
| Consensus Rules | Configurable: 3/3 vs 2/3 vs union, PEP thresholds |
| Versioned Snapshots | San José v1.0, v1.1… frozen, reproducible datasets |
| Release Gates | Bias checks, cross-lab validation, human sign-off |

## Quick Start

### Prerequisites

1. **FragPipe** (v24+) - https://fragpipe.nesvilab.org/ (academic license)
2. **DIA-NN** (v2.3+) - https://github.com/vdemichev/DiaNN (academic license)
3. **Sage** - https://github.com/lazear/sage (open source Rust)
4. **Snakemake** (v8+) - Workflow orchestration
5. **Singularity** (optional) - For containerized HPC execution

### Configuration

Edit `config/config.yaml`:

```yaml
# Search engine paths
fragpipe:
  path: "/path/to/fragpipe-24.0"
  workflow: "LFQ-MBR"

diann:
  path: "/path/to/diann-linux"
  threads: 16
  qvalue: 1.0  # San José: report ALL results

sage:
  path: "/path/to/sage"
  config: "config/sage_config.json"

# San José mode (default)
san_jose:
  report_all: true

# Organism FASTA mapping
organisms:
  human:
    proteome_id: "UP000005640"
    local_fasta: "/path/to/human.fasta"

# Dataset to organism mapping
dataset_metadata:
  PXD019086:
    organism: "human"
```

### Running the Pipeline

```bash
# Phase 1: Triple-engine search
snakemake process_fragpipe --cores 16  # Spectrum-centric
snakemake process_diann --cores 16     # Peptide-centric
snakemake process_sage --cores 16      # Fast Rust engine
snakemake process_all --cores 16       # All three (recommended)

# Build precursor index (bridges raw data to search results)
snakemake build_precursor_index --config accession=PXD019086

# Create versioned snapshot
snakemake create_snapshot --config version=v1.0

# Phase 2: Train models
snakemake train_model --config snapshot=v1.0 model=ccs_v1
```

### Interactive Dashboard

Browse precursors with full 4D visualization:

```bash
# Backend (FastAPI)
.venv/bin/python dashboard/backend/main.py \
    --store data/processed/PXD019086/precursor_index_v3.parquet \
    --port 8000

# Frontend (React + TypeScript)
cd dashboard/frontend && npm run dev  # http://localhost:5173
```

Features: fragment spectrum, IM vs m/z scatter, XIC, mobilogram, isotope envelope.

## Project Structure

```
claudius-proteomics/
├── Snakefile                   # Main workflow
├── config/
│   ├── config.yaml             # Pipeline configuration
│   └── sage_config.json        # Sage search parameters
├── profiles/mogon2/            # HPC cluster profile
│
├── rules/                      # Snakemake rule modules
│   ├── fasta.smk               # FASTA resource management
│   ├── download.smk            # PRIDE download
│   ├── fragpipe.smk            # FragPipe search
│   ├── diann.smk               # DIA-NN search
│   ├── sage.smk                # Sage search
│   ├── extract.smk             # Raw feature extraction
│   ├── database.smk            # Database insertion
│   ├── snapshot.smk            # Snapshot creation
│   └── train.smk               # Model training
│
├── scripts/
│   ├── run_fragpipe.py         # FragPipe wrapper (San José mode)
│   ├── build_precursor_index.py # Unified precursor index
│   ├── sequence_utils.py       # UNIMOD normalization
│   ├── fragment_matching.py    # Theoretical fragment matching
│   ├── precursor_store.py      # Parquet store utilities
│   └── ...
│
├── dashboard/
│   ├── backend/main.py         # FastAPI server
│   └── frontend/               # React + TypeScript + Deck.gl
│
├── docs/
│   ├── SAN_JOSE_PITCH.md       # Vision and architecture
│   ├── RAW_DATA_ARCHITECTURE.md # Parquet schema, dashboard API
│   ├── RUNNER_OUTPUT_SCHEMA.md # Per-dataset output format
│   └── BUGS_FOUND.md           # Bug documentation
│
├── data/
│   ├── raw/                    # Downloaded .d files
│   └── processed/{accession}/  # Search results + precursor store
│
├── CLAUDIUS-PROTEOMICS.md      # Full project plan
└── CLAUDE.md                   # Claude Code instructions
```

## Data Storage

### Precursor Parquet Schema

Each precursor observation contains:

```
# Identity
precursor_id, raw_file, mz, charge, rt_seconds, mobility

# Search engine results (NULL if not identified)
fragpipe_peptide, sage_peptide, diann_peptide
consensus_peptide, n_engines

# Quality metrics
fragpipe_probability, diann_qvalue, sage_qvalue

# Fragment spectrum (variable length lists)
fragment_mz, fragment_intensity, fragment_mobility

# MS1 projections
xic_rt, xic_intensity
mobilogram_im, mobilogram_intensity
isotope_mz, isotope_intensity

# Raw 4D MS1 data
raw_rt, raw_mz, raw_mobility, raw_intensity
```

See [docs/RAW_DATA_ARCHITECTURE.md](docs/RAW_DATA_ARCHITECTURE.md) for full specification.

## Technical Stack

- **Orchestration**: Snakemake v8+
- **Raw processing**: rustims/imspy (Rust backend + Python bindings)
- **Search engines**: FragPipe + DIA-NN + Sage
- **Database**: DuckDB/Parquet (columnar, versioned snapshots)
- **Dashboard**: FastAPI + React + TypeScript + Deck.gl (WebGL)
- **ML**: TensorFlow/keras via imspy-predictors
- **Containers**: Singularity (HPC) / Docker (local)
- **HPC**: SLURM (Mogon2 @ JGU Mainz)

## Data Sources

### POC: Meier et al. 2021

- **Paper**: [Deep learning the collisional cross sections of the peptide universe](https://doi.org/10.1038/s41467-021-21352-8)
- **PRIDE**: PXD019086 (primary), PXD010012, PXD017703
- **Data**: >1 million CCS measurements from timsTOF
- **Target**: R > 0.99, <1.4% median error

## Why "San José"?

Named after the *San José*, a sunken ship whose rediscovery revealed immense, carefully preserved value beneath the surface.

PRIDE is similar: a vast repository where the real value is hidden in raw experimental data, waiting to be systematically recovered, catalogued, and understood.

## License

[Add license]

## Citation

[Add citation when published]
