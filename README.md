# CLAUDIUS-PROTEOMICS

Building the world's largest peptide property prediction resource through systematic reanalysis of public proteomics data.

## Overview

This project provides a two-phase pipeline:

**Phase 1: Database Building** - Download, process, and extract peptide properties from public proteomics data (PRIDE). The database is the primary output - an extensible resource of peptide properties.

**Phase 2: Model Training** - Train prediction models on database snapshots. Models are versioned and reproducible.

```
PRIDE Archive                    Database                      Models
┌─────────────┐                 ┌─────────────┐               ┌─────────────┐
│ PXD019086   │ ──┐             │ peptides/   │               │ ccs_v1/     │
│ PXD010012   │   │  download   │   PXD019086 │   snapshot    │   model.h5  │
│ PXD...      │   ├──────────►  │   PXD010012 │ ───────────►  │   metrics   │
│             │   │  process    │   ...       │   train       │             │
└─────────────┘ ──┘  extract    │ snapshots/  │               └─────────────┘
                                │   v1.0/     │
                                └─────────────┘
```

See [CLAUDIUS-PROTEOMICS.md](CLAUDIUS-PROTEOMICS.md) for the full project plan.

## Quick Start

### Prerequisites

1. **FragPipe** (v24+) - Download from https://fragpipe.nesvilab.org/
2. **Singularity** - For containerized HPC execution
3. **Snakemake** (v8+) - Workflow orchestration

### Configuration

Edit `config/config.yaml`:
```yaml
fragpipe:
  path: "/path/to/your/fragpipe"  # Required

database:
  fasta: "/path/to/database.fasta"  # Required
```

### Build Containers

```bash
# Build on a machine with sudo access
singularity build containers/imspy.sif containers/imspy/Singularity.def
singularity build containers/fragpipe-base.sif containers/fragpipe-base/Singularity.def
singularity build containers/download.sif containers/download/Singularity.def
```

### Phase 1: Build Database

```bash
# Add a dataset to the database
snakemake --profile profiles/mogon2 add_dataset --config accession=PXD019086

# Add all configured datasets
snakemake --profile profiles/mogon2 add_all_datasets

# Create a snapshot for training
snakemake --profile profiles/mogon2 create_snapshot --config version=v1.0
```

### Phase 2: Train Models

```bash
# Train CCS model on snapshot
snakemake --profile profiles/mogon2 train_model --config snapshot=v1.0 model=ccs_v1

# Evaluate model
snakemake --profile profiles/mogon2 evaluate_model --config snapshot=v1.0 model=ccs_v1
```

### Local Development

```bash
# Create conda environment
conda env create -f envs/imspy.yaml
conda activate claudius-proteomics

# Install imspy from source
./scripts/setup_imspy.sh

# Run locally
snakemake --cores 8 add_dataset --config accession=PXD019086
```

## Project Structure

```
claudius-proteomics/
├── Snakefile                   # Main workflow
├── config/config.yaml          # Pipeline configuration
├── profiles/mogon2/            # HPC cluster profile
│
├── rules/                      # Snakemake rule modules
│   ├── download.smk            # PRIDE download
│   ├── fragpipe.smk            # PSM identification
│   ├── extract.smk             # Raw feature extraction (imspy)
│   ├── database.smk            # Database insertion
│   ├── snapshot.smk            # Snapshot creation
│   └── train.smk               # Model training
│
├── containers/                 # Singularity definitions
│   ├── fragpipe-base/          # Java runtime
│   ├── download/               # PRIDE tools
│   └── imspy/                  # imspy + TensorFlow
│
├── scripts/                    # Python/Bash scripts
│   ├── download_pride.py
│   ├── run_fragpipe.sh
│   ├── extract_raw_features.py # imspy raw extraction
│   ├── insert_to_database.py
│   ├── create_snapshot.py
│   ├── train_ccs.py
│   ├── evaluate_model.py
│   └── setup_imspy.sh          # Local imspy setup
│
├── database/                   # Peptide property database
│   ├── peptides/               # Per-accession data
│   └── snapshots/              # Versioned training sets
│
├── models/                     # Trained models
│
├── data/                       # Temporary processing data
│   ├── raw/                    # Downloaded .d files
│   ├── processed/              # FragPipe output
│   └── extracted/              # Raw features
│
├── envs/imspy.yaml             # Conda environment
├── CLAUDIUS-PROTEOMICS.md      # Full project plan
└── README.md
```

## Data Sources

### POC: Meier et al. 2021

- **Paper**: [Deep learning the collisional cross sections of the peptide universe](https://doi.org/10.1038/s41467-021-21352-8)
- **PRIDE**: PXD019086 (primary), PXD010012, PXD017703
- **Data**: >1 million CCS measurements from timsTOF
- **Target**: R > 0.99, <1.4% median error

## Technical Stack

- **Data processing**: rustims/imspy (Rust + Python)
- **PSM identification**: FragPipe
- **ML framework**: TensorFlow, imspy-predictors
- **Workflow**: Snakemake
- **Containers**: Singularity
- **HPC**: SLURM (Mogon2)

## License

[Add license]

## Citation

[Add citation when published]
