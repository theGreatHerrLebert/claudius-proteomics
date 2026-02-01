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
2. **DIA-NN** (v2.3+) - Download from https://github.com/vdemichev/DiaNN (required for DDA support)
3. **Snakemake** (v8+) - Workflow orchestration
4. **Singularity** (optional) - For containerized HPC execution

### Configuration

Edit `config/config.yaml`:
```yaml
# FragPipe installation
fragpipe:
  path: "/path/to/fragpipe-24.0"
  workflow: "LFQ-MBR"  # Good for DDA timsTOF

# DIA-NN installation (v2.3+ for DDA support)
diann:
  path: "/path/to/diann-linux"
  threads: 16

# Organism FASTA mapping
organisms:
  human:
    proteome_id: "UP000005640"
    local_fasta: "/path/to/human.fasta"  # Optional
    includes_contaminants: true  # If local FASTA already has cRAP

# Dataset to organism mapping
dataset_metadata:
  PXD019086:
    organism: "human"
```

### Search Methods

The pipeline supports two orthogonal search methods for robust peptide identification:

| Method | Tool | Approach | Strengths |
|--------|------|----------|-----------|
| **Spectrum-centric** | FragPipe/MSFragger | Match spectra to peptides | Fast, established, open modification support |
| **Peptide-centric** | DIA-NN | Match peptides to spectra | Deep learning rescoring, excellent for ion mobility |

Running both provides orthogonal validation - consensus identifications have higher confidence.

```bash
# Spectrum-centric only (FragPipe)
snakemake process_only --cores 16

# Peptide-centric only (DIA-NN)
snakemake process_diann --cores 16

# Both methods (recommended for production)
snakemake process_both --cores 16
```

### FASTA Resource Management

The pipeline automatically manages FASTA databases:

1. **Organism mapping**: Each dataset maps to an organism (human, yeast, mouse, etc.)
2. **Local or download**: Uses local FASTA if available, otherwise downloads from UniProt
3. **Contaminants**: Adds cRAP database (unless already included)
4. **Decoys**: Generates reversed decoy sequences for target-decoy analysis

```bash
# Prepare FASTA for a dataset
snakemake prepare_fasta --config accession=PXD019086

# List available databases
snakemake list_fasta_databases
```

### Phase 1: Build Database

```bash
# Add a dataset to the database
snakemake add_dataset --config accession=PXD019086

# Add all configured datasets
snakemake add_all_datasets

# Create a snapshot for training
snakemake create_snapshot --config version=v1.0
```

### Phase 2: Train Models

```bash
# Train CCS model on snapshot
snakemake train_model --config snapshot=v1.0 model=ccs_v1

# Evaluate model
snakemake evaluate_model --config snapshot=v1.0 model=ccs_v1
```

### Local Development

```bash
# Create Python venv
python3.12 -m venv .venv
source .venv/bin/activate
pip install snakemake

# Install imspy from source
./scripts/setup_imspy.sh

# Run locally (test mode: 3 files)
snakemake process_both --cores 8
```

## Project Structure

```
claudius-proteomics/
├── Snakefile                   # Main workflow
├── config/config.yaml          # Pipeline configuration
├── profiles/mogon2/            # HPC cluster profile
│
├── rules/                      # Snakemake rule modules
│   ├── fasta.smk               # FASTA resource management
│   ├── download.smk            # PRIDE download
│   ├── fragpipe.smk            # MSFragger search (spectrum-centric)
│   ├── diann.smk               # DIA-NN search (peptide-centric)
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
│   ├── download_pride.py       # PRIDE download
│   ├── run_fragpipe.py         # FragPipe wrapper
│   ├── extract_raw_features.py # imspy raw extraction
│   ├── insert_to_database.py
│   ├── create_snapshot.py
│   ├── train_ccs.py
│   ├── evaluate_model.py
│   └── setup_imspy.sh          # Local imspy setup
│
├── resources/                  # Resource files
│   └── fasta/                  # FASTA databases
│       ├── {organism}.fasta    # Per-organism databases
│       └── search_db/          # Per-dataset search databases
│           ├── {accession}.fasta
│           └── {accession}_decoys.fasta
│
├── database/                   # Peptide property database
│   ├── peptides/               # Per-accession data
│   └── snapshots/              # Versioned training sets
│
├── models/                     # Trained models
│
├── data/                       # Temporary processing data
│   ├── raw/                    # Downloaded .d files
│   ├── processed/              # Search results
│   │   └── {accession}/
│   │       ├── psm.tsv         # FragPipe output
│   │       └── diann/
│   │           └── report.parquet  # DIA-NN output
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

- **Raw data processing**: rustims/imspy (Rust + Python)
- **Search engines**:
  - FragPipe/MSFragger (spectrum-centric)
  - DIA-NN v2.3+ (peptide-centric, DDA mode)
- **ML framework**: PyTorch, imspy-predictors
- **Workflow**: Snakemake
- **Containers**: Singularity
- **HPC**: SLURM (Mogon2)

## License

[Add license]

## Citation

[Add citation when published]
