# Dataset Definition

## Definition

**DATASET = PRIDE Accession**

One dataset corresponds to one PRIDE accession (e.g., `PXD019086`). This is the atomic unit for:
- Job distribution to workers
- Database partitioning
- Progress tracking
- Result aggregation

## Hierarchy

```
PRIDE Accession (DATASET)
└── Raw Files (N per dataset)
    └── PSMs (M per raw file)
        └── Peptides (deduplicated across PSMs)
```

## Assumptions

Each dataset is assumed to have:
- **Single organism** (or dominant organism)
- **Compatible instrument type** (timsTOF for now)
- **Consistent acquisition method** (DDA or DIA, not mixed)

## Handling Edge Cases

### Multi-organism Studies

If a PRIDE accession contains multiple organisms:
1. Add to `config/blacklist.yaml` with reason
2. Or: manually split into sub-jobs with explicit organism override

```yaml
# config/blacklist.yaml
blacklist:
  PXD012345:
    reason: "Mixed human/mouse samples"
    action: "manual_split"
```

### Multi-instrument Studies

If a PRIDE accession has incompatible instruments (e.g., timsTOF + Orbitrap):
1. Filter raw files by instrument type in job manifest
2. Or: blacklist and process separately

### Problematic Datasets

Maintain a blacklist for datasets that:
- Have mixed organisms
- Have incompatible file formats
- Are known to have data quality issues
- Require special handling

## Configuration

```yaml
# config/config.yaml
datasets:
  - PXD019086
  - PXD010012
  - PXD017703

# Dataset-specific overrides (optional)
dataset_metadata:
  PXD019086:
    organism: "human"
    instrument: "timsTOF Pro"
    notes: "Meier et al. 2021 CCS benchmark"
```

## Job Manifest

Each worker job processes exactly one dataset:

```yaml
job_id: "PXD019086_20240127"
accession: "PXD019086"      # The dataset
organism: "human"           # Auto-detected or overridden
raw_files: [...]            # All files in this accession
```

## Database Partitioning

Data is partitioned by dataset (accession):

```
database/peptides/
├── accession=PXD019086/
│   └── data.parquet
├── accession=PXD010012/
│   └── data.parquet
```

This enables:
- Easy addition/removal of datasets
- Parallel processing
- Dataset-level quality filtering
