# San José Database Design

**Owner:** TBD
**Priority:** Critical
**Status:** Not Started

## Overview

Define the complete schema for the San José peptide property database. Design tables, metadata structures, and partitioning strategy for scalability.

## Database Name

**San José** - The central peptide property database for the CLAUDIUS-PROTEOMICS project.

## Core Tables

### 1. `peptides/` - Main Peptide Observations

Primary table storing all peptide identifications with measured properties.

| Column | Type | Description |
|--------|------|-------------|
| sequence | string | Unmodified peptide sequence |
| modified_sequence | string | Sequence with modification annotations |
| charge | int8 | Precursor charge state |
| mass | float64 | Monoisotopic mass (Da) |
| mz | float64 | Precursor m/z |
| retention_time | float64 | Retention time (seconds) |
| mobility | float64 | Ion mobility (1/K0) |
| ccs | float64 | Collisional cross section (Å²) |
| accession | string | PRIDE dataset accession |
| raw_file | string | Source raw file name |
| search_engine | string | 'fragpipe' or 'diann' |

**Partitioning:** By `accession` (PRIDE ID)

### 2. `psm_scores/` - Search Engine Scores

Detailed scoring information for each PSM.

| Column | Type | Description |
|--------|------|-------------|
| psm_id | string | Unique ID (hash of seq+charge+file+scan) |
| sequence | string | Peptide sequence |
| charge | int8 | Charge state |
| hyperscore | float64 | MSFragger primary score |
| expectation | float64 | E-value |
| q_value | float64 | FDR-controlled q-value |
| pep | float64 | Posterior error probability |
| delta_score | float64 | Score difference to 2nd hit |
| matched_ions | int16 | Number of matched fragments |
| total_ions | int16 | Total theoretical fragments |
| accession | string | PRIDE accession |

**Partitioning:** By `accession`

### 3. `raw_features/` - Extracted Signal Features

Raw signal data extracted from instrument files.

| Column | Type | Description |
|--------|------|-------------|
| psm_id | string | Links to psm_scores |
| xic_rt | list[float64] | Chromatogram retention times |
| xic_intensity | list[float64] | Chromatogram intensities |
| mobilogram_mobility | list[float64] | Ion mobility values |
| mobilogram_intensity | list[float64] | Mobility intensities |
| isotope_mz | list[float64] | Isotope envelope m/z |
| isotope_intensity | list[float64] | Isotope intensities |
| ms2_mz | list[float64] | MS2 fragment m/z |
| ms2_intensity | list[float64] | MS2 fragment intensities |

**Partitioning:** By `accession`

### 4. `datasets/` - Dataset Metadata

Metadata for each processed PRIDE dataset.

| Column | Type | Description |
|--------|------|-------------|
| accession | string | PRIDE accession (PK) |
| title | string | Dataset title |
| description | string | Dataset description |
| organism | string | Source organism |
| instrument | string | Mass spectrometer model |
| acquisition_type | string | DDA, DIA, or diaPASEF |
| n_raw_files | int32 | Number of raw files |
| n_psms | int64 | Total PSMs identified |
| n_unique_peptides | int64 | Unique peptides |
| date_added | timestamp | When added to San José |
| processing_status | string | queued/processing/completed/failed |
| pride_url | string | Link to PRIDE entry |

**Storage:** Single parquet file (not partitioned)

### 5. `consensus/` - Deduplicated Peptide Properties

Aggregated peptide properties across all observations.

| Column | Type | Description |
|--------|------|-------------|
| sequence | string | Peptide sequence |
| charge | int8 | Charge state |
| median_ccs | float64 | Median CCS across observations |
| std_ccs | float64 | Standard deviation of CCS |
| median_rt | float64 | Median retention time |
| std_rt | float64 | RT standard deviation |
| n_observations | int32 | Number of PSMs |
| n_datasets | int16 | Number of source datasets |
| sources | list[string] | List of accessions |
| confidence | string | high/medium/low |

**Storage:** Single parquet file

### 6. `snapshots/` - Training Dataset Versions

Versioned, immutable exports for model training.

| Column | Type | Description |
|--------|------|-------------|
| version | string | Snapshot version (v1.0) |
| created_at | timestamp | Creation timestamp |
| n_peptides | int64 | Unique peptides |
| n_observations | int64 | Total observations |
| consensus_type | string | union/overlap/weighted |
| split_seed | int32 | Random seed for splits |
| train_fraction | float32 | Training set fraction |
| val_fraction | float32 | Validation set fraction |
| test_fraction | float32 | Test set fraction |
| filter_criteria | json | Applied filters |
| ccs_range | tuple[float, float] | CCS bounds |

**Storage:** Directory per version with parquet + manifest

## Directory Structure

```
database/
├── peptides/                    # Partitioned by accession
│   ├── accession=PXD019086/
│   │   └── data.parquet
│   ├── accession=PXD010012/
│   │   └── data.parquet
│   └── ...
├── psm_scores/                  # Partitioned by accession
│   └── accession=PXD019086/
│       └── data.parquet
├── raw_features/                # Partitioned by accession
│   └── accession=PXD019086/
│       └── data.parquet
├── datasets.parquet             # Dataset metadata
├── consensus.parquet            # Deduplicated properties
├── snapshots/                   # Training snapshots
│   ├── v1.0_union/
│   │   ├── train.parquet
│   │   ├── val.parquet
│   │   ├── test.parquet
│   │   └── manifest.json
│   └── v1.0_overlap/
│       └── ...
├── metrics.parquet              # Time-series metrics
├── schema.json                  # Column definitions
└── metadata.json                # Database-level stats
```

## Metadata Files

### `schema.json`
```json
{
  "version": "1.0",
  "tables": {
    "peptides": {
      "columns": {
        "sequence": {"type": "string", "nullable": false, "description": "..."},
        ...
      },
      "partition_by": "accession"
    }
  }
}
```

### `metadata.json`
```json
{
  "database_name": "San José",
  "version": "1.0",
  "created_at": "2024-01-27T...",
  "last_updated": "2024-01-27T...",
  "n_datasets": 15,
  "n_unique_peptides": 250000,
  "n_total_psms": 1500000,
  "schema_version": "1.0"
}
```

## Implementation Tasks

- [ ] Create `database/schema.json` with full definitions
- [ ] Update `scripts/insert_to_database.py` for new schema
- [ ] Create `scripts/build_consensus.py` for aggregation
- [ ] Update Snakemake rules for new table structure
- [ ] Add data validation on insert
- [ ] Update documentation with "San José" branding

## Acceptance Criteria

- [ ] All 6 tables implemented and documented
- [ ] Schema validation on every insert
- [ ] Efficient queries via partitioning
- [ ] Manifest files for provenance tracking
- [ ] Migration path from current structure
