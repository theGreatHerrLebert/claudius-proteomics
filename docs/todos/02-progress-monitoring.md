# Progress Monitoring

**Owner:** TBD
**Priority:** High
**Status:** Not Started

## Overview

Build a comprehensive tracking system to monitor pipeline progress, peptide identification rates, and dataset quality metrics. Enable visibility into the San José database growth over time.

## Current State

- Basic manifest.json files per dataset with counts
- Database-level metadata.json with aggregate stats
- Snakemake logs in `logs/` directory

## TODOs

### 1. Peptide/PSM Statistics

- [ ] Total unique peptides identified (across all datasets)
- [ ] Total PSMs (with breakdown by dataset)
- [ ] Peptide ID rate: #peptides / #datasets mined
- [ ] New peptides per dataset (marginal contribution)
- [ ] Charge state distribution (overall and per dataset)
- [ ] Modification coverage (PTM types identified)
- [ ] Sequence length distribution

### 2. Dataset Progress Tracking

- [ ] Datasets queued / processing / completed / failed
- [ ] Raw files processed per dataset
- [ ] Processing time per dataset (download, search, extract, insert)
- [ ] Storage usage (raw data, processed, database)
- [ ] Failure reasons and retry counts

### 3. Quality Metrics

- [ ] Average PSMs per raw file
- [ ] FDR distribution across datasets
- [ ] CCS coverage (min/max/median per charge state)
- [ ] RT range coverage
- [ ] Outlier detection (unusual CCS values, etc.)

### 4. Search Engine Comparison

- [ ] FragPipe vs DIA-NN overlap rate per dataset
- [ ] Unique IDs per method
- [ ] Concordance metrics (same peptide, similar scores)
- [ ] Score correlation when found by both

### 5. Implementation Tasks

- [ ] Create `scripts/generate_stats.py` for metrics computation
- [ ] Add `stats` rule to Snakefile
- [ ] Store time-series metrics in `database/metrics.parquet`
- [ ] Create Streamlit dashboard for visualization (optional)
- [ ] Slack/email notifications for milestones (10k, 100k, 1M peptides)

## Proposed Metrics Schema

```python
# database/metrics.parquet
{
    "timestamp": datetime,
    "n_datasets": int,
    "n_datasets_completed": int,
    "n_unique_peptides": int,
    "n_total_psms": int,
    "n_peptides_fragpipe_only": int,
    "n_peptides_diann_only": int,
    "n_peptides_overlap": int,
    "storage_gb": float,
    "avg_psms_per_file": float,
}
```

## Files to Create/Modify

- `scripts/generate_stats.py` - NEW: Metrics computation
- `scripts/dashboard.py` - NEW: Optional Streamlit app
- `rules/stats.smk` - NEW: Statistics rules
- `database/metrics.parquet` - NEW: Time-series metrics

## Acceptance Criteria

- [ ] `snakemake stats` generates comprehensive report
- [ ] Metrics tracked over time in parquet file
- [ ] Clear visibility into pipeline health
- [ ] Automated alerts for failures or milestones
