# San José Runner Architecture

## Overview

The San José Runner is a self-contained 6-step pipeline designed for SLURM cluster submission. It searches PRIDE datasets with independent engines, extracts raw signal features, and produces dashboard-ready datasets.

> **Engines actually run.** The production path is **two** engines, **Sage +
> FragPipe**. A DIA-NN harness exists (`runner/engines/diann_job.py`) and the
> schema below carries `diann_*` columns and `diann_only` strata, but DIA-NN is
> **off by default** — measured as a third vote it added +8.5% cross-engine
> corroboration and 0.7% unique precursors at q ≤ 0.01, which did not justify the
> per-dataset cost. Everything described below is engine-agnostic; with the
> default configuration the `diann_*` fields are NULL and `n_engines` is 0–2.

```
RUNNER (Steps 1-5) - Submittable to SLURM
┌──────────────────────────────────────────────────────────────────┐
│  Step 1      Step 2       Step 3          Step 4       Step 5   │
│  Download → Search → Stratify/Merge → Extract → Final Merge     │
└──────────────────────────────────────────────────────────────────┘
                                                          ↓
                                                    Step 6: Dashboard
```

---

## Directory Structure

```
runner/
├── __init__.py
├── run_dataset.py           # Main entry point
├── slurm_wrapper.sh         # SLURM submission script
├── steps/
│   ├── __init__.py
│   ├── step1_download.py    # Download from PRIDE or link local data
│   ├── step2_search.py      # Run FragPipe, DIA-NN, Sage
│   ├── step3_stratify.py    # UNIMOD standardization, consensus
│   ├── step4_extract.py     # 4D raw signal extraction
│   └── step5_merge.py       # Final merge for dashboard
├── state.py                 # Checkpointing for resumability
└── summary.py               # Summary generation utilities

lib/                         # Shared utilities
├── __init__.py
├── sequence_utils.py        # UNIMOD standardization
├── precursor_matching.py    # Tolerance matching logic
├── quality_metrics.py       # Gaussian fits, isotope cosine
└── parquet_utils.py         # Parquet read/write helpers
```

---

## Step-by-Step Specification

### Step 1: Download Raw Data

**Purpose**: Fetch timsTOF .d files from PRIDE or symlink local data

**Input**: PRIDE accession (e.g., `PXD019086`)

**Outputs**:
```
data/raw/{accession}/*.d/
data/metadata/{accession}/metadata.yaml
data/processed/{accession}/step1_summary.json
```

**Summary file** (`step1_summary.json`):
- `n_raw_files`, `total_size_gb`, `raw_files[]`
- `metadata`: organism, instrument, lab_id, completeness
- `status`: success/error

---

### Step 2: Execute Third-Party Software

**Purpose**: Run the configured engines (default Sage + FragPipe; DIA-NN optional) — all with FDR=1.0 so nothing is filtered before the store

**Input**: `data/raw/{accession}/*.d`, FASTA database

**Outputs**:
```
data/processed/{accession}/
├── fragpipe/combined_ion.tsv
├── diann/report.parquet
├── sage/results.sage.parquet
└── step2_summary.json
```

**Summary file** (`step2_summary.json`):
- Per-engine: version, n_psms, n_peptides, duration
- FASTA: path, n_proteins

---

### Step 3: Stratify and Merge Third-Party Output

**Purpose**: UNIMOD-standardize sequences, compute consensus, stratify by engine agreement

**Input**: FragPipe/DIA-NN/Sage outputs from Step 2

**Outputs**:
```
data/processed/{accession}/
├── precursor_index.parquet    # Unified index: precursor_id → all engine IDs
├── consensus/
│   ├── overlap_stats.json     # 3-way Venn statistics
│   ├── overlap_report.html    # Visual for human review
│   └── stratified/
│       ├── all_three.parquet
│       ├── two_plus.parquet
│       ├── fragpipe_only.parquet
│       ├── diann_only.parquet
│       └── sage_only.parquet
└── step3_summary.json
```

**Summary file** (`step3_summary.json`):
- `n_all_three`, `n_two_plus`, engine-specific counts
- `validation_rates`: per-engine % confirmed by another engine
- `match_tiers`: distribution of match quality

**Implementation Notes**:
- Uses efficient DataFrame merges instead of row-by-row iteration
- Processes ~1M PSMs in ~13 seconds
- Groups by (sequence_normalized, charge) for unique precursors

---

### Step 4: Raw Data Extraction with Quality Metrics

**Purpose**: Extract 4D raw signals with Gaussian fits and quality metrics

**Input**: `data/raw/{accession}/*.d`

**Outputs**:
```
data/extracted/{accession}/
├── {raw_file}/
│   ├── index.parquet          # Per-file precursor metadata
│   └── blobs.bin              # Serialized fragment/MS1 data
├── raw_features.parquet       # Merged across all files
└── step4_summary.json
```

**Quality metrics computed**:
- `ms1_rt_sigma`, `ms1_rt_r2` - RT Gaussian fit
- `ms1_im_sigma`, `ms1_im_r2` - IM Gaussian fit
- `isotope_cosim` - Cosine similarity vs theoretical envelope

**Summary file** (`step4_summary.json`):
- `n_precursors_extracted`
- Quality metric distributions: mean, median, q10, q90
- `blob_size_gb`, `duration_minutes`

---

### Step 5: Merge Search Results with Extracted Raw Data

**Purpose**: Join engine IDs with raw features into final dashboard-ready dataset

**Input**:
- `data/processed/{accession}/precursor_index.parquet` (Step 3)
- `data/extracted/{accession}/raw_features.parquet` (Step 4)

**Outputs**:
```
data/merged/{accession}/
├── precursor_store.parquet    # Final unified dataset
├── manifest.json              # Runner completion manifest
└── step5_summary.json
```

**Final schema** (dashboard-ready):
- Identifiers: `precursor_id`, `raw_file`
- Properties: `mz`, `charge`, `rt_seconds`, `mobility`
- Engine IDs (UNIMOD): `fragpipe_modified`, `diann_modified`, `sage_modified`
- Consensus: `n_engines`, `consensus_peptide`, `confidence_weight`
- Quality: `ms1_rt_r2`, `ms1_im_r2`, `isotope_cosim`

**Manifest file** (`manifest.json`):
```json
{
  "accession": "PXD019086",
  "pipeline_version": "1.0",
  "generated_at": "2026-02-03T15:45:23",
  "n_total_precursors": 1249006,
  "n_per_engine": {
    "fragpipe": 51124,
    "diann": 81315,
    "sage": 543109
  },
  "quality_summary": {
    "n_all_three": 28927,
    "pct_high_quality": 1.3
  }
}
```

---

### Step 6: Dashboard Visualization

**Purpose**: Serve and visualize merged data

**Input**: `data/merged/{accession}/precursor_store.parquet`

**Components**:
- FastAPI backend (`dashboard/backend/main.py`)
- React frontend (`dashboard/frontend/`)

---

## Usage

### Local Execution

```bash
# Full run
python runner/run_dataset.py PXD019086

# Test mode (limited files)
python runner/run_dataset.py PXD019086 --test-mode --max-files 3

# Resume from checkpoint
python runner/run_dataset.py PXD019086 --resume

# Run specific steps
python runner/run_dataset.py PXD019086 --steps 1 2 3

# Use local data (skip download)
python runner/run_dataset.py PXD019086 --local-data /path/to/data
```

### SLURM Submission

```bash
# Basic submission
sbatch runner/slurm_wrapper.sh PXD019086

# Test mode on cluster
sbatch runner/slurm_wrapper.sh PXD019086 --test-mode

# Custom resources
sbatch --cpus-per-task=32 --mem=128G runner/slurm_wrapper.sh PXD019086

# Resume failed job
sbatch runner/slurm_wrapper.sh PXD019086 --resume
```

### Command-Line Options

| Option | Description |
|--------|-------------|
| `--config, -c` | Config file (default: config/config.yaml) |
| `--output-dir, -o` | Output base directory (default: data) |
| `--test-mode, -t` | Process limited files for testing |
| `--max-files` | Maximum files in test mode (default: 3) |
| `--resume, -r` | Resume from checkpoint if exists |
| `--steps` | Specific steps to run (1-5) |
| `--threads` | Number of threads (default: 16) |
| `--local-data` | Path to local data (skip download) |

---

## Checkpointing

The runner uses file-based checkpointing for resumability:

```
checkpoints/{accession}/
├── state.json         # Full runner state
├── step1.done         # Step completion markers
├── step2.done
├── step3.done
├── step4.done
└── step5.done
```

**State file** (`state.json`):
```json
{
  "accession": "PXD019086",
  "current_step": "step3_stratify",
  "started_at": "2026-02-03T15:30:00",
  "test_mode": true,
  "max_files": 3,
  "steps": {
    "step1_download": {"status": "completed", "duration_seconds": 0.5},
    "step2_search": {"status": "completed", "duration_seconds": 1.2},
    "step3_stratify": {"status": "running", "started_at": "2026-02-03T15:30:02"}
  }
}
```

---

## Shared Libraries (lib/)

### sequence_utils.py

UNIMOD standardization across search engines:

```python
from lib.sequence_utils import (
    standardize_fragpipe_modified_peptide,
    standardize_diann_sequence,
    standardize_sage_sequence,
    normalize_sequence_il,
)

# FragPipe: n[42.0106]PEPTIDEK → [UNIMOD:1]-PEPTIDEK-[]
# DIA-NN: (UniMod:1)PEPTIDEK → [UNIMOD:1]-PEPTIDEK-[]
# Sage: [+42.011]PEPTIDEK → [UNIMOD:1]-PEPTIDEK-[]
```

### precursor_matching.py

Tiered matching strategy:

```python
from lib.precursor_matching import PrecursorMatcher, MatchConfig, MatchTier

config = MatchConfig(mz_tol_ppm=20, rt_tol_sec=30, im_tol=0.05)
matcher = PrecursorMatcher(reference_df, config)

tier = matcher.match(query_row)
# MatchTier.SEQUENCE_EXACT, SEQUENCE_IL_NORM, COORDINATE_FULL, etc.
```

### quality_metrics.py

Gaussian fitting and isotope similarity:

```python
from lib.quality_metrics import fit_gaussian, compute_isotope_cosine_similarity

# Fit Gaussian to RT/IM profiles
mu, sigma, r2 = fit_gaussian(x_values, intensities)

# Compare observed vs theoretical isotope envelope
cosim = compute_isotope_cosine_similarity(observed, theoretical)
```

---

## Performance

Tested on PXD019086 (24 .d files, ~1M PSMs):

| Step | Duration | Notes |
|------|----------|-------|
| 1. Download | 0.0s | Used local data |
| 2. Search | 0.4s | Used existing outputs |
| 3. Stratify | 13.1s | 585,982 unique precursors |
| 4. Extract | 0.9s | Fallback mode |
| 5. Merge | 34.1s | 1.25M final rows |
| **Total** | **~48s** | |

**Optimization**: Step 3 uses DataFrame merges instead of row iteration, reducing time from 12+ minutes to ~13 seconds for 1M PSMs.

---

## Output Statistics Example

> ⚠️ **Historical.** These numbers come from an early three-engine POC run on
> PXD019086 and are kept to show the shape of the output, not as representative
> figures. They predate the current search configuration, and the production path
> no longer runs DIA-NN, so a "3 engines" row does not occur. For current corpus
> figures see the [README](../README.md).

From that test run on PXD019086:

**Engine Agreement**:
- FragPipe: 51,124 unique precursors
- DIA-NN: 81,315 unique precursors
- Sage: 543,109 unique precursors
- Union: 585,982 unique precursors
- All 3 engines: 28,927 (4.9%)
- At least 2: 60,639 (10.3%)

**Quality Metrics**:
- RT R² median: 0.59
- IM R² median: 0.69
- Isotope cosim median: 0.69
- High quality (R² ≥ 0.8, cosim ≥ 0.9): 1.3%

---

## Error Handling

The runner handles errors gracefully:

1. **Step failure**: State is saved, can resume with `--resume`
2. **Missing dependencies**: Falls back to alternative methods (e.g., imspy fallback)
3. **Missing data**: Clear error messages with file paths
4. **SLURM timeout**: Job can be resubmitted with `--resume`

---

## Future Enhancements

- [ ] Add unit tests for each step
- [ ] Parallel file processing in Step 4
- [ ] GPU support for extraction
- [ ] Snakemake integration for dependency management
- [ ] Cloud storage support (S3, GCS)
