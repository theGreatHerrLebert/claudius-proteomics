# Runner Output Schema

Runner processes a single PRIDE accession and outputs standardized results.

## Architecture

```
Map phase:  N runners → N datasets (parallel, independent)
Reduce phase: Aggregator → unified database + snapshots
```

## Runner Input

```yaml
accession: "PXD019086"
software_config:
  fragpipe_workflow: "LFQ-MBR"
  fragpipe_threads: 16
  diann_qvalue: 0.01
  extraction_tolerances:
    mz_ppm: 20
    rt_sec: 30
    mobility: 0.05
fasta:
  path: "/path/to/human.fasta"   # Direct path
  # OR
  organism: "human"               # Downloads from UniProt
```

## Runner Output

```
{accession}/
├── manifest.json
├── fragpipe.parquet
├── diann.parquet
├── features.parquet
├── consensus/
│   ├── overlap_stats.json    # FragPipe vs DIA-NN comparison metrics
│   ├── overlap.parquet       # High-confidence: found by BOTH engines
│   └── union.parquet         # Maximum coverage: found by EITHER engine
└── spectra/                  # Optional: raw signal blobs
    ├── index.parquet
    └── blobs.bin
```

---

## Modified Sequence Format (UniMod)

Canonical format for all outputs:

```
[N-term]-SEQUENCE[UNIMOD:ID]-[C-term]
```

**Examples:**
```
[]-PEPTIDEK-[]                      # Unmodified
[]-AAC[UNIMOD:4]PEPTIDEK-[]         # Carbamidomethyl (C)
[]-M[UNIMOD:35]PEPTIDEK-[]          # Oxidation (M)
[UNIMOD:1]-PEPTIDEK-[]              # N-term Acetyl
[]-AAC[UNIMOD:4]M[UNIMOD:35]PK-[]   # Multiple mods
```

**Common UNIMOD IDs:**
| ID | Name | Residue |
|----|------|---------|
| 1 | Acetyl | N-term |
| 4 | Carbamidomethyl | C |
| 35 | Oxidation | M |
| 21 | Phospho | S, T, Y |
| 7 | Deamidated | N, Q |

**Join key:** `(modified_sequence, charge)` → unique precursor identity

---

## manifest.json

```json
{
  "accession": "PXD019086",
  "runner_version": "1.0.0",
  "started_at": "2024-01-15T09:00:00Z",
  "completed_at": "2024-01-15T10:30:00Z",

  "input": {
    "n_raw_files": 48,
    "raw_files": ["file1.d", "file2.d"],
    "fasta_source": "human",
    "fasta_sha256": "abc123..."
  },

  "software": {
    "fragpipe": {"version": "24.0", "workflow": "LFQ-MBR"},
    "diann": {"version": "2.3.2", "qvalue": 0.01},
    "imspy": {"version": "0.4.0"}
  },

  "counts": {
    "fragpipe_psms": 1250000,
    "diann_precursors": 1180000,
    "features_extracted": 1100000
  },

  "qc": {
    "fragpipe_median_mass_error_ppm": 2.1,
    "diann_median_mass_error_ppm": 1.8
  }
}
```

---

## fragpipe.parquet

PSM-level results from FragPipe/MSFragger (spectrum-centric search).

| Column | Type | Description |
|--------|------|-------------|
| `precursor_id` | int64 | Unique ID within dataset (for joins) |
| `raw_file` | string | Source .d file name |
| `scan_id` | int64 | Spectrum scan number |
| `peptide_sequence` | string | Stripped sequence (no mods) |
| `modified_sequence` | string | UniMod format: `[]-AAC[UNIMOD:4]PY-[]` |
| `charge` | int8 | Precursor charge |
| `precursor_mz` | float64 | Observed m/z |
| `retention_time` | float32 | RT in seconds |
| `ion_mobility` | float32 | 1/K0 (timsTOF) |
| `hyperscore` | float32 | MSFragger score |
| `expectation` | float64 | E-value |
| `pep` | float64 | Posterior error probability |
| `protein_ids` | string | Semicolon-separated |
| `is_decoy` | bool | Target/decoy flag |

---

## diann.parquet

Precursor-level results from DIA-NN (peptide-centric search).

| Column | Type | Description |
|--------|------|-------------|
| `precursor_id` | int64 | Unique ID within dataset (for joins) |
| `raw_file` | string | Source .d file name |
| `peptide_sequence` | string | Stripped sequence (no mods) |
| `modified_sequence` | string | UniMod format: `[]-AAC[UNIMOD:4]PY-[]` |
| `charge` | int8 | Precursor charge |
| `precursor_mz` | float64 | Observed m/z |
| `retention_time` | float32 | RT in seconds |
| `ion_mobility` | float32 | 1/K0 |
| `ccs` | float32 | Calculated CCS (Å²) |
| `qvalue` | float64 | FDR q-value |
| `pep` | float64 | Posterior error probability |
| `ms1_intensity` | float64 | Precursor intensity |
| `protein_ids` | string | Semicolon-separated |

---

## features.parquet

Raw signal features extracted via imspy from timsTOF .d files.

| Column | Type | Description |
|--------|------|-------------|
| `precursor_id` | int64 | Unique ID within dataset (FK to fragpipe/diann) |
| `raw_file` | string | Source .d file name |
| `scan_id` | int64 | Reference scan (MS2 event) |
| `precursor_mz` | float64 | Target m/z |
| `charge` | int8 | Target charge |
| `rt_apex` | float32 | Chromatographic apex (s) |
| `rt_fwhm` | float32 | Peak width (s) |
| `rt_skew` | float32 | Peak asymmetry |
| `mobility_apex` | float32 | 1/K0 apex |
| `mobility_fwhm` | float32 | Mobilogram width |
| `ccs_measured` | float32 | CCS from mobility (Å²) |
| `isotope_mz` | list[float32] | Isotope m/z values |
| `isotope_intensity` | list[float32] | Isotope intensities |
| `total_intensity` | float64 | Summed signal |

---

## overlap_stats.json

Consensus metrics comparing FragPipe and DIA-NN identifications.

```json
{
  "accession": "PXD019086",
  "timestamp": "2024-01-15T10:30:00Z",
  "match_type": "sequence+charge",
  "stats": {
    "n_fragpipe": 51064,
    "n_diann": 33758,
    "n_overlap": 29780,
    "n_union": 55042,
    "n_fragpipe_only": 21284,
    "n_diann_only": 3978,
    "jaccard_similarity": 0.541,
    "fragpipe_overlap_rate": 0.583,
    "diann_overlap_rate": 0.882
  },
  "analysis": {
    "fragpipe_seq_length": {"mean": 15.7, "median": 15, "min": 7, "max": 50},
    "diann_seq_length": {"mean": 14.1, "median": 13, "min": 7, "max": 30},
    "fragpipe_charges": {"2": 44735, "3": 12650, "4": 2295},
    "diann_charges": {"2": 33811, "3": 7048, "4": 838}
  }
}
```

---

## consensus/overlap.parquet

High-confidence peptides found by BOTH search engines.

| Column | Type | Description |
|--------|------|-------------|
| `peptide_sequence` | string | Stripped sequence |
| `modified_sequence` | string | UniMod format |
| `charge` | int8 | Precursor charge |
| `precursor_mz` | float64 | Average m/z |
| `retention_time` | float32 | Average RT (seconds) |
| `ion_mobility` | float32 | Average 1/K0 |
| `ccs` | float32 | CCS (Å²) |
| `fragpipe_score` | float32 | Best hyperscore |
| `fragpipe_pep` | float64 | FragPipe PEP |
| `diann_qvalue` | float64 | DIA-NN q-value |
| `diann_pep` | float64 | DIA-NN PEP |
| `n_observations` | int16 | Total PSMs across both engines |
| `protein_ids` | string | Semicolon-separated |

---

## consensus/union.parquet

Maximum coverage: all peptides found by EITHER search engine.

| Column | Type | Description |
|--------|------|-------------|
| `peptide_sequence` | string | Stripped sequence |
| `modified_sequence` | string | UniMod format |
| `charge` | int8 | Precursor charge |
| `precursor_mz` | float64 | Observed m/z |
| `retention_time` | float32 | RT (seconds) |
| `ion_mobility` | float32 | 1/K0 |
| `ccs` | float32 | CCS (Å²) |
| `source` | string | 'fragpipe', 'diann', or 'both' |
| `best_score` | float32 | Best score from available engine |
| `best_pep` | float64 | Best PEP from available engine |
| `confidence_weight` | float32 | Weight for training (1.0 if both, 0.7 if single) |
| `protein_ids` | string | Semicolon-separated |

---

## Join Strategy

**Within a runner output:**
- FragPipe ↔ Features: `precursor_id` (assigned during extraction)
- DIA-NN ↔ Features: `precursor_id` (assigned during extraction)

**Note:** Avoid tolerance-based joins on (mz, rt, charge) - they're fragile. Instead, assign `precursor_id` during the extraction phase by matching each extracted feature to the nearest MS2 scan event.

**Across runners (aggregator):**
- `(modified_sequence, charge)` → unique precursor identity
- Allows merging same peptide observed in different datasets
