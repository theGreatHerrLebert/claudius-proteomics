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
└── features.parquet
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
| `raw_file` | string | Source .d file name |
| `scan_id` | int64 | Reference scan |
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

## Join Strategy

**Within a runner output:**
- FragPipe ↔ Features: `(raw_file, scan_id)`
- DIA-NN ↔ Features: `(raw_file, precursor_mz, charge, retention_time)` with tolerance

**Across runners (aggregator):**
- `(modified_sequence, charge)` → unique precursor identity
- Allows merging same peptide observed in different datasets
