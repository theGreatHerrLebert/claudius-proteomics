# Third-Party ID Fields Extraction

**Owner:** TBD
**Priority:** High
**Status:** Not Started

## Overview

Extract comprehensive scoring and quality metrics from FragPipe and DIA-NN search results. These fields are critical for filtering, confidence assessment, and downstream analysis.

## Current State

- Basic PSM fields extracted (sequence, charge, m/z, RT, mobility)
- Limited score extraction
- No systematic documentation of available fields

## Fields to Extract

### 1. From FragPipe/MSFragger (psm.tsv)

| Field | Type | Description | Use Case |
|-------|------|-------------|----------|
| `Hyperscore` | float | MSFragger primary score | Ranking, filtering |
| `Nextscore` | float | Second-best match score | Delta score calculation |
| `Expectation` | float | E-value (lower = better) | Statistical significance |
| `PeptideProphet Probability` | float | Probability of correct ID | Confidence filtering |
| `Protein` | string | Primary protein accession | Protein inference |
| `Mapped Proteins` | string | All matching proteins | Proteotypic check |
| `Is Unique` | bool | Proteotypic peptide flag | Protein quantification |
| `Purity` | float | Precursor isolation purity | Quality filtering |
| `Assigned Modifications` | string | PTM annotations | Modification analysis |
| `Delta Mass` | float | Mass error (Da) | Calibration QC |
| `Spectrum` | string | Scan reference | Raw data linking |
| `Matched Peaks` | int | # matched fragment ions | Spectral quality |
| `Longest Sequence Match` | int | Consecutive matched ions | Spectral quality |

### 2. From FragPipe ion.tsv (Ion-Level)

| Field | Type | Description | Use Case |
|-------|------|-------------|----------|
| `Ion` | string | Fragment annotation (b3+, y7++) | MS2 feature extraction |
| `Charge` | int | Fragment charge state | Ion series analysis |
| `MZ` | float | Observed fragment m/z | Mass accuracy |
| `Intensity` | float | Fragment intensity | Intensity prediction |
| `Correlation` | float | XIC correlation | Peak quality |

### 3. From DIA-NN (report.parquet)

| Field | Type | Description | Use Case |
|-------|------|-------------|----------|
| `Q.Value` | float | Global Q-value (FDR) | Primary filtering |
| `PEP` | float | Posterior Error Probability | Per-PSM confidence |
| `Global.Q.Value` | float | Run-level FDR | Cross-run consistency |
| `Lib.Q.Value` | float | Library-level FDR | Library matching |
| `CScore` | float | Confidence score | Ranking |
| `Decoy` | bool | Target/decoy flag | FDR calculation |
| `Ms1.Area` | float | MS1 quantification | Abundance |
| `Proteotypic` | int | Unique to protein (0/1) | Protein inference |
| `Precursor.Quantity` | float | LFQ value | Quantification |
| `RT.Start` | float | Peak start RT | Peak boundaries |
| `RT.Stop` | float | Peak end RT | Peak boundaries |
| `IM` | float | Ion mobility (1/K0) | CCS calculation |
| `iIM` | float | Indexed ion mobility | Internal reference |
| `Predicted.RT` | float | DIA-NN RT prediction | Prediction comparison |
| `Predicted.IM` | float | DIA-NN mobility prediction | Prediction comparison |

### 4. Computed/Derived Fields

| Field | Formula | Description |
|-------|---------|-------------|
| `delta_score` | Hyperscore - Nextscore | Score discrimination |
| `ppm_error` | (obs - theo) / theo × 1e6 | Mass accuracy |
| `spectral_angle` | cosine(obs, library) | Library similarity |
| `log10_intensity` | log10(precursor_intensity) | Normalized intensity |
| `n_missed_cleavages` | count(K/R not at C-term) | Digestion efficiency |

## Implementation Tasks

### Phase 1: FragPipe Extraction
- [ ] Parse psm.tsv for all score columns
- [ ] Parse ion.tsv for fragment-level data
- [ ] Handle missing values appropriately
- [ ] Validate against expected ranges

### Phase 2: DIA-NN Extraction
- [ ] Read report.parquet native format
- [ ] Map column names to schema
- [ ] Extract prediction columns for comparison
- [ ] Handle v2.3+ format changes

### Phase 3: Derived Fields
- [ ] Implement delta_score calculation
- [ ] Implement ppm_error calculation
- [ ] Add missed cleavage counting
- [ ] Normalize scores to 0-1 scale (optional)

### Phase 4: Storage & Documentation
- [ ] Store in `psm_scores/` table
- [ ] Update `database/schema.json`
- [ ] Document field meanings and units
- [ ] Add data validation checks

## Files to Create/Modify

- `scripts/extract_fragpipe_scores.py` - NEW: FragPipe parser
- `scripts/extract_diann_scores.py` - NEW: DIA-NN parser
- `scripts/insert_to_database.py` - Update for new fields
- `database/schema.json` - Add field definitions
- `docs/field_definitions.md` - NEW: Field documentation

## Score Normalization (Optional)

For cross-engine comparisons, normalize scores to [0, 1]:

```python
# Hyperscore: higher is better, log-transform
normalized_hyperscore = (log(score) - min) / (max - min)

# E-value: lower is better, negative log-transform
normalized_evalue = -log10(evalue) / max_neglog

# Q-value: already in [0, 1], but invert
normalized_qvalue = 1 - qvalue
```

## Acceptance Criteria

- [ ] All listed fields extracted from both engines
- [ ] Fields documented in schema.json with types and descriptions
- [ ] Derived fields computed correctly
- [ ] Unit tests for parsing logic
- [ ] No data loss during extraction
