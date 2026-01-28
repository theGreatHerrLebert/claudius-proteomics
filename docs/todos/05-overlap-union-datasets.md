# Overlap vs Union Dataset Design

**Owner:** TBD
**Priority:** High
**Status:** Not Started

## Overview

Design and implement a consensus strategy for combining peptide identifications from FragPipe and DIA-NN. Define UNION (maximum coverage), OVERLAP (high confidence), and WEIGHTED (balanced) dataset types.

## Problem Statement

FragPipe (spectrum-centric) and DIA-NN (peptide-centric) identify different but overlapping sets of peptides. We need a principled approach to:

1. Maximize training data (UNION)
2. Ensure high confidence labels (OVERLAP)
3. Balance coverage and confidence (WEIGHTED)

## Dataset Types

### 1. UNION Dataset (Maximum Coverage)

**Definition:** Include ALL peptides identified by EITHER search engine.

**Characteristics:**
- Highest peptide count
- Maximum sequence diversity
- May include more false positives
- Some peptides have scores from only one engine

**Use Cases:**
- Maximize training data coverage
- Explore rare peptides
- Initial model training

**Filter Criteria:**
```python
# Peptide included if:
(fragpipe_qvalue <= 0.01) OR (diann_qvalue <= 0.01)
```

### 2. OVERLAP Dataset (High Confidence)

**Definition:** Include ONLY peptides found by BOTH search engines.

**Characteristics:**
- Lower peptide count
- Higher confidence (orthogonal validation)
- Scores available from both engines
- Conservative approach

**Use Cases:**
- Benchmark model evaluation
- High-precision applications
- Publication-quality results

**Filter Criteria:**
```python
# Peptide included if:
(fragpipe_qvalue <= 0.01) AND (diann_qvalue <= 0.01)
AND (same_sequence) AND (same_charge)
AND (abs(rt_diff) <= 30)  # seconds
AND (abs(mobility_diff) <= 0.02)  # 1/K0
```

### 3. WEIGHTED Dataset (Balanced)

**Definition:** Include all peptides with confidence weights for training.

**Characteristics:**
- All peptides included
- Weight reflects confidence level
- Enables sample-weighted training
- Balances coverage and confidence

**Weight Calculation:**
```python
def calculate_weight(peptide):
    if found_by_both:
        # High confidence: found by both engines
        base_weight = 1.0
        score_bonus = (1 - fragpipe_pep) * (1 - diann_pep)
    elif found_by_fragpipe:
        base_weight = 0.7
        score_bonus = (1 - fragpipe_pep)
    else:  # found by diann
        base_weight = 0.7
        score_bonus = (1 - diann_pep)

    # Observation bonus
    obs_bonus = min(n_observations / 10, 0.2)

    return base_weight * score_bonus + obs_bonus
```

## Consensus Matching Criteria

### Defining "Same Peptide"

Two PSMs are considered the same peptide if:

| Criterion | Threshold | Rationale |
|-----------|-----------|-----------|
| Sequence | Exact match | Must be identical |
| Charge | Exact match | Different charge = different precursor |
| Mass error | < 10 ppm | Allows for calibration differences |
| RT difference | < 30 sec | Or < 2% of gradient length |
| Mobility difference | < 0.02 1/K0 | ~2% relative tolerance |

### Handling Modifications

- Unmodified sequences: Direct comparison
- Modified sequences: Normalize notation first
  - FragPipe: `M[147]` or `M(ox)`
  - DIA-NN: `M(UniMod:35)`
  - Normalize to: `M[+15.9949]` (delta mass)

### Handling Isobaric Peptides

**Problem:** Same precursor mass can match different peptide sequences. FragPipe and DIA-NN may assign different sequences to the same MS2 spectrum.

**Detection:**
```python
# Flag potential isobaric conflicts
conflicts = merged_df.groupby(['raw_file', 'scan']).filter(
    lambda x: x['sequence'].nunique() > 1
)
```

**Resolution Strategies:**
1. **Trust higher score:** Use peptide with better combined score
2. **Require agreement:** Exclude conflicting assignments from OVERLAP
3. **Flag for review:** Keep both with `isobaric_conflict=True` flag
4. **Mass difference:** If Δmass > 0.01 Da, treat as different precursors

**Recommendation:** For OVERLAP dataset, require sequence agreement. For UNION, include both with conflict flag.

## Implementation

### New Scripts

```python
# scripts/compute_consensus.py

def compute_consensus(fragpipe_df, diann_df, consensus_type='union'):
    """
    Merge FragPipe and DIA-NN results.

    Args:
        fragpipe_df: FragPipe PSMs
        diann_df: DIA-NN PSMs
        consensus_type: 'union', 'overlap', or 'weighted'

    Returns:
        Consensus DataFrame with source annotations
    """
    pass

def match_peptides(fp_peptides, dn_peptides, tolerances):
    """
    Find matching peptides between engines.
    """
    pass

def calculate_overlap_stats(matched_df):
    """
    Compute Jaccard similarity, unique counts, etc.
    """
    pass
```

### Snakemake Integration

```python
# rules/consensus.smk

rule compute_consensus:
    input:
        fragpipe="data/{accession}/fragpipe/psm.tsv",
        diann="data/{accession}/diann/report.parquet"
    output:
        "data/{accession}/consensus.parquet"
    params:
        consensus_type=config.get("consensus_type", "union")
    script:
        "../scripts/compute_consensus.py"

rule create_consensus_snapshot:
    input:
        expand("data/{acc}/consensus.parquet", acc=config["datasets"])
    output:
        directory("database/snapshots/{version}_{consensus_type}")
    params:
        consensus_type="{consensus_type}"
    script:
        "../scripts/create_snapshot.py"
```

## Statistics to Track

### Per-Dataset Metrics

```json
{
  "accession": "PXD019086",
  "n_fragpipe_only": 15000,
  "n_diann_only": 12000,
  "n_overlap": 45000,
  "n_union": 72000,
  "jaccard_similarity": 0.625,
  "overlap_rate": 0.625
}
```

### Aggregate Metrics

- Overall Jaccard similarity across all datasets
- Peptide-level vs PSM-level overlap
- Overlap by:
  - Sequence length
  - Charge state
  - Modification status
  - CCS range

## Analysis Questions

After implementation, analyze:

1. **Which peptide types differ?**
   - Length distribution of engine-specific peptides
   - Charge state bias
   - Modification prevalence

2. **Are differences systematic or random?**
   - Consistent bias across datasets?
   - Instrument/acquisition-dependent?

3. **Score correlation for overlap peptides**
   - Hyperscore vs CScore correlation
   - PEP agreement

4. **Impact on model training**
   - Train on UNION vs OVERLAP
   - Compare model performance

## Snapshot Directory Structure

```
database/snapshots/
├── v1.0_union/
│   ├── train.parquet
│   ├── val.parquet
│   ├── test.parquet
│   └── manifest.json
├── v1.0_overlap/
│   ├── train.parquet
│   ├── val.parquet
│   ├── test.parquet
│   └── manifest.json
├── v1.0_weighted/
│   ├── train.parquet      # includes 'weight' column
│   ├── val.parquet
│   ├── test.parquet
│   └── manifest.json
└── v1.0_fragpipe_only/    # single-engine baseline
    └── ...
```

## Worker Architecture

Each worker processes exactly **one PRIDE accession** (PID) atomically:

```
┌─────────────────────────────────────────────────────────┐
│  Worker Job: PXD019086                                  │
├─────────────────────────────────────────────────────────┤
│  1. Download raw files                                  │
│  2. Run FragPipe (all files in accession)              │
│  3. Run DIA-NN (all files in accession)                │
│  4. Compute consensus (per-accession)                  │
│  5. Extract scores + raw features                       │
│  6. Insert to San José database                        │
│  7. Report metrics                                      │
└─────────────────────────────────────────────────────────┘
```

**Key Design Decisions:**
- Consensus computed **per-accession** (not globally) - allows incremental database growth
- Global consensus table rebuilt periodically from per-accession results
- Worker is stateless - can be killed and restarted
- Output is idempotent - re-running same PID overwrites cleanly

**Worker Input:**
```yaml
# job.yaml
accession: PXD019086
organism: human           # from PRIDE metadata or override
fasta_source: uniprot     # or custom path
consensus_type: union     # default for initial processing
```

**Worker Output:**
```
database/peptides/accession=PXD019086/
├── peptides.parquet      # All peptide observations
├── consensus.parquet     # Per-accession consensus
├── scores.parquet        # Search engine scores
├── overlap_stats.json    # FragPipe vs DIA-NN metrics
└── manifest.json         # Job metadata + provenance
```

## Acceptance Criteria

- [ ] `compute_consensus.py` handles all three types
- [ ] Matching criteria configurable via config.yaml
- [ ] Overlap statistics computed and stored
- [ ] Venn diagram visualization generated
- [ ] Snapshots created for each consensus type
- [ ] Documentation of which type to use when
- [ ] Worker processes single PID end-to-end
- [ ] Output is idempotent (safe to re-run)
