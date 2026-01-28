# San José: A Large-Scale Peptide Property Database for Ion Mobility Mass Spectrometry

## The Problem

Deep learning models for predicting peptide properties (CCS, retention time, MS2 spectra) are limited by training data. Current public datasets are:

- **Small** (10-100k peptides per study)
- **Inconsistent** (different instruments, methods, quality)
- **Incomplete** (missing ion mobility, limited charge states)

Building better models requires a **unified, high-quality database** with millions of validated peptide observations.

## The Solution: San José

San José is an automated pipeline that:

1. **Mines public repositories** (PRIDE) for timsTOF datasets
2. **Searches with three orthogonal methods** (FragPipe + DIA-NN + Sage) for validation
3. **Extracts raw signal features** (XICs, mobilograms, isotope envelopes)
4. **Builds a consensus database** with confidence-weighted peptide properties

### Key Innovation: Triple Orthogonal Validation

By running **three independent search engines** on each dataset, we achieve:

| Metric | PXD019086 (HeLa) |
|--------|------------------|
| FragPipe (spectrum-centric) | 51,064 |
| DIA-NN (peptide-centric) | 33,747 |
| Sage (fast Rust-based) | 42,787 |
| **All three engines** | **25,382** |
| At least two engines | 44,382 |
| UNION (any engine) | 57,834 |

| Validation Rates | |
|------------------|---|
| FragPipe (confirmed by others) | 86.3% |
| DIA-NN (confirmed by others) | 89.2% |
| Sage (confirmed by others) | **93.4%** |

**76.7% of identifications are confirmed by at least two engines.** Triple orthogonal agreement provides publication-quality confidence without manual validation.

**Runtime comparison (3 test files):**
| Engine | Time |
|--------|------|
| FragPipe | 73 min |
| DIA-NN | 260 min |
| Sage | **4 min** |

## Architecture

```
PRIDE Repository
       │
       ▼
┌─────────────────────────────────┐
│  Worker (1 per PRIDE accession) │
├─────────────────────────────────┤
│  1. Download raw files          │
│  2. FragPipe search (73 min)    │
│  3. DIA-NN search (260 min)     │
│  4. Sage search (4 min)         │
│  5. Compute 3-way consensus     │
│  6. Extract raw features        │
│  7. Upload to San José DB       │
└─────────────────────────────────┘
       │
       ▼
  San José Database
  (Parquet/DuckDB)
       │
       ▼
  Model Training
  (CCS, RT, MS2)
```

Workers are **self-contained** and can run on HPC clusters or collaborator machines. Each processes exactly one PRIDE accession atomically.

## Database Schema

| Table | Contents | Est. Size |
|-------|----------|-----------|
| `peptides` | Sequence, charge, RT, CCS, intensity | ~100 bytes/peptide |
| `psm_scores` | Hyperscore, PEP, q-value per engine | ~50 bytes/PSM |
| `raw_features` | XIC, mobilogram, isotope envelope | ~5 KB/precursor |
| `consensus` | OVERLAP/UNION flags, confidence weights | ~20 bytes/peptide |

**Estimated scale:** 1000 datasets × 50k peptides = **50M peptide observations**

## Current Status

- Pipeline running on PXD019086 (HeLa benchmark)
- FragPipe + DIA-NN integration complete
- Overlap analysis with I/L normalization and UNIMOD standardization
- HTML report generation for per-dataset QC

## Roadmap

| Phase | Deliverable | Timeline |
|-------|-------------|----------|
| 1 | Database schema + single dataset processing | Done |
| 2 | Score extraction + consensus tables | 2 weeks |
| 3 | Raw feature extraction (imspy) | 4 weeks |
| 4 | Worker containerization + HPC deployment | 6 weeks |
| 5 | Scale to 100+ datasets | Ongoing |

## Why "San José"?

Named after the San José scale insect - a prolific species that spreads across orchards. Like the insect, this database aims to systematically spread across public proteomics data, extracting value from every dataset.

## Resource Requirements

| Resource | Requirement |
|----------|-------------|
| Storage | ~2 TB (raw: 1TB, processed: 500GB, database: 500GB) |
| Compute | 32 cores, 64GB RAM per worker |
| Time | ~4 hours per dataset (search + extraction) |

## Collaboration Opportunity

The worker design enables **distributed processing**:

1. We provide the container + job queue
2. Collaborators run workers on their infrastructure
3. Results automatically upload to central database
4. Everyone benefits from the aggregated data

This could accelerate database growth significantly while distributing compute costs.

## Contact

[Your contact information]

---

*Generated from the San José Pipeline - claudius-proteomics*
