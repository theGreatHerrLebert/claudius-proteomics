# San José: Let's Build the Largest Peptide Property Database Ever

## The Vision

We're going to **mine every single timsTOF dataset on PRIDE** and build a peptide property database that dwarfs anything that exists today. Not 100k peptides. Not 1M. We're talking **50M+ validated peptide observations** with CCS, retention time, and raw signal features.

And we're going to do it *fast*.

## Why Now?

Three things just clicked into place:

1. **Sage** - A Rust-based search engine that's 65x faster than DIA-NN. We can now run triple orthogonal validation in the time it used to take to run one engine.

2. **Claude Code** - AI-assisted development means we can iterate on pipeline code at unprecedented speed. The entire 3-way validation system was built and tested in one session.

3. **Public data explosion** - PRIDE has thousands of timsTOF datasets sitting there, waiting to be systematically processed.

The window is open. Let's go.

## The Approach: Automation WITH Human Checkpoints

This is **not** about removing humans from the loop. That's the wrong way to think about it.

This is about **putting humans at the RIGHT checkpoints** while automating everything in between:

```
┌─────────────────────────────────────────────────────────────────┐
│                    HUMAN CHECKPOINTS                            │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  [1] DATASET SELECTION                                          │
│      - Curated list of PRIDE accessions                         │
│      - Organism/instrument filters                              │
│      - Blacklist for known problematic datasets                 │
│                                                                 │
│  [2] QC DASHBOARD                                               │
│      - Per-dataset HTML reports (auto-generated)                │
│      - Overlap statistics, charge distributions                 │
│      - Flag anomalies for human review                          │
│                                                                 │
│  [3] CONSENSUS THRESHOLDS                                       │
│      - Configurable: "all 3 engines" vs "at least 2"            │
│      - PEP thresholds (we use 0.05 = 95% confidence)            │
│      - Human decides quality/quantity tradeoff                  │
│                                                                 │
│  [4] MODEL TRAINING SNAPSHOTS                                   │
│      - Versioned database snapshots                             │
│      - Human approves before training                           │
│      - Reproducible: snapshot v1.0 always = same data           │
│                                                                 │
│  [5] RELEASE GATES                                              │
│      - Model performance benchmarks                             │
│      - Human sign-off before public release                     │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

**Between checkpoints?** Full automation. Workers churn through PRIDE accessions 24/7. No babysitting required.

## Triple Orthogonal Validation: The Numbers

We ran PXD019086 (HeLa benchmark) through three independent search engines:

| Engine | Philosophy | Precursors | Runtime | Validation Rate |
|--------|------------|------------|---------|-----------------|
| FragPipe | Spectrum-centric | 48,799 | 73 min | 81.1% |
| DIA-NN | Peptide-centric | 31,690 | 260 min | 90.2% |
| **Sage** | Fast Rust-based | 32,788 | **4 min** | **97.1%** |

**Consensus (PEP ≤ 0.05, 95% confidence):**

| Level | Precursors | % of Union |
|-------|------------|------------|
| All 3 engines agree | 20,559 | 38.8% |
| At least 2 engines | 39,720 | **75.0%** |
| Union (any engine) | 52,998 | 100% |

**Key insight:** Sage at 95% confidence has only 965 unique identifications (2.9%). Almost everything Sage finds is confirmed by another engine. That's the power of orthogonal validation.

## Architecture

```
PRIDE Repository (1000s of datasets)
            │
            ▼
    ┌───────────────┐
    │  Job Queue    │◄─── [HUMAN] Dataset selection + blacklist
    └───────────────┘
            │
            ▼
┌─────────────────────────────────────┐
│  Worker Pool (HPC / Collaborators)  │
├─────────────────────────────────────┤
│  Each worker:                       │
│  • Pulls 1 PRIDE accession          │
│  • Downloads raw .d files           │
│  • Runs FragPipe → DIA-NN → Sage    │
│  • Computes 3-way consensus         │
│  • Generates QC report              │
│  • Uploads to central DB            │
│  • Moves to next accession          │
└─────────────────────────────────────┘
            │
            ▼
    ┌───────────────┐
    │  QC Dashboard │◄─── [HUMAN] Review flagged datasets
    └───────────────┘
            │
            ▼
    ┌───────────────┐
    │  San José DB  │◄─── [HUMAN] Approve snapshot for training
    │  (DuckDB)     │
    └───────────────┘
            │
            ▼
    ┌───────────────┐
    │  Model Train  │◄─── [HUMAN] Validate before release
    │  (CCS/RT/MS2) │
    └───────────────┘
```

## What We Have Working Today

- **Triple search engine pipeline** (FragPipe + DIA-NN + Sage)
- **3-way overlap analysis** with I/L normalization, UNIMOD standardization
- **PEP filtering** (Probability ≥ 0.95 for FragPipe, PEP ≤ 0.05 for others)
- **HTML QC reports** - self-contained, ready for human review
- **Consensus parquet export** - overlap.parquet, union.parquet

## What's Next

| Phase | What | Human Checkpoint |
|-------|------|------------------|
| **Now** | Process 10 more HeLa datasets | Review QC reports |
| **Week 2** | Raw feature extraction (XICs, mobilograms) | Validate extraction quality |
| **Week 4** | Worker containerization | Test on HPC |
| **Week 6** | Scale to 100 datasets | Approve v1.0 snapshot |
| **Ongoing** | Process all timsTOF on PRIDE | Periodic QC review |

## The Ask

1. **Compute time** - HPC allocation for worker pool
2. **Storage** - ~2TB for raw + processed + database
3. **Occasional human review** - QC dashboard checks, snapshot approval

In return: **The largest validated peptide property database for ion mobility MS.**

## Why "San José"?

Named after the San José scale insect - a species so prolific it spreads across entire orchards. That's exactly what we're doing: systematically spreading across PRIDE, extracting value from every dataset.

Let's scale this thing.

---

*Built with the San José Pipeline (claudius-proteomics) + Claude Code*
