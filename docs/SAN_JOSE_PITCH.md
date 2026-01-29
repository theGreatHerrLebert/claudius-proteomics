# San José: A Reprocessed, Bias-Aware Reference Layer for timsTOF Data on PRIDE

## The Vision

PRIDE is not just an archive.

It is the **largest unmined experimental corpus of LC-IMS-MS/MS data ever created**.

We are going to systematically reprocess every timsTOF dataset on PRIDE and build:

> **A reproducible, metadata-rich, bias-aware reference layer of validated peptide observations**

This is not “a big peptide database”.

This is:

> **PRIDE-IMS-Processed** — a standardized, re-interpreted layer on top of PRIDE for ion mobility proteomics.

The outcome is a versioned dataset that can be cited, reused, and used as a training corpus for CCS, RT, and MS2 models — grounded in real experimental physics across thousands of labs.

---

## Why This Represents a Paradigm Shift (for PIs and Reviewers)

Until now, we had no practical way to systematically compare peptide behavior across labs, instruments, gradients, and experiments at scale.

PRIDE has contained the raw material for this for years, but the engineering cost of reprocessing thousands of raw datasets in a consistent way made such questions infeasible.

That cost barrier is now gone.

For the first time, we can treat PRIDE not as an archive of independent studies, but as a **single, reprocessable experimental corpus**.

This enables questions that previously could only be asked within a single study, for example:

* How stable is CCS across labs and instrument tunings?
* How does gradient length affect RT reproducibility globally?
* Where do search engines systematically disagree, and why?
* What constitutes a truly "validated peptide observation" when seen across hundreds of experiments?

A useful analogy is:

> PRIDE today is like having thousands of microscopy images without a common coordinate system.
> San José provides the coordinate system.

This project changes what we can consider *data* in proteomics: not isolated datasets, but a unified, comparable experimental space.

---

## Why Now?

Three things changed recently:

1. **Sage** makes high-quality searching cheap and fast enough to run multi-engine validation at scale.
2. **rustims + timsTOF raw data expertise** makes feature extraction (XICs, mobilograms, peak shapes) feasible and reproducible.
3. **LLM-assisted pipeline engineering** makes it realistic to build and iterate on complex processing infrastructure quickly.

For the first time, the bottleneck is not engineering. It is **designing the data layer correctly**.

---

## Core Principle

> We are not collecting peptides.
> We are collecting **peptide observations in experimental context**.

Every stored observation is inseparable from:

* Lab / dataset provenance
* Instrument configuration
* Gradient / LC conditions
* Acquisition method
* Engine agreement / disagreement profile
* Raw signal features

This turns the project from “big data” into:

> **A research platform for understanding peptide behavior across experimental space.**

---

## Triple Orthogonal Validation as a Curation Primitive

Each dataset is processed with:

* FragPipe (spectrum-centric)
* DIA-NN (peptide-centric)
* Sage (fast, independent Rust engine)

We do **not** treat consensus as “ground truth”.

Instead, we store:

* Full union
* 2/3 agreement
* 3/3 agreement
* Engine-specific unique IDs

Because:

> The most scientifically valuable data are often where engines disagree.

Consensus reduces random FPs.
Disagreement reveals systematic assumptions.

Both are first-class citizens in San José.

---

## The Hidden Problem We Address Explicitly: Lab Bias

PRIDE datasets are not independent samples. They cluster heavily by:

* Lab
* Sample prep
* Column chemistry
* Gradient length
* Instrument tuning
* Organism

If aggregated naively, a model trained on this data will learn:

> “How 20 labs run timsTOF” instead of “how peptides behave”.

Therefore, from day one, San José tracks and exposes:

```
lab_id
dataset_id
organism
gradient_length
column_type
acquisition_mode
instrument_metadata
```

This allows:

* Stratified sampling
* Bias analysis
* Weighted model training
* Cross-lab validation studies

This is a **design requirement**, not an afterthought.

---

## What an Entry in the San José Database Looks Like

A single peptide observation contains:

* Sequence (+ UNIMOD PTMs)
* Charge
* CCS
* RT
* MS2 features
* XIC / mobilogram characteristics
* Dataset / lab metadata
* Engine scores and PEPs
* Engine agreement profile (FragPipe / DIA-NN / Sage)

This makes every row a **fully contextualized physical observation**, not just an identification.

---

## Human Checkpoints in the Right Places

Automation runs the pipeline. Humans guard the scientific integrity.

### [1] Dataset Selection

* Curated PRIDE accessions
* Blacklist for problematic datasets
* Instrument / organism filters

### [2] QC Dashboard

* Auto-generated per-dataset reports
* Charge distributions, overlap stats, anomaly detection
* Human review only when flagged

### [3] Consensus & Inclusion Rules

* Configurable: 3/3 vs 2/3 vs union
* PEP thresholds
* Explicit quality vs quantity decision

### [4] Versioned Snapshots

* San José v1.0, v1.1, …
* Frozen, reproducible datasets for model training
* Snapshot always reproducible from PRIDE + pipeline version

### [5] Release Gates

* Bias checks
* Cross-lab validation tests
* Human sign-off before public release

---

## Architecture Overview

```
PRIDE (timsTOF datasets)
            │
            ▼
      Job Queue  ◄── [Human dataset curation]
            │
            ▼
  Worker Pool (HPC / collaborators)
            │
            ├─ FragPipe
            ├─ DIA-NN
            ├─ Sage
            ├─ Feature extraction (rustims)
            ├─ QC report
            ▼
        San José DB (DuckDB/Parquet)
            │
            ▼
     Versioned Snapshots (v1.0, v1.1…)
            │
            ▼
   Model training / community reuse
```

---

## What This Enables (Beyond “a big DB”)

Researchers can ask:

* How does CCS drift across labs?
* How does gradient length affect RT reproducibility?
* Where do search engines systematically disagree?
* How does instrument tuning affect mobility distributions?
* How reproducible are peptide properties across thousands of experiments?

This becomes a **scientific instrument**, not just a dataset.

---

## What Exists Today

* Triple-engine pipeline
* Overlap analysis with UNIMOD/I-L normalization
* PEP-filtered consensus export
* Automated HTML QC reports
* Parquet exports for union and consensus sets

---

## Roadmap

| Phase   | Goal                   | Human Checkpoint               |
| ------- | ---------------------- | ------------------------------ |
| Now     | 10 HeLa datasets       | QC report review               |
| Week 2  | Raw feature extraction | Validate feature quality       |
| Week 4  | Containerized workers  | HPC test run                   |
| Week 6  | 100 datasets           | Approve San José v1.0 snapshot |
| Ongoing | All timsTOF on PRIDE   | Periodic bias & QC review      |

---

## The Ask

* HPC compute allocation
* ~2TB storage
* Occasional QC review

---

## The Outcome

San José is not:

> “The largest peptide database”

It is:

> **The reference, reproducible, bias-aware reprocessing layer for timsTOF data on PRIDE**

A dataset that can be cited like ProteomeTools — but based on real experimental diversity across thousands of labs.

---

## Why “San José”?

Named after the *San José*, a sunken ship whose rediscovery revealed immense, carefully preserved value beneath the surface.

PRIDE is similar: a vast repository where the real value is hidden in raw experimental data, waiting to be systematically recovered, catalogued, and understood.

This project is about diving down, recovering that value methodically, and bringing it back in a structured, usable form for the community.
