# San José

> A reproducible, metadata-rich, bias-aware reference layer for timsTOF data on PRIDE.

## Vision

**We are not collecting peptides. We are collecting peptide observations in experimental context.**

San José systematically reprocesses every timsTOF dataset on PRIDE to build a versioned, citable dataset for CCS, RT, and MS2 model training — grounded in real experimental physics across thousands of labs.

```
PRIDE Archive                    San José Database              Models
┌─────────────┐                 ┌─────────────────┐           ┌─────────────┐
│ PXD019086   │ ──┐             │ precursors/     │           │ ccs_v1/     │
│ PXD010012   │   │  download   │   PXD019086     │ snapshot  │   model.h5  │
│ PXD...      │   ├──────────►  │   PXD010012     │ ────────► │   metrics   │
│             │   │  search     │ snapshots/      │  train    │             │
└─────────────┘ ──┘  extract    │   v1.0/         │           └─────────────┘
                                └─────────────────┘
```

See [docs/DESIGN.md](docs/DESIGN.md) for the full project plan and [docs/SAN_JOSE_PITCH.md](docs/SAN_JOSE_PITCH.md) for the vision.

## 🤗 The corpus is public

**[huggingface.co/datasets/theGreatHerrLebert/timstof-dda-pasef-cc0](https://huggingface.co/datasets/theGreatHerrLebert/timstof-dda-pasef-cc0)** — CC0-1.0, no gate, no request form.

Release **v0.2**: **58 datasets, 33.7M precursors, 504M b/y fragments**, every one
of them reprocessed from raw `.d` by this pipeline rather than harvested from
submitter search results.

| Config | Rows | What it is |
|---|---|---|
| `tier1_psms` | 33.7M | one row per precursor: sequence, charge, CCS, aligned RT, m/z |
| `tier3_fragments` | 504M | matched b/y fragments with intensities |

Split 27.2M / 3.27M / 3.27M train/validation/test by **peptide hash** on
`sequence_normalized` (seed 0), so no peptide appears in more than one split.
Dataset-level splitting is deferred to a later version — a model scored on the
test split is held out from the *peptides*, not from the labs or instruments.
Schema in [docs/CORPUS_SCHEMA.md](docs/CORPUS_SCHEMA.md), card in
[docs/DATASET_CARD.md](docs/DATASET_CARD.md).

```python
from datasets import load_dataset
psms = load_dataset("theGreatHerrLebert/timstof-dda-pasef-cc0", "tier1_psms", split="train")
```

**Known composition limits, stated up front** — 97.6% of rows are human, there is
no timsTOF Ultra data in v0.2, and the top 8 labs contribute about half the rows.
Broadening this is what v0.3+ is for; do not read v0.2 as a species-balanced
sample.

## Key Design Principles

### Orthogonal Validation

Every dataset is searched by **two independent engines**, both reading Bruker
`.d` directly:

| Engine | Approach | Role | Status |
|--------|----------|------|--------|
| **Sage** | Spectrum-centric | Fast Rust engine, PTM-aware | production |
| **FragPipe** | Spectrum-centric | Gold standard for DDA | production |
| **DIA-NN** | Peptide-centric | Library-free, NN predictions | evaluated, deferred |

We store the full union, 2-engine agreement, and engine-specific unique IDs.
Disagreements are scientifically valuable and are kept, not resolved away.

**Why DIA-NN is not in the production path.** A harness exists
(`runner/engines/diann_job.py`, `rules/diann.smk`) and it was benchmarked as a
third vote rather than assumed useful. At the FDR the corpus publishes
(q ≤ 0.01), Sage ∩ FragPipe already agree on 1.08M precursors; DIA-NN adds
+91,611 cross-engine corroborations (+8.5%) and 10,500 unique precursors (0.7%
of the union). That is a confidence tier, not a coverage win, and it does not
justify its cost per dataset — DDA search is the slow CPU path, and the Linux
builds are CPU-only. If it is enabled, use **2.5.1**: measured on DDA-PASEF,
2.6.0 is 2.0× slower for identical IDs. Revisit for v0.4.

A caveat worth stating because it shapes the schema: engine agreement is
recorded at the level of the peptidoform each engine reported, and the engines
are not configured identically for PTM datasets. See the conformance note under
*Configuration conformance* below.

### Bias-Awareness

PRIDE datasets cluster heavily by lab, sample prep, column chemistry, and gradient length. San José tracks:

```
lab_id, dataset_id, organism, gradient_length, column_type, acquisition_mode
```

This enables stratified sampling, bias analysis, and cross-lab validation.

### Configuration Conformance

Asking two engines the same question is harder than it looks. A dataset's
`mod_profile` is the intended search configuration; each engine renders it into
its own dialect, and the rendering can silently disagree with the intent.

Measured on a PTM cohort: from one identical `mod_profile`, the rendered Sage
config was **97% conformant** with the requested modifications while the rendered
FragPipe workflow was **77.5%** — FragPipe was running effectively PTM-blind for
part of the cohort. Downstream that shows up as a peptidoform that has had its
PTM stripped, which is a corpus correctness issue, not a cosmetic one.

The pipeline therefore emits a **TrustReport** per rendered engine config,
recording what was asked for versus what the engine was actually told, so
non-conformance is visible in provenance instead of being discovered later in the
data. Reports land in `data/provenance/<accession>/trustreport-*.json`.

### San José Mode (Report ALL Results)

All search engines are configured to report ALL PSMs without FDR filtering. San José stores the full union with per-engine scores; thresholds are applied at query time for downstream analysis.

### Human Checkpoints

| Checkpoint | Purpose |
|------------|---------|
| Dataset Selection | Curated PRIDE accessions, blacklist for problematic datasets |
| QC Dashboard | Auto-generated reports, anomaly detection, human review when flagged |
| Consensus Rules | Configurable: 2/2 agreement vs union, q-value thresholds (3/3 only if DIA-NN is enabled) |
| Versioned Snapshots | San José v1.0, v1.1… frozen, reproducible datasets |
| Release Gates | Bias checks, cross-lab validation, human sign-off |

## Quick Start

### Prerequisites

1. **FragPipe** (v24+) - https://fragpipe.nesvilab.org/ (academic license)
2. **DIA-NN** (2.5.1 recommended; *optional* — not in the production path, see above) - https://github.com/vdemichev/DiaNN (academic license)
3. **Sage** - https://github.com/lazear/sage (open source Rust)
4. **Snakemake** (v8+) - Workflow orchestration
5. **Singularity** (optional) - For containerized HPC execution

### Configuration

Edit `config/config.yaml`:

```yaml
# Search engine paths
fragpipe:
  path: "/path/to/fragpipe-24.0"
  workflow: "LFQ-MBR"

diann:
  path: "/path/to/diann-linux"
  threads: 16
  qvalue: 1.0  # San José: report ALL results

sage:
  path: "/path/to/sage"
  config: "config/sage_config.json"

# San José mode (default)
san_jose:
  report_all: true

# Organism FASTA mapping
organisms:
  human:
    proteome_id: "UP000005640"
    local_fasta: "/path/to/human.fasta"

# Dataset to organism mapping
dataset_metadata:
  PXD019086:
    organism: "human"
```

### Running the Pipeline

The 5-step runner processes a dataset with checkpointing:

```bash
# Full pipeline (download → search → stratify → extract → merge)
.venv/bin/python runner/run_dataset.py PXD019086 --config config/config.yaml

# Test mode (1 file for quick validation)
.venv/bin/python runner/run_dataset.py PXD019086 --test-mode --max-files 1

# Use existing local data
.venv/bin/python runner/run_dataset.py PXD019086 --local-data data/raw/PXD019086

# Resume from checkpoint after failure
.venv/bin/python runner/run_dataset.py PXD019086 --resume
```

Pipeline steps:
1. **Download** - Fetch .d files from PRIDE (or link local data)
2. **Search** - Run Sage + FragPipe on .d files (DIA-NN optional, off by default)
3. **Stratify** - Unify results, compute engine agreement
4. **Extract** - Extract raw 4D signal via imspy
5. **Merge** - Join IDs with raw features → `precursor_store.parquet`
6. **Package** (optional) - Create distributable archive

```bash
# Create self-contained archive for sharing/upload
.venv/bin/python runner/run_dataset.py PXD019086 --package --package-version 1.0
# Output: data/packages/PXD019086_v1.0.zip (~2.8 GB)
```

```bash
# Phase 2: Train models
snakemake train_model --config snapshot=v1.0 model=ccs_v1
```

### Interactive Dashboard

Browse precursors with full 4D visualization:

```bash
# Single dataset mode
.venv/bin/python dashboard/backend/main.py \
    --store data/merged/PXD019086/precursor_store.parquet \
    --port 8000

# Collection mode (browse multiple datasets)
.venv/bin/python dashboard/backend/main.py \
    --collection /path/to/san_jose_collection/ \
    --port 8000

# Frontend (React + TypeScript)
cd dashboard/frontend && npm run dev  # http://localhost:5173
```

Features: fragment spectrum, IM vs m/z scatter, XIC, mobilogram, isotope envelope.

### Collection Management

Organize processed datasets into study-based collections:

```bash
# Create collection with studies.yaml
mkdir my_collection
cp config/studies.yaml.example my_collection/studies.yaml
# Edit studies.yaml to define your studies

# Add dataset archive to collection
python scripts/add_to_collection.py \
    --archive data/packages/PXD019086_v1.0.zip \
    --collection my_collection/ \
    --study meier_2021_ccs

# Rebuild manifest from existing files
python scripts/add_to_collection.py --collection my_collection/ --rebuild
```

Collection structure:
```
my_collection/
├── studies.yaml                # Study definitions
├── collection_manifest.json    # Auto-generated index
├── meier_2021_ccs/            # Study folder
│   └── PXD019086_v1.0/        # Extracted dataset
└── archives/                   # Original zip backups
```

## Project Structure

```
claudius-proteomics/
├── Snakefile                   # Main workflow
├── config/
│   ├── config.yaml             # Pipeline configuration
│   └── sage_config.json        # Sage search parameters
├── profiles/mogon2/            # SLURM profile (name predates the move to MOGON-NHR)
│
├── rules/                      # Snakemake rule modules
│   ├── fasta.smk               # FASTA resource management
│   ├── download.smk            # PRIDE download
│   ├── fragpipe.smk            # FragPipe search
│   ├── diann.smk               # DIA-NN search
│   ├── sage.smk                # Sage search
│   ├── extract.smk             # Raw feature extraction
│   ├── database.smk            # Database insertion
│   ├── snapshot.smk            # Snapshot creation
│   └── train.smk               # Model training
│
├── scripts/
│   ├── run_fragpipe.py         # FragPipe wrapper (San José mode)
│   ├── build_precursor_index.py # Unified precursor index
│   ├── add_to_collection.py    # Collection management
│   ├── sequence_utils.py       # UNIMOD normalization
│   ├── fragment_matching.py    # Theoretical fragment matching
│   ├── precursor_store.py      # Parquet store utilities
│   └── ...
│
├── dashboard/
│   ├── backend/main.py         # FastAPI server
│   └── frontend/               # React + TypeScript + Deck.gl
│
├── docs/
│   ├── SAN_JOSE_PITCH.md       # Vision and architecture
│   ├── RAW_DATA_ARCHITECTURE.md # Parquet schema, dashboard API
│   ├── RUNNER_OUTPUT_SCHEMA.md # Per-dataset output format
│   ├── CORPUS_SCHEMA.md        # Published HF corpus schema
│   └── DATASET_CARD.md         # HF dataset card
│
├── runner/                    # 6-step pipeline runner
│   ├── run_dataset.py         # Main entry point
│   ├── state.py               # Checkpoint management
│   └── steps/                 # Step implementations (step1-6)
│
├── data/
│   ├── raw/{accession}/       # Step 1: Bruker .d files
│   ├── processed/{accession}/ # Steps 2-3: Search results + consensus
│   ├── extracted/{accession}/ # Step 4: Raw 4D features
│   ├── merged/{accession}/    # Step 5: Final precursor_store.parquet
│   ├── packages/              # Step 6: Distributable archives
│   └── checkpoints/           # Pipeline state (state.json)
│
├── paper_inbox/               # Per-dataset publication pointers (SOURCE.md)
├── LICENSE
└── README.md
```

## Data Storage

### Precursor Parquet Schema

Each precursor observation contains:

```
# Identity
precursor_id, raw_file, mz, charge, rt_seconds, mobility

# Search engine results (NULL if not identified)
fragpipe_peptide, sage_peptide, diann_peptide   # diann_* stay NULL unless DIA-NN is enabled
consensus_peptide, n_engines                    # n_engines is 0-2 in the production path

# Quality metrics
fragpipe_probability, diann_qvalue, sage_qvalue

# Fragment spectrum (variable length lists)
fragment_mz, fragment_intensity, fragment_mobility

# MS1 projections
xic_rt, xic_intensity
mobilogram_im, mobilogram_intensity
isotope_mz, isotope_intensity

# Raw 4D MS1 data
raw_rt, raw_mz, raw_mobility, raw_intensity
```

See [docs/RAW_DATA_ARCHITECTURE.md](docs/RAW_DATA_ARCHITECTURE.md) for full specification.

## How the Models Are Trained

The corpus exists to train peptide-property predictors. All three are
**fine-tuned from public base models rather than trained from scratch**, and each
was evaluated on held-out data before being deployed.

| Property | Model | Held-out result |
|---|---|---|
| Retention time | [Chronologer](https://github.com/searlelab/chronologer) (Searle Lab residual CNN, Apache-2.0) | ~4× tighter median residual than the imspy transformer on timsTOF data |
| CCS / ion mobility | fine-tuned `DeepPeptideIonMobilityApex` | median abs. residual **38.0 → 17.6 mK0** (−54%), Pearson r **0.84 → 0.96** |
| Fragment intensity | fine-tuned `DeepPeptideIntensityPredictor` | spectral angle **0.689 → 0.833** |

These ship as rustims release `models-v0.6.0`, so a bare
`DeepPeptideIonMobilityApex()` or `DeepPeptideIntensityPredictor()` resolves the
fine-tuned weights.

**Method notes.** The CCS and intensity fine-tunes use three modern deposits
(pig `PXD046675`, HIV `PXD046777`, Chlamydomonas `PXD068782`), chosen because
they post-date the early-timsTOF era the base models were trained on — headline
metrics on training-era data over-count noise that modern instruments no longer
produce. Intensity is fine-tuned at **per-PSM optimal collision energy** rather
than one NCE per deposit; the per-PSM optimum has std ≈ 10 NCE across 24k PSMs,
so a uniform value averages over real variation (+0.005 SA for the per-PSM
paradigm). Generalisation was checked out-of-distribution on Meier `PXD019086`
(early timsTOF): SA 0.68 → 0.80.

**Where the gain actually lands.** Same weights, two consumers, on the same
diaPASEF bundle:

| Stack | stock | fine-tuned | Δ |
|---|---|---|---|
| spectrum-centric rescoring | 28,369 | 28,822 | +1.6% |
| peptide-centric, top-k=8 | 12,772 | 23,262 | **+82%** |

The asymmetry is the finding, not noise. Peptide-centric extraction is
detection-by-predicted-spectrum, so predictor accuracy is a *prerequisite* for
detection; spectrum-centric search has observed spectra to score against, and the
predictor is one feature among many for the discriminator.

**Training directly on the corpus** (`peptide-property-ng`, a unified
intensity/CCS/RT/charge model) currently reaches SA **0.798**, CCS MAE **0.0148**,
RT r **0.918** on a 56-dataset fine-tune.

**The current ceiling is the target encoding, not the data.** b/y ions account
for only ~17.5% of observed MS2 intensity, which bounds achievable spectral
angle regardless of how many datasets are added. Extending the encoding —
neutral losses and immonium ions are implemented and validated, internal ions
next — is the lever that moves it.

## Technical Stack

- **Orchestration**: Snakemake v8+
- **Raw processing**: rustims/imspy (Rust backend + Python bindings)
- **Search engines**: Sage + FragPipe (DIA-NN available, not in the production path)
- **Database**: DuckDB/Parquet (columnar, versioned snapshots)
- **Dashboard**: FastAPI + React + TypeScript + Deck.gl (WebGL)
- **ML**: PyTorch via imspy-predictors (`torch>=2.0`); Chronologer for RT
- **Containers**: Singularity (HPC) / Docker (local)
- **HPC**: SLURM (MOGON-NHR @ JGU Mainz, cluster `mogonki`)

## Data Sources

### POC: Meier et al. 2021

- **Paper**: [Deep learning the collisional cross sections of the peptide universe](https://doi.org/10.1038/s41467-021-21352-8)
- **PRIDE**: PXD019086 (primary), PXD010012, PXD017703
- **Data**: >1 million CCS measurements from timsTOF
- **Target**: R > 0.99, <1.4% median error

## Why "San José"?

Named after the *San José*, a sunken ship whose rediscovery revealed immense, carefully preserved value beneath the surface.

PRIDE is similar: a vast repository where the real value is hidden in raw experimental data, waiting to be systematically recovered, catalogued, and understood.

## License

Code: **MIT** — see [LICENSE](LICENSE).

Corpus: **CC0-1.0**. Every dataset in the release was license-verified against
the PRIDE API before inclusion; nothing non-CC0 is redistributed.

## Citation

A manuscript is in preparation. In the meantime, cite the corpus by its
Hugging Face repository:

```
timstof-dda-pasef-cc0 (v0.2). https://huggingface.co/datasets/theGreatHerrLebert/timstof-dda-pasef-cc0
```

The approach was presented as poster #1015 at ECCB 2026.
