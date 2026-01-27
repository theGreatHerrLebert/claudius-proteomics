# CLAUDIUS-PROTEOMICS

> Building the world's largest peptide property prediction resource through systematic reanalysis of public proteomics data.

---

## 1. Vision (Full Scope)

### The Problem
Public proteomics repositories (PRIDE, MassIVE, PeptideAtlas) contain petabytes of raw mass spectrometry data. This data is:
- Processed with inconsistent, often outdated pipelines
- Underutilized for machine learning applications
- Fragmented across studies with incompatible formats

### The Solution
CLAUDIUS-PROTEOMICS will:
1. **Systematically download** raw data from public repositories
2. **Reprocess uniformly** with state-of-the-art tools
3. **Extract peptide properties** (CCS, retention time, fragmentation spectra, charge states)
4. **Build prediction models** trained on this unprecedented scale of consistent data

### Target Properties (Long-term)
| Property | Instrument Requirement | Use Case |
|----------|------------------------|----------|
| Collisional Cross Section (CCS) | Ion mobility (timsTOF, FAIMS) | Ion mobility prediction, validation |
| Retention Time (RT) | Any LC-MS | Peptide identification confidence |
| MS2 Fragmentation | Any MS/MS | Spectral library generation |
| Charge State Distribution | Any MS | Improved precursor selection |
| Detectability | Any MS | Proteome coverage optimization |

---

## 2. Proof of Concept: CCS Prediction from timsTOF Data

### Source Study

**Meier F, Köhler ND, Brunner AD, Wanka JMH, Voytik E, Strauss MT, Theis FJ, Mann M.**
*Deep learning the collisional cross sections of the peptide universe from a million experimental values.*
Nature Communications 12, 1185 (2021). https://doi.org/10.1038/s41467-021-21352-8

#### Dataset Summary
| Metric | Value |
|--------|-------|
| Data points | >1 million CCS measurements |
| Organisms | 5 (whole-proteome digests) |
| Instrument | timsTOF Pro (TIMS-PASEF) |
| CCS precision | CV < 1% |
| Published model performance | R > 0.99, median relative error 1.4% |

#### PRIDE Accessions
| Accession | Description |
|-----------|-------------|
| **PXD019086** | Primary dataset - MS raw files + MaxQuant output |
| **PXD010012** | Previously acquired HeLa data |
| **PXD017703** | diaPASEF raw files |

#### Code Repositories (Meier et al.)
- Data analysis: https://github.com/MannLabs/DeepCollisionalCrossSection
- Deep learning model: https://github.com/theislab/DeepCollisionalCrossSection

### Scope Constraints
| Dimension | POC Choice | Rationale |
|-----------|------------|-----------|
| Instrument | timsTOF only | Native CCS measurement via TIMS |
| Processing | FragPipe | SOTA, handles timsTOF, open source |
| Data source | **PXD019086** (Meier et al.) | Landmark CCS study, >1M data points, CV <1% |
| Target property | CCS | Well-defined, measurable ground truth |
| Model | ionmob (GRU-based) | Existing architecture, proven performance |
| ML Framework | imspy-predictors | Integrated with rustims, TensorFlow-based |

### Existing Stack (rustims ecosystem)

We leverage an existing, battle-tested infrastructure:

| Component | Repository | Purpose |
|-----------|------------|---------|
| **rustims** | [theGreatHerrLebert/rustims](https://github.com/theGreatHerrLebert/rustims) | Rust backend for timsTOF data processing |
| **mscore** | rustims/mscore | In-memory data structures, algorithms |
| **rustdf** | rustims/rustdf | Read/write Bruker TDF files (.d folders) |
| **imspy** | rustims/imspy | Python interface, TensorFlow integration |
| **imspy-predictors** | rustims/packages/imspy-predictors | CCS, RT, intensity prediction |
| **ionmob** | [theGreatHerrLebert/ionmob](https://github.com/theGreatHerrLebert/ionmob) | GRU-based CCS predictor |

#### ionmob Architecture (in imspy-predictors)
- **Two-component model**: Initial projection (mass + charge → baseline CCS) + GRU refinement
- **Input**: Peptide sequence, charge state, PTMs
- **Performance**: Competitive with Meier et al. on tryptic peptides
- **Extensions**: Non-tryptic, phosphorylated peptides, MHC ligands

#### imspy-predictors Capabilities
- CCS/ion mobility prediction (ionmob-based)
- Retention time prediction (GRU-based)
- Fragment intensity prediction (Prosit 2023 wrapper)
- Charge state distribution (binomial + deep learning)
- Koina integration for remote model access
- TensorFlow backend

### Dual Extraction Strategy

**Why not just use FragPipe output?**

FragPipe/MaxQuant export summarized results (PSMs, quantification). We go deeper by extracting **raw signal features** directly from timsTOF .d files via imspy/rustdf:

| Data Type | FragPipe | imspy (raw extraction) |
|-----------|----------|------------------------|
| PSM identifications | ✓ | — |
| CCS values | ✓ | ✓ (full distribution) |
| Retention time | ✓ (apex) | ✓ (full chromatogram) |
| Ion mobility | summarized | full mobilogram |
| Peak shape | — | ✓ |
| Isotope distribution | — | ✓ |
| Signal intensity | summarized | raw trace |

**This enables prediction targets no one else has:**
- Peak shape prediction
- Isotope pattern prediction
- Ion mobilogram shape
- Chromatographic behavior

### The Minimal Loop
```
┌──────────────────┐
│   PRIDE Archive  │
│   (PXD019086)    │
└────────┬─────────┘
         │ download .d files
         ▼
┌──────────────────┐
│    FragPipe      │────────────► PSMs (identifications)
└────────┬─────────┘                      │
         │ .d files                       │
         ▼                                │
┌──────────────────┐                      │
│  imspy/rustdf    │────────────► Raw extraction:
│  (raw extraction)│              - ion distributions
└────────┬─────────┘              - chromatograms
         │                        - mobilograms
         │                                │
         ▼                                ▼
┌──────────────────────────────────────────┐
│              Merge & Align               │
│  (PSM identities + raw signal features)  │
└────────────────────┬─────────────────────┘
                     │
                     ▼
┌──────────────────────────────────────────┐
│           imspy-predictors               │
│            train / eval                  │
└──────────────────────────────────────────┘
```

### Success Criteria (POC)
- [ ] Automated download of raw data from PRIDE
- [ ] Singularity containers for processing pipeline
- [ ] FragPipe processing of timsTOF .d files (user-provided FragPipe)
- [ ] imspy raw extraction (ion distributions, chromatograms, mobilograms)
- [ ] Merge PSMs with raw features
- [ ] Train/test split with held-out evaluation
- [ ] CCS prediction achieving R² > 0.9 (target: match Meier et al. R > 0.99)
- [ ] Documentation for users to set up their own FragPipe

---

## 3. Architecture Decisions

### 3.1 Orchestration Framework

**Decision: Snakemake**

Rationale:
- Python-native, integrates naturally with rustims/imspy stack
- Readable rules, easy to maintain
- Strong in genomics/proteomics communities
- Good Docker/Singularity support
- Scales to HPC via cluster profiles

Key features we'll use:
- `container:` directive for Docker integration
- `rule` definitions for each pipeline step
- Dependency-based execution (only re-run what's needed)
- Config files for dataset parameters

### 3.2 HPC Environment (Mogon2)

**Cluster details:**
- Location: Johannes Gutenberg University Mainz
- Scheduler: SLURM
- Container runtime: Singularity (not Docker)
- GPU nodes: Available for training

**Snakemake SLURM integration:**
```bash
# Run with SLURM executor
snakemake --slurm --jobs 100

# Or use a cluster profile
snakemake --profile mogon2/
```

**Cluster profile structure:**
```
profiles/mogon2/
├── config.yaml          # SLURM defaults (partition, time, memory)
└── slurm-status.py      # Job status script (optional)
```

### 3.3 External Dependencies & Licensing

**FragPipe/MSFragger** requires academic license agreement - cannot be redistributed.

Users must:
1. Download FragPipe from https://fragpipe.nesvilab.org/
2. Accept academic license terms
3. Provide path via config

**Config approach:**
```yaml
# config/config.yaml
fragpipe:
  path: "/path/to/fragpipe"      # User provides this
  workflow: "LFQ-MBR"            # Or timsTOF-specific workflow
  threads: 16
```

**Snakefile injection:**
```python
# FragPipe path from config, bound into container
rule fragpipe:
    input: "data/raw/{accession}"
    output: "data/processed/{accession}/psm.tsv"
    params:
        fragpipe_path=config["fragpipe"]["path"]
    singularity: "containers/fragpipe-base.sif"
    shell:
        """
        singularity exec --bind {params.fragpipe_path}:/opt/fragpipe \
            {input.container} /opt/fragpipe/bin/fragpipe ...
        """
```

### 3.4 Containerization Strategy (Singularity)

```
claudius-proteomics/
├── containers/
│   ├── fragpipe-base/
│   │   └── Singularity.def     # Java + dependencies (NO FragPipe binaries)
│   ├── download/
│   │   └── Singularity.def     # PRIDE API tools
│   └── imspy/
│       └── Singularity.def     # rustims + imspy + imspy-predictors
```

**Build workflow:**
```bash
# Build locally (requires sudo) or use remote builder
singularity build fragpipe-base.sif containers/fragpipe-base/Singularity.def

# Or convert from Docker Hub
singularity pull docker://ghcr.io/claudius-proteomics/fragpipe-base:latest
```

**Key points:**
- `fragpipe-base.sif` contains Java + dependencies, NOT FragPipe itself
- FragPipe binaries mounted at runtime via `--bind`
- Users responsible for obtaining their own FragPipe installation
- imspy container: Python 3.8+, Rust toolchain, TensorFlow
- GPU training: use `--nv` flag with Singularity for CUDA passthrough

### 3.5 Snakemake Workflow Structure

```
claudius-proteomics/
├── Snakefile                   # Main workflow
├── config/
│   └── config.yaml             # Dataset parameters (accessions, paths)
├── profiles/
│   └── mogon2/
│       └── config.yaml         # SLURM cluster profile
├── rules/
│   ├── download.smk            # PRIDE download rules
│   ├── fragpipe.smk            # FragPipe processing rules
│   ├── extract.smk             # Data extraction/cleaning rules
│   └── train.smk               # Model training rules
├── containers/
│   ├── fragpipe/Singularity.def
│   ├── download/Singularity.def
│   └── imspy/Singularity.def
├── scripts/
│   ├── download_pride.py
│   ├── run_fragpipe.sh
│   ├── extract_ccs.py
│   └── train_ccs.py
└── envs/
    └── imspy.yaml              # Conda environment (alternative to Singularity)
```

#### Core Rules
```python
# Simplified Snakefile structure
rule all:
    input: "models/ccs_v1/metrics.json"

rule download:
    output: directory("data/raw/{accession}")
    singularity: "containers/download.sif"
    shell: "python scripts/download_pride.py {wildcards.accession}"

rule fragpipe:
    input: "data/raw/{accession}"
    output: "data/processed/{accession}/psm.tsv"
    params:
        fragpipe_path=config["fragpipe"]["path"]
    singularity: "containers/fragpipe-base.sif"
    resources:
        mem_mb=64000,
        time="4:00:00"
    shell: "bash scripts/run_fragpipe.sh {params.fragpipe_path} {input} {output}"

rule extract_raw:
    input: "data/raw/{accession}"
    output: "data/extracted/{accession}/raw_features.parquet"
    singularity: "containers/imspy.sif"
    resources:
        mem_mb=32000,
        time="2:00:00"
    script: "scripts/extract_raw_features.py"
    # Extracts: ion distributions, chromatograms, mobilograms

rule merge:
    input:
        psm="data/processed/{accession}/psm.tsv",
        raw="data/extracted/{accession}/raw_features.parquet"
    output: "data/merged/{accession}/peptides_full.parquet"
    singularity: "containers/imspy.sif"
    script: "scripts/merge_psm_raw.py"

rule train:
    input: "data/merged/{accession}/peptides_full.parquet"
    output: "models/ccs_v1/metrics.json"
    singularity: "containers/imspy.sif"
    resources:
        mem_mb=32000,
        gpu=1,
        time="2:00:00"
    script: "scripts/train_ccs.py"
```

### 3.6 Data Flow

```
data/
├── raw/                        # Downloaded .d files (large, can delete after extraction)
│   └── PXD019086/
│       └── *.d/
│
├── processed/                  # FragPipe output (PSM identifications)
│   └── PXD019086/
│       ├── psm.tsv
│       ├── ion.tsv
│       └── ...
│
├── extracted/                  # imspy raw extraction
│   └── PXD019086/
│       └── raw_features.parquet
│           # - ion distributions
│           # - chromatographic peaks
│           # - mobilograms
│           # - isotope patterns
│
├── merged/                     # PSMs + raw features aligned
│   └── PXD019086/
│       └── peptides_full.parquet
│
└── datasets/                   # Combined training sets
    └── ccs_training_v1.parquet

models/                         # Trained models
└── ccs_v1/
    ├── model.h5                # TensorFlow model
    └── metrics.json
```

---

## 4. POC Implementation Plan

### Two-Phase Architecture

The project is split into two decoupled phases:

```
┌─────────────────────────────────────────────────────────────┐
│                    PHASE 1: DATABASE                        │
│                 (Continuous, Extensible)                    │
├─────────────────────────────────────────────────────────────┤
│  PRIDE → FragPipe → imspy extraction → Peptide Database    │
│                                                             │
│  - Can add new datasets incrementally                       │
│  - Versioned snapshots                                      │
│  - The "resource" itself                                    │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                   PHASE 2: MODELING                         │
│                    (On-demand)                              │
├─────────────────────────────────────────────────────────────┤
│  Database snapshot → Train/Val/Test split → Model training │
│                                                             │
│  - Multiple models from same data                           │
│  - Reproducible experiments                                 │
│  - Separate versioning                                      │
└─────────────────────────────────────────────────────────────┘
```

---

### Phase 1: Database Building Pipeline

#### 1.1 Infrastructure Setup
1. Set up Snakemake workflow for data ingestion
2. Create Singularity containers (FragPipe-base, download, imspy)
3. Create Mogon2 SLURM cluster profile
4. Design database schema (parquet-based, partitioned by accession)

#### 1.2 Per-Dataset Pipeline
For each PRIDE accession:
1. Download raw timsTOF .d files
2. Run FragPipe → PSM identifications
3. Run imspy raw extraction (XICs, mobilograms, isotopes)
4. Merge PSMs + raw features
5. Validate and insert into database
6. Update database metadata/index

#### 1.3 Database Structure
```
database/
├── metadata.json           # Database version, stats, index
├── peptides/               # Partitioned by accession
│   ├── PXD019086/
│   │   ├── peptides.parquet
│   │   └── manifest.json   # Processing info, QC metrics
│   ├── PXD010012/
│   │   └── ...
│   └── ...
├── snapshots/              # Versioned exports for training
│   ├── v1.0/
│   │   ├── training_set.parquet
│   │   └── manifest.json
│   └── ...
└── schema.json             # Column definitions, types
```

#### 1.4 Database API
- `add_dataset(accession)` - Process and add a new dataset
- `list_datasets()` - Show all datasets in database
- `export_snapshot(version)` - Create training-ready export
- `get_stats()` - Database statistics

---

### Phase 2: Model Training Pipeline

#### 2.1 Training Infrastructure
1. Load database snapshot
2. Apply train/val/test split (by peptide, not PSM)
3. Configure model architecture
4. Train with early stopping
5. Evaluate and save metrics

#### 2.2 Model Versioning
```
models/
├── ccs/
│   ├── v1.0/
│   │   ├── model.h5
│   │   ├── config.yaml
│   │   ├── metrics.json
│   │   └── training_data_version: "database/snapshots/v1.0"
│   └── v1.1/
│       └── ...
└── rt/                     # Future: retention time
    └── ...
```

---

### POC Milestones

**Milestone 1: Database MVP**
- [ ] Process PXD019086 (Meier et al.)
- [ ] Store in database format
- [ ] Export first snapshot

**Milestone 2: Model MVP**
- [ ] Train CCS model on snapshot
- [ ] Achieve R² > 0.9
- [ ] Document results

**Milestone 3: Validation**
- [ ] Compare to Meier et al. (target: R > 0.99, <1.4% error)
- [ ] Add second dataset to database
- [ ] Retrain and compare

---

## 5. Open Questions

1. ~~**PRIDE accession number** for Meier et al. CCS dataset?~~ → **PXD019086** (resolved)
2. ~~**Orchestration choice:**~~ → **Snakemake** (resolved)
3. ~~**Compute environment:**~~ → **Mogon2 HPC** (JGU Mainz) (resolved)
4. ~~**FragPipe version:**~~ → **v24 (latest)**, user-provided (resolved)
5. ~~**Model framework:**~~ → **imspy-predictors** (resolved)

**All major decisions resolved.**

---

## 6. Future Scaling (Post-POC)

### Immediate Extensions (raw feature predictions)
Unique to this project - enabled by imspy raw extraction:
- **Peak shape prediction** - chromatographic elution profiles
- **Mobilogram shape prediction** - ion mobility distributions
- **Isotope pattern prediction** - full isotope envelopes
- **Ion arrival time distribution** - beyond single CCS values

### Standard Extensions
- Add retention time prediction (any LC-MS data)
- Add MS2 fragmentation prediction (Prosit-style)
- Charge state distribution

### Scaling
- Ingest multiple PRIDE studies automatically
- Build dataset versioning and provenance tracking
- Publish open dataset + models
- Koina integration for online prediction serving

---

## References

### Primary Data Source
- Meier F, Köhler ND, Brunner AD, Wanka JMH, Voytik E, Strauss MT, Theis FJ, Mann M. (2021). Deep learning the collisional cross sections of the peptide universe from a million experimental values. *Nature Communications* 12, 1185. https://doi.org/10.1038/s41467-021-21352-8

### PRIDE Datasets
- PXD019086: https://www.ebi.ac.uk/pride/archive/projects/PXD019086
- PXD010012: https://www.ebi.ac.uk/pride/archive/projects/PXD010012
- PXD017703: https://www.ebi.ac.uk/pride/archive/projects/PXD017703

### Code (Meier et al.)
- MannLabs/DeepCollisionalCrossSection: https://github.com/MannLabs/DeepCollisionalCrossSection
- theislab/DeepCollisionalCrossSection: https://github.com/theislab/DeepCollisionalCrossSection

### Our Stack
- rustims: https://github.com/theGreatHerrLebert/rustims
- imspy-predictors: https://github.com/theGreatHerrLebert/rustims/tree/feature/koina-online/packages/imspy-predictors
- ionmob: https://github.com/theGreatHerrLebert/ionmob

### Tools
- FragPipe: https://fragpipe.nesvilab.org/ (academic license required - not redistributable)
- PRIDE Archive: https://www.ebi.ac.uk/pride/
- Snakemake: https://snakemake.readthedocs.io/
- Singularity: https://sylabs.io/singularity/
