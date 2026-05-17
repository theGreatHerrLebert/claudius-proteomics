# Lab-Drift Latent Model

**Status:** draft / design notes
**Created:** 2026-05-17
**Scope:** Phase 2 modeling — predicting how raw measurements drift as a function
of lab/instrument/acquisition metadata.

---

## 1. Purpose

San José collects *peptide observations in experimental context*. The first
scientifically interesting result is a **coordinate system of lab drift**:
quantify how technical factors (lab, timsTOF model, individual instrument, TIMS
calibration, LC column, gradient, acquisition date) shift measured properties —
**retention time (RT)**, **ion mobility / CCS (1/K0)**, and **fragment
intensity** — when the biology is held constant.

The data substrate is the **human timsTOF DDA/PASEF corpus** measured across
many labs and devices. Because RT and CCS are properties of the peptide *ion*
rather than the sample (see §1.1), any peptide that recurs across runs is an
anchor: the spread in its measured values *is* the drift signal, and therefore
the supervision target for a drift model. HeLa — the universal proteomics QC
standard and Bruker's shipped timsTOF QC sample — is the densest such anchor
block, but not a requirement.

The end goal is a **drift-conditioned predictor**: a base sequence→property
model plus a learned drift embedding that conditions it, analogous to how
collision energy (CE) is fed as a conditioning input into intensity models.

### 1.1 What is held constant — and what is not

The invariant we anchor on is **peptide identity**, not the sample. RT and CCS
are physicochemical properties of the peptide *ion*: a given
`modified_sequence + charge` has the same CCS (gas-phase, sample-independent)
and the same intrinsic RT whether it came from HeLa, HEK293, or a tissue biopsy.
So whenever a peptide recurs across runs, any difference in its measured RT/CCS
*is* drift — and the sample of origin is irrelevant.

Consequences:

- The supervision substrate for RT/CCS drift is the **whole human DDA/PASEF
  corpus** (~573 timsTOF datasets), not a HeLa-only subset. Anchor peptides
  emerge *after* search, from peptide recurrence across runs — not from
  dataset-level sample labels. In principle this extends cross-organism: a
  conserved tryptic peptide shared by human and mouse has the same CCS.
- HeLa is therefore a **convenience, not a requirement** — the densest, cleanest
  shared-anchor block (deep standard proteome; commercial digest standards as
  the purest tier).
- The exception is **fragment intensity**, which *is* abundance-dependent and
  therefore sample-dependent. Same-sample (HeLa) matters much more there.

**Tissue/sample type is a covariate, not a constraint.** It does not change a
peptide's CCS or RT, but it changes (a) which peptides are observed — the
anchor-overlap structure — and (b) matrix complexity (co-elution, dynamic range,
ionization suppression), a second-order effect on RT and a first-order effect on
intensity. Left unmodelled, matrix variance leaks into the drift embedding — a
single-matrix lab acquires a spurious "drift" signature. Cluster by
**matrix-complexity tier** — cell line → primary cells → solid tissue → body
fluid (plasma the extreme outlier) — for which tissue is the readable proxy.
Housekeeping / core-proteome peptides, observed in every tissue, are the
**bridge anchors** that keep the global coordinate system identifiable across
tissue blocks.

---

## 2. Core reframe: three properties, three drift geometries

RT, CCS, and intensity are **not one problem**. Their drift geometry differs,
and that should drive model design and effort allocation.

| Property | Transferability | Drift geometry | Drift dimensionality | Difficulty |
|----------|----------------|----------------|---------------------|------------|
| **CCS / 1/K0** | High — physical observable | Near-affine in 1/K0 (TIMS calibration: tune-mix ions, gas T/P) | ~1–2 params/run | Low — cheap win |
| **RT** | Low, but structured | Monotonic warp (gradient, column, flow, temp, column age) | Smooth high-dim curve, monotone-constrained | Medium |
| **Fragment intensity** | Low, highest-dim | Per-(ion type, position, charge) modulation; dominated by CE | High-dim | High |

**Implications:**

- **CCS** is the cheap, high-certainty win. The existing CCS-reproducibility
  analysis already captured it with a per-file scalar offset (median CV 1.31% →
  0.67%). Cross-lab the raw spread will be larger and the gain bigger.
- **RT** should be modeled as a *warp*, not an offset: predict an intrinsic
  coordinate (iRT / Chronologer space), then fit a per-run monotone spline
  iRT→RT.
- **Intensity** is the frontier. CE is the dominant knob — but on timsTOF CE is
  *ramped* as a function of mobility, so it is not even a scalar per run.

Stage the effort: **CCS → RT → intensity**, increasing difficulty, decreasing
certainty.

---

## 3. Modeling approaches

Ordered by increasing ambition. The intended build path is residual + FiLM on
known metadata, then a latent `z` for the residual-of-the-residual.

### 3.1 Mixed-effects / ComBat baseline
Empirical-Bayes batch model. Fixed effects = sequence-driven prediction; random
effects = lab/instrument/run nested hierarchy. Worth running first as a
yardstick — partial pooling handles the ragged design for free.

### 3.2 Residual modeling *(non-negotiable)*
Do not predict the absolute property with the drift model. Predict the
**residual** of a strong base model (Chronologer for RT, ionmob / imspy
predictor for CCS). Base model = biology; drift model = technical. Far less
data-hungry, and it cleanly factorizes biological from technical variance.

### 3.3 Explicit covariate conditioning (FiLM)
Embed categoricals (instrument model, column, lab, tissue / matrix-complexity
tier), scale continuous covariates (gradient length, temperature), produce
feature-wise (γ, β) that modulate the base model's hidden state. This is the CE-input idea, generalized. Interpretable;
fails on unseen categories.

### 3.4 Latent run-embedding
Each run / instrument / lab gets a learned embedding vector (scVI-style batch
embedding). Captures *unrecorded* factors. Problem: a new run at inference time
has no embedding.

### 3.5 VAE / amortized latent drift
An encoder maps a run's observed **anchor peptides** (consensus-identified,
known property values) → a posterior over a latent drift code `z`. The
predictor is conditioned on `z`. This generalizes at inference time: new run →
measure known peptides → infer `z` → predict everything else.

> The existing per-file 1/K0 offset (inferred from reference peptides seen in
> ≥40% of files) is the **1-D linear special case** of this. The latent drift
> code is that, generalized: higher-dimensional, nonlinear, and covering RT and
> intensity as well as CCS.

---

## 4. Prior art

- **RT:** DeepLC — most directly relevant; has an explicit *calibration mode*
  that recalibrates to a new run from a handful of peptides. Also Chronologer,
  Prosit-RT, AlphaPeptDeep.
- **CCS:** ionmob (explicitly models CCS; PXD026463 is ionmob training data),
  AlphaPeptDeep, the Meier et al. 2021 model (benchmark target, R > 0.99).
- **Intensity:** Prosit, MS²PIP, AlphaPeptDeep, the local
  DeepPeptideIntensityPredictor — all take CE as conditioning input.
- **Architecture template:** scVI / scANVI (single-cell — latent + batch
  covariate is structurally identical to this problem).
- **Classical precedent:** chemometrics "calibration transfer" (piecewise direct
  standardization, PDS) — the same problem, solved for NIR spectroscopy in the
  1990s.

---

## 5. Missing data — three distinct kinds

1. **Missing metadata.** PRIDE metadata is sparse (datasets routinely resolve to
   `organism: unknown`, completeness 0). Do not impute-and-pretend. Give each
   categorical a learnable "unknown" token; let the latent `z` absorb the rest —
   that is the latent's job. Some metadata is recoverable *from the data*:
   gradient length from RT span, instrument family from mobility/resolution
   fingerprints.
2. **Missing observations (MNAR).** Peptides are missing *because* of their
   properties. Mild for RT/CCS (observed → measured). **Severe for intensity** —
   only intense fragments are observed, so training intensity drift on observed
   peaks alone is selection-biased. Needs a censored/truncated likelihood or an
   explicit detection model.
3. **Ragged experimental design.** instrument × lab cells are sparse and
   unbalanced. Pure per-cell estimates fail; hierarchical partial pooling is
   mandatory.

---

## 6. Pitfalls — what to watch out for

- **Confounding / design rank.** Lab, instrument, date, and sample prep are
  correlated. If instrument X only ever appears in lab A, their effects are
  *not separable*. Check the design-matrix overlap before modeling: you need
  instruments that cross labs and labs that cross instruments. **This feeds back
  into dataset selection — prefer datasets that create crossings.**
- **Gauge / reference frame.** "Drift" is undefined without a reference. The
  model has a gauge freedom; pin it — a canonical reference run, a constraint
  that mean drift = 0, or anchoring to trusted physical CCS values. Get this
  wrong and the embedding is uninterpretable.
- **Matrix / tissue confound.** Tissue does not change a peptide's RT/CCS, but
  matrix complexity (dynamic range, co-elution, suppression) shifts apparent
  intensity and — mildly — RT. Model tissue / matrix-tier as a covariate, or a
  single-matrix lab acquires a spurious "drift" signature (see §1.1).
- **Hold the search side constant.** Search engine / FASTA / FDR / modification
  settings change the peptide set and its reported values — anchor on the
  triple-engine consensus subset so engine choice is not itself a drift axis.
  Commercial HeLa digest standards (Pierce/Promega) remain the cleanest tier
  when same-sample purity is needed (intensity).
- **Evaluation must be leave-one-*group*-out.** Held-out peptides will fool you.
  Use leave-one-instrument-out and leave-one-lab-out. Cold-start on an unseen
  instrument model (e.g. timsTOF Ultra) is the real test — and the case where
  latent-`z` beats explicit covariates.
- **Overfitting run identity.** Free per-run embeddings can memorize run ID and
  learn nothing generalizable. Use a small `z`, KL pressure, anchor-based
  inference.
- **CE is dirty metadata.** Reported inconsistently across labs and
  mobility-ramped on timsTOF. Treat as noisy.
- **Intra-batch / temporal drift.** Instruments drift *within* a batch (cf. the
  22-Prot vs 23-Prot ±0.02 offsets in PXD069059). Acquisition time is a
  covariate, not just run identity.

---

## 7. Expected gains

- **CCS cross-lab alignment** — the certain win. Affine per-run correction makes
  the corpus poolable; benchmark against Meier R > 0.99.
- **Poolability itself** — once everything aligns to one reference frame, the
  whole San José corpus becomes a single training set instead of N incompatible
  ones. This multiplies effective training data for *every* downstream model —
  the core strategic payoff.
- **Tighter identification windows** — drift-aware RT/CCS prediction narrows
  search tolerances → more IDs, fewer false matches, especially library-free
  DIA.
- **A drift atlas** — the `z`-space becomes a QC and instrument-fingerprinting
  artifact: flag outlier runs, characterize labs and instruments.
- **Uncertainty-aware prediction** — a posterior over `z` yields run-specific
  error bars on every prediction.
- **Intensity** — high-risk frontier: CE- and instrument-aware fragment
  intensity enables cross-instrument spectral libraries.

---

## 8. Proposed roadmap

1. **v0.5 baseline.** On the human DDA/PASEF corpus (HeLa-dense subset first),
   restrict to consensus-identified peptides seen across many runs. Fit (a) a per-run affine CCS correction and
   (b) a per-run monotone RT warp. Report variance explained **and the
   design-matrix overlap**. This establishes the reference frame and bounds the
   headroom a learned model can claim.
2. **Residual + FiLM.** Add explicit metadata conditioning on top of a base
   model (Chronologer / ionmob), predicting residuals.
3. **Latent `z`.** Introduce an amortized latent drift code inferred from anchor
   peptides; evaluate cold-start on held-out instruments/labs.
4. **Extend to intensity** with a censored likelihood for the MNAR problem.

---

## 9. Open questions

- Reference-frame choice: canonical run vs zero-mean constraint vs physical-CCS
  anchor — which gives the most stable, interpretable embedding?
- How much metadata can be reliably recovered from the data itself, and should
  that be a separate model feeding this one?
- Does the design-matrix overlap from the available PRIDE human DDA/PASEF
  datasets actually support separating lab from instrument effects? (Drives
  dataset selection.)
- How granular should the matrix-complexity tiering be, and how reliably can it
  be derived from PRIDE's noisy `tissue` / `organismsPart` / `sampleAttributes`
  fields?
- Single multi-task model over RT/CCS/intensity with a shared `z`, or separate
  models per property with property-specific `z`?
