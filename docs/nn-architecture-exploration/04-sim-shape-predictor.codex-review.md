Reading additional input from stdin...
OpenAI Codex v0.136.0
--------
workdir: /home/administrator/Documents/promotion/claudius-proteomics
model: gpt-5.5
provider: openai
approval: never
sandbox: read-only
reasoning effort: medium
reasoning summaries: none
session id: 019ea71b-6137-7111-9c5c-a5d9ea6df5ff
--------
user
Review this ML design doc as an independent ML/proteomics engineer. It plans extending a peptide property predictor (Depthcharge-based, multi-task: intensity/CCS/RT/charge) to emit per-peptide peak-SHAPE parameters that configure a timsTOF simulator (timsim): RT as an exponentially-modified Gaussian (sigma, lambda) and IM as a Gaussian (sigma). Labels come from already-extracted per-precursor profile fits (ms1_rt_sigma/skew/r2, ms1_im_sigma/r2) in precursor_store.parquet, joined to peptide IDs, R2-masked, with RT shape gradient-normalized.

Focus, concrete and specific only:
1) CORRECTNESS: Does training on per-observation precursor peak widths actually yield what timsim needs (a generative per-peptide shape)? Is regressing to the conditional mean of sigma/lambda the right target, or does it bias the simulator toward too-narrow peaks (Jensen / averaging-of-widths effect)?
2) The gradient-normalization claim rests on corr=0.973 across only 4 datasets (n=4 points). Is that statistically meaningful? What confounds gradient length here (instrument, peptide population, column)? Is normalize-and-rescale safe across gradient types (linear vs multi-step)?
3) EMG parameter identifiability: is lambda recoverable from noisy XICs, and is the moment->EMG conversion (from a Gaussian LSQ fit's sigma+skew) valid at all? Risk that ms1_rt_sigma is a Gaussian-fit width, not an EMG sigma.
4) R2 masking selection bias: masking to R2>=0.8 keeps high-abundance/clean peptides — does the predictor then generalize to the low-abundance peptides timsim must also simulate? 
5) The 4/99 dataset coverage: any way to de-risk before MOGON re-extraction? Which of the 5 open questions in the doc are answerable from the doc itself.
Cap at ~800 words.

<stdin>
# Design — Per-Peptide Peak-Shape Predictor to Configure timsim

## Goal

Extend the `peptide-property-ng` predictor to emit the **per-peptide peak-shape
parameters that configure the timsim simulator**, replacing timsim's current
non-peptide-specific (sampled/fixed) widths with predicted, sequence-aware
shapes.

This is **not** apex/center prediction (already covered: RT via Chronologer, IM
via the CCS √(m/z) head). It is the *width and asymmetry of the real ion signal
distribution* — the elution peak shape in RT and the mobility peak shape in IM.

### Key reframe (and a dead end ruled out)

The target is the **physical ion distribution per observation**, captured from
the 4D raw data, **not** the run-to-run scatter of the apex. We verified the
latter is degenerate: `aligned_rt` is Sage's per-peptide aligned value, so the
within-peptide std of `aligned_rt` is ≈0 (median and p90 = 0.000 on every
dataset). Apex jiggle is the wrong signal; peak FWHM/σ is the right one.

## timsim contract (verified in imspy-simulation/timsim)

What the simulator actually ingests per peptide/ion:

- **RT = Exponentially-Modified Gaussian (EMG)**. Params: `rt_sigma` (Gaussian
  width) + `rt_lambda` (exponential tail rate); `rt_mu` is derived downstream
  from the apex. **Currently Beta-sampled** (`sigma_*_rt`, `k_*_rt` priors) —
  not peptide-specific.
- **IM = Gaussian**. Params: `inv_mobility_gru_predictor` (μ, 1/K0) +
  `inv_mobility_gru_predictor_std` (σ). **σ currently fixed** (~2% relative, or
  scaled to target mean 0.009).

A predictor that configures one peptide must therefore emit: RT `(sigma,
lambda)` and IM `sigma` (apex/μ come from the existing apex predictors).

## Labels — already extracted, local, unblocked

`precursor_store.parquet` (per precursor) carries fitted profile shapes:

```
ms1_rt_apex, ms1_rt_sigma, ms1_rt_skew, ms1_rt_fwhm, ms1_rt_r2
ms1_im_apex, ms1_im_sigma, ms1_im_skew, ms1_im_fwhm, ms1_im_r2
```

Peptide identity per precursor comes from `precursor_index.parquet`
(`sage_peptide`/`fragpipe_*`/`diann_*` + q-values), joined on
`precursor_id` + `raw_file`. Filter to confident IDs (q ≤ 0.01) and good fits.

**Coverage constraint:** only **4 of 99 datasets** have `precursor_store`
(PXD019086 + 3 trypsin sets — the POC batch). Scaling labels to the rest
requires re-running the 4D extraction on the `.d` files, which lives on
MOGON-NHR (currently IP-blocked). The prototype runs on the 4 local datasets.

## Data-validated design decisions (measured on the 4 datasets)

1. **IM is a direct hit.** The IM/CCS head already outputs `(mean, std)`; today
   `std` is computed but unsupervised (loss is L1 on mean only). Add a **masked
   L1 on `std`** against `ms1_im_sigma`, masked by `ms1_im_r2 ≥ 0.8`. Gradient-
   independent. IM Gaussian R² median ≈ 0.74 (peaks genuinely Gaussian).

2. **RT needs reparameterization + is gradient-dependent.**
   - Extraction gives Gaussian `ms1_rt_sigma` + `ms1_rt_skew`; timsim wants EMG
     `(sigma, lambda)`. Prefer an **EMG refit** of the raw XIC (stored in
     `_extracted_poc/*/blobs.bin`) to get `(sigma, lambda)` + an EMG-R² mask;
     fallback is analytic moment→EMG conversion. RT Gaussian R² median ≈ 0.62 —
     lower than IM precisely because RT peaks are skewed (confirms EMG/lognormal
     over Gaussian).
   - **RT σ scales ~linearly with gradient length**: corr(gradient, median
     rt_σ) = 0.973 across the 4 datasets; normalizing σ by gradient drops the
     cross-dataset CV from 0.373 → 0.181. → **Predict gradient-normalized shape**
     (`σ/G`, tail `1/λ` likewise `/G`); timsim rescales by target gradient `G'`.
     `gradient_length` metadata field exists but is null → derive `G` per run
     from p99 of `rt_seconds`.
   - RT **apex stays Chronologer**; we add shape only.

3. **Masking is essential.** Many precursors have no usable distribution (poor
   fit). R²-gating: `R²≥0.8` keeps ~33% RT / ~39% IM (~1.6–1.9M precursors per
   dataset — ample); `R²≥0.9` keeps ~17–19%. Use **L1 on the variance/shape**
   terms (not Gaussian-NLL, which was prototype-unstable), each with its own
   R²-gated validity mask.

## Implementation sketch (`peptide-property-ng`)

- **IM head**: no structural change — add `im_sigma` target + `im_r2` mask in
  `sage_dataset.py`/collate; add masked L1 term in `losses.py`.
- **RT head** (`heads/scalar.py`): extend to emit `(sigma_norm, lambda_norm)`
  (shape only). Optionally condition on a gradient-length scalar, but the data
  favors normalize-and-rescale over conditioning.
- **Dataset**: new joined loader (`precursor_store` ⋈ `precursor_index`),
  q-filtered, R²-masked, gradient-normalized RT targets.
- **Losses**: `+ masked_L1(im_std, im_sigma)` `+ masked_L1(rt_shape, emg_params)`.
- **Prototype** on the 4 POC datasets; gate full corpus on MOGON access to
  extract the remaining ~95.

## Open questions for review

1. **EMG refit vs moment-conversion** for RT `(sigma, lambda)` — is refitting
   the XIC worth it, or is converting `ms1_rt_sigma`+`ms1_rt_skew` adequate?
   (Depends on whether the extraction's sigma/skew are LSQ-Gaussian-fit values
   or empirical moments — to confirm.)
2. **Normalize-and-rescale vs gradient-conditioning** for RT σ — normalization
   wins on cross-dataset CV, but does it lose peptide×gradient interactions a
   conditioned head could capture?
3. **Single shared σ per peptide vs per-charge/per-observation** — precursor_store
   is per-observation; do we regress to the conditional mean, or model the
   spread of shapes too?
4. **Only 4 datasets** — is a 4-dataset prototype enough to trust the head before
   investing in MOGON re-extraction of the other 95?
5. **EMG λ identifiability** — is `rt_lambda` well-constrained by the data, or
   should we predict only σ and keep λ sampled?
</stdin>
codex
## Findings

**1. The proposed labels do not yet establish a generative per-peptide shape model.**

Each row is a *per-observation* estimate combining peptide effects with charge, abundance, interference, run, chromatography, acquisition, and fitting error. A deterministic sequence model trained with L1 estimates the conditional median, not a peptide’s shape distribution. This is sufficient only if timsim wants a typical clean observation.

Using one mean/median width can distort aggregate simulated peaks. A mixture of Gaussians/EMGs with variable widths is not equivalent to one profile at the average width. The central region is generally sharper and the tails heavier than the corresponding single-width profile. This is a mixture-distribution issue, not simply Jensen implying “too narrow.”

Recommended target:

- Decompose log-width into peptide, charge, and run effects using repeated observations.
- Predict a distribution, e.g. `log σ ~ Normal(μ, τ²)`, rather than only σ.
- Sample observation-level shape in timsim.
- At minimum, validate deterministic predictions by comparing simulated versus real FWHM, central height, and tail quantiles, not parameter MAE alone.

A single shared value per peptide is premature because the current labels are precursor/observation-level.

**2. The gradient-normalization evidence is insufficient.**

`r=0.973` from four dataset-level medians is not persuasive. With `n=4`, there are only two degrees of freedom; the approximate two-sided p-value is about 0.027, but the confidence interval is extremely wide and one dataset can determine the result. More importantly, these are four datasets, not independent gradient perturbations under controlled conditions.

Gradient length is confounded with:

- LC method, column dimensions, particle size, flow rate, and temperature
- Instrument/acquisition method and extraction behavior
- Sample complexity and peptide composition
- Loading amount and saturation
- Definition of the effective separation window

Using `p99(rt_seconds)` is especially unsafe: it depends on identification density and peptide population, not just the programmed gradient.

Linear rescaling is physically plausible only for closely related linear gradients and otherwise matched chromatography. It should not be assumed for multi-step, nonlinear, shortened, or extended gradients. Local gradient slope around the peptide, ideally from the programmed LC method, is more relevant than total duration.

Treat normalization as a hypothesis. Include run/method effects and test leave-one-dataset-out extrapolation. Do not use `p99 RT` as production gradient metadata.

**3. The moment-to-EMG fallback is not valid without confirming label semantics.**

EMG moments can be inverted only when the supplied variance and skewness are empirical moments of an EMG-like distribution. A Gaussian least-squares fit’s `sigma` is the width of the fitted symmetric approximation; it is not the Gaussian component σ of an EMG. Likewise, a separately reported fit-skew parameter may not equal standardized third moment. Combining them can produce meaningless or unstable `(σ, λ)`.

Even direct EMG fitting may poorly identify λ when:

- The tail is truncated by extraction windows
- Sampling is sparse
- Baseline/interference is present
- Signal-to-noise is low
- The peak is nearly symmetric, where `λ → ∞`

σ and λ will also be strongly correlated. Fit in a stable parameterization such as log Gaussian width plus log tail timescale, impose bounds/priors, and retain uncertainty or identifiability diagnostics. Compare EMG against simpler Gaussian and possibly alternative asymmetric models on held-out XICs. Until then, predicting σ while sampling or hierarchically modeling λ is safer.

**4. R² masking creates target-population selection bias.**

The retained examples are likely abundant, isolated, well-sampled precursors. The model then estimates shape conditional on being cleanly measurable, while timsim must generate low-abundance and interfered observations too.

R² is also not a pure label-quality score: asymmetric real peaks can receive low Gaussian R². Thus the RT mask may preferentially remove exactly the tailing behavior the new head should learn.

Use:

- Soft reliability weighting or uncertainty-aware fitting rather than only hard masking
- Label-quality analyses versus intensity, charge, RT, m/z, sequence, and dataset
- Validation stratified by abundance and fit quality
- Repeated-observation aggregation to recover peptide shape from noisy individual fits
- Separate physical shape from simulated measurement noise/interference

**5. Four datasets can validate plumbing, not generalization.**

Before re-extraction, de-risk by:

- Inspecting stored XICs and confirming every existing field’s exact fitting definition
- Refitting a stratified subset with EMG and evaluating parameter recovery via synthetic XIC injection
- Measuring within-peptide versus between-peptide variance across runs and charges
- Comparing sequence-only, charge/run baselines, and peptide random-effects models
- Leave-one-dataset-out testing
- Testing whether predicted shapes improve timsim observables, not merely fit-label loss

## Open Questions Resolvable Now

From the document alone:

1. **EMG conversion:** answerable provisionally: conversion is not defensible until `sigma/skew` semantics are confirmed; refitting is required for a valid target.
2. **Normalization:** not answerable from four datasets; current evidence is inadequate.
3. **Shared versus observation-level:** partly answerable: shared deterministic σ is not justified without variance decomposition.
4. **Prototype sufficiency:** yes: sufficient for implementation and identifiability studies, insufficient for trusting deployment/generalization.
5. **λ identifiability:** not answerable from metadata or R²; requires profile-level fits, uncertainty/profile-likelihood analysis, and synthetic recovery experiments.
tokens used
5,019
## Findings

**1. The proposed labels do not yet establish a generative per-peptide shape model.**

Each row is a *per-observation* estimate combining peptide effects with charge, abundance, interference, run, chromatography, acquisition, and fitting error. A deterministic sequence model trained with L1 estimates the conditional median, not a peptide’s shape distribution. This is sufficient only if timsim wants a typical clean observation.

Using one mean/median width can distort aggregate simulated peaks. A mixture of Gaussians/EMGs with variable widths is not equivalent to one profile at the average width. The central region is generally sharper and the tails heavier than the corresponding single-width profile. This is a mixture-distribution issue, not simply Jensen implying “too narrow.”

Recommended target:

- Decompose log-width into peptide, charge, and run effects using repeated observations.
- Predict a distribution, e.g. `log σ ~ Normal(μ, τ²)`, rather than only σ.
- Sample observation-level shape in timsim.
- At minimum, validate deterministic predictions by comparing simulated versus real FWHM, central height, and tail quantiles, not parameter MAE alone.

A single shared value per peptide is premature because the current labels are precursor/observation-level.

**2. The gradient-normalization evidence is insufficient.**

`r=0.973` from four dataset-level medians is not persuasive. With `n=4`, there are only two degrees of freedom; the approximate two-sided p-value is about 0.027, but the confidence interval is extremely wide and one dataset can determine the result. More importantly, these are four datasets, not independent gradient perturbations under controlled conditions.

Gradient length is confounded with:

- LC method, column dimensions, particle size, flow rate, and temperature
- Instrument/acquisition method and extraction behavior
- Sample complexity and peptide composition
- Loading amount and saturation
- Definition of the effective separation window

Using `p99(rt_seconds)` is especially unsafe: it depends on identification density and peptide population, not just the programmed gradient.

Linear rescaling is physically plausible only for closely related linear gradients and otherwise matched chromatography. It should not be assumed for multi-step, nonlinear, shortened, or extended gradients. Local gradient slope around the peptide, ideally from the programmed LC method, is more relevant than total duration.

Treat normalization as a hypothesis. Include run/method effects and test leave-one-dataset-out extrapolation. Do not use `p99 RT` as production gradient metadata.

**3. The moment-to-EMG fallback is not valid without confirming label semantics.**

EMG moments can be inverted only when the supplied variance and skewness are empirical moments of an EMG-like distribution. A Gaussian least-squares fit’s `sigma` is the width of the fitted symmetric approximation; it is not the Gaussian component σ of an EMG. Likewise, a separately reported fit-skew parameter may not equal standardized third moment. Combining them can produce meaningless or unstable `(σ, λ)`.

Even direct EMG fitting may poorly identify λ when:

- The tail is truncated by extraction windows
- Sampling is sparse
- Baseline/interference is present
- Signal-to-noise is low
- The peak is nearly symmetric, where `λ → ∞`

σ and λ will also be strongly correlated. Fit in a stable parameterization such as log Gaussian width plus log tail timescale, impose bounds/priors, and retain uncertainty or identifiability diagnostics. Compare EMG against simpler Gaussian and possibly alternative asymmetric models on held-out XICs. Until then, predicting σ while sampling or hierarchically modeling λ is safer.

**4. R² masking creates target-population selection bias.**

The retained examples are likely abundant, isolated, well-sampled precursors. The model then estimates shape conditional on being cleanly measurable, while timsim must generate low-abundance and interfered observations too.

R² is also not a pure label-quality score: asymmetric real peaks can receive low Gaussian R². Thus the RT mask may preferentially remove exactly the tailing behavior the new head should learn.

Use:

- Soft reliability weighting or uncertainty-aware fitting rather than only hard masking
- Label-quality analyses versus intensity, charge, RT, m/z, sequence, and dataset
- Validation stratified by abundance and fit quality
- Repeated-observation aggregation to recover peptide shape from noisy individual fits
- Separate physical shape from simulated measurement noise/interference

**5. Four datasets can validate plumbing, not generalization.**

Before re-extraction, de-risk by:

- Inspecting stored XICs and confirming every existing field’s exact fitting definition
- Refitting a stratified subset with EMG and evaluating parameter recovery via synthetic XIC injection
- Measuring within-peptide versus between-peptide variance across runs and charges
- Comparing sequence-only, charge/run baselines, and peptide random-effects models
- Leave-one-dataset-out testing
- Testing whether predicted shapes improve timsim observables, not merely fit-label loss

## Open Questions Resolvable Now

From the document alone:

1. **EMG conversion:** answerable provisionally: conversion is not defensible until `sigma/skew` semantics are confirmed; refitting is required for a valid target.
2. **Normalization:** not answerable from four datasets; current evidence is inadequate.
3. **Shared versus observation-level:** partly answerable: shared deterministic σ is not justified without variance decomposition.
4. **Prototype sufficiency:** yes: sufficient for implementation and identifiability studies, insufficient for trusting deployment/generalization.
5. **λ identifiability:** not answerable from metadata or R²; requires profile-level fits, uncertainty/profile-likelihood analysis, and synthetic recovery experiments.
