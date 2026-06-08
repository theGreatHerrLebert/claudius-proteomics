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
session id: 019ea74b-2a49-7a22-990a-44bf4ce5bd05
--------
user
This is v2 of an ML design doc (revised after your first review) for a predictor that emits per-peptide peak-shape params to configure the timsTOF simulator timsim: RT as exponentially-modified Gaussian (sigma, lambda), IM as Gaussian (sigma). Since v1 we VERIFIED in the timsim source that the shape params deterministically generate the integrated signal (normalized frame/scan abundance = analytic EMG/Gaussian density binned over frames/scans; timsim does NOT sum stochastic observations and already injects its own variability + already scales sigma with gradient length). We also verified the label semantics in the extraction code: ms1_rt_sigma is a Gaussian curve_fit width, ms1_rt_skew is an empirical moment, ms1_rt_r2 is the Gaussian fit R²; raw XICs are stored in blobs.bin so an EMG refit is feasible.

You previously raised 4 points: (1) point vs distribution target, (2) weak n=4 gradient evidence, (3) EMG identifiability / no moment-conversion, (4) R2 masking selection bias. Assess SPECIFICALLY:
A) Given the verified integrated-signal mechanism, is the v2 decision to predict a SINGLE per-peptide shape (dropping the sample-a-distribution requirement) now correct, or is there a residual case where it still biases timsim output?
B) Is the validation plan (synthetic-XIC lambda recovery, within/between-peptide variance decomposition, LODO gradient test, end-to-end timsim-observable comparison) sufficient, or what's missing?
C) Any NEW risks introduced by the v2 changes (EMG refit step, gradient-normalized sigma + timsim rescale double-counting since timsim ALSO has calculate_rt_defaults(gradient_length))? The double-application of a gradient factor is a specific worry — flag if real.
D) Anything in v2 that is now over-engineered or under-specified.
Concrete only, cap ~700 words.

<stdin>
# Design v2 — Per-Peptide Peak-Shape Predictor to Configure timsim

> v2 incorporates an independent Codex review + verification of the timsim
> source. Key change from v1: we predict a **single per-peptide analytic shape**
> (not a sampled distribution), because the timsim code shows the shape
> parameters *deterministically generate the integrated signal*. The EMG-refit,
> EMG-R² masking, and gradient-validation upgrades from the review are kept.

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

What the simulator ingests per peptide/ion, and **how it uses it**:

- **RT = Exponentially-Modified Gaussian (EMG)**. Params: `rt_sigma` (Gaussian
  width) + `rt_lambda` (tail rate); `rt_mu` derived from the apex
  (`estimate_mu_from_mode_emg`, `simulate_frame_distributions_emg.py:250`).
  These feed `calculate_frame_occurrences_emg_par` /
  `calculate_frame_abundances_emg_par` (lines 254, 265) → a **normalized
  `frame_abundance`** vector (line 306). **Currently Beta-sampled**
  (`sigma_*_rt`, `k_*_rt` priors).
- **IM = Gaussian**. Params: `inv_mobility_gru_predictor` (μ, 1/K0) +
  `inv_mobility_gru_predictor_std` (σ) → `calculate_scan_*_gaussian_par`
  (`simulate_scan_distributions_with_variance.py:42-61`) → normalized
  `scan_abundance`. **σ currently fixed** (~2% relative, scaled to 0.009).

**The parameters deterministically create the integrated signal.** One peptide →
one analytic integrated peak (the EMG/Gaussian density binned across
frames/scans, normalized to sum 1). timsim does **not** sum many stochastic
observations. Consequences:

1. The consumer wants **one shape per peptide**; predicting the *central*
   per-peptide shape is the correct target type (no mixture-of-widths bias —
   that critique assumes stochastic accumulation, which timsim does not do).
2. A per-precursor profile fit (`ms1_rt_sigma`, …) is a **direct measurement of
   exactly the integrated shape timsim consumes**.
3. Run-to-run variability is timsim's own job (`sample_sigma_lambda_emg` +
   `add_uniform_noise`, lines 239, 303). So predicting a per-peptide *spread* τ
   is **optional realism, off the critical path**.

A predictor that configures one peptide therefore emits: RT `(sigma, lambda)`
and IM `sigma` (apex/μ come from the existing apex predictors).

**timsim already scales σ with gradient length** —
`sigma_middle_rt = gradient_length/3600 * 0.75 + 1.125`
(`calculate_rt_defaults`, `simulate_frame_distributions_emg.py:23`). So the
σ∝gradient dependence is a *pre-existing simulator assumption*, not just our
correlation. Our predictor refines this gradient-only default into a
**sequence-aware** σ.

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

## Design decisions (data-validated on the 4 datasets; v2-adjusted)

1. **IM is a direct hit.** The IM/CCS head already outputs `(mean, std)`; today
   `std` is computed but unsupervised. Add a **masked L1 on `std`** against
   `ms1_im_sigma`, masked by `ms1_im_r2 ≥ 0.8`. Gradient-independent (1713 s and
   5017 s datasets have ~identical IM σ). IM Gaussian R² median ≈ 0.74 (peaks
   genuinely Gaussian, so `ms1_im_r2` is a valid mask here).

2. **RT needs an EMG refit and a corrected mask.**
   - **Refit, do not convert** — *confirmed from the extractor code*
     (`scripts/extract_precursors.py`): `ms1_rt_sigma` ← `fit_gaussian()["sigma"]`
     (line 606) is a **Gaussian `curve_fit` width**, while `ms1_rt_skew` ←
     `rt_moments.skewness` (line 586) is an **empirical standardized 3rd moment**
     — *different estimators*, so backing `(σ, λ)` out of them is doubly invalid.
     Refit the raw XIC (stored in `_extracted_poc/*/blobs.bin` as
     `ms1_rt_coords`/`ms1_rt_intensities`, confirmed via `sanjose/blob.py:86`)
     with an EMG in a **stable log-parameterization** (log σ, log tail-timescale
     1/λ), with bounds/priors.
   - **Mask on EMG-R², not `ms1_rt_r2`.** *Confirmed:* `ms1_rt_r2` ←
     `fit_gaussian()["r2"]` (`extract_precursors.py:607`) is a **Gaussian**
     goodness-of-fit, which is *low* precisely for tailing peaks — masking on it
     would discard the EMG tail behavior we want to learn (selection bias against
     the signal). Use the EMG-fit R² from the refit; prefer soft reliability
     weighting over a hard cut.
   - **λ identifiability is a real risk** (truncated tails, low SNR, near-
     symmetric → λ→∞; σ,λ strongly correlated). **Start by predicting σ only and
     leave λ to timsim's sampling**; promote λ to a predicted output only after
     synthetic-XIC recovery shows it is identifiable.
   - RT **apex stays Chronologer**; we add shape only.

3. **Gradient handling.** RT σ scales with gradient (timsim already assumes this;
   our 4-dataset corr=0.973 is *consistent but weak* — n=4, and the `p99(rt)`
   proxy tracks ID density, not the programmed gradient). Plan: predict a
   **gradient-normalized σ** and let timsim rescale (it already knows the
   gradient via `calculate_rt_defaults`). But **treat the linear law as a
   hypothesis**: source the real LC gradient from the programmed method (ideally
   *local slope* at the peptide's RT, not total duration), and **leave-one-
   dataset-out test** before trusting cross-gradient extrapolation. Do **not**
   ship `p99(rt_seconds)` as production gradient metadata.

4. **Single per-peptide shape, not a sampled distribution.** Justified by the
   timsim integrated-signal mechanism above. Regress the central shape (L1).
   Predicting a peptide-level spread τ is deferred (timsim injects variability
   itself).

## Implementation sketch (`peptide-property-ng`)

- **IM head**: no structural change — add `im_sigma` target + `im_r2` mask in
  `sage_dataset.py`/collate; add masked L1 term in `losses.py`.
- **RT head** (`heads/scalar.py`): extend to emit `sigma_norm` (and later
  `lambda_norm`), shape only.
- **Dataset**: new joined loader (`precursor_store` ⋈ `precursor_index`),
  q-filtered, EMG-R²-masked, gradient-normalized RT targets.
- **Losses**: `+ masked_L1(im_std, im_sigma)` `+ masked_L1(rt_sigma_norm, emg_sigma)`.
- **Prototype** on the 4 POC datasets; gate full corpus on MOGON access.

## Validation plan (before trusting the head / investing in re-extraction)

- ~~**Confirm label semantics** / `blobs.bin` holds XICs~~ — **DONE**:
  `ms1_rt_sigma`=Gaussian `curve_fit` width, `ms1_rt_skew`=empirical moment,
  `ms1_rt_r2`=Gaussian fit R² (`extract_precursors.py:586,606,607`); XICs are in
  `blobs.bin` (`ms1_rt_coords/intensities`, `sanjose/blob.py:86`). Refit is local-feasible.
- **EMG parameter recovery**: inject synthetic XICs with known `(σ, λ)`, refit,
  measure recovery and the σ–λ correlation / λ identifiability envelope.
- **Variance decomposition**: within-peptide vs between-peptide (and per-charge,
  per-run) variance of the shape labels — quantifies how much a sequence model
  can explain vs irreducible per-observation noise.
- **Baselines**: sequence-only vs charge/run-aware vs peptide random-effects.
- **Gradient law**: leave-one-dataset-out extrapolation of normalized σ.
- **End-to-end**: does feeding predicted shapes improve *timsim observables*
  (simulated vs real FWHM, central height, tail quantiles) — not just label MAE.

## Open questions (v2 status)

1. **EMG refit vs moment-conversion** — *resolved*: refit (conversion not
   defensible until `sigma/skew` semantics confirmed).
2. **Normalize-and-rescale vs gradient-conditioning** — lean normalize (timsim
   already rescales), but validate with real gradient metadata + LODO; n=4 is
   not yet conclusive.
3. **Shared vs per-observation shape** — *resolved*: single per-peptide shape
   (timsim consumes one analytic profile); spread τ deferred.
4. **4-dataset prototype** — sufficient for plumbing + identifiability studies,
   **not** for trusting generalization; full corpus needs MOGON re-extraction.
5. **λ identifiability** — *open, de-risked*: predict σ first, keep λ sampled;
   promote λ only after synthetic-recovery evidence.
</stdin>
codex
**A) Single-shape target**

Mostly correct. Since timsim renders one normalized analytic profile per simulated ion, predicting one conditional shape is the right primary interface. A distribution-valued model is not required.

Residual bias remains if “per-peptide” collapses systematic conditioning variables:

- Charge state, modifications, run/instrument settings, intensity/SNR, or local gradient slope may alter shape.
- The average of EMG profiles is not generally the EMG at the average parameters. If timsim’s existing variability is not calibrated conditionally around the prediction, population-level FWHM/tail distributions can still be wrong.
- Clarify whether the prediction is peptide-, precursor-, or peptide/run-specific. Your labels are precursor observations, while the proposed sequence-only output implies peptide-level invariance.

Thus: one point per simulated precursor is correct; blindly using one sequence-only constant across all contexts may not be.

**B) Validation**

The plan covers the main questions, but add:

1. **Real-XIC fit validation.** Synthetic recovery is necessary but insufficient. Manually inspect or stratify real fits by truncation, SNR, interference, number of frames, and tail coverage.
2. **Leakage-safe splits.** Hold out peptide sequences and raw files. Random precursor splits will leak repeated peptides and overstate sequence generalization.
3. **Default and simple baselines.** Compare against timsim defaults, gradient-only regression, peptide median, and charge-aware gradient regression. End-to-end improvement must exceed these.
4. **Fit uncertainty.** Obtain Hessian/bootstrap uncertainty for `(σ, λ)` and use uncertainty-aware weights. EMG R² alone does not detect parameter non-identifiability.
5. **Profile-level metric.** Compare binned normalized XICs directly, using likelihood/integrated error or Wasserstein distance. FWHM and tail quantiles can miss localized profile errors.
6. **Resolution/truncation tests.** Synthetic XICs must use actual frame spacing, extraction windows, censoring, baseline, noise, and saturation/interference regimes.
7. **Ablation.** Run end-to-end comparisons for default, predicted σ only, predicted σ plus sampled λ, and eventually predicted λ.

LODO with four datasets is a useful failure test, not evidence of a general gradient law.

**C) New risks**

The gradient double-application concern is real unless the interface is defined more precisely.

`calculate_rt_defaults(G)` appears to produce an **absolute default σ location**, not necessarily a dimensionless rescaling factor. Safe normalization is:

\[
\sigma_{\text{norm}}=\sigma_{\text{observed}} /
[f(G)/f(G_{\text{ref}})]
\]

and reconstruction is:

\[
\sigma_{\text{absolute}}=\hat{\sigma}_{\text{norm}}
[f(G)/f(G_{\text{ref}})]
\]

Apply this exactly once. If training labels are divided by a gradient factor and timsim subsequently applies both its default calculation and another scaling path, gradient dependence is doubled. Conversely, if the predicted value replaces the sampled absolute `rt_sigma`, timsim may not rescale it at all. Trace and document the precise injection point with unit tests at two gradient lengths.

Further EMG-refit risks:

- Baseline subtraction, truncated tails, overlapping peaks, and extraction-window censoring can bias λ.
- Joint location/shape fitting can couple apex error to σ and λ. Fit location freely, then verify that mode-recentering matches `estimate_mu_from_mode_emg`.
- R² is scale-dependent and weak for comparing tail quality. Reliability should incorporate parameter uncertainty and boundary hits.
- Soft weighting by fit quality can still create selection bias if fit quality correlates with peptide properties.

**D) Specification issues**

Under-specified:

- Units and exact definition of `sigma_norm`.
- Whether outputs condition on charge/modification/run.
- How predictions replace or modify `sample_sigma_lambda_emg`.
- Constraints and parameterization ensuring positive σ and numerically safe λ.
- EMG objective, baseline model, fitting window, minimum observations, boundary-hit policy, and uncertainty calculation.
- Whether IM labels are absolute or relative σ; timsim’s “2% scaled to 0.009” needs an explicit unit conversion contract.

Potentially over-engineered: peptide random-effects are useful as an explainability ceiling, but not required for the first plumbing prototype. Keep λ prediction deferred as proposed. The highest-priority correction is specifying and testing the single gradient-scaling path.
tokens used
6,213
**A) Single-shape target**

Mostly correct. Since timsim renders one normalized analytic profile per simulated ion, predicting one conditional shape is the right primary interface. A distribution-valued model is not required.

Residual bias remains if “per-peptide” collapses systematic conditioning variables:

- Charge state, modifications, run/instrument settings, intensity/SNR, or local gradient slope may alter shape.
- The average of EMG profiles is not generally the EMG at the average parameters. If timsim’s existing variability is not calibrated conditionally around the prediction, population-level FWHM/tail distributions can still be wrong.
- Clarify whether the prediction is peptide-, precursor-, or peptide/run-specific. Your labels are precursor observations, while the proposed sequence-only output implies peptide-level invariance.

Thus: one point per simulated precursor is correct; blindly using one sequence-only constant across all contexts may not be.

**B) Validation**

The plan covers the main questions, but add:

1. **Real-XIC fit validation.** Synthetic recovery is necessary but insufficient. Manually inspect or stratify real fits by truncation, SNR, interference, number of frames, and tail coverage.
2. **Leakage-safe splits.** Hold out peptide sequences and raw files. Random precursor splits will leak repeated peptides and overstate sequence generalization.
3. **Default and simple baselines.** Compare against timsim defaults, gradient-only regression, peptide median, and charge-aware gradient regression. End-to-end improvement must exceed these.
4. **Fit uncertainty.** Obtain Hessian/bootstrap uncertainty for `(σ, λ)` and use uncertainty-aware weights. EMG R² alone does not detect parameter non-identifiability.
5. **Profile-level metric.** Compare binned normalized XICs directly, using likelihood/integrated error or Wasserstein distance. FWHM and tail quantiles can miss localized profile errors.
6. **Resolution/truncation tests.** Synthetic XICs must use actual frame spacing, extraction windows, censoring, baseline, noise, and saturation/interference regimes.
7. **Ablation.** Run end-to-end comparisons for default, predicted σ only, predicted σ plus sampled λ, and eventually predicted λ.

LODO with four datasets is a useful failure test, not evidence of a general gradient law.

**C) New risks**

The gradient double-application concern is real unless the interface is defined more precisely.

`calculate_rt_defaults(G)` appears to produce an **absolute default σ location**, not necessarily a dimensionless rescaling factor. Safe normalization is:

\[
\sigma_{\text{norm}}=\sigma_{\text{observed}} /
[f(G)/f(G_{\text{ref}})]
\]

and reconstruction is:

\[
\sigma_{\text{absolute}}=\hat{\sigma}_{\text{norm}}
[f(G)/f(G_{\text{ref}})]
\]

Apply this exactly once. If training labels are divided by a gradient factor and timsim subsequently applies both its default calculation and another scaling path, gradient dependence is doubled. Conversely, if the predicted value replaces the sampled absolute `rt_sigma`, timsim may not rescale it at all. Trace and document the precise injection point with unit tests at two gradient lengths.

Further EMG-refit risks:

- Baseline subtraction, truncated tails, overlapping peaks, and extraction-window censoring can bias λ.
- Joint location/shape fitting can couple apex error to σ and λ. Fit location freely, then verify that mode-recentering matches `estimate_mu_from_mode_emg`.
- R² is scale-dependent and weak for comparing tail quality. Reliability should incorporate parameter uncertainty and boundary hits.
- Soft weighting by fit quality can still create selection bias if fit quality correlates with peptide properties.

**D) Specification issues**

Under-specified:

- Units and exact definition of `sigma_norm`.
- Whether outputs condition on charge/modification/run.
- How predictions replace or modify `sample_sigma_lambda_emg`.
- Constraints and parameterization ensuring positive σ and numerically safe λ.
- EMG objective, baseline model, fitting window, minimum observations, boundary-hit policy, and uncertainty calculation.
- Whether IM labels are absolute or relative σ; timsim’s “2% scaled to 0.009” needs an explicit unit conversion contract.

Potentially over-engineered: peptide random-effects are useful as an explainability ceiling, but not required for the first plumbing prototype. Keep λ prediction deferred as proposed. The highest-priority correction is specifying and testing the single gradient-scaling path.
