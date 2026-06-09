# Plan v2 — Recalibrate timsim's gradient→peak-width prior on the corpus

> v2 (post-Codex review, `05-gradient-sigma-prior.codex-review.md`). Pivots from
> "fit Gaussian/EMG σ then feed timsim" to **calibrate against the model-free
> observable FWHM via simulation-based quantile matching**; uses a cluster-robust
> fit + dataset-LODO; fixes the wiggle (hierarchical decomposition, no 20%
> double-count); flags the λ-interaction risk with per-metric acceptance.

## Why (evidence)

timsim configures each peptide's RT peak as an EMG; the Gaussian width σ is set by
`calculate_rt_defaults(gradient_length)`:
`sigma_middle = G/3600 * 0.75 + 1.125` (linear), bounds ±25%, then σ is **sampled**
from a Beta in `[lower, upper]` (the random "wiggle"). λ (tail) similar.

Empirical study on 102 datasets synced from MOGON (`scripts/validate_gradient_law.py`,
`study_width_properties.py`, `study_width_replication.py`):
- σ scales **sub-linearly** with gradient: corr(G, σ)=0.57, log-log slope≈**0.49 (~√G)**,
  *not* linear. timsim's current linear default is mis-shaped.
- Per-peptide σ replicates only **modestly**: cross-run η²≈0.36/0.25 vs df-inflation
  baseline 0.18/0.08 → genuine peptide-explained variance ≈ **18–22%**; within-peptide
  cross-run **CV≈0.31** (the wiggle). Peptide **length does not** predict σ.
- ⇒ **σ ≈ scale(G) [dominant] + ~20% non-length peptide term + random wiggle.**

So the high-value, cheap lever is the **gradient→σ scale**, not a sequence→σ head
(parked as secondary). This plan recalibrates the scale (and the wiggle width) on the
corpus and feeds it back into timsim.

## Goal

Replace `calculate_rt_defaults`'s linear `sigma_middle(G)` with an empirically-fit
sub-linear function, and set the sampling spread from the observed population spread.
Keep the Beta sampling (the wiggle). Per-peptide and λ modeling are out of scope.

## Method (v2 — revised after Codex review)

Pivot: calibrate against a **model-free observable width** and match it by simulation,
rather than trusting Gaussian/EMG fit parameters as ground truth.

1. **Target = empirical FWHM, not σ.** Use `ms1_rt_fwhm` (half-max crossing width,
   computed model-free in the extractor — *not* a Gaussian/EMG fit), per LC run, all
   102 datasets, masked by a fit-independent quality gate. This removes the entire
   Gaussian→EMG correction problem (former step 3 — dropped).
2. **Scale law via a proper model, not naive OLS.** Fit
   `log FWHM_run = α + b·log G_run + u_dataset + ε_run` (random dataset intercept;
   runs are not independent). Uncertainty via **dataset-cluster bootstrap**. Do all
   model selection (power-law vs monotone/piecewise spline) **inside dataset-clustered
   LODO**. `G = p99(rt_seconds)` is an endogenous proxy → `b` is attenuated; report it
   as predictive-for-this-proxy, and run **sensitivity over p95/p97.5/p99/max** + the
   feature-count / quality thresholds. Compare against dataset-median fit as robustness.
3. **Wiggle via hierarchical decomposition (the main statistical task).**
   `log width_ijk = f(G_j) + run_j + peptide_i + residual_ijk`. The per-peptide
   sampling width timsim needs = **cross-sectional variance among peptides conditional
   on a run = peptide effect + within-run residual** — NOT the cross-run CV (0.31),
   which contaminates with run instability + fit noise. Since we are **not** modeling
   peptides, the wiggle must **absorb the ~20% peptide heterogeneity** — so do **not**
   add the 20% term separately (double-count). Calibrate the Beta from conditional
   **log-scale quantiles**, not symmetric ±bounds. Estimate the decomposition where
   we have cross-run peptides (PXD061039, PXD040716).
4. **Feed back by simulation-based quantile matching.** Set timsim's σ-scale and Beta
   wiggle, simulate, measure *simulated* FWHM, and tune until the **simulated FWHM
   distribution matches the real one** (per dataset, dataset-clustered LODO). This
   jointly absorbs the σ↔λ interaction and avoids inverting uncertain EMG params.
   Rewrite `calculate_rt_defaults` with the matched scale + spread; keep the sampler;
   unit-test at several G (30/60/120/197 min).

## Validation / acceptance

- Dataset-clustered LODO: matched prior's **FWHM-distribution divergence** (e.g.
  Wasserstein/KS on real vs simulated) **lower** than the current linear default.
- **Per-metric** acceptance, reported separately: core width, total FWHM, tail/asymmetry,
  peak area/height. Claim improvement **only** on metrics actually improved (a σ-only
  change can leave — or worsen — tail/FWHM if λ is mismatched).
- Proxy-choice + threshold sensitivity reported.

## Risks / open questions

- **G = `p99(rt_seconds)` proxy** — endogenous (ID depth / elution), attenuates `b`,
  per-dataset bias. Mitigation: proxy sensitivity sweep; interpret `b` as proxy-predictive,
  not a physical exponent. Real programmed gradient unavailable.
- **λ interaction (Codex #4):** changing σ-scale with a **stale λ prior can worsen the
  joint EMG**. The simulation-based FWHM match (step 4) partly absorbs this, but if the
  per-metric tail/asymmetry check fails, λ must be co-calibrated — pulling more EMG
  refits (blobs) then becomes in-scope. Decision deferred to the tail-metric result.
- **Confounding**: one global G→width law averages over column/flow/instrument/sample;
  fine for a global prior, per-setup refinement later.
- **Cross-run peptide coverage** for the decomposition is thin outside PXD061039 — the
  wiggle estimate may be dataset-specific; sanity-check on ≥2 datasets.

## Out of scope (unchanged)

- Peptide σ-head (~20% ceiling, non-length features) — its variance is folded into the
  wiggle here, not modeled explicitly.
- Explicit λ modeling — unless the tail-metric acceptance (above) forces it in.
