# RT & IM Peak-Shape Distributions — Empirical Findings

What governs the **shape** (width, asymmetry, modality) of chromatographic (RT) and
ion-mobility (IM) peaks in timsTOF DDA data, and what is predictable. Derived from
**102 datasets** (`raw_features` extracted on MOGON; model-free observables + EMG
refits on the local subset). Companion engineering docs:
`docs/nn-architecture-exploration/05-timsim-rt-width-patch.md`,
`05-timsim-im-width-note.md`, `05-gradient-sigma-prior.md` (plan, 2× Codex-reviewed).

## TL;DR

> **Peak *centers* are predictable; peak *widths* are not sequence-predictable.**
> RT apex → Chronologer, IM apex/CCS → √(m/z) physics — both good. But peak **width**
> is governed by the **acquisition setup**, not the peptide: RT width by the LC
> gradient (sub-linear, ~√G), IM width by **charge**. Per-peptide sequence explains
> only ~6–20% of width, with no usable features. **Don't build a sequence→width
> predictor.** Set timsim's width priors from the physical covariate + keep its
> random sampler.

---

## Retention-time (RT) peak width

**Width scales sub-linearly with gradient length** (model-free `ms1_rt_fwhm`, 102
datasets, dataset-clustered fit):

```
FWHM(G) ≈ 0.200 · G^0.424        (G = gradient length proxy = p99 retention_time, s)
                                  log-log r = 0.58; exponent stable 0.41–0.43 across p95→max
EMG core σ ≈ 0.441 · FWHM         (empirical ratio, paired EMG fits, n=111k)
```

- **It is NOT gradient-independent.** Short gradients (~11–47 min) have FWHM ~5 s;
  long (~74–197 min) ~8 s (1.6× tercile-median; **~3–5× at the true extremes**:
  ~2–4 s at ~12 min vs ~9–13 s at ~2–3 h). A constant width is unfaithful — it is
  **45% off at the short end**.
- **But gradient explains only the gross scale, not dataset-to-dataset scatter.**
  Out-of-sample (leave-one-dataset-out), the power law barely beats a constant in
  the *median* (0.232 vs 0.226) — because most datasets cluster mid-range where
  setup confounds dominate. Stratified, the gradient law clearly wins at the short
  extreme (0.314 vs 0.446) and ties mid/long. ⚠️ *A global-median LODO hid this —
  always stratify by the covariate.*
- **Per-peptide signal is low.** Cross-run replication (same peptide, different runs,
  per-run scale removed): genuine peptide-explained variance ≈ **18–22%**;
  within-peptide cross-run CV ≈ **0.31** (large run-to-run wiggle). The naive
  "within-run 62%" is a **shared-elution confound** (repeated PSMs of a peptide in
  one run share the same XIC).
- **Length does not predict width** (corr ≈ 0, inconsistent sign across datasets;
  partial corr controlling RT-position = 0.003). No usable sequence feature found.
- **Width reliability is intensity-dependent.** σ recovery is reliable above ~SNR 20
  (≤6% error) and degrades fast below SNR 10 (27% error at SNR≈4). A fit-quality (R²)
  mask alone is **insufficient** at low intensity — weight by intensity/SNR.
- **The EMG tail (λ) is unidentifiable** at timsTOF MS1 frame sampling (7–13-point
  XICs): predict/sample σ, leave λ to the simulator.

## Ion-mobility (IM) peak width

- **IM width is charge-driven, NOT sequence-determined.** Cross-run replication on
  the robust model-free `ms1_im_fwhm` (high-intensity): genuine peptide-explained
  ≈ 0.16 *pooled* but only **~0.06 within a charge state** — the pooled signal is the
  charge effect. High-intensity gating did not raise it. (The IM *apex*/CCS *is*
  structure-determined — only the width isn't.)
- **IM width is gradient-independent** (datasets at 1713 s vs 5017 s gradients have
  ~identical IM σ).
- **Per-charge width** (PXD061039, high-intensity, 1/K0 units):

  | charge | σ (≈FWHM/2.355) | relative σ/apex |
  |---|---|---|
  | 1 | 0.017 | 1.8% |
  | 2–4 | ~0.008 | ~0.8% |

  Multiply-charged cluster tightly; singly-charged ~2× wider.

- **Multi-modality is real but a modest, 3+-enriched minority.** Single-Gaussian fits
  are poor (im_r2 ~0.65; 76–80% <0.8), but mostly because **mobilograms are noisy**,
  not bimodal. A noise-robust, intensity-stratified mode test shows apparent
  bimodality **falls with SNR** (z3 53%→36% low→high). At high SNR ~30–36% is flagged,
  but **charge-1 is already ~33%** — a charge-independent false-positive floor
  (chimeric/co-isolated ions, shoulders). The *genuine* conformer excess for 3+ over
  2+ is small (~6 pts). Real, but not worth explicit bimodal modeling.

## What this means

| quantity | predictable? | by what | status |
|---|---|---|---|
| RT apex | ✅ | sequence | Chronologer (production) |
| IM apex / CCS | ✅ | sequence/structure | √(m/z) physics head |
| **RT width** | ✗ (sequence) | LC gradient (~√G) + noise | use gradient prior + sampler |
| **IM width** | ✗ (sequence) | charge + noise | use charge prior + sampler |
| IM multi-modality | — | charge (3+, modest) | documented limitation, not modeled |

**For timsim:** set each width from its physical covariate (RT: sub-linear gradient
law or per-dataset empirical FWHM; IM: charge-dependent σ) and keep the existing
per-peptide/per-ion random sampling for the within-setup wiggle. **No sequence→width
head; no bimodal model.** See the two `05-timsim-*` docs for the concrete changes.

## Methodological lessons (why earlier reads were wrong)

- **Small-n correlation is treacherous:** the gradient↔width law looked airtight at
  n=4 (r=0.97, "normalize by G") and collapsed at n=102 (r=0.57, ~√G, normalization
  *hurts*). Validate scaling laws at corpus scale with LODO.
- **Global-median metrics mask covariate-dependent failure:** "constant beats gradient
  OOS" was true only because mid-range datasets dominate the median; stratifying
  exposed the short-gradient failure. Always stratify by the covariate.
- **η² has a degrees-of-freedom floor:** with many small groups, random grouping gives
  η² ≈ (G−1)/(N−1) (~0.18 here). Subtract the shuffled baseline before reading a
  "peptide explains X%" number.
- **Don't trust model-fit parameters as ground truth:** the Gaussian `ms1_rt_sigma`
  over-widens (absorbs the EMG tail) and `ms1_rt_r2`/`ms1_im_r2` penalize skew/noise.
  Calibrate against **model-free observables** (FWHM, half-height width).
- **IM mobilograms are noisy:** gate shape/modality analysis on **high intensity** and
  use **noise-aware** detectors (MAD floor + significant-dip test), or "multi-modality"
  is just noise.
- **Within-run repeats are confounded:** repeated PSMs of a peptide in one run share
  the elution event → inflated per-peptide consistency. Use **cross-run** replication.
