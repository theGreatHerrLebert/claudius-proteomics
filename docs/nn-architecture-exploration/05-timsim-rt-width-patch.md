# timsim RT peak-width prior — recalibration patch (corpus-derived)

Deliverable from plan `05-gradient-sigma-prior.md`. **Proposed change** to
`packages/imspy-simulation/.../timsim/jobs/simulate_frame_distributions_emg.py`
`calculate_rt_defaults`. Not applied here — timsim is on an active feature branch;
apply when convenient.

## Evidence (102-dataset corpus, model-free `ms1_rt_fwhm`)

- RT peak FWHM scales **sub-linearly** with gradient length G:
  **FWHM ≈ 0.200 · G^0.424** (G in seconds; log-log r=0.58; exponent stable
  0.41–0.43 across p95→max proxies; dataset-clustered LODO).
- EMG Gaussian-core σ ≈ **0.441 · FWHM** (empirical ratio from paired EMG fits on
  frac01, n=111k; ~Gaussian's 0.425 but measured, not assumed).
- The current **linear** form `σ = G/3600·0.75 + 1.125` is the wrong shape —
  ≈right at 10 min but **too narrow at long gradients** (σ≈3.6 vs faithful ~4.7 s
  at 197 min).
- Width is **not** gradient-independent: short gradients (11–47 min) have FWHM
  ~5 s, long (74–197 min) ~8 s (1.6× tercile-median; ~3–5× at the true extremes).
  A constant default is 45% off at the short end — rejected.
- Gradient sets the **gross scale**; it does **not** explain dataset-to-dataset
  scatter within a gradient range (setup confounds) — that stays in the wiggle.

## Patch — replace `calculate_rt_defaults`

```python
def calculate_rt_defaults(gradient_length: float,
                          empirical_fwhm: float | None = None) -> dict:
    """Corpus-calibrated EMG Gaussian-core sigma (seconds) for the RT peak prior.

    Calibrated on 102 timsTOF datasets (claudius-proteomics): the model-free RT
    peak FWHM scales sub-linearly with gradient, FWHM ~= 0.200 * G**0.424 (G in
    seconds), and EMG core sigma ~= 0.441 * FWHM. Replaces the prior linear form,
    which was too narrow at long gradients.

    Pass `empirical_fwhm` (median observed RT FWHM in seconds for the acquisition
    being simulated, e.g. from that dataset's extracted features) to bypass the
    gradient model entirely -- most accurate when simulating a known setup.
    """
    SIGMA_PER_FWHM = 0.441          # EMG core sigma / model-free FWHM (measured)
    if empirical_fwhm is not None and empirical_fwhm > 0:
        fwhm = empirical_fwhm                      # known setup: use it directly
    else:
        fwhm = 0.200 * gradient_length ** 0.424    # corpus sub-linear law
    sigma_middle_rt = SIGMA_PER_FWHM * fwhm
    # +/-25% = within-setup scatter the gradient cannot predict (the wiggle).
    # Corpus between-dataset IQR/median ~= 0.45 supports ~+/-22%; keep 25%.
    sigma_lower_rt = sigma_middle_rt * 0.75
    sigma_upper_rt = sigma_middle_rt * 1.25
    return {'sigma_lower_rt': sigma_lower_rt, 'sigma_upper_rt': sigma_upper_rt}
```

σ(G) sanity (seconds): 10 min ≈ 1.3 · 30 min ≈ 2.1 · 60 min ≈ 2.85 · 120 min ≈ 3.8
· 197 min ≈ 4.7. (Old linear: 1.25 / 1.5 / 1.875 / 2.625 / 3.59.)

## Per-dataset empirical lookup

`/scratch/claudius-proteomics/_mogon_raw_features/_fwhm_gradient.csv`
(`dataset, G, fwhm`) — median observed RT FWHM per dataset. When simulating one of
these accessions, pass its `fwhm` as `empirical_fwhm`. For a new dataset, compute
`median(ms1_rt_fwhm | ms1_rt_r2>=0.8)` from its extracted features.

## Caveat to verify before landing (Codex review #4)

timsim samples its **own** λ (tail). Setting σ from FWHM while λ is sampled
independently can make the *output* FWHM drift from target. **Validate with one
simulated run**: compare simulated peak FWHM to the target above; if it drifts,
nudge `SIGMA_PER_FWHM` or co-set λ. This is a 1-run check, not a blocker.

## Not doing (evidence-based)

- No sequence→σ head (per-peptide width has a ~20% ceiling, no usable features).
- No explicit λ model (EMG λ unidentifiable at timsTOF MS1 sampling) — unless the
  output-FWHM check above forces it.
