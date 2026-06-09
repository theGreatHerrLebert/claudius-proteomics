# timsim IM peak-width prior — recommendation (corpus-derived)

Companion to `05-timsim-rt-width-patch.md`. Same conclusion shape as RT: the IM
peak **width** is not sequence-predictable, so don't model it per-peptide — set it
from the physical covariate (here **charge**) + keep timsim's sampler.

## Evidence (`study_im_width*.py`, `study_im_multimodality.py`)

- **IM width is NOT sequence-determined.** Cross-run peptide replication on the
  robust model-free `ms1_im_fwhm` (high-intensity, deduped): genuine peptide-
  explained ≈ 0.16 *pooled* but only **~0.06 within a charge state** — the pooled
  signal is the charge effect, not sequence. High-intensity gating didn't raise it.
  (The IM *apex*/CCS *is* structure-determined — that's the √(m/z) physics head;
  only the width isn't.)
- **Single-Gaussian IM fit is poor (im_r2 ~0.65, 76–80% <0.8)** — but mostly because
  mobilograms are **noisy**, not because most are truly bimodal.
- **Multi-modality is a modest, 3+-enriched minority.** Noise-robust, intensity-
  stratified bimodality (frac01): apparent bimodality falls with SNR (z3 53%→36%
  low→high); at high SNR ~30–36% flagged, but charge-1 is already ~33% (a charge-
  independent chimera/shoulder false-positive floor). Genuine conformer excess for
  3+ over 2+ is small (~6 pts). Real, but a minority — **not worth explicit modeling.**

## Per-charge IM width (PXD061039, high-intensity; units = 1/K0)

| charge | σ (≈FWHM/2.355) | relative σ/apex |
|---|---|---|
| 1 | **0.017** | **1.8%** |
| 2 | 0.0076 | 0.8% |
| 3 | 0.0085 | 0.9% |
| 4 | 0.0081 | 0.8% |

Multiply-charged (2–4+) cluster at **σ ≈ 0.008 (1/K0) ≈ 0.8% relative**; singly-
charged ~2× wider. **timsim's current `inv_mobility_..._std` default (~0.009 abs /
~2% relative) is ~2× too wide for the dominant 2–4+ population.**

## Recommendation

- Set `inv_mobility_gru_predictor_std` as a **charge-dependent** value:
  `~0.017 (1/K0)` for z=1, `~0.008` for z≥2 (or relative: 1.8% / 0.8%). Keep
  timsim's existing per-ion sampling for the wiggle.
- **Single Gaussian** IM peak — fine for the bulk. Document the modest high-SNR-3+
  multi-modality as a known limitation; do **not** model it (low yield, hard to
  identify, large false-positive baseline).
- **No sequence→IM-width head** (~6% within-charge ceiling).

## Caveats

- Anchors from one dataset (PXD061039); IM width can vary by instrument/setup —
  corpus-aggregate the per-charge medians across the 102 datasets to finalize
  (cheap: one pass over `_mogon_raw_features`). Per-dataset empirical override is
  available the same way as RT.
- FWHM→σ uses the Gaussian 2.355 (IM peaks are more Gaussian than RT, so this is
  defensible); confirm against a simulated mobilogram if precision matters.
