# THE timsTOF Model — Synthesis

What it takes to model **every aspect of a timsTOF DDA experiment that you
cannot get from first principles**: which quantities are sequence-determined
(worth a neural net), which are setup-determined (worth a calibrated prior, not
a net), and what the empirical investigations found about the ceiling of each.

This is the capstone over the architecture docs (`00`–`05`), the shape findings
(`RT_IM_DISTRIBUTION.md`), and the intensity-ceiling investigation. It serves
two consumers at once:

- **Rescoring** — predicted properties (intensity, RT, IM, charge) feed PSM
  scoring / peptide-centric extraction.
- **Simulation** — the same model + a calibrated shape layer configure `timsim`
  to generate faithful synthetic runs.

---

## The one organising idea: two regimes

A timsTOF observable is either **identity-determined** (set by the peptide
sequence + its modifications + charge — a function a network can learn) or
**setup-determined** (set by the LC/instrument method, with sequence
contributing little — a covariate prior + a stochastic sampler, *not* a net).

The whole investigation was, in effect, sorting each observable into one of
these two buckets — and the surprise was that the *centers* are identity-
determined (learnable) while the *widths* are setup-determined (not).

| Observable | First-principles? | Governed by | Sequence-predictable? | Empirical ceiling | Decision |
|---|---|---|---|---|---|
| **MS2 b/y intensity** | ✗ | sequence + CE + instrument | ✅ strongly | SA ~0.73–0.75 oracle; bounded by the b/y *encoding* (17.5% of MS2), not data | **Learn** (intensity head); lever = PTM coverage, not more data |
| **RT apex** | ✗ | sequence (hydrophobicity) | ✅ | Chronologer, production | **Learn** (Chronologer) |
| **IM apex / CCS** | partial (∝√(m/z)) | sequence/structure + m/z | ✅ | √(m/z) physics-prior head | **Learn** (physics-anchored head) |
| **Precursor charge** | ✗ | sequence + pH/source | ✅ | charge head ~0.80 acc | **Learn** (charge head) |
| **RT peak width** | ✗ | **LC gradient** (~√G) + noise | ✗ (~20%, no usable feature) | low ceiling; OOS no lift over a stratified prior | **Prior** (gradient law) + sampler |
| **IM peak width** | ✗ | **charge** + noise | ✗ (~6% within charge) | charge-driven | **Prior** (per-charge σ) + sampler |
| **IM multi-modality** | ✗ | charge (3+, modest) | — | real but a minority over a noise floor | **Don't model** (documented limitation) |

The boundary in the table — the line between the four "Learn" rows and the
three "Prior" rows — is the central empirical result. It was not assumed; each
row below was stress-tested to land it on the correct side.

---

## Regime 1 — identity-determined (the neural net: `peptide-property-ng`)

One shared Depthcharge encoder, hybrid UNIMOD modification encoding (learnable
token **+** atomic-composition feature, open-vocab), task heads for intensity /
RT / IM / charge. This is the part of "THE model" that genuinely needs to be
learned, because these quantities are functions of sequence that no closed form
recovers.

### MS2 fragment intensity — the headline, and its real ceiling
- The model reaches **SA ~0.73–0.75 (oracle CE)** on phospho held-out — **+13 pp
  over Prosit-PTM, +4 pp over PeptDeep** — once you feed it *normalized* CE.
- **The "phospho gap" was a benchmarking artifact**, not a model deficit: raw
  Bruker CE (~25 V) fed to a model trained on normalized CE (~0.26) silently
  invented a 10-pp gap. Rule: match the *input* distribution, not just the
  output format, when benchmarking a model you didn't train.
- **The true ceiling is the b/y target encoding.** Only **17.5% of MS2
  intensity** sits at theoretical b/y positions (Sage's "19% matched" ≈ the same
  number from the other side). The other ~80% is immonium / internal / neutral-
  loss / chimera / noise.
- **But extending the encoding does NOT help** on tryptic data: a reverse-decoy
  test showed NL + immonium + internal ions carry **~0 information beyond b/y**
  (the +5.4 pp raw match is random matching against a denser peak set — the decoy
  explains the same). Do *not* add NL/internal as predictor targets blindly.
- **Catastrophic forgetting did NOT happen** (also a CE-bug artifact). Adding
  PTM-enriched datasets to the fine-tune was a wash; ~any sensible fine-tune
  clusters at 0.73–0.75.
- **The actual lever:** PTM *training coverage* + capacity, and — still untested
  — **residue/mod-conditional** neutral losses on PTM-rich data (−H₃PO₄ only on
  phospho-bearing fragments), which *should* be informative where the
  unconditional version was null.

### RT apex, IM apex/CCS, charge
- RT apex → **Chronologer** (production); IM apex → **√(m/z) physics-prior head**;
  charge → charge head. All sequence/structure-determined and working. The IM
  *apex* is structure-determined even though the IM *width* is not — keep that
  distinction sharp.

---

## Regime 2 — setup-determined (the timsim shape layer, NOT a net)

The decisive negative result of the investigation: **peak widths are not
sequence-predictable, so do not build a sequence→width head.** Configure timsim
from the physical covariate + keep its random sampler for the within-setup
wiggle. This killed an over-engineered predictor *before any model code shipped*.

### RT peak width — gradient-driven, sub-linear
- Model-free FWHM across **102 datasets**: `FWHM(G) ≈ 0.200 · G^0.424` (≈√G, NOT
  linear; exponent stable 0.41–0.43). EMG core `σ ≈ 0.441 · FWHM`.
- **It is NOT gradient-independent** — ~2–4 s FWHM at ~12 min vs ~9–13 s at
  ~2–3 h (3–5× at the extremes). A constant width is **45% off at the short end**.
- Genuine per-peptide signal is only **~18–22%** (cross-run, after subtracting the
  η² df-floor), with **no usable feature** (length corr ≈ 0). Within-peptide
  cross-run CV ≈ 0.31 — large stochastic wiggle.
- The EMG **tail λ is unidentifiable** at timsTOF MS1 sampling (7–13-pt XICs) →
  set σ, let the simulator sample λ.
- **Decision:** recalibrate timsim's `calculate_rt_defaults` from its current
  linear σ(G) to the sub-linear law (per-dataset empirical FWHM when simulating a
  known run); keep the random sampler. See `05-timsim-rt-width-patch.md`.

### IM peak width — charge-driven
- Genuine per-peptide signal **~0.06 within a charge state** (pooled ~0.16 is the
  charge effect, not sequence). High-intensity gating didn't raise it.
- **Gradient-independent.** Per-charge σ (1/K0): z1 ~0.017 (~1.8%), z2–4 ~0.008
  (~0.8%) — timsim's 2%-relative default is ~2× too wide for multiply-charged.
- **Decision:** charge-dependent σ prior + sampler. See `05-timsim-im-width-note.md`.

### IM multi-modality — real but not worth modeling
- A modest, 3+-enriched minority over a charge-independent false-positive floor
  (~33% even at z1, from chimera/shoulders). Genuine 3+ excess ~6 pts. Document
  as a limitation; don't add a bimodal model.

---

## What "THE model" actually is

```
                    ┌─────────────────────────────────────────┐
   peptide          │  peptide-property-ng (shared encoder)    │
   + mods    ───────▶  intensity · RT apex · IM apex · charge  │──▶ rescoring
   + charge         └─────────────────────────────────────────┘        │
                                                                        │
                    ┌─────────────────────────────────────────┐        ▼
   LC gradient ─────▶  timsim shape layer (calibrated priors)  │     timsim
   charge      ─────▶  RT σ = f(G, ~√G)  ·  IM σ = f(charge)   │──▶  synthetic
                    │  + existing random sampler (the wiggle)   │      runs
                    └─────────────────────────────────────────┘
```

- **One learned model** for the identity-determined properties (the genuinely
  hard-from-first-principles part).
- **One thin calibrated layer** for the setup-determined distributions — *not* a
  network, because the data says sequence doesn't carry the signal.
- The learned model serves rescoring directly and feeds timsim its centers; the
  shape layer adds faithful widths. Together they generate a complete, faithful
  synthetic timsTOF experiment.

---

## Highest-leverage open work (ordered)

1. **Conditional neutral losses on PTM-rich data** — the one intensity lever not
   yet tested; the unconditional NL extension was a decoy-null, but residue/mod-
   conditional losses on phospho data should differ. Retest before any encoding
   change.
2. **Recalibrate timsim's gradient→σ prior** to the sub-linear law (the RT patch).
   Concrete, validated, ready.
3. **Corpus-wide per-charge IM σ anchors** — finalize beyond the single
   PXD061039 numbers across the 102 datasets (raw_features sync, no blobs needed).
4. **PTM training coverage** for intensity — capacity + data on modifications,
   not more unmodified tryptic data (scaling sweep shows ~0 lift from volume).

---

## Methodological lessons (why the easy reads were wrong)

These recur and are worth internalizing — most of the wrong turns came from one
of these, and each was caught by stress-testing a confident claim:

- **Small-n correlation is treacherous** — the gradient↔width law looked airtight
  at n=4 (r=0.97) and collapsed at n=102 (~√G; linear normalization *hurts*).
  Validate scaling laws at corpus scale with leave-one-dataset-out.
- **Global-median metrics mask covariate-dependent failure** — "constant beats
  gradient OOS" was true only because mid-range datasets dominate the median;
  stratifying exposed the short-gradient failure. Always stratify by the covariate.
- **η² has a degrees-of-freedom floor** ≈ (G−1)/(N−1) — subtract the shuffled
  baseline before reading a "peptide explains X%" number.
- **Don't trust model-fit parameters as ground truth** — the Gaussian σ over-
  widens by absorbing the EMG tail; fit-R² penalizes legitimate skew/noise.
  Calibrate against model-free observables (FWHM).
- **Within-run repeats are confounded** — repeated PSMs of a peptide in one run
  share the elution event → inflated per-peptide consistency. Use cross-run.
- **Verify the input distribution when benchmarking** — the normalized-CE bug
  invented a 10-pp phospho "gap" and made every downstream ablation interpret
  noise. Match what the model was trained on, not just the output format.
- **"Learn the distribution first"** — empirically characterizing each observable
  before building a head killed an over-engineered width predictor with zero model
  code written. The cheapest model is the one you prove you don't need.
