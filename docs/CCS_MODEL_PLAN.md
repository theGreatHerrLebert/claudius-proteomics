# CCS Prediction Model — Design Plan

Status: draft, revised after Codex review (see `CCS_MODEL_PLAN.codex-review.md`)
Target: San José Phase 2 (prediction models)

## 1. Motivation

San José's three prediction heads are CCS, RT and MS2 intensity. An exploration of
the current predictors on PXD046777 (see `notebook/rt_im_exploration.ipynb`,
`notebook/intensity_exploration.ipynb`) found:

- **RT** has a strong external model — Chronologer (Searle Lab) — Pearson r 0.990,
  ~33 s median residual, ~5× tighter than the imspy GRU/transformer baseline.
- **IM/CCS has no equivalent.** The imspy `DeepPeptideIonMobilityApex` is good
  *within* a charge state (within-charge r ≈ 0.96) but mis-registers the per-charge
  1/K0 bands, so its pooled correlation collapses to ~0.87 (Simpson's paradox).
  Sage's own `predicted_mobility` is the better-registered baseline.

So the missing piece is a strong, charge-correct CCS model. This plan applies the
**Chronologer recipe** to CCS.

### What makes Chronologer work

1. **A universal target scale.** Chronologer predicts C18 retention coefficient in
   `% ACN`; every source dataset is aligned onto that one axis, so 11 heterogeneous
   LC setups become one coherent training problem.
2. **Scale + pooling.** >2.6 M retention observations / 2.25 M peptides, 11 datasets.
3. **A small dilated residual CNN** (~500 K params): residue embedding → 3 dilated
   residual conv blocks (kernel 7) → flatten → linear. Small + lots of data, well
   regularized (dropout 0.1), growing batch size 64→1024, 100 epochs, Adam 1e-3.

### Why CCS is an even better fit than RT

CCS (Å²) is a **physical, instrument-independent molecular property** — it *is* the
universal scale, for free. RT needed the empirical `% ACN` harmonization; CCS does
not. Two caveats remain and are *not* free: per-run TIMS **calibration** drift
(systematic) and per-measurement **imprecision** (random). The corpus design (§5)
treats these as two separate problems — conflating them is the main risk this plan
guards against.

## 2. Goals & non-goals

**Goals**
- A CCS predictor that is correct *across* charge states, not just within them.
- Beat `DeepPeptideIonMobilityApex` and Sage `predicted_mobility` on per-charge and
  pooled metrics, on held-out peptides **and** a held-out dataset.
- Per-dataset fine-tune from a checkpoint (Chronologer-style), cheaply.

**Non-goals (for v1)**
- RT — Chronologer already covers it; reuse it.
- Multi-conformer / mixture CCS modelling (v2 — see §3, §6).
- Non-tryptic / exotic modification coverage beyond the San José corpus's profile.
- Replacing the search-engine IM predictors used inside search.

## 3. Prediction target

Predict **CCS (Å²)**, not 1/K0 directly.

- Training labels: observed `1/K0 → CCS` via Mason–Schamp (needs m/z + charge, both
  known). Inference: predict CCS → convert to 1/K0 per dataset's m/z + charge.
- CCS is instrument-independent, so observations from different labs/instruments are
  poolable — *after* calibration (§5), which corrects systematic per-run drift.
- **Train on a consensus CCS per `(peptide, charge)`**, built fold-aware (§5) so the
  label never sees held-out data.
- **Single-conformer only for v1.** A `(peptide, charge)` whose pooled 1/K0
  observations are multi-modal is excluded from v1 training (§5) — averaging it into
  one consensus value is a silent error. Multimodal precursors are detected, kept
  aside, and reported separately at evaluation; mixture/uncertainty modelling is v2.

## 4. Architecture — physics-hybrid

A pure sequence→CCS CNN is what `DeepPeptideIonMobilityApex` effectively is, and it
mis-registers charges. We add an explicit per-charge physics baseline so the bulk of
the charge structure is fixed by construction and the CNN learns only the
sequence-specific residual.

```
inputs:  tokenized modified sequence (UNIMOD) | charge z | precursor m/z

  CCS_baseline(z, m/z)        per-charge smooth size trend, e.g.
                              CCS ≈ a_z · sqrt(m/z) + b_z, or a per-charge
                              monotone spline if sqrt-linear is too rigid

  ΔCCS(sequence, z)           Chronologer backbone:
                              residue embedding → 3 dilated residual conv
                              blocks (kernel 7) → flatten → charge-embedding
                              concat → MLP head → residual

  CCS_pred = CCS_baseline(z, m/z) + ΔCCS(sequence, z)
```

**Honest framing of what the baseline buys (Codex catch).** The per-charge baseline
does **not** make charge mis-registration *impossible* — it anchors the per-charge
mean trend, but the residual head (which sees a charge embedding) retains the
freedom to re-introduce a charge-dependent offset. The baseline *reduces the degrees
of freedom* for that error; it does not remove them. So it must be enforced and
checked, not assumed:

- Fit `CCS_baseline` on **training folds only**.
- After training, **report residual mean/median per charge** — the registration
  check; mean signed residual per charge should be ≈ 0.
- Add a soft penalty driving the per-charge mean residual on calibration anchors
  toward 0 (or constrain the residual head to be charge-mean-free).
- Decide empirically whether to freeze the baseline or fine-tune it jointly.

**Baseline-form choice.** `sqrt(m/z)` is a reasonable size proxy but may be too
rigid (charge state affects compactness and protonation, not just size). Evaluate a
per-charge monotone spline / small monotonic MLP over `(mass, m/z)` as the baseline.

**The CNN must earn its place.** Benchmark it against simple baselines —
`CCS ~ f_z(m/z)`, `CCS ~ f_z(mass, length)`, per-charge splines. The CNN ships only
if it beats these meaningfully (§7).

Size: ~0.5–1 M params. Optional v2: a heteroscedastic head predicting per-peptide
CCS uncertainty.

## 5. Training corpus & quality weighting

**The Meier finding.** PXD019086 (Meier et al. 2021) is the obvious large CCS
reference (~2 M peptide-charge observations) and is already in the dataset list.
But initial exploration shows its CCS values are **imprecise** — early timsTOF
hardware — so it is *diverse but noisy*.

**Precision is not accuracy (Codex catch).** Reproducibility weighting fixes
*random* error (imprecision). It does **nothing** about *systematic* error: a
dataset can be internally reproducible and still biased relative to modern
calibrated measurements. The two are handled separately:

- **Systematic bias → calibration.** Per-dataset (and where needed per-run) CCS
  correction fit on anchor peptides shared across datasets. A *scalar* offset is
  insufficient if drift is nonlinear; the calibration is fit **per charge** and, if
  anchor diagnostics show structure, as a function of mobility / m/z. Validate with
  residual plots vs CCS, m/z, charge, length, modification class and run date.
- **Random error → precision weighting.** San José computes cross-run CCS
  reproducibility per precursor (the stability analysis in
  `rt_im_exploration.ipynb`). Each observation carries a measured precision used as
  a weight.

**Pooling.** Pool every processed timsTOF dataset's `(peptide, charge, 1/K0)` → CCS.
Meier contributes sequence/charge **diversity**; modern datasets contribute
**precision**. Guards:

- **Cap per-dataset (and per-run) weight** so a single large dataset cannot dominate
  the consensus even after down-weighting.
- **Watch weighting bias:** high reproducibility correlates with abundant, canonical
  tryptic z2/z3 peptides. Track and report per-class coverage so the model is not
  silently excellent on easy peptides and poor on rare ones.

**Consensus CCS labels — fold-aware (Codex catch).** Aggregate observations of a
`(peptide, charge)` into a precision-weighted consensus CCS + uncertainty. Train on
the consensus (denoised), weighting each example by its uncertainty.

> **Leakage guard.** The consensus value, the calibration parameters, the precision
> weights and the per-charge baseline are all *learned from data*. For any held-out
> evaluation (held-out peptide or held-out dataset), they must be **rebuilt with the
> held-out observations masked out** — otherwise the training *target itself* leaks
> test information. P0c produces consensus-generation code that takes a fold spec.

**Multi-conformer detection.** Before forming a consensus, test each
`(peptide, charge)` 1/K0 distribution for multimodality. Multimodal precursors are
**excluded from v1 training**, kept in a separate set, and evaluated separately
(§7). Do not average them into a single consensus value.

**The consensus table is a deliverable.** A versioned, fold-aware CCS reference,
independent of the model.

**Meier: in or out — decided in P0, not later.** Include Meier in *training* only
after P0b shows its post-calibration residuals are unbiased against cleaner anchors.
If residual bias survives calibration, Meier is used as low-weight augmentation or
as a held-out stress-test set only — not a major label source.

## 6. Training procedure

- Chronologer hyperparameters as the starting point: ~100 epochs, Adam lr 1e-3,
  dropout 0.1, batch size 64 → ×2 every 30 epochs → 1024.
- Loss: Huber on the CCS residual, **weighted by consensus label precision** (or a
  heteroscedastic Gaussian NLL if the uncertainty head is built in v2).
- Split **by stripped (unmodified) sequence** — no stripped sequence in both train
  and test, so modified variants of the same peptide cannot leak. Additionally hold
  out **whole datasets** for the generalization test.
- Fit the per-charge baseline on training folds only; freeze vs joint fine-tune
  decided empirically in P2.
- Charge imbalance: z2/z3 are plentiful, z4+ sparse and noisier. Use a shared
  sequence encoder with charge-specific residual heads, enforce minimum sample
  thresholds per charge, and validate per charge.

## 7. Evaluation

- **Per charge, always** (the IM exploration lesson): per-charge Pearson r, median
  absolute CCS error (Å² and %), **mean signed error** (the registration check —
  median absolute alone hides a charge offset), and the within- vs pooled-r gap.
- **Baselines & ablations:** baseline-only performance; ablations of no-baseline /
  frozen-baseline / jointly-trained-baseline / spline-baseline. The CNN ships only
  if it beats the simple baselines.
- **Leakage-safe generalization:** leave-one-dataset-out, with the consensus,
  calibration and baseline rebuilt with that dataset masked.
- **Stratified error:** by charge, peptide length, m/z, CCS, modification class,
  missed cleavages, identification confidence, and dataset.
- **Benchmarks:** `DeepPeptideIonMobilityApex`, Sage `predicted_mobility`, and
  Meier's own deep model (noting Meier's reported R is against its own noisy labels).
- **Search-relevance metric:** fraction of predicted 1/K0 within practical
  extraction/search tolerance.
- **Multi-conformer subset** reported separately, never folded into the headline.
- **Noise floor, stated carefully:** cross-run reproducibility of the *observed*
  values is not a universal lower bound — against a denoised consensus the model can
  beat single-run noise; it cannot beat irreducible biological/instrumental
  ambiguity.
- Uncertainty head (v2): reliability curves, NLL, prediction-interval coverage.

## 8. Fine-tuning / extensibility

Chronologer-style: ship a base checkpoint; allow cheap per-dataset fine-tuning
(re-train from checkpoint, reduced LR, on a dataset's confident anchors). The
intensity and RT fine-tune tests already show this works and helps; the CCS model
exposes the same `fine_tune` entry point.

## 9. Phased plan

P0 holds the core scientific risk and is split. P1 does not start until P0 produces
frozen train/validation/test manifests and fold-aware consensus-generation code.

| Phase | Deliverable |
|-------|-------------|
| **P0a** | Raw harmonization: pool timsTOF `(peptide, charge, 1/K0)`; validate the Mason–Schamp `1/K0 → CCS` conversion. |
| **P0b** | Calibration: per-charge (and if needed mobility-dependent) per-dataset correction on shared anchors; anchor diagnostics; **the Meier in/out decision** (residual-bias test). |
| **P0c** | Consensus policy: precision weighting + per-dataset weight caps, multimodality detection, fold-aware consensus-generation code, frozen train/val/test manifests. Deliverable: a versioned CCS consensus reference. |
| **P1** | Physics-hybrid model in `imspy-predictors` (`ccs/`): per-charge baseline + dilated-resnet residual + charge embedding + per-charge residual-mean enforcement. |
| **P2** | Train on the consensus corpus; baseline & ablation benchmarks; per-charge, held-out-peptide and leave-one-dataset-out evaluation; decide freeze vs joint baseline. |
| **P3** | Per-dataset fine-tune hook + base checkpoint release. |

## 10. Open questions & risks

- **Charge coverage.** z4+ is sparse and noisier (per-charge r fell to ~0.81 for z4
  in the exploration). Mitigations in §6; may still need z4+ pooling.
- **Modifications.** Scope the supported UNIMOD set (1/4/35 dominate; HLA datasets
  add others), define unknown-modification behavior, and benchmark per modification
  class. Tokenizer must match.
- **Baseline coupling.** Freeze the per-charge baseline or fine-tune it jointly —
  decided empirically in P2.
- **Consensus aggregation.** Weighting scheme, minimum observations per peptide,
  multimodality test and threshold.
- **v2 directions.** Heteroscedastic uncertainty head; mixture/multi-conformer CCS.

## 11. Where it lives

- Model + training code: `imspy-predictors` (`ccs/`), alongside the existing
  `DeepPeptideIonMobilityApex`, which becomes the baseline-to-beat.
- Consensus CCS reference table: a versioned San José database artifact.
- Chronologer stays the RT predictor; this plan is CCS-only.
