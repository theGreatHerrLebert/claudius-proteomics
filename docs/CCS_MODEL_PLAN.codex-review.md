# Codex review — CCS_MODEL_PLAN.md

Independent second-opinion review of `docs/CCS_MODEL_PLAN.md` via `codex exec`.

> Overall: the model idea is sensible, but the plan currently overclaims
> architectural protection and under-specifies target construction. The success
> of this project will be determined less by the CNN and more by whether the
> consensus CCS table is leakage-free, bias-aware, and honest about conformers.

## 1. Physics-hybrid architecture
- The per-charge `sqrt(m/z)` baseline does **not** make charge misregistration
  "structurally impossible" — it only anchors the per-charge mean trend. The
  residual CNN (with a charge embedding) can still learn charge-dependent offsets
  that undo the baseline. Weaken the claim to "reduces the degrees of freedom for
  that error", then enforce empirically: report residual mean/median by charge,
  regularize residual to ~zero mean per charge on anchors, fit baseline on training
  folds only.
- Compare against simple baselines (`f_z(m/z)`, `f_z(mass,length)`, per-charge
  splines) — the CNN must beat them meaningfully.
- `sqrt(m/z)` may be too rigid; a per-charge spline / GAM / small monotonic MLP may
  be a better baseline.

## 2. Data strategy
- **Reproducibility weighting fixes precision, not bias.** A dataset can be
  internally reproducible and still systematically wrong. Meier may be reproducibly
  biased; weighting will not catch that.
- Failure modes: weighting bias (high reproducibility ~ easy canonical peptides →
  model good on z2/z3 tryptic, bad on rare classes); dataset dominance (cap weights
  per dataset/run); anchor bias (scalar offset insufficient if drift is nonlinear in
  mobility/charge/m/z/time); consensus collapse (a single weighted mean erases real
  multi-conformer populations); **circularity** (if the consensus table includes a
  held-out dataset's observations, the target leaks test info).
- Require **fold-aware consensus**: rebuild/mask the consensus, calibration params,
  uncertainty estimates and baseline fit so test observations never influence train
  labels. Calibrate per-charge (and possibly mobility-dependent), not a scalar
  per-dataset offset; validate with residual plots vs CCS/m-z/charge/length/mod/date.

## 3. Hidden assumptions / edge cases
- **Multi-conformer peptides need an explicit v1 policy** — exclude / dominant-apex
  only / predict mean+uncertainty / mixture. Do not silently average multimodal CCS.
- Charge imbalance underdeveloped: charge-stratified validation, min sample
  thresholds, partial parameter sharing (shared encoder + charge-specific heads).
- Modification handling underspecified: explicit supported vocab, unknown-mod
  behavior, mod/unmod split tests, benchmark by modification class, guard leakage
  where a modified peptide's unmodified counterpart is in train.
- Split by **stripped sequence** (or protein/family for stricter tests), not exact
  modified peptide — near-duplicate modified variants otherwise leak.

## 4. Evaluation gaps
- Add: baseline-only performance; ablations (no baseline / frozen / joint / spline);
  leave-one-dataset-out with no held-out contribution to consensus or calibration;
  error stratified by charge/length/m-z/CCS/mod/missed-cleavages/confidence/dataset;
  **mean signed error** by charge & dataset (not only median absolute); uncertainty
  calibration (reliability curves, NLL, interval coverage) if the uncertainty head is
  built; a search-relevance metric (predicted 1/K0 within extraction tolerance);
  multi-conformer subset reported separately.
- The reproducibility "noise floor" is not a universal lower bound — the model can
  beat individual-run noise relative to consensus.

## 5. Phasing
- P0 does too much and holds the core scientific risk. Split: **P0a** raw
  harmonization + Mason-Schamp validation; **P0b** calibration model + anchor
  diagnostics; **P0c** consensus policy + multimodality detection + leakage-safe fold
  generation.
- Do not start P1 until P0 yields frozen train/val/test manifests + consensus code
  that can exclude held-out data.
- The Meier in/out decision belongs in P0, not open into P2: include Meier in
  training only after its post-calibration residuals are shown unbiased vs cleaner
  anchors — otherwise low-weight augmentation or stress-test only.
