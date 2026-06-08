# Recommended Architecture — and the package built from it

Implemented as **`rustims/packages/peptide-property-ng/`** — a clean-slate,
Depthcharge-based unified peptide property predictor. A research-track sibling
of `imspy-predictors`, not a production replacement.

## Stack

- **Encoder** — subclass of Depthcharge's `AnalyteTransformerEncoder`. Presets:
  `SMALL` (d_model 256, 6 layers, 8 heads, FF 1024 — the prototype) and
  `RESEARCH` (384 / 8 / 8 / 1536). The Rust ProForma tokenizer is kept as the
  token-id source.
- **Hybrid residue embedding** — per residue: a learnable token embedding **+**
  an atomic-composition feature. The composition encoder is bias-free, so a
  bare residue (zero composition) maps to an exact zero chemistry term and
  falls back to the pure token embedding; a modified residue gets token +
  chemistry. Additive fusion (a gated variant is available).
- **Composition table** — built offline from `sagepy.core.unimod`
  (`modification_atomic_composition()`): 1568 modifications resolved, a
  31-element alphabet (incl. stable isotopes and metals), *signed* counts
  (neutral losses are negative — encoded with a sign-aware log compression).
  Using sagepy's UNIMOD view keeps the chemistry consistent with the Sage
  search labels. A new UNIMOD release only needs the table rebuilt.
- **Conditioning** — instrument model + acquisition mode via the encoder's
  `global_token_hook`; charge / m/z / collision energy at the heads (so the
  charge head never sees its label — verified by a test).

## Heads (unified multi-task)

| Head | Design |
|---|---|
| **Intensity** | Fragment-indexed — one prediction per cleavage site from the two flanking residue outputs → `(L-1, 6)` (b/y × charge 1-3). No 30-residue cap. FiLM-conditioned on charge + CE. |
| **CCS / ion mobility** | Ports the imspy `SquareRootProjectionLayer` √(m/z) per-charge physics prior + a learned NN correction. Predicts 1/K0 directly (avoids a Mason-Schamp conversion); outputs `(mean, std)`. |
| **RT** | MLP on the encoder global token; trained on Sage `aligned_rt`. Research target only — production RT stays Chronologer. |
| **Charge** | MLP on the global token → charge logits. Leak-free by construction. |

## Training

- Sage `results.sage.parquet` + `matched_fragments.sage.parquet` →
  filter (decoy-free, rank 1, q ≤ 0.01, ≥6 matched peaks) → UNIMOD conversion
  via `sagepy_rescore._parse_sage_peptide` (residue-specificity-aware) →
  tokenise → per-task targets.
- Multi-task loss: spectral-angle (intensity) + L1 (CCS, prototype-stable —
  Gaussian-NLL revisited later) + L1 (RT) + cross-entropy (charge), fixed
  weights.
- **Peptide-level split** — a peptide is hashed to a fixed bucket, so it never
  crosses train/val/test regardless of dataset, charge or run.

## Module layout

```
src/peptide_property_ng/
  modifications/  composition.py, build_table.py, data/mod_composition_table.npz
  model/          config.py, embedding.py, encoder.py, multitask.py,
                  heads/{intensity,ccs,scalar}.py
  data/           sage_dataset.py, fragment_targets.py, collate.py, splits.py
  losses.py
  train/train.py    eval/{metrics,evaluate}.py
tests/            test_composition.py, test_encoder.py, test_multitask.py
```

## First-prototype result

Trained on the ~25 processed datasets in `/scratch/claudius-proteomics`
(~76 k peptide-level-split training PSMs), `SMALL` preset, ~6 M parameters.
Runnable end-to-end:

```bash
python -m peptide_property_ng.modifications.build_table   # one-off
python -m peptide_property_ng.train.train --datasets-glob '/scratch/claudius-proteomics/*'
```

Reports per-task held-out metrics — intensity spectral angle, CCS / RT MAE,
charge accuracy. (See the run's `metrics.json`.)

## What is deferred (flag-gated upgrades)

Self-supervised MLM pretraining; uncertainty-based loss weighting; gated
fusion; extra ion types; the `RESEARCH` preset; the **unseen-modification
evaluation split** — train with a modification held out entirely, compare the
full hybrid model against a composition-zeroed ablation. That split is the
decisive, falsifiable test that the chemistry encoding delivers the
generalisation it promises.
