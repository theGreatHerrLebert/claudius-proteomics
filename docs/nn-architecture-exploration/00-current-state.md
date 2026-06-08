# Current State — the `imspy_predictors` baseline

The production model being improved on. `UnifiedPeptideModel`
(`rustims/packages/imspy-predictors/.../models/`):

- **Encoder** (`transformer.py`): token embedding → sinusoidal positional
  encoding → pre-LN transformer. BASE preset = d_model 256, 6 layers, 8 heads,
  feed-forward 1024 (also small / large presets).
- **Heads** (`heads.py`): intensity, ccs, rt, charge — each a small head on a
  pooled representation. The CCS head has a physics prior worth keeping —
  `SquareRootProjectionLayer`, CCS ≈ slope·√(m/z) + intercept per charge, plus
  a learned NN correction.
- **Tokenizer**: Rust-backed ProForma tokenizer, ~2200-token vocabulary — 20
  amino acids + extended + ~30 composite AA+common-mod tokens + ~2100
  `[UNIMOD:N]` slots + specials.

## Three structural gaps

| Gap | Detail |
|---|---|
| **Modifications are opaque IDs** | A `[UNIMOD:N]` token has a *learnable* embedding, but only ~10 common modifications ever appear in training data — the other ~1900 vocabulary slots stay at random initialisation. The model sees an *ID*, not chemistry, so it cannot generalise across modifications. (Prosit-style intensity prediction is even stricter — effectively ~2 modifications.) |
| **Instrument conditioning is dormant** | A `use_instrument` flag plus a 15-type instrument embedding is wired into every head — but the production predictors never pass it; `instrument` defaults to "unknown". The metadata exists (`timstof_catalog.tsv`, Bruker `.d/analysis.tdf` `GlobalMetadata`) and is simply unused. |
| **Output length ceiling** | The intensity head emits a fixed 174-value vector (29 positions × 6 ion types) → peptides longer than 30 residues overflow it (a Rust `index 175` panic was hit in practice). RT / CCS scalar heads are unaffected. |

## What is worth carrying forward

- The **CCS √(m/z) physics prior** — a sound, interpretable inductive bias.
- The **shared-encoder + per-task-head** shape — already foundation-model-like.
- The **Rust ProForma tokenizer** — fast, and its `[UNIMOD:N]` vocabulary is a
  fine *token* layer to build the hybrid encoding on top of.

The gap is not vocabulary size — it is that a vocabulary slot is not a trained
embedding, and an opaque ID carries no chemistry. That is what the next-gen
design targets (see `02-design-space.md`).
