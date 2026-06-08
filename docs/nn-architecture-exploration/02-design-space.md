# Design Space — axes and decisions

## A. Input modalities (the "SMILES → UNIMOD → metadata" stack)

| Level | Options | Decision |
|---|---|---|
| Residue | amino-acid identity + position | keep |
| **Modification** | opaque token / atomic composition / SMILES-graph | **hybrid: token + composition** |
| Precursor | charge, m/z, collision energy | head-level inputs |
| Acquisition | instrument model, mode (DDA/diaPASEF), gas, source | encoder global token |
| Sample | organism, lab | skipped — likely noise |

## B. Modification representation — the crux

The instinct "just put all ~2K UNIMOD modifications in the vocabulary" is
reasonable — and imspy already does (the tokenizer has ~2100 `[UNIMOD:N]`
slots). The catch: **a vocabulary slot is not a trained embedding**. ~99 % of
*PSMs* are a handful of common modifications (those train well); the other
~1900 slots see ≈0 training examples, so their embeddings never leave random
initialisation. "Covers 99 %" is true for *throughput*, false for *modification
diversity* — and the rare tail is often the scientifically interesting part.

The fix is not a bigger vocabulary, it is giving each modification a
**chemistry** representation. The **hybrid** keeps a learnable per-mod token
(common mods train well from data) *and* folds in an atomic-composition feature
(rare/unseen mods inherit meaning from chemically similar ones). This follows
the atom-count PTM-encoding line of **pDeep2** (originator), **DeepLC** (RT) and
**AlphaPeptDeep** (MS2/RT/CCS), with an added learnable per-mod token branch
alongside it.

Open-vocab is *partial*, to be precise: the composition feature is keyed by
tokenizer token id, so a new UNIMOD entry is free only if the tokenizer already
has a stable id for it (it carries ~2100 `[UNIMOD:N]` slots, so most are
covered); a change in vocabulary *size* still needs a config bump and a table
rebuild. The win over opaque tokens is real, but it is "rebuild the table, not
retrain the vocabulary", not literally unbounded.

## C. Output formulation

- RT / CCS / charge — scalar / small heads.
- **Intensity — fragment-indexed**, one prediction per cleavage site, replacing
  Prosit's fixed 174-vector. Removes the 30-residue cap.

## D. Architecture families considered

1. Extend the production `UnifiedPeptideModel` — incremental, low-risk.
2. **Clean-slate on Depthcharge** — new encoder from Depthcharge primitives;
   room for a structure-aware modification branch and self-supervised
   pretraining. Higher ceiling.
3. Structure-aware modification branch (GNN / SMILES sub-encoder) — fused into
   (1) or (2); the frontier option.

## E. Conditioning

Instrument model + acquisition mode are conditioned **inside** the encoder (the
global token) — they cannot leak any task label. Charge / m/z / collision
energy are conditioned at the **heads** that use them, so the charge head —
which reads the shared encoder output — never sees its own label. One encoder
pass, no leakage.

## Locked decisions (with the user)

| Question | Decision |
|---|---|
| Modification representation | **Hybrid** — token + atomic composition |
| Starting point | **Clean-slate on Depthcharge** |
| First-prototype target | **Unified** — one encoder, all four heads |

See [03-recommended-architecture.md](03-recommended-architecture.md) for the
resulting architecture.
