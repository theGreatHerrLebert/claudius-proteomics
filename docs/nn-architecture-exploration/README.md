# NN Architecture Exploration — Next-Gen Peptide Property Predictor

Captured research and design from the 2026-05-22 deep dive into building a
next-generation neural network for **peptide property prediction** — fragment
intensity, ion mobility / CCS, retention time and precursor charge — on timsTOF
proteomics data.

## Contents

| File | What |
|---|---|
| [00-current-state.md](00-current-state.md) | The production `imspy_predictors` baseline and its structural limits |
| [01-landscape.md](01-landscape.md) | What the field does — Prosit, AlphaPeptDeep, MS2PIP, DeepLC, Depthcharge/Casanovo, ChemBERTa/PeptideCLM, GNN frontier |
| [02-design-space.md](02-design-space.md) | The design axes, the modification-representation question, the locked decisions |
| [03-recommended-architecture.md](03-recommended-architecture.md) | The chosen architecture — and the package built from it |
| [sources.md](sources.md) | References |

## Status

The recommended architecture has been **implemented** as a new research-track
package: `rustims/packages/peptide-property-ng/`. The first prototype trains a
unified multi-task model on the processed Sage results in
`/scratch/claudius-proteomics`. See `03-recommended-architecture.md`.

## One-paragraph summary

The production predictors work but are bounded by their *architecture*:
modifications are opaque tokens (rare/unseen PTMs are unpredictable), instrument
conditioning is built but dormant, and the intensity output is locked to
Prosit's 174-vector (a hard 30-residue cap). The next-gen design is a
clean-slate, Depthcharge-based **unified** model (one encoder → intensity + CCS
+ RT + charge) with a **hybrid modification encoding** — every residue carries
both a learnable per-mod token *and* an atomic-composition chemistry feature, so
rare/unseen modifications generalise instead of failing, with no vocabulary
re-mint when UNIMOD grows.
