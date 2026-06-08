# Landscape — what the field does

Surveyed 2026-05-22: Prosit, AlphaPeptDeep, MS2PIP, DeepLC, Chronologer;
Depthcharge + Casanovo + InstaNovo; ChemBERTa/ChemBERTa-2, MolFormer,
PeptideCLM, SMILES-BERT; the GNN frontier (MassFormer, SigmaCCS2, MoMS-Net).
References in [sources.md](sources.md).

## Fragment-intensity / property predictors

| Tool | Core model | Modifications | Output | Instrument cond. |
|---|---|---|---|---|
| **Prosit** | bi-GRU + attention (later: transformer) | ~2 fixed (Cam-C, Ox-M) | fixed 174-vector, ≤30 aa | NCE + charge |
| **AlphaPeptDeep** | transformer (MS2), CNN+BiLSTM (RT/CCS) | **atomic composition** | per-fragment | instrument + NCE meta-embedding |
| **MS2PIP** | gradient-boosted trees (XGBoost) | data-dependent | per-fragment | per-instrument models |
| **DeepLC** | multi-path CNN | **atomic composition** | scalar RT | global features |
| **Chronologer** | residual CNN | mass-bracket UNIMOD | scalar RT | — |

## The central axis — modification encoding

| Paradigm | Who | Unseen-modification behaviour |
|---|---|---|
| Opaque token | Prosit, **imspy today** | fails — embedding stays ~random |
| **Atomic composition** | pDeep2 (originator), DeepLC, AlphaPeptDeep | **generalises** — chemically similar mods share a representation |
| Full structure (SMILES / graph) | PeptideCLM, MassFormer | best generalisation (+3D shape, relevant to CCS); heaviest |

Field consensus (AlphaPeptDeep, DeepLC): **atomic composition is the cheap,
proven unlock** for arbitrary-PTM generalisation — represent a modification by
its element-count delta (C, H, N, O, S, P, …) so an unseen mod inherits meaning
from chemically similar seen ones. Full structure (SMILES/graph) is the
frontier, needed for 3D-shape awareness (ion mobility) and non-canonical /
cyclic peptides.

## Output formulation

Prosit's fixed whole-spectrum vector caps peptide length (29 positions × 6 ions
= 174 → ≤30 residues). MS2PIP and AlphaPeptDeep predict **per fragment** —
fragment-indexed, no length cap, natural for variable-length peptides. UniSpec
(NIST, 2024) goes further still — full-range fragment-ion-series prediction
beyond the b/y charge-1–3 channels (internal ions, neutral losses, immonium);
a later extension axis if 6 b/y channels prove limiting.

## Reusable infrastructure — Depthcharge

Depthcharge (Fondrie / Noble Lab) is a PyTorch library of mass-spectrometry
transformer primitives: an analyte (peptide) encoder/decoder, sinusoidal float
encoders, a ProForma tokenizer, and a `global_token_hook` for injecting
batch-level metadata. Casanovo (de novo sequencing) and InstaNovo build on it.
A clean-slate model starts here rather than from scratch.

## Chemical / peptide language models

- **ChemBERTa-2** — RoBERTa on SMILES; chemistry-aware tokenisation + multi-task
  regression beat generic masked-LM.
- **MolFormer** — SMILES transformer scaled to ~1.1 B molecules (rotary
  embeddings, linear attention).
- **PeptideCLM** — represents a whole modified/cyclic peptide as one SMILES
  string; handles non-canonical residues that fixed amino-acid alphabets cannot.
- **SELFIES** — a string representation with 100 % chemical validity.

Takeaway: structure-aware encodings are powerful but heavier; for *modified
standard peptides* the atomic-composition middle ground captures most of the
generalisation benefit at a fraction of the cost.

## Pretraining trend

The field is moving to self-supervised pretraining on large unlabelled corpora
(ChemBERTa: masked-LM on 77 M SMILES; PRISM: 1.2 B spectra; the 2025
MS-proteomics foundation model). The MOGON campaign's ~28 M-PSM corpus is large
enough to pretrain a peptide encoder before per-task fine-tuning.
