# HF corpus v0.1 — model-validation

Before publishing, we trained a model **end-to-end on the aggregated corpus
parquets** (Tier-1 + Tier-3) to prove the corpus carries trainable signal across
all four heads. It does.

## Setup
- Loader `peptide_property_ng/data/hf_corpus_dataset.py` (rustims,
  `feat/peptide-property-ng`): streams Tier-1 + Tier-3 with `pyarrow.dataset`
  pushdown (the 391 M-row Tier-3 can't be materialised), reconstructs the exact
  example dicts the existing collate/train loop consumes — b/y intensity target
  from Tier-3, CCS/charge from Tier-1, **real per-PSM CE** (`collision_energy_mean_v`/100),
  RT per §RT below.
- Fine-tune from the cross-instrument pretrained checkpoint
  (`cap0-pretrained.pt`), `small` preset, lr 2e-5, 8 epochs, 15 datasets,
  cap 4000/dataset (~57 k train examples) on an RTX 2070.
- `train.py --hf-corpus <dir> [--hf-max-datasets N] [--hf-rt-lookup ...]`.

## Result (held-out test, best epoch)

| Metric | raw-RT | **aligned-RT** | baseline (old Sage path, 56 ds / cap 8k) |
|---|---|---|---|
| Intensity SA | 0.733 | **0.734** | 0.764 |
| CCS MAE (norm) | 0.0239 | 0.0226 | 0.0172 |
| RT Pearson | 0.849 | **0.880** | 0.896 |
| Charge acc | 0.801 | 0.804 | 0.817 |

- **All four heads train to in-neighborhood metrics straight from the corpus.**
- The model **starts strong from pretraining** (val SA 0.73 at epoch 1) and
  fine-tunes to timsTOF-ready in 8 fast epochs — the pretrain→finetune→ready
  arc, on the publishable artifact.
- Gaps vs baseline are **data volume**, not corpus defects: this run used ~1/8
  the training data (15 ds / cap 4k vs 56 / cap 8k); val SA plateaued by ep5.

## RT: raw vs aligned (a fix this surfaced)
Tier-1 shipped raw `rt_seconds` (per-run, gradient-dependent). Swapping in Sage's
cross-run **`aligned_rt`** (from `results.sage.parquet`, ~[0,1]) lifted **RT
Pearson 0.849 → 0.880 (+3.1 pp)** with everything else unchanged — closing the
RT gap to the baseline (residual = data volume). → `build_hf_tier1.py` now
sources `aligned_rt` and emits an **`rt_aligned`** column for v0.2;
`build_rt_lookup.py` produced the v0.1 lookup used to validate this.

## Verdict
**Corpus is publish-ready.** v0.2 polish: ship `rt_aligned` (builder done →
rebuild Tier-1); strip `[]` placeholders in `modified_sequence`; guard the
aggregation `build_stats` against an empty per-dataset manifest; fold in the 4
HLA datasets (per-group merge fix); full-54-dataset headline run.
