---
license: cc0-1.0
pretty_name: Claudius timsTOF DDA-PASEF PSM corpus
tags:
  - proteomics
  - mass-spectrometry
  - timsTOF
  - DDA-PASEF
  - peptide-property-prediction
  - fragment-intensity
  - retention-time
  - ion-mobility
  - collision-energy
size_categories:
  - 10M<n<100M
configs:
  - config_name: tier1_psms
    data_files:
      - split: train
        path: tier1/train.parquet
      - split: validation
        path: tier1/val.parquet
      - split: test
        path: tier1/test.parquet
  - config_name: tier3_fragments
    data_files:
      - split: train
        path: tier3/train.parquet
      - split: validation
        path: tier3/val.parquet
      - split: test
        path: tier3/test.parquet
---

# Claudius timsTOF DDA-PASEF PSM corpus

A large, **CC0**, multi-workflow corpus of peptide-spectrum matches from
**timsTOF DDA-PASEF** runs, built to train and benchmark peptide-property
predictors (fragment intensity, ion mobility / CCS, retention time, charge) and
to configure timsTOF simulators.

> **Design principle — a reference corpus, not a pre-baked training set.**
> We include broadly, filter minimally at build, and **label richly** so you can
> reproduce any filtering decision yourself. Every confidence/quality signal is
> a column; nothing you might want to filter on was thrown away.

## What's novel

To our knowledge the first public timsTOF DDA-PASEF PSM resource that ships:

1. **Dual-engine (Sage + FragPipe) agreement** labels per PSM and per fragment;
2. per-precursor **RT and ion-mobility peak shapes** (FWHM, σ, skew, fit R²,
   trace SNR) — not just apex point estimates;
3. per-precursor, **per-PASEF-event collision energies** (raw volts);
4. **matcher-consistent dual-engine b/y fragment** intensities.

## Cohort

**58 datasets**, 100% **CC0** (PRIDE-verified, per-dataset evidence + retrieval
dates retained), **timsTOF family** instruments (Pro / Pro2 / HT / SCP / fleX /
Ultra), DDA / PASEF acquisition. Tryptic, HLA / immunopeptidomics, and several
PTM-enriched workflows.

| Tier | rows | train / validation / test |
|---|---:|---|
| **Tier 1** — precursors | **33,717,207** | 27,128,597 / 3,205,642 / 3,382,968 |
| **Tier 3** — b/y fragments | **504,097,214** | 404,005,449 / 48,374,948 / 51,716,817 |

(Per-dataset counts + provenance: `manifest.json`.)

## Tiers

- **`tier1_psms`** — one row per accepted precursor identification
  `(raw_file, precursor_id)`: sequence + modifications, charge, m/z, RT apex
  (raw + Sage cross-run-aligned `rt_aligned`), ion mobility, per-PASEF-event
  collision energy, **RT/IM peak-shape labels** (`ms1_{rt,im}_{fwhm,sigma,skew,r2,snr}`),
  isotope envelope, per-engine scores, engine-agreement, and derived
  `width_reliable` / `strict` presets.
- **`tier3_fragments`** — one row per matched **b/y** fragment (imspy-rematched
  for both engines under identical settings): calculated + experimental m/z,
  intensity, ppm error, per-engine match flags, Sage-native provenance.

The full field-by-field schema (units, null semantics) is in `SCHEMA.md`.

## Filtering & certainty policy (read this before filtering)

- **Inclusion floor:** rank-1 **target** PSMs accepted by **at least one engine
  at that engine's reported q ≤ 0.01**. ⚠️ This is a *union of two
  engine-specific 1%-FDR sets* — **not** a guaranteed corpus-level 1% FDR. Treat
  the per-engine `q`/`pep` as the confidence signals.
- **Filter up, freely:** every score is exposed (`sage_qvalue`, `sage_pep`,
  `fragpipe_qvalue`, `fragpipe_pep`, `fragpipe_probability`, hyperscores,
  `n_engines`, `sage_cosine`, `isotope_cosim`). Applying a stricter reported-q
  threshold *selects* higher-confidence PSMs — it does **not** recompute FDR;
  post-hoc combinations don't retain a nominal FDR.
- **Engine agreement is a label, not a gate** — single-engine PSMs are kept and
  flagged (`n_engines == 1`). Want the high-confidence subset? Use the
  documented **`strict`** boolean column.
- **Peak-shape reliability:** width labels carry `ms1_*_r2` + `ms1_*_snr` + the
  `*_width_reliable` flags. ⚠️ The SNR-reliability band is **provisional**; gate
  width labels on `width_reliable` for peak-shape work.

## Splits

**Peptide-hash** on `sequence_normalized` → ~80 / 10 / 10 train / val / test
(seed 0) — **no peptide leakage** across splits. Dataset/group-level (sibling
PXD) splitting is deferred to a later version. The `split` column is in every row.

## Validation — a model trained on this corpus

Training a unified peptide-property model **directly on these parquets** (from a
cross-instrument pretrained checkpoint, all 58 datasets) reaches, on the
held-out test split:

| Head | Metric |
|---|---|
| MS2 b/y intensity | spectral angle **0.744** |
| Ion mobility (CCS) | normalized MAE **0.0184** |
| Retention time | Pearson **0.876** (native `rt_aligned`) |
| Charge | accuracy **0.809** |

The model **fine-tunes to timsTOF-ready in a few epochs** from pretraining — the
corpus carries trainable signal across all four heads, at parity with an
equivalent model trained from the raw Sage outputs. (Using raw `rt_seconds`
instead of `rt_aligned` costs ~3 pp RT Pearson — prefer `rt_aligned`.)

## Usage

```python
from datasets import load_dataset

# precursor properties (RT / IM / charge / CE / peak shapes)
t1 = load_dataset("claudius-proteomics/timstof-psms-cc0", "tier1_psms", split="train")

# high-confidence subset (one boolean)
strict = t1.filter(lambda r: r["strict"])

# b/y fragment intensities (spectrum prediction)
t3 = load_dataset("claudius-proteomics/timstof-fragments-cc0", "tier3_fragments", split="train")
```

Or read the parquet directly with `pyarrow` / `polars` (the files are plain
zstd parquet; use predicate pushdown on `accession` / `split` for subsets).

## Limitations

- The inclusion floor is a **per-engine q-union, not a corpus FDR** (see above).
- **b/y fragments only** (no neutral-loss / immonium / internal ions).
- Fragment matching uses a fixed mass→UNIMOD converter; **fragments for
  unmapped modifications are skipped** — PTM-heavy datasets have reduced fragment
  coverage (tracked, not silently dropped).
- The **SNR-reliability band is provisional**; `rt_aligned` is Sage's per-dataset
  cross-run alignment (not a universal iRT).
- Decoys are not shipped in this release (a decoy companion for FDR
  recomputation / rescorer training is planned).

## Versioning

Semantic, **immutable** revision tags on this repo (`v0.1`, `v0.2`, …); a
published revision is never edited in place. Each revision ships a `manifest.json`
(dataset list, build-pipeline commit, filter parameters, schema version) so the
release is reproducible from the source artifacts.

## License & citation

All datasets are **CC0** (PRIDE-verified; per-dataset license + retrieval date
in `manifest.json`). Please cite the originating PRIDE accessions (DOIs listed in
the dataset metadata) **and** the Claudius processing pipeline (citation TBD).
