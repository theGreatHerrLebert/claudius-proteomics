# San Jose Dataset Priority List — DDA/PASEF Only

Generated: 2026-02-11
Source: PRIDE catalog (1,773 timsTOF datasets, 1,029 DDA/PASEF)

## Selection Criteria

1. **DDA only** — DDA-PASEF or generic PASEF (excludes diaPASEF, DIA, MALDI, PRM)
2. **Modern instruments first** — Ultra > HT > SCP > Pro 2 > fleX > Pro
3. **Organism diversity** — Expand beyond current 7 organisms
4. **Benchmark value** — Published benchmarks with known ground truth
5. **Recency** — Prefer 2024+ submissions

## Current Instrument Coverage on PRIDE (DDA)

| Model | DDA Datasets | Status |
|-------|-------------|--------|
| timsTOF Pro | 695 | Covered (PXD019086, PXD046675) |
| timsTOF Pro 2 | 118 | **Need** |
| timsTOF | 130 | Low priority (older/unspecified) |
| timsTOF HT | 31 | **Need** |
| timsTOF SCP | 27 | **Need** |
| timsTOF fleX | 26 | **Need** |
| timsTOF Ultra | 2 | **Need** (only 2 exist!) |

## Current Organism Coverage

| Organism | Supported | DDA Datasets |
|----------|-----------|-------------|
| Human | Yes | 474 |
| Mouse | Yes | 172 |
| E. coli | Yes | 15 |
| Yeast | Yes | 10 |
| Pig | Yes | 12 |
| Drosophila | Yes | 5 |
| C. elegans | Yes | 5 |
| **Rat** | **No** | 25 |
| **Arabidopsis** | **No** | 16 |
| **Zebrafish** | **No** | 6 |
| **Bovine** | **No** | 7 |

---

## TIER 1: ACTIVE (already processing)

| # | Accession | Instrument | Organism(s) | Description |
|---|-----------|-----------|-------------|-------------|
| 1 | PXD019086 | timsTOF Pro | human, yeast, ecoli, drosophila, c_elegans | Meier 2021 CCS benchmark — 1M+ CCS values |
| 2 | PXD046675 | timsTOF Pro | pig | Pig muscle 4D-LFQ (6 runs, running on monster2) |

---

## TIER 2: HIGH VALUE — Process next

These are benchmark or strategically important datasets.

| # | Accession | Instrument | Organism | Why |
|---|-----------|-----------|----------|-----|
| 3 | **PXD014777** | timsTOF Pro | human+yeast+ecoli | MaxQuant benchmark. 3-species mixture + blood plasma. Widely cited. |
| 4 | **PXD026463** | timsTOF Pro | mouse | Ionmob CCS training. 21K phosphopeptides + 17K MHC ligands. Ground truth. |
| 5 | **PXD051790** | timsTOF Pro 2 | human | Time-segmented DDA-PASEF. HeLa dilution + 85K phosphosites. **First Pro 2 dataset.** |
| 6 | **PXD060544** | timsTOF Ultra | human | Only 2 Ultra DDA datasets on PRIDE. Molecular glue study. **First Ultra.** |
| 7 | **PXD064497** | timsTOF Ultra | mouse | The other Ultra DDA dataset. Mouse brain. **Must take both Ultra.** |
| 8 | **PXD069027** | timsTOF HT | human | Most recent HT DDA (2026-01). Mitochondrial oxidative phosphorylation. **First HT.** |

---

## TIER 3: INSTRUMENT DIVERSITY — Fill hardware gaps

| # | Accession | Instrument | Organism | Why |
|---|-----------|-----------|----------|-----|
| 9 | **PXD073076** | timsTOF HT | mouse | HT + mouse (2026-01). Sperm flagellar proteome. |
| 10 | **PXD067573** | timsTOF HT | c_elegans | HT + C. elegans (2025-10). Phosphoproteomics. |
| 11 | **PXD059519** | timsTOF SCP | human | SCP + human (2026-02). Stem cell screening. **First SCP.** |
| 12 | **PXD070363** | timsTOF SCP | mouse | SCP + mouse (2025-12). Endosomal proteomics. |
| 13 | **PXD063326** | timsTOF fleX | mouse | fleX + mouse (2026-01). Serum IgG analysis. **First fleX.** |
| 14 | **PXD050818** | timsTOF fleX | human | fleX + human (2025-11). HEK293 ATP6AP1 knockout. |

---

## TIER 4: ORGANISM EXPANSION — Add new species

| # | Accession | Instrument | Organism | Why |
|---|-----------|-----------|----------|-----|
| 15 | **PXD050342** | timsTOF Pro | rat | Rat lysine acetylation (2026-02). Most recent rat DDA. |
| 16 | **PXD066353** | timsTOF Pro 2 | rat | Rat microglia (2026-02). Pro 2 + rat. |
| 17 | **PXD068118** | timsTOF Pro 2 | arabidopsis | Arabidopsis PCNA (2025-12). Pro 2 + plant. |
| 18 | **PXD065109** | timsTOF Pro | zebrafish | Zebrafish endothelium (2025-06). Best recent zebrafish DDA. |
| 19 | **PXD068410** | timsTOF Pro | bovine | Bovine milk peptides (2026-01). Most recent bovine. |

---

## TIER 5: DEPTH — Large human/mouse datasets for training mass

| # | Accession | Instrument | Organism | Why |
|---|-----------|-----------|----------|-----|
| 20 | **PXD069989** | timsTOF Pro 2 | human | 32 leukemia patient proteomes. Pro 2, clinical, large. |
| 21 | **PXD070354** | timsTOF HT | human | Chromatin architecture. HT, epigenetics. |
| 22 | **PXD056515** | timsTOF HT | human | VHL mutations. HT, renal cell carcinoma. |
| 23 | **PXD071487** | timsTOF Pro 2 | mouse | Brain SIRT1 metabolomics. Pro 2, neuroscience. |
| 24 | **PXD065562** | timsTOF SCP | mouse | Neutrophil histone clipping. SCP, hypoxia. |
| 25 | **PXD062684** | timsTOF fleX | mouse | Spatial tissue proteomics (laser capture). fleX, spatial. |

---

## Recommended Processing Order

**Phase A — Instruments + Benchmarks (datasets 3-8):**
Adds timsTOF Pro 2, Ultra, HT. Gets benchmark ground truth from MaxQuant/Ionmob.

**Phase B — Full instrument coverage (datasets 9-14):**
Adds timsTOF SCP, fleX. Multiple organisms per instrument.

**Phase C — Organism expansion (datasets 15-19):**
Adds rat, arabidopsis, zebrafish, bovine.

**Phase D — Training depth (datasets 20-25):**
Large clinical/tissue datasets for model training mass.

---

## After adding all tiers, coverage would be:

**Instruments:** timsTOF Pro, Pro 2, HT, SCP, fleX, Ultra (6/7 models)
**Organisms:** human, mouse, yeast, ecoli, drosophila, c_elegans, pig, rat, arabidopsis, zebrafish, bovine (11 species)
**Total datasets:** ~25 (manageable for Phase 1)

---

## Notes

- File counts and sizes not yet validated (run `discover_pride.py --stage3` for Bruker file check)
- Some datasets may contain DIA runs mixed with DDA — verify per-file after download
- PXD019086 already covers human/yeast/ecoli/drosophila/c_elegans on timsTOF Pro
- timsTOF Ultra has only 2 DDA datasets on all of PRIDE — both are must-haves
- The `discover_pride.py` script can auto-update `config.yaml` with `--update-config`
