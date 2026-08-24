# A2 — Configuration-conformance sweep (full)

**Date:** 2026-08-24. **Scope:** every rendered `fragpipe.workflow` under
`data/processed/<ACC>/*/fragpipe_output/` on MOGON (stale `.bak-*` run dirs
excluded). Compares the *rendered* MSFragger configuration against the
`mod_profile` that `config/config.mogon.yaml` resolves for each accession.
Read-only. Script: `a2_conformance.py` (parses `msfragger.table.var-mods`, the
enzyme/termini/digest keys). This is the poster's **figure 2** — the full sweep
that extends the 3-dataset preview in `TO_ECCB_POSTER.md` §4.1.

## Headline

| metric | value |
|---|---|
| rendered workflows | 143 (128 datasets); `.bak*`/`.failed*` excluded |
| **conforming workflows** | **99 / 143 = 69.2%** |
| **conforming datasets** | **99 / 128 = 77.3%** |
| non-conformance among no-profile (generic tryptic) datasets | **0 / 90** |
| non-conformance among PTM/HLA-profile datasets | **29 / 38** |

Reproduce: `python scripts/analysis/a2_conformance.py --config config/config.mogon.yaml --root <data>/processed`

**All non-conformance is concentrated in PTM/HLA-profile datasets.** The 91
generic-tryptic corpus datasets render exactly as expected.

## Two distinct mechanisms

### 1. PTM variable mod never enabled in FragPipe (25 datasets, the dominant class)
The rendered `table.var-mods` for every PTM dataset enables only the FragPipe
LFQ-MBR default — `Oxidation/M` and `N-term Acetyl` — and the profile's target
PTM is absent or set `false`:

| profile | non-conf / total | missing/disabled PTM |
|---|---|---|
| acyl_lactyl | 11/11 | Lactyl 72.0211 / K (absent) |
| acyl_succinyl | 1/1 | Succinyl 100.016 / K (absent) |
| acyl_crotonyl | 1/1 | Crotonyl 68.0262 / K (absent) |
| acyl_malonyl | 1/1 | Malonyl 86.0004 / K (absent) |
| phospho | 5/5 | Phospho 79.9663 / STY (present but `false`) |
| ubiquitin (+ubiq) | 5/5 | GlyGly 114.0429 / K (absent) |
| multi_ptm (PXD042416) | 1/1 (15 subgroup workflows) | all target PTMs absent — even the synthetic kit's per-PTM subgroups only have Ox+N-term-Acetyl |

Directly verified in-file for PXD050906/PXD040514 (lactyl), PXD048158 (succinyl),
PXD042416/Kmod_Succinyl: the target mass appears nowhere as an enabled entry.

### 2. HLA rendered generic-tryptic instead of nonspecific (4 datasets — heterogeneity)
`hla` profile sets `enzyme_override: nonspecific`, `fragpipe_workflow:
Nonspecific-HLA`, digest 7–25. **8 of 12 HLA datasets render correctly**
(nonspecific, termini 0); **4 render generic tryptic** (termini 2):
PXD054248, PXD055547, PXD058376, PXD064980. Same pipeline, same class, two
different rendered configurations — the handoff's "heterogeneity, not uniform
breakage", now measured across the full set.

## Caveats (state these on the poster)
- **Mass-only PTM matching** conflates isobaric sites: the `acetyl` profile
  (Acetyl 42.0106 on **K**) is scored conforming because N-term Acetyl 42.0106
  is default-enabled — the K-site is not actually enabled. Affects 1 dataset
  (PXD050342). A site-aware check would reclassify it non-conforming, i.e. the
  true PTM-class non-conformance is marginally *higher*, not lower.
- **Intent unresolved (the key open question):** this is the *FragPipe* side
  only. The pipeline is dual-engine — the PTM may be carried on the **Sage**
  side (`mod_profile.sage_overrides`) while FragPipe runs generic. If so the
  finding is "the two engines search different modifications for the same
  dataset", not "the PTM was never searched". Settling this needs the rendered
  Sage configs — the natural §4.1-style next verification (an O5 analogue).
- Dates/counts are from the current cluster state; `.bak-*` run dirs excluded.

## Sage-side check — the intent question, SETTLED (2026-08-24)

Same comparison against rendered `data/processed/<ACC>/*/sage/sage_config.json`
(`database.variable_mods`, `database.enzyme`). 153 configs / 132 datasets.

| | FragPipe | Sage |
|---|---|---|
| conforming datasets | 99/128 = **77.3%** | 128/132 = **97.0%** |
| acyl (lactyl/succinyl/crotonyl/malonyl) | **0/15** (PTM absent) | **15/15** (PTM present) |
| HLA | 8/12 (4 render generic tryptic) | 13/13 |
| ubiquitin | 0/5 (GlyGly absent) | 4/4 |
| phospho | 0/5 (Phospho disabled) | 5/7 |

**Directly verified:** PXD059168 (acyl_lactyl) Sage config has
`variable_mods: {"M":[15.994915], "K":[72.021129]}` — Sage searched Lactyl-K;
the same dataset's FragPipe workflow enabled only Ox/M + N-term-Acetyl.

**Conclusion (reframes figure 2):** the PTM is *not* globally missing — the
profile's `variable_modifications` reach **Sage** but never render into the
**FragPipe** workflow. So the finding is not "the PTM was never searched" but:

> **The same `mod_profile` produced two engines searching different modification
> spaces** — Sage PTM-aware, FragPipe PTM-blind — for ~25 PTM datasets, and the
> pipeline's cross-engine agreement step therefore compared a PTM-aware search
> against a PTM-blind one. A real defect in the FragPipe config path, not a
> design choice.

Residual Sage gaps (small, genuine): PXD021789 + PXD037501 (phospho PTM absent
on Sage too), PXD031371 (ubiq GlyGly absent). PXD042416 (multi_ptm synthetic
kit) is mis-scored on **both** engines by the aggregate-profile comparison — its
per-PTM subgroups each target one mod; exclude it or score site/subgroup-aware.
