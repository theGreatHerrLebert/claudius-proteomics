# TO: HF corpus v0.3 expansion — PTM/HLA gap-fillers

**Status:** plan, 2026-06-12. Follows the v0.2 publish
(`theGreatHerrLebert/timstof-dda-pasef-cc0`). Targeted, CC0-only expansion of
the *thin* workflows — deepen, don't brute-force ([[curated-cohort-strategy]]).

## Scope

From the 2026-06-12 discovery (`corpus_expansion_cc0_2026-06-12.tsv`): **115 new
CC0 timsTOF DDA(-PASEF) datasets**, all license-verified, in the workflows the
corpus is thinnest on. Skip the ~485 redundant human/mouse tryptic.

| workflow | new CC0 | corpus today | mod_profile |
|---|---:|---|---|
| phospho | 35 | ~8 | `phospho` ✅ exists |
| **acyl** (lactyl/succinyl/crotonyl/malonyl) | **24** | ~2 | **NEW — §2** |
| acetyl | 18 | ~5 | `acetyl` ✅ |
| methyl | 11 | ~4 | `methyl` ✅ (mono/di/tri) |
| ubiq (K-GG) | 11 | ~4 | `ubiquitin` ✅ |
| HLA | 9 | ~6 | no-enzyme ✅ (`TO_CORPUS_DESIGN` §2.1) |
| glyco | 7 | ~2 | §3 (deglyco vs intact — decide per dataset) |

**Headline value: the 24 acyl datasets** — lysine lactylation / succinylation /
crotonylation across diverse organisms (barley, turnip, virus, *Staph*,
*Magnaporthe*, guinea pig, mouse, human). A near-absent mod class + organism
diversity in one stroke. The phospho/acetyl/ubiq/methyl reuse existing profiles.

## 2. NEW acyl mod_profiles (add to `config/config.mogon.yaml` `mod_profiles:`)

Each dataset is enriched for ONE acyl type. Masses are canonical
(monoisotopic, on K); Sage searches on mass+residue (the unimod_id is metadata).
Crotonyl id 1363 verified via sagepy; succinyl/malonyl/lactyl use canonical
UniMod ids — confirm against the lab UniMod build before the run.

```yaml
  acyl_lactyl:
    description: "K-lactylation (Kla) enrichment"
    fixed_modifications:
      - {unimod_id: 4, name: "Carbamidomethyl", mass: 57.021464, residues: ["C"]}
    variable_modifications:
      - {unimod_id: 35, name: "Oxidation", mass: 15.994915, residues: ["M"]}
      - {unimod_id: 1926, name: "Lactyl", mass: 72.021129, residues: ["K"]}
    max_variable_mods: 3
  acyl_succinyl:
    description: "K-succinylation (Ksucc) enrichment"
    fixed_modifications:
      - {unimod_id: 4, name: "Carbamidomethyl", mass: 57.021464, residues: ["C"]}
    variable_modifications:
      - {unimod_id: 35, name: "Oxidation", mass: 15.994915, residues: ["M"]}
      - {unimod_id: 64, name: "Succinyl", mass: 100.016044, residues: ["K"]}
    max_variable_mods: 3
  acyl_crotonyl:
    description: "K-crotonylation (Kcr) enrichment"
    fixed_modifications:
      - {unimod_id: 4, name: "Carbamidomethyl", mass: 57.021464, residues: ["C"]}
    variable_modifications:
      - {unimod_id: 35, name: "Oxidation", mass: 15.994915, residues: ["M"]}
      - {unimod_id: 1363, name: "Crotonyl", mass: 68.026215, residues: ["K"]}
    max_variable_mods: 3
  acyl_malonyl:
    description: "K-malonylation (Kmal) enrichment"
    fixed_modifications:
      - {unimod_id: 4, name: "Carbamidomethyl", mass: 57.021464, residues: ["C"]}
    variable_modifications:
      - {unimod_id: 35, name: "Oxidation", mass: 15.994915, residues: ["M"]}
      - {unimod_id: 747, name: "Malonyl", mass: 86.000394, residues: ["K"]}
    max_variable_mods: 3
```

Per-dataset assignment of the four acyl types comes from the title (the
discovery already tagged `acyl`; sub-type is in the title: "lactyl"/"succinyl"/
"crotonyl"/"malonyl"). Add a `dataset_metadata.<acc>.mod_profile` entry each.

## 3. Glyco (7) — decide per dataset

Per `TO_CORPUS_DESIGN` §2.1: deglyco/PNGase-F → `N:0.9840` deamidation marker +
N-X-S/T motif filter; intact glycopeptide → needs glycan-mass search (the
modification encoder only handles glycans with a UniMod id). **Recommend: admit
only deglyco evidence in v0.3; defer intact glyco.**

## 4. Pipeline (reuses the existing machinery)

Per candidate, the standard flow (gated on MOGON):
1. **Download** `.d` from PRIDE (`scripts/download_pride.py`).
2. **Search** Sage + FragPipe with the workflow's `mod_profile`
   (`runner/steps/step2_search.py`); HLA → no-enzyme.
3. **Extract** `raw_features` + blobs (the tier1/tier3 inputs).
4. **Audit** against `TO_CORPUS_DESIGN` gates (volume ≥10k, mod-localization
   ≥floor, CE `[0,1]`, dual-engine agreement where both ran).
5. **Build** tier1 + tier3 (the v0.2 builders — `_pick`, `rt_aligned`,
   `_clean_modseq` all already in).
6. **Re-aggregate** to `hf_corpus/v0.3/` (the parameterized aggregator) over the
   union of v0.2's 58 + the new passers → bump the HF repo to tag `v0.3`.

## 5. Risks / gates

- **Mod-localization quality** is the make-or-break for acyl (low-stoichiometry
  PTMs) — enforce the `delta_best ≥ floor` audit gate (`TO_CORPUS_DESIGN` §4.6);
  drop datasets where site ID is unreliable (degrades the intensity signal).
- **Compression-type-1** / readability checks as usual (`[[sage-compression-type1]]`).
- **Cluster-bound**: 115 datasets × (download + dual-engine search + extract) is
  a real campaign; prioritize the **24 acyl + 9 HLA** first (highest diversity
  gain), then phospho/acetyl/ubiq/methyl.
- **Sage build** with the new mods needs the [[sage-build-memory-optimization]]
  fork; verify the acyl masses don't blow the fragment index for K-heavy search.

## 6. Suggested first wave (highest diversity / dataset)

The 24 acyl (4 sub-types, multi-organism) + 9 HLA = 33 datasets — the biggest
per-dataset diversity gain. Run these first as v0.3-wave-1; phospho/acetyl/ubiq/
methyl/glyco as wave-2. Full candidate list:
`corpus_expansion_cc0_2026-06-12.tsv`.
