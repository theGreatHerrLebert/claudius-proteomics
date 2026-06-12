# TO: HF corpus v0.3 expansion — PTM/HLA gap-fillers

**Status:** plan v2, 2026-06-12 (Codex-reviewed, `TO_V03_EXPANSION.codex-review.md` — 6 catches folded in: lactyl id 2114, register acyl masses, no localization claim, MC 3-4, FragPipe-HLA bug, title-vs-protocol). Follows the v0.2 publish
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

Each dataset is enriched for ONE acyl type. Masses are canonical (monoisotopic,
on K). **Correct UniMod ids (Codex-verified): Lactyl 2114, Succinyl 64,
Crotonyl 1363, Malonyl 747.** ⚠️ "Sage searches mass+residue so the id is just
metadata" is **only true inside Sage** — it is unsafe corpus-wide:
`runner/engines/sage_job.py` discards the id, and sequence canonicalization
(`scripts/sequence_utils.py`) infers UniMod from a coarse integer-mass table
that **currently lacks all four acyl masses**, so Sage/FragPipe output would
keep raw `[+72.021129]` (or encode inconsistently) and DIA-NN consumes the id.
→ **§2.1 is a hard prerequisite, not optional.**

### 2.1 PREREQUISITE — register the acyl masses before any search

Add exact (mass, residue) → UniMod mappings for all four acyl mods to the
canonicalization table in `scripts/sequence_utils.py` (and any DIA-NN id map),
with a round-trip test (`[+72.021129]K → [UNIMOD:2114]K`). Without this the
`modified_sequence` / `sequence_normalized` columns are wrong for the whole
acyl cohort — and the tier3 fragment matcher's mass→UNIMOD converter (already
small, `[[imspy-fragment-ion-extension]]`) must know them too, or acyl peptides
get skipped at Tier-3.

```yaml
  acyl_lactyl:
    description: "K-lactylation (Kla) enrichment"
    fixed_modifications:
      - {unimod_id: 4, name: "Carbamidomethyl", mass: 57.021464, residues: ["C"]}
    variable_modifications:
      - {unimod_id: 35, name: "Oxidation", mass: 15.994915, residues: ["M"]}
      - {unimod_id: 2114, name: "Lactyl", mass: 72.021129, residues: ["K"]}  # 2114 = L-lactyl (1926 is a methylglyoxal artefact — Codex)
    max_variable_mods: 4   # 4 (not 3) so multi-acyl peptides aren't truncated (Codex)
  acyl_succinyl:
    description: "K-succinylation (Ksucc) enrichment"
    fixed_modifications:
      - {unimod_id: 4, name: "Carbamidomethyl", mass: 57.021464, residues: ["C"]}
    variable_modifications:
      - {unimod_id: 35, name: "Oxidation", mass: 15.994915, residues: ["M"]}
      - {unimod_id: 64, name: "Succinyl", mass: 100.016044, residues: ["K"]}
    max_variable_mods: 4   # 4 (not 3) so multi-acyl peptides aren't truncated (Codex)
  acyl_crotonyl:
    description: "K-crotonylation (Kcr) enrichment"
    fixed_modifications:
      - {unimod_id: 4, name: "Carbamidomethyl", mass: 57.021464, residues: ["C"]}
    variable_modifications:
      - {unimod_id: 35, name: "Oxidation", mass: 15.994915, residues: ["M"]}
      - {unimod_id: 1363, name: "Crotonyl", mass: 68.026215, residues: ["K"]}
    max_variable_mods: 4   # 4 (not 3) so multi-acyl peptides aren't truncated (Codex)
  acyl_malonyl:
    description: "K-malonylation (Kmal) enrichment"
    fixed_modifications:
      - {unimod_id: 4, name: "Carbamidomethyl", mass: 57.021464, residues: ["C"]}
    variable_modifications:
      - {unimod_id: 35, name: "Oxidation", mass: 15.994915, residues: ["M"]}
      - {unimod_id: 747, name: "Malonyl", mass: 86.000394, residues: ["K"]}
    max_variable_mods: 4   # 4 (not 3) so multi-acyl peptides aren't truncated (Codex)
```

Per-dataset assignment of the four acyl types comes from the title (the
discovery already tagged `acyl`; sub-type is in the title: "lactyl"/"succinyl"/
"crotonyl"/"malonyl"). Add a `dataset_metadata.<acc>.mod_profile` entry each.

## 3. Glyco (7) — decide per dataset

Per `TO_CORPUS_DESIGN` §2.1: deglyco/PNGase-F → `Deamidated(N)` marker
(UniMod:7, +0.984016, **N only — not unrestricted N/Q**) + `N-X-S/T, X≠P` motif
filter, AND require evidence of actual **PNGase-F treatment + glyco enrichment**
(ideally ¹⁸O) — chemical deamidation alone passes the motif filter, so the
marker is not glyco proof (Codex). **Recommend: admit only verified-deglyco
evidence in v0.3; defer intact glyco** (the encoder only handles UniMod-id'd
glycans).

## 4. Pipeline (reuses the existing machinery)

Per candidate, the standard flow (gated on MOGON):
1. **Download** `.d` from PRIDE (`scripts/download_pride.py`).
2. **Search** Sage + FragPipe with the workflow's `mod_profile`
   (`runner/steps/step2_search.py`). **Acyl: missed_cleavages 3–4** (not 2) —
   acyl-K blocks tryptic cleavage and consumes the MC allowance; pilot the
   yield/cost per acyl type. (Semi-tryptic only if the original study used it —
   blocked K alone doesn't create a non-tryptic terminus.) ⚠️ **HLA dual-engine
   is currently broken on the FragPipe side**: `fragpipe_workflow: Nonspecific-HLA`
   in `config.mogon.yaml` is **never consumed** — `run_fragpipe.py` always runs
   `LFQ-MBR` (tryptic), and only the enzyme *name* (not MC/terminus) is passed.
   FIX `run_fragpipe.py` to consume the workflow, **or** treat HLA as Sage-only
   in v0.3 (drop the dual-engine claim for HLA).
3. **Extract** `raw_features` + blobs (the tier1/tier3 inputs).
4. **Audit** against `TO_CORPUS_DESIGN` gates (volume ≥10k, mod-localization
   ≥floor, CE `[0,1]`, dual-engine agreement where both ran).
5. **Build** tier1 + tier3 (the v0.2 builders — `_pick`, `rt_aligned`,
   `_clean_modseq` all already in).
6. **Re-aggregate** to `hf_corpus/v0.3/` (the parameterized aggregator) over the
   union of v0.2's 58 + the new passers → bump the HF repo to tag `v0.3`.

## 5. Risks / gates

- **Site localization is NOT established by the `delta_best` gate** (Codex):
  that's a top-vs-next *peptide-score* gap, not a positional-isomer test, and
  the canonical parsers don't retain PTMProphet site probabilities. So either
  add real localization evidence (PTMProphet / site-determining ions + a
  PTM-specific PSM FDR), or **label acyl records as modification-bearing but
  NOT confidently site-localized** in the corpus (honest, and fine for intensity
  training). Don't claim localized acyl-K.
- ⚠️ **Title-only classification has false positives** (Codex examples:
  PXD036950 "glycolytic"≠glyco; PXD053580 Chlamydomonas≠HLA; PXD057704 "MHCC"
  cell-line≠MHC; PXD041207 "acylated segment" w/o subtype; **PXD040481 is
  malonyl but my tagger called it acetyl**). → **Confirm workflow + acyl subtype
  from the PRIDE protocol/sample-metadata (not the title), generate a per-dataset
  manifest, and human-review it before submission.**
- **Pilot ONE dataset per acyl type end-to-end first** (search → extract →
  tier1/tier3 → check canonical `modified_sequence` + PTM yield) before
  launching all 115. Confirms §2.1 encoding + the MC setting actually work.
- Malonyl has a documented **CO₂ neutral loss** — may reduce localization
  sensitivity if the score ignores it.
- **Compression-type-1** / readability checks as usual (`[[sage-compression-type1]]`).
- **Cluster-bound**: prioritize the verified **acyl + HLA** subset first
  (highest diversity gain), then phospho/acetyl/ubiq/methyl.
- **Sage build** with the new mods needs the [[sage-build-memory-optimization]]
  fork; verify the acyl masses don't blow the fragment index for K-heavy search.

## 6. Suggested first wave (highest diversity / dataset)

The 24 acyl (4 sub-types, multi-organism) + 9 HLA = 33 datasets — the biggest
per-dataset diversity gain. Run these first as v0.3-wave-1; phospho/acetyl/ubiq/
methyl/glyco as wave-2. Full candidate list:
`corpus_expansion_cc0_2026-06-12.tsv`.

## 7. Existing-HLA reprocess (the FragPipe no-enzyme bug side-effect)

The `fragpipe_job.py` fix (now passes `--workflow Nonspecific-HLA`) means every
HLA dataset previously FragPipe-searched ran **tryptic LFQ-MBR** on no-enzyme
HLA peptides — so their FragPipe IDs / `n_engines==2` / fragpipe_* columns are
unreliable. **Verified 2026-06-12: `Nonspecific-HLA.workflow` ships** in
FragPipe 24.0 (`engines/fragpipe-run/workflows/`), so the fix resolves.

**Affected in the v0.2 corpus (mod_profile=hla): PXD026463, PXD038273** (only 2 —
the other recovered per-group datasets PXD035301/PXD050342 are phospho/acetyl,
not HLA). Broader already-processed HLA pool that benefits if pulled into v0.3:
PXD042316, PXD058376, PXD046370, PXD051880 (+ the 9 new HLA candidates process
fresh = already correct).

Reprocess (folds into the v0.3 rebuild):
1. **Audit first**: for the 2 in-corpus HLA, report current `n_engines==2` rate
   + FragPipe HLA peptide-length distribution. Expect FragPipe contributed
   little (tryptic on ~9-mer no-enzyme peptides) → low blast radius, but confirm.
2. **Re-run FragPipe** (`Nonspecific-HLA`, now wired) on PXD026463 + PXD038273
   (+ any HLA added to v0.3). Per-dataset check: **Cam(C) only if the sample was
   alkylated** — many HLA preps are not (`TO_CORPUS_DESIGN` §2.1; FragPipe ships
   a `Nonspecific-HLA-C57` variant for the alkylated case).
3. **Re-merge** dual-engine (step5) → corrected `precursor_index` (real HLA
   `n_engines`/agreement).
4. **Rebuild** tier1+tier3 (v0.2 builders) → re-aggregate to v0.3.

Net: v0.3 HLA = 2 corrected + (optionally) 4 already-processed re-FragPipe'd + 9
new = up to ~15 properly no-enzyme dual-engine HLA datasets.
