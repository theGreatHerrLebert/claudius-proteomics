# TO_V03_EXPANSION.md — Codex review (gpt-5.5, 2026-06-12, code-grounded)

**Review Findings**

1. **Critical: lactyl UniMod ID is wrong.**  
   `+72.021129` on K is correct, but `UNIMOD:1926` is methylglyoxal-derived carboxyethyllysine, classified as an artefact. Lysine lactylation is **UNIMOD:2114**.  
   The other entries are correct:

   | Modification | Mass | Correct ID |
   |---|---:|---:|
   | Lactyl-K | 72.021129 | **2114**, not 1926 |
   | Succinyl-K | 100.016044 | 64 |
   | Crotonyl-K | 68.026215 | 1363 |
   | Malonyl-K | 86.000394 | 747 |

   Sources: [UniMod 1926](https://www.unimod.org/modifications_view.php?editid1=1926), [2114](https://www.unimod.org/modifications_view.php?editid1=2114), [64](https://www.unimod.org/modifications_view.php?editid1=64), [1363](https://www.unimod.org/modifications_view.php?editid1=1363), [747](https://www.unimod.org/modifications_view.php?editid1=747).

   “Sage searches mass+residue, so ID is metadata” is only true inside Sage. It is unsafe corpus-wide. Sage config generation discards IDs ([sage_job.py](/home/administrator/Documents/promotion/claudius-proteomics/runner/engines/sage_job.py:117)), while sequence canonicalization infers UniMod from a coarse integer-mass table that currently lacks all four acyl masses ([sequence_utils.py](/home/administrator/Documents/promotion/claudius-proteomics/scripts/sequence_utils.py:26)). Consequently Sage/FragPipe output will retain raw `[+72.021129]`, etc., or be inconsistently encoded, while DIA-NN consumes the supplied ID. Fix the ID and add exact mass-plus-residue mappings and tests before searching.

2. **High: the proposed localization gate does not test localization.**  
   `delta_next ≥ 5` is explicitly described as a “mod-localization proxy” in [corpus_v2_audit_phase1.md](/home/administrator/Documents/promotion/claudius-proteomics/corpus_v2_audit_phase1.md:5). A top-versus-next peptide score gap does not distinguish two positional isomers of the same modified sequence. The canonical parsers also do not retain PTMProphet site probabilities. Do not claim localized acyl-K from this gate. Require PTMProphet/site-determining-ion evidence or an explicit localization score, plus PTM-specific PSM FDR. Otherwise label records as modification-bearing but not confidently site-localized.

3. **High: digestion settings must be PTM-aware, but semi-tryptic is not automatically required.**  
   Acyl-K blocks tryptic cleavage, so a modified internal K consumes a missed-cleavage allowance under an ordinary tryptic digest. Two missed cleavages permit at most two such blocked sites plus any incomplete digestion. That is likely restrictive for enriched multi-acyl peptides. Use **3–4 missed cleavages**, ideally matching each publication, and test identification yield/search cost on one dataset per acyl type. Semi-tryptic searching is justified only where sample preparation or the original analysis used semi-specific digestion; blocked K alone does not create a non-tryptic terminus.

   This does not affect HLA: nonspecific digestion has no meaningful missed-cleavage concept. However, Sage’s enzyme override is implemented, while FragPipe only receives the enzyme name; missed-cleavage and terminus parameters are not propagated. Furthermore, `fragpipe_workflow: Nonspecific-HLA` in [config.mogon.yaml](/home/administrator/Documents/promotion/claudius-proteomics/config/config.mogon.yaml:811) is never consumed; FragPipe always defaults to `LFQ-MBR` ([run_fragpipe.py](/home/administrator/Documents/promotion/claudius-proteomics/scripts/run_fragpipe.py:317)). Therefore “HLA → no-enzyme dual-engine search” is not established.

4. **Medium: the minimal acyl profile is reasonable, but must be dataset-specific.**  
   Cam(C) fixed + Ox(M) + one K-acyl is a good narrow baseline. `max_variable_mods=3` is adequate for one acyl plus oxidation, but potentially truncates multi-acyl peptides; use 4 for enriched acyl searches after a pilot. Protein N-terminal acetylation is optional and probably not worth the search-space cost unless reported. More important: verify alkylation reagent, labeling, digestion, and original modifications per dataset. Malonyl also has a documented CO₂ neutral loss, which may reduce localization sensitivity if the scoring model ignores it.

5. **Medium: deglyco-only is defensible, but deamidation is not glycosylation proof.**  
   Search `Deamidated(N)` rather than unrestricted N/Q for the glyco marker, enforce `N-X-S/T, X≠P`, and require actual PNGase-F treatment plus enrichment or ^18O evidence. Ordinary chemical deamidation can pass the motif filter. Intact datasets should be deferred, not searched as deglyco.

6. **Campaign-stopping metadata risks.**  
   Title-only classification has obvious false positives: `PXD036950` is “glycolytic,” not glycoproteomics; `PXD053580` is Chlamydomonas CTI1 biology, not HLA; `PXD057704` uses “MHCC” in a cell-line context; `PXD041207` says “acylated segment” but provides no acyl PTM subtype. `PXD040481` is malonylation but tagged acetyl. Require protocol/file-level confirmation and a generated per-dataset manifest reviewed before submission. Also pilot one dataset end-to-end, including canonical modified-sequence output and PTM yield, before launching 115 datasets.
