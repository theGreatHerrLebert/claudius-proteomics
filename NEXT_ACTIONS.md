# NEXT_ACTIONS — resume point, 2026-08-24

Handoff written at end of session 2026-08-24. Everything below is either
verified today or explicitly marked as unverified. Related memory:
`v03-campaign-plan`, `quota-wall-2026-07`, `nonhuman-diversity-gap`,
`pride-resweep-2026-07`, `no-poll-mogon-nhr`.

---

## 0. TL;DR — where we stand

| thread | state |
|---|---|
| acyl-23 resume | ✅ **diagnosed + preflight deployed**; 3 jobs running, 2 datasets need a decision (§1) |
| repo ↔ cluster reconciliation | ✅ done; **2 commits unpushed**, 1 fix **not deployed** (§2) |
| ECCB poster | ✅ plan v2 + Codex review landed; **A2 is the next work** (§3) |
| non-human corpus expansion | ✅ audited, 544/544 CC0, wave-1 shortlist ready (§4) |

---

## 1. Array status — DIAGNOSED AND FIXED (2026-08-24, end of session)

**Running when the session closed** (these survive independently, just check them):

```
569993_4    PXD046631   RUNNING 53min   extracted=yes, in feature fitting — healthy
569993_20   PXD075276   RUNNING 53min   still downloading — healthy
570106_21   PXD078736   RUNNING          re-fired after repair (job cp-acyl23r2)
```

`ssh mogon-nhr 'squeue -u dateschn'` — one shot, never in a loop.

### What went wrong: three unrelated causes, all now handled

| task | dataset | root cause | status |
|---|---|---|---|
| 21 | PXD078736 | zip double-nested (`X.d/X.d/`) + `__MACOSX`; `validate_d.py` saw no `analysis.tdf` and quarantined all 6 | ✅ **repaired + re-running** |
| 19 | PXD058676 | organism `platycodon_grandiflo` has **no UniProt reference proteome** (0 proteomes, 155 UniProtKB entries, 1 reviewed) and the deposit ships no FASTA. Raw data itself is fine (8 valid `.d`) | ⛔ **needs a decision** |
| 22 | PXD045802 | deposit's RAW zips contain **no timsTOF `.d`** — they unpack to `data.ser` + `p0/*.apl` (MaxQuant output). No `analysis.tdf` anywhere | ⛔ **drop from cohort** |

### The systemic fix (this is what stops future waste)

New `scripts/cluster/preflight_raw.py`, wired into `_process_runner.sh` between
download and search (`--steps 1` → preflight → full run). **Deployed to MOGON**
(backup `_process_runner.sh.bak-pre-preflight-20260824`).

It *repairs* the packaging variants — flattens nested `.d`, drops `__MACOSX`,
adds a missing `.d` suffix to dirs holding an `analysis.tdf`, restores anything
in `_quarantine/` that becomes valid — then *gates*: at least one valid `.d`,
and a resolvable FASTA for the organism. Exit `2` = no usable raw, `3` = no
FASTA source.

Verified against all three real failures on the cluster:

```
PXD078736 -> exit=0   6 __MACOSX removed, 6 .d restored from quarantine, FASTA UP000005640
PXD058676 -> exit=3   "organism 'platycodon_grandiflo' has no FASTA source"
PXD045802 -> exit=2   "no usable .d ... deposit may not ship timsTOF .d"
```

Those three failures each burned ~35 min of a 128-core node. They now fail in
seconds, before the engines start.

### Decisions waiting for you

- **PXD058676** — search against a related Campanulaceae proteome, or drop it?
  Using a wrong-species database has real FDR consequences, so this is a
  scientific call, not a config tweak. Until decided it is gated out
  automatically. Cohort would go 23 → 22.
- **PXD045802** — recommend dropping from `acyl_cohort23.txt`; the deposit has
  no raw timsTOF data to process. Cohort → 21.

## 2. Repo ↔ cluster state

**Local repo** `/home/administrator/Documents/promotion/claudius-proteomics`,
branch `feat/hf-corpus`:

```
4539533  fragpipe: make the enzyme override apply, or fail loudly   ← unpushed
53b538c  fragpipe: pass the profile's base workflow and skip_msbooster through  ← unpushed
c11079f  v0.3: land the MOGON campaign config, acyl-23 cohort and array driver
936643b  sage: honour per-profile database.* overrides
e1c898b  HF corpus: resolve blobs relative to each raw_features source dir
```

- **2 commits unpushed** (`53b538c`, `4539533`). `git push origin feat/hf-corpus`
  when ready.
- The main checkout sits on `fix/fragpipe-prot-fdr` with **62 dirty entries**
  (long-lived operational edits, pre-existing).
- **MOGON checkout** is on `main` @ `cfa0d6d` (2026-05-26) with 11 modified
  tracked files. Only `runner/engines/fragpipe_job.py` was deployed today
  (backup `.bak-pre-workflow-20260824`).

**⚠️ NOT deployed to the cluster, deliberately:** `scripts/run_fragpipe.py`
(commit `4539533`). It adds a `ValueError` when the requested enzyme
contradicts the base workflow. Only the `hla` mod_profile sets
`enzyme_override` (config line 975), so the acyl array cannot trip it — but
deploying a behaviour change mid-campaign was left as a deliberate decision.
Deploy when the array is idle.

---

## 3. ECCB poster — A2 is the next real work

- `TO_ECCB_POSTER.md` — **v2**, rewritten after the Codex review.
- `TO_ECCB_POSTER.codex-review.md` — the raw review (the actual content is at
  the end of the file, after the echoed prompt).

**A2 — configuration-conformance sweep.** This is both poster figure 2 and the
analysis the submitted abstract promised. Read-only, safe to run while jobs are
going.

> For every `data/processed/<ACC>/*/fragpipe_output/fragpipe.workflow` on
> MOGON, parse `msfragger.search_enzyme_name_1`, `search_enzyme_cut_1`,
> `num_enzyme_termini`, `digest_min_length`, `digest_max_length`, and the
> variable/fixed mods. Compare against the `mod_profile` that
> `config/config.mogon.yaml` resolves for that accession. Report the
> **fraction of runs whose rendered configuration matches the profile**, and
> list every non-conforming dataset.

Preview from 3 datasets (already verified today, §4.1 of the plan):

| dataset | enzyme | termini | digest | verdict |
|---|---|---|---|---|
| PXD046535 | nonspecific | 0 | 7–25 | conforms |
| PXD042316 | nonspecific | 0 | 7–25 | conforms |
| PXD046755 | trypsin | 2 | 7–50 | **non-conforming** |

Remaining actions from the plan: **A5** (does `DATASET_CARD.md` record a
FragPipe version? if not it's a reproducibility gap — the same workflow *name*
means different digestion in FragPipe 23.1 vs 24), **O1** (is the submitted
abstract still editable? decides the poster spine), **O5** (which invocation
path selected `Nonspecific-HLA.workflow` for the two conforming datasets?).

A1, A3, A4 are **done** (commit `4539533`).

---

## 4. Corpus expansion — ready to go, nothing blocking

Artifacts created today, all in `data/discovery/`:

- `nonhuman_pool_sized_2026-08-24.tsv` — **544 non-human DDA candidates,
  544/544 CC0**, total 41.3 TB, median 21 GB.
- `nonhuman_wave1_candidates_2026-08-24.txt` — **22-dataset breadth wave,
  660 GB**, 2 each across mouse/rat/arabidopsis/E. coli/S. aureus/pig/bovine/
  rice/yeast/S. pombe/*Flavobacterium*. Includes a timsTOF **Ultra** set (corpus
  has none), the Kbhb rice histone dataset, and the K48-polyUb rat pair.
- `acyl23_license_check_2026-08-24.tsv` — acyl cohort **23/23 CC0**, subtypes
  match config `mod_profile` 23/23.
- `new_since_2026-07-23.txt` — 70 new datasets from today's re-sweep, only
  22 DDA-compatible, 9 non-human.
- Catalog backups: `*.bak-pre-resweep-20260824`.

**Why this matters:** the corpus is **97.6% human by row**. 457 of the 544
non-human candidates were already catalogued back in February. A month of new
PRIDE yields ~22 usable datasets vs 544 sitting unused — **diversity is a
selection problem, not a discovery problem.**

⚠️ `timstof_catalog.yaml` is still the Feb-11 baseline (1773 entries) — the
July and August runs only regenerated the CSV. `--incremental` therefore
reports "new since Feb"; month deltas must be computed by differencing CSVs.

---

## 5. Gotchas worth not rediscovering

- **PRIDE v2 `files/byProject` is dead** — HTTP 200 with an empty body, which
  reads as "no files" rather than as an error. Use **v3**
  `/projects/{acc}/files` (`fileSizeBytes`, paginated 1000/page) and size on
  `fileCategory == RAW`, not filename globs (globs miss `raw.zip`, `*.tar.gz`,
  `Topology.rar` — they silently zeroed 191/544 datasets).
- **`BATCH_QC_HALT` is scoped to the array job id**, so a stale flag never
  blocks a new submission. Nothing to clear.
- **The runner always passes `--resume`** — a partial dataset resumes just by
  being re-fired.
- **The QC gate reports on *other* recent datasets**, not the one that just ran:
  a task can fail while the log's last line says `qc_gate: batch OK`. Do not
  read that as success.
- `du -sh $B/*` on the cluster takes ~47 s. A previous attempt timed out at
  900 s because it did three traversals piped through `sort` (which buffers,
  so it printed nothing at all). Disk: dateschn ≈ 20.2 TB, group out of grace.
- **Two claims I made this session were wrong** and are corrected in the plan
  and in memory: HLA runs are *not* uniformly tryptic (they are heterogeneous),
  and the 0 GB audit rows were *not* converted-only deposits (they were
  archive-wrapped raw). Both were caught by checking artifacts rather than
  reading code.

---

## 6. Suggested order next session

1. `squeue` once — did tasks 4/20 finish? (§1)
2. Diagnose the step-2 failure (§1 diagnostics 1–3), re-fire the 3 failed tasks.
3. **A2 conformance sweep** (§3) — the poster's figure 2.
4. Push the 2 commits; deploy `run_fragpipe.py` once the cluster is idle (§2).
5. Then: tier1/tier3 build array for the acyl increment → v0.3, and/or the
   22-dataset non-human wave (§4).
