# ECCB 2026 poster — LLM-assisted scientific infrastructure PLAN

**Status:** v2, 2026-08-24. Revises v1 (same day) after an independent Codex
review (`TO_ECCB_POSTER.codex-review.md`, gpt-5.6-terra) and — more
importantly — after **verifying the review's central objection against actual
run artifacts**, which showed both v1 and the review to be wrong. For ECCB 2026
submission #1015 (Posters track, submitted 2026-04-20).

**What changed from v1:**

1. **Thesis reframed** — the causal claim ("LLM assistance relocates the work")
   is not supported: there is no non-LLM baseline, no tool-use log, no commit
   attribution. Replaced with an amplification claim (§1).
2. **"Verification debt" renamed and demoted** — it is deployment–VCS
   divergence, an illustration, not a metric. Replaced by two auditable
   measures (§5).
3. **The HLA instance was materially wrong in v1 and is now measured, not
   asserted** (§4.1). It is also now the poster's strongest panel, because the
   ground truth turned out to be *heterogeneity*, which is measurable.
4. **A second, independent defect found while checking** — the FragPipe enzyme
   override is a silent no-op (§4.1b).
5. **Novelty positioned against prior work** rather than claimed outright (§6).
6. **cu-ims demoted** from co-equal system to a single inset (§3).
7. All sentences the review flagged as overclaiming are removed or softened
   (§10 records them, so the edit is auditable).

---

## 0. Repos, paths, artifacts

| system | local path | GitHub | state |
|---|---|---|---|
| **claudius-proteomics** | `/home/administrator/Documents/promotion/claudius-proteomics` | `theGreatHerrLebert/claudius-proteomics` | `feat/hf-corpus` @ `53b538c`; working tree ~58 uncommitted entries |
| **evident** | `/scratch/TMAlign/evident` | `theGreatHerrLebert/evident` | `docs/overview-examples` @ `93092b4` |
| **cu-ims-primitives** | `/scratch/ims/cu-ims-primitives` | `theGreatHerrLebert/cu-ims-primitives` | `feat/mbi-ms1-only` @ `9c86d15` |

Key files for a reviewer:

- `evident/OVERVIEW.md`, `evident/concepts/typed-trust.md` — typed trust.
- `evident/cases/cu-ims-primitives.md` — EVIDENT already treats cu-ims as a case.
- `cu-ims-primitives/docs/primitives/TEST_MATRIX.md` — four-tier parity
  taxonomy (note: parts are marked "to be created" — prospective, not applied).
- `runner/engines/fragpipe_job.py`, `scripts/run_fragpipe.py`,
  `runner/steps/step2_search.py`, `config/config.mogon.yaml` (`hla` profile) —
  the §4.1 configuration path.
- `scripts/cluster/_process_runner.sh` — the provenance record (§5).
- `scripts/audit_pool_size_license.py` — §4.3, including its own defect.
- Cluster (not reachable off-site): MOGON-NHR
  `/lustre/project/ki-proanagi/dateschn/data/processed/<ACC>/*/fragpipe_output/fragpipe.workflow`
  — **the rendered-configuration artifacts that settle §4.1.**

---

## 1. Thesis

v1 claimed LLM assistance *relocates* work from implementation to verification.
That is a causal claim about a counterfactual, and nothing in these repos
supports it: there is no non-LLM control, no time accounting, and no
attribution of commits to model assistance. Dropped.

**v2 thesis:**

> Fast, single-author, AI-assisted development amplifies verification and
> provenance gaps at **configuration, deployment, and external-service
> boundaries** — gaps that unit tests and successful-looking runs do not
> expose, and that are invisible to code review by humans *and* by models.

The final clause is new in v2 and is earned by §4.4: during the review of this
very plan, an independent model reviewer read the relevant code, cited line
numbers, and reached a confident wrong conclusion — as had the plan's author.
Only the run artifacts settled it.

---

## 2. What exists, quantitatively

Measured 2026-08-24 from the three git repos.

| | claudius-proteomics | cu-ims-primitives | evident |
|---|---|---|---|
| commits | 180 (2026-01-27 → 05-20) | 393 in 21 days (04-14 → 05-05) | 154 (04-24 → 06-06) |
| tracked files | 736 | 441 | 234 |
| code LOC | 42,292 (109 py) | 40,649 (45 CUDA/C++, 109 py) | 42,187 (33 rs, 69 py) |
| prose LOC | 7,273 (166 md) | 7,047 (36 md) | 15,794 (54 md) |
| test files | 15 | 29 | 90 |

**727 commits, ~125k LOC, three systems, ~4 months, one author.** Report these
as *scale context*, not as evidence for the thesis — file counts do not
establish causality (v1 claimed they did; removed).

Delivered output, to establish this is operating infrastructure:

- **HF corpus v0.2** — 58 datasets, 33,717,207 precursor rows, 504,097,214 b/y
  fragment rows, 100% CC0.
- **cu-ims** — 17,814 peptides @1% FDR vs DiaNN 2.3's 50,482; reported as behind.
- **EVIDENT** — `typed-trust` engine + READ MCP server, `evident-agent` + EXEC
  MCP server, Claude/Codex driver. Orchestrator **sketched only**.

---

## 3. Spine and scope

**Spine: claudius-proteomics.** It matches the submitted abstract and it is
where the measurable findings are. **EVIDENT is the proposed remedy** panel.
**cu-ims is cut to a single inset** inside the oracle figure — it is off-abstract,
it enlarges the claim surface, and its parity taxonomy is partly prospective.

**O1 (still open):** is the submitted abstract editable? If yes, EVIDENT can be
named in it. If not, the poster stays recognisably within the submitted text.

---

## 4. The empirical core

### 4.1 Configuration non-conformance — the HLA case (measured)

**This is the poster's strongest panel and v1 got it wrong.**

v1 asserted "every HLA dataset was searched tryptic on the FragPipe side while
Sage ran nonspecific". The Codex review challenged it and asserted the opposite
(that digestion *was* nonspecific and only the base workflow differed). Both are
wrong. The rendered `fragpipe.workflow` artifacts from real runs show:

| dataset | `search_enzyme_name_1` | `num_enzyme_termini` | digest length | verdict |
|---|---|---|---|---|
| PXD046535 | `nonspecific` | 0 | 7–25 | conforms to `hla` profile |
| PXD042316 | `nonspecific` | 0 | 7–25 | conforms |
| PXD046755 | `trypsin` | 2 | 7–50 | **generic tryptic; mod_profile not applied at all** |

**The finding is heterogeneity, not uniform breakage:** the same pipeline
produced different search configurations for datasets of the same class, and
nothing in the system detected it. That is a better result than v1's claim
because it is *measurable* rather than anecdotal.

Mechanism, as far as the code supports:

- Until `53b538c` (2026-08-24), `FragPipeJob` never passed the profile's
  `fragpipe_workflow`, so a run through that path could not select
  `Nonspecific-HLA.workflow`.
- Yet two of three runs *did* use it (`num_enzyme_termini=0`,
  `search_enzyme_cut_1=@` exist only in that base workflow). So those runs
  reached it by a **different invocation path** — direct `run_fragpipe.py
  --workflow`, or a per-dataset script under `scripts/cluster/process/`.
- **Multiple invocation paths with different effective configuration is the
  actual defect.** Not a wrong constant.

**Open question O5:** identify the second path and confirm it. Until then, the
mechanism is *Judged*; the artifact table above is *Verified*.

### 4.1b The silent no-op override (found while checking 4.1)

`scripts/run_fragpipe.py:132` rewrites workflow lines matching
`msfragger.search_enzyme_name=`. FragPipe workflow files key the enzyme as
**`msfragger.search_enzyme_name_1=`** (multi-enzyme format, FragPipe ≥20). The
prefix never matches, the loop replaces nothing, and **nothing checks that the
replacement applied**. The `--enzyme` override is therefore inert.

By contrast `msfragger.digest_{min,max}_length` (lines 215–216) use exact keys
and *do* apply — which is why the conforming runs show 7–25 and PXD046755 shows
the LFQ-MBR default 7–50.

This is a clean, self-contained example of the failure class: an override that
silently does nothing, in a system where "the run succeeded" carries no
information about whether the intended configuration was used.

### 4.1c A workflow name is not a stable configuration identifier

Found while fixing 4.1b. The *same* named base workflow ships different enzyme
semantics in different FragPipe versions:

| | `search_enzyme_name_1` | `cut_1` | `nocut_1` | digestion |
|---|---|---|---|---|
| FragPipe 23.1 `LFQ-MBR` (local) | `stricttrypsin` | `KR` | *(empty)* | cleaves K/R **including** before P |
| FragPipe 24 `LFQ-MBR` (cluster) | `trypsin` | `KR` | `P` | does **not** cleave before P |

So "we used LFQ-MBR" does not identify the search that ran. Reproducing the
corpus with a locally installed FragPipe 23.1 would silently apply a different
digestion rule than the cluster used. This strengthens §5.1: conformance has to
be measured on the **rendered** workflow, and the corpus provenance should
record the engine version alongside it.

**A5** — check whether `DATASET_CARD.md` / `SCHEMA.md` state a FragPipe version;
if not, that is a reproducibility gap in a shipped CC0 resource.

### 4.2 Deployment–VCS divergence (illustration, **not** a metric)

v1 called this "verification debt" and presented four rows as if they were a
result. Renamed and demoted: "verification debt" already denotes unreviewed
generated code, and deploy-to-commit lag does not measure verification at all —
a long lag can be deliberate batching or untrusted experimental code, and a
zero-lag commit can be entirely unverified.

Retained as **illustration** of how far a deployed system can drift from any
committed state:

| change | live on cluster | committed | lag |
|---|---|---|---|
| `sage_overrides` | 2026-05-27 | 2026-08-24 | 89 d |
| `skip_msbooster` | 2026-06-01 | 2026-08-24 | 84 d |
| blob-resolver fix | 2026-07-23 | 2026-08-24 | 32 d |
| `fragpipe_workflow` pass-through | not deployed until 2026-08-24 | 2026-08-24 | — |

**Stated confounds** (say these on the board): dates come from file mtimes and
`.bak-*` names, which establish neither first execution nor which outputs used
which version; n=4 is not a distribution; the author's branch strategy is a
confound.

### 4.3 Silent nulls at external-service boundaries

PRIDE's **v2 `files/byProject` returns HTTP 200 with an empty body**. An audit
read that as "no raw files" and zeroed **191 of 544 datasets (35%)**. Caught
only by 20 controls where two independent sizing methods must agree (19/20 to
the byte, 1 differing by 8.6%); v3 + `fileCategory == RAW` recovered 188 of 191
and 24.4 TB.

**Admission that belongs on the poster:** `scripts/audit_pool_size_license.py`
— the tool written to *detect* this — has the same defect. `get()` returns
`None` on failure and the pagination loop treats that as end-of-pages, so a
fetch failure is indistinguishable from an empty listing. Flagged by the Codex
review. **Action A1: fix, and commit the raw control table and archived
responses** so the instance is Verified rather than asserted.

### 4.4 The review instance (meta, and the reason for the thesis's last clause)

Reviewing this plan, an independent model read `fragpipe_job.py`,
`run_fragpipe.py`, `step2_search.py` and `config.mogon.yaml`, cited line
numbers, and concluded the enzyme override reached MSFragger. It does not
(§4.1b). The plan's author had independently concluded the opposite error from
a code comment — a comment that is itself inaccurate and still in the file.

Neither careful human reading nor careful model reading resolved it. The
rendered artifact did, in one command. **This is the poster's argument
demonstrated on its own production**, and it costs nothing to include because
the evidence is already in hand.

---

## 5. Two sound metrics (replacing "verification debt")

Both are computable from what already exists on the cluster.

1. **Configuration conformance** — fraction of runs whose *rendered* engine
   configuration matches the profile the config resolves to. Directly
   computable by parsing every
   `data/processed/<ACC>/*/fragpipe_output/fragpipe.workflow` and comparing
   against the `mod_profile` in `config.mogon.yaml`. §4.1 is a 3-dataset
   preview; **the full sweep is the replacement for v1's figure 3.**
2. **Execution provenance completeness** — fraction of runs carrying immutable
   code/config/input hashes. Current provenance
   (`scripts/cluster/_process_runner.sh`) records git SHA, a `working_tree_dirty`
   boolean, and a config SHA-256. The dirty flag proves ambiguity but does not
   identify the diff, so affected runs are **not reproducibly attributable from
   the recorded provenance alone** (v1 said "unanswerable" — too absolute).

---

## 6. Novelty, positioned honestly

Prior art covers more than v1 admitted: **W3C PROV** and **RO-Crate** /
**Workflow Run RO-Crate** model workflow-run provenance; **FAIR / FAIR4RS**
cover research-software reuse; **MLflow / DVC** cover parameter, artifact and
lineage tracking. The poster must cite these rather than imply a vacuum.

What plausibly remains:

- A **concrete, executable case study** of configuration non-conformance in
  operating scientific infrastructure, with rendered-artifact evidence — as
  opposed to a provenance *model* or *standard*.
- **Typed trust**: Verified / Judged / Absent as distinct types with **Absent as
  a first-class result**, plus the binding rule that **synthesis is
  deterministic and calls no model**. This is a mechanism, not a checklist.
- The **four-tier GPU parity taxonomy** — described as a project contract, not
  an established contribution (it is partly prospective and unevaluated).

**O6:** does typed trust survive contact with the differential-testing and
metamorphic-testing literature, or is it a repackaging? Worth a lit check before
printing.

---

## 7. Figures (five panels)

1. **The problem** — claudius at scale: corpus numbers, dual-engine design,
   what the pipeline promises.
2. **Configuration conformance** ⭐ — the §4.1 artifact table, then the full
   sweep across all processed datasets. Replaces v1's "verification debt
   timeline". This is the money figure *because it is measured*.
3. **The failure classes** — silent no-op override (§4.1b), silent null at the
   API boundary (§4.3), deployment drift (§4.2, labelled as illustration).
4. **The oracle question** — what do you check against when there is no ground
   truth? claudius (cross-engine agreement — correlated evidence, *not* an
   oracle), cu-ims inset (simulator truth + tiered CPU/GPU parity), EVIDENT
   (the manifest as oracle).
5. **Typed trust, worked** — one real TrustReport: a claim, Verified evidence,
   a Judged line, an Absent line.

**Device:** tag each claim on the poster with Verified / Judged / Absent. §4.1
is the natural showcase — artifact table Verified, mechanism Judged, second
invocation path Absent until O5 closes.

---

## 8. Explicitly exploratory

- **n = 1, no control**; author = developer = analyst = subject.
- **No LLM attribution** in the commit record — the "AI-assisted" framing is
  contextual, not measured.
- The full conformance sweep (§5.1) **is not done yet**; §4.1 is 3 datasets.
- **EVIDENT's orchestrator is sketched, not built.**
- cu-ims is behind SOTA and its parity matrix is partly prospective.
- Generalisation beyond one developer, one domain, one toolchain is unestablished.

---

## 9. Hook

Lead with §4.1: *three HLA datasets, same pipeline, same config file — two
searched nonspecifically, one tryptic, and nothing noticed.* Then §4.4: neither
the author nor an independent model reviewer could settle it by reading the
code; the rendered artifact settled it immediately. Then hand the audience typed
trust as the proposed answer.

---

## 10. Claims removed or softened from v1 (audit trail)

- "The bottleneck moves, and the new bottleneck has no established tooling" —
  removed; PROV/RO-Crate/MLflow/DVC contradict it.
- "The thesis is visible in the file counts" — removed.
- "All are real, dated, and verifiable" — softened; cluster artifacts are not
  third-party reachable.
- "FragPipe silently fell back to its default tryptic workflow" / "wrong enzyme
  model for months" / "every HLA dataset" — **replaced with the measured table**
  in §4.1.
- "which code produced this data is unanswerable for the whole corpus" →
  "not reproducibly attributable from the recorded provenance alone".
- "cross-engine agreement as a mutual oracle" → "correlated evidence, not an
  oracle".
- "The four-tier parity taxonomy is a contribution in its own right" → "a
  project contract".
- "No equivalent published" — removed entirely.

---

## 11. Actions and open questions

- **A1** — fix the fetch-failure/empty-listing conflation in
  `scripts/audit_pool_size_license.py`; commit it plus the control table and
  archived responses (§4.3).
- **A2** — run the full configuration-conformance sweep over every rendered
  `fragpipe.workflow` on the cluster (§5.1). This is figure 2 and the abstract's
  promised analysis.
- **A3** — fix the inaccurate comment in `runner/engines/fragpipe_job.py:90-95`
  (it asserts the tryptic-fallback story that §4.1 disproves). It is the
  proximate cause of v1's error and should not survive into the poster period.
- **A4** — fix or remove the inert `--enzyme` override in `run_fragpipe.py`
  (§4.1b), with a check that the replacement applied.
- **O1** — is the submitted abstract editable? (§3)
- **O2** — how much cluster-side deploy history is recoverable at all? (§4.2)
- **O5** — which invocation path selected `Nonspecific-HLA.workflow` for the two
  conforming datasets? (§4.1)
- **O6** — does typed trust survive a differential/metamorphic-testing lit
  check? (§6)
- **O7** — publishing §4.1 admits a shipped corpus has non-conforming runs.
  Right call for credibility, or reprocess first and show the fix?
