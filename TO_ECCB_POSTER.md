# ECCB 2026 poster — LLM-assisted scientific infrastructure PLAN

**Status:** v4, 2026-08-24. Revises v3 after an independent Codex plan review
(`TO_ECCB_POSTER.codex-review.md`). v4 is a **de-overclaiming pass**: it reframes
the work as a **single-case study** (not a quasi-comparative "empirical study"),
replaces the mismatched A2 denominators with a **paired** engine comparison,
softens "defect"/"amplifies"/"remedy" to what the evidence supports, and reframes
EVIDENT from a claimed remedy into a **conformance gate + evidence-status
reporting** demonstrated prospectively. For ECCB 2026 submission #1015 (Posters
track, submitted 2026-04-20).

**What changed from v3 (all from the Codex review):**

1. **Thesis de-causalised** (§1). "Amplifies" implied a non-LLM baseline we don't
   have. Reframed as: in *this one* project, rendered-artifact validation exposed
   failures the development-time checks did not. Lead with **case study**.
2. **A2 paired** (§4.1). v3 compared FragPipe 99/128 vs Sage 128/132 — different
   sets. v4 reports the **128 datasets run through both engines**: FragPipe 77.3%,
   Sage 96.9%, **25/128 (19.5%) divergent, fully asymmetric** (FragPipe never the
   conformant one).
3. **"Real defect" → "conformance defect vs the (reconstructed) `mod_profile`
   contract"** (§4.1). Claim scoped to what's shown: FragPipe omitted the
   profile-specified target variable PTM; Sage included it.
4. **New caveat: rendered ≠ executed config** (§4.1) — proving execution needs
   command line + tool versions + logs; supplement to archive them (A9).
5. **Abstract de-scaled** (§A): drops "rare setting / team-scale / previously
   required teams"; LOC demoted; structural biology trimmed to "contextual
   breadth" (portfolio stays on the poster, lighter in the abstract).
6. **EVIDENT reframed** (§6, §7 fig 5): from "systematic validation strategy /
   remedy" to a **small executable conformance gate + Verified/Judged/Absent
   reporting**, demonstrated **prospectively** (the gate failing *before* a search
   runs), not by labeling our own past finding.
7. **Novelty** (§6): explicitly cites **differential + metamorphic testing** and
   concedes PROV / Workflow Run RO-Crate already model run-provenance-for-trust.

**v4.1 update (A9 folded in, 2026-08-24):** the rendered→executed linkage is now
**recovered**, closing Codex round-2's top residual. A committed supplement
(`A9_provenance_supplement.tsv`, 150 (accession,group) rows) hashes each rendered
config and links it to the run's git commit, SLURM job, UTC time, and engine
versions (FragPipe 24.0 / MSFragger 4.4 / Sage 0.15.0-beta.2). Two findings are
now quantified, not hedged: **100% of runs had `working_tree_dirty: true`**
(§5 metric #2 = 0% code-hash-complete) and only **12/150** recorded config SHAs
match the current config (drift). Figure 2 is a **run-level** claim; code-pinning
is the stated, measured gap. §4.1, §5, §8, A9 updated below.

**v4.2 update (A6/A10 gate built, 2026-08-24):** the figure-5 conformance gate is
now **built + demonstrated** (`scripts/analysis/conformance_gate.py`) — it exits 2
and **blocks on non-conformance**, emitting an EVIDENT TrustReport
(Verified/Judged/Absent). Shown failing divergent FragPipe configs *pre-search*
(PXD059168 Lactyl absent; PXD058376 HLA→tryptic) and passing the conforming Sage
config. This closes Codex's fig-5 circularity concern: EVIDENT now *prevents*
something, it does not just label a past finding. §7, §8, A6/A10 updated.

---

## A. Revised abstract (editable — draft for submission)

**Title options:**
- (keep, safest) the submitted title — *LLM-Assisted Development of Large-Scale
  Proteomics Infrastructure: An Empirical Study from PRIDE Reanalysis*.
- (sharper, if body-and-title editable) *…: A Case Study of Failure Modes and a
  Conformance Gate*. Avoid "remedy" — the gate is a partial, demonstrated check,
  not a validated remedy.

> In data-intensive fields such as proteomics, the bottleneck is often not
> algorithm development but the construction and operation of scalable data
> infrastructure. We present a **case study** of such infrastructure built by a
> single researcher with LLM assistance: an operational PRIDE timsTOF reanalysis
> pipeline — multi-engine search, raw-signal extraction, and a versioned,
> bias-documented reference layer — that has produced a 58-dataset,
> 33.7-million-precursor, 100% CC0 reference corpus. Two further research-software
> projects, a GPU-native diaPASEF stack and a published Rust
> structural-bioinformatics toolkit, provide contextual breadth.
>
> Analyzing the pipeline through version history and *rendered run artifacts*, we
> identify failure modes not exposed by the development-time checks in use — unit
> tests, successful-looking runs, and (for the configuration case) code review by
> human and model: configuration non-conformance — across 128 datasets searched
> by two engines under one modification profile, the engines' rendered configs
> diverged on 25 (21 a target PTM, 4 the nonspecific-enzyme setting), one engine
> applying the profile and the other a generic default — plus a silent null at an
> external-service boundary and deployment–version-control divergence. These cluster at configuration, deployment, and external-service
> boundaries.
>
> As a partial answer we present a small, executable **conformance gate** that
> checks each engine's rendered configuration against the resolved profile,
> reported through EVIDENT, a typed-trust layer separating Verified, Judged, and
> Absent evidence deterministically. We argue that in this case validation was the
> binding constraint, and that lightweight rendered-artifact checks address a
> concrete part of it.

Open: does ECCB permit editing the **title**? If not, the body edit stands within
the existing title.

---

## 0. Repos, paths, artifacts

| system | domain | local path | GitHub | state |
|---|---|---|---|---|
| **claudius-proteomics** | proteomics reanalysis (spine) | `~/Documents/promotion/claudius-proteomics` | `theGreatHerrLebert/claudius-proteomics` | `feat/hf-corpus`; poster docs + fig-2 committed |
| **proteon** | structural biology | `/scratch/TMAlign/proteon` | `theGreatHerrLebert/proteon` | `main` @ `c66b975`, v0.4.0, pip-published |
| **cu-ims-primitives** | GPU diaPASEF | `/scratch/ims/cu-ims-primitives` | `theGreatHerrLebert/cu-ims-primitives` | `feat/mbi-ms1-only` @ `9c86d15` |
| **evident** | validation meta / reporting | `/scratch/TMAlign/evident` | `theGreatHerrLebert/evident` | `docs/overview-examples` @ `93092b4` |

Key files for a reviewer:

- `scripts/analysis/a2_conformance.py`, `A2_CONFORMANCE_RESULTS.md` — figure 2,
  committed (`065cc82`). The paired analysis is A8 (to fold into the report).
- claudius: `runner/engines/fragpipe_job.py`, `scripts/run_fragpipe.py`,
  `runner/engines/sage_job.py`, `config/config.mogon.yaml` — the §4.1 path.
- `evident/OVERVIEW.md`, `evident/concepts/typed-trust.md` — Verified/Judged/
  Absent; deterministic synthesis. `evident/cases/{proteon,cu-ims-primitives}.md`.
- `proteon/docs/ORACLE_SETUP.md`, `proteon/STABILITY.md` — validation-first
  structure (context for §5).
- Cluster (not off-site reachable): MOGON-NHR
  `.../data/processed/<ACC>/*/{fragpipe_output/fragpipe.workflow,sage/sage_config.json}`
  — the rendered configs. **A9: archive a hashed extraction + exact
  FragPipe/MSFragger/Sage versions as a supplement** (reproducibility).

---

## 1. Thesis (a case study, stated as one)

> In one large, single-author, AI-assisted infrastructure project, validation
> against **rendered run artifacts** exposed configuration, deployment, and
> external-service-boundary failures that the development-time checks — unit
> tests, successful-looking runs, and the code reviews performed here (human and
> model) — did not surface.

**What this is and is not.** This is a **case study / incident analysis** with a
reproducible artifact-level audit. n = 1, no non-LLM control, author =
developer = analyst = subject. "AI-assisted" is the **setting**, not a measured
explanatory variable — we do **not** claim LLM assistance *causes* or *amplifies*
these failures relative to a baseline (that would need a control we don't have).
The defensible, and still interesting, claim is that these failure classes exist,
are artifact-detectable, and were missed by the checks in use.

The "not resolved by human or model review" clause is earned narrowly by §4.4:
one author misreading + one model review, both wrong, both settled by the
rendered artifact. It is a demonstration, not a prevalence claim.

---

## 2. What exists (scale context only)

Measured 2026-08-24 from the git repos. **Scale context, not evidence** — LOC is
sensitive to language, generated/vendor code, tests, and counting tool, so it is
reported with a fixed policy (`git ls-files '*.rs|*.py'` tracked lines, no vendor
exclusion) and **not** used as a headline.

| | claudius | cu-ims | proteon | EVIDENT |
|---|---|---|---|---|
| domain | proteomics | GPU diaPASEF | structural bio | validation/reporting |
| commits | 180 | 393 | 517 | 154 |
| code LOC (policy above) | 42,292 | 40,649 | 167,340 | 42,187 |
| maturity | operational | behind SOTA (research) | v0.4 pip-published | early, tested |

**The headline is operational output, not LOC:**

- **claudius / HF corpus v0.2** — 58 datasets, 33,717,207 precursor rows,
  504,097,214 b/y fragment rows, 100% CC0.
- **proteon** — `pip install proteon`; SASA/DSSP/H-bonds, prep/minimization,
  TM-align, CHARMM19+BALL electrostatics, Vina docking.
- **cu-ims** — 17,814 peptides @1% FDR vs DiaNN 2.3's 50,482 (reported as behind).
- **EVIDENT** — typed-trust engine + READ MCP server, evident-agent + EXEC MCP
  server, Claude/Codex driver. Deterministic orchestrator sketched only.

---

## 3. Spine and scope

**Spine: claudius-proteomics** — the abstract's subject and the *only* system with
a quantified audit. **proteon and cu-ims are contextual breadth**, not evidence of
generalization: the validation conclusion rests on claudius alone (§8). They
appear as a portfolio strip (§7 fig 1) and one oracle inset (fig 4); they do not
carry the argument.

Discipline: if space forces cuts, cu-ims goes first, then proteon detail (keep it
in the strip only). **In the abstract, structural biology is one phrase**
("contextual breadth") so it does not dilute a proteomics-track submission —
while the full portfolio still shows on the poster (author's call).

---

## 4. The empirical core

### 4.1 Configuration non-conformance — paired engine comparison (measured)

**The poster's strongest panel (figure 2).** For datasets searched by *both*
engines, the *rendered* configuration of each is compared against the
`mod_profile` that `config.mogon.yaml` resolves. Reproducible:
`scripts/analysis/a2_conformance.py`.

**Paired set: 128 datasets run through both FragPipe and Sage** (every FragPipe
dataset is paired; 4 Sage-only datasets excluded for the comparison):

| engine | conforming (same 128) | |
|---|---|---|
| **FragPipe** | **99/128 = 77.3%** | all non-conformance is PTM/HLA profiles; 0/90 generic-tryptic fail |
| **Sage** | **124/128 = 96.9%** | |
| **divergent** (Sage conforms, FragPipe does not) | **25/128 = 19.5%** | acyl-lactyl 11, HLA 4, ubiquitin 4, phospho 3, crotonyl/malonyl/succinyl 1 each |
| divergent the other way (FragPipe conforms, Sage not) | **0** | one-sided discordance (25 vs 0) |
| both non-conformant | **4/128** | Sage is not universally conformant |

**The scoped claim (what the artifacts show, and only that):**

> For the tested profiles, on 25 of 128 paired datasets **FragPipe omitted the
> profile-specified target variable PTM (or, for 4 HLA sets, the nonspecific
> enzyme) while Sage included it** — one modification profile, two engines, two
> different rendered search configurations. The pairwise discordance is one-sided
> (25 vs 0) — Sage was never the less-conformant engine on a shared dataset —
> though Sage is itself non-conformant on 4/128 (both marginal rates retained
> above); this is directional discordance, not "Sage always conformed."

Mechanism, verified in-file: PXD059168 Sage `variable_mods {M:[15.99],
K:[72.021129]}` searched Lactyl-K; the same dataset's FragPipe workflow enabled
only Oxidation/M + N-term-Acetyl. FragPipe stores nothing for the target mass.

**Framed as conformance, not "defect."** There was no pre-existing, versioned
contract stating both engines must implement the resolved profile identically;
that contract is **reconstructed post hoc** from the config's intent. So the
claim is *non-conformance relative to the (reconstructed) `mod_profile`
contract*, not a proven bug in intended behavior. (It makes an undocumented
intentional divergence unlikely, but does not prove intent.)

**Rendered→executed linkage — recovered (A9, ✅).** A per-(accession,group)
supplement (`A9_provenance_supplement.tsv`, 150/150 rows linked to a provenance
run) hashes each rendered `fragpipe.workflow` / `sage_config.json` and joins it to
the run's `git_commit`, SLURM `job_id`, UTC timestamp, and engine versions —
**FragPipe 24.0 / MSFragger 4.4** (workflow header) and **Sage 0.15.0-beta.2**
(run json). So figure 2 is a **run-level** conformance claim, not merely a config
audit. **Quantified honest residual:** 100% of the 150 runs carried
`working_tree_dirty: true` (code not bit-pinned) and only 12/150 recorded config
SHAs match the current config (drift); the literal command line is not in the run
json and is not claimed. The conformance verdict itself stands on the rendered
artifacts regardless of code-pinning.

**Other caveats (on the board):**
- **PTM matching is engine-semantic**, not string/table equality — mass +
  tolerance + residue/terminus site; mass-only matching conflates isobaric sites
  (N-term vs K Acetyl → PXD050342 false-conforms; true non-conformance is
  marginally *higher*).
- **Search space is broader than a single variable-mod mass** (fixed mods, max
  var-mods/peptide, length/mass windows, enzyme specificity, localization). The
  claim is scoped to the target PTM and enzyme, not to full search-space identity.
- PXD042416 (synthetic multi-PTM kit) needs subgroup-specific expected sets;
  currently mis-scored by the aggregate comparison.

### 4.1b The silent no-op override (found while checking 4.1)

`scripts/run_fragpipe.py:132` rewrites `msfragger.search_enzyme_name=` lines;
FragPipe keys the enzyme as `msfragger.search_enzyme_name_1=` (multi-enzyme
format, ≥20). The prefix never matches, nothing checks the replacement applied —
the `--enzyme` override is inert. Same class as §4.1: an intended configuration
silently not applied, where "the run succeeded" says nothing about whether the
intended config was used.

### 4.1c A workflow name is not a stable configuration identifier

The same named base workflow ships different enzyme semantics across FragPipe
versions (23.1 `LFQ-MBR` `stricttrypsin`, cleaves before P; 24 `trypsin`, does
not). **A5** — record the FragPipe version in `DATASET_CARD.md`/`SCHEMA.md`; its
absence is a reproducibility gap in a shipped CC0 resource.

### 4.2 Deployment–VCS divergence (illustration, not a metric)

Illustration only, confounds stated (mtimes/`.bak-*` names establish neither
first execution nor which outputs used which version; n=4):

| change | live | committed | lag |
|---|---|---|---|
| `sage_overrides` | 2026-05-27 | 2026-08-24 | 89 d |
| `skip_msbooster` | 2026-06-01 | 2026-08-24 | 84 d |
| blob-resolver fix | 2026-07-23 | 2026-08-24 | 32 d |

### 4.3 Silent nulls at external-service boundaries

PRIDE v2 `files/byProject` returns HTTP 200 with an empty body; an audit read that
as "no raw files" and zeroed **191/544 datasets (35%)**. Caught only by 20
controls where two sizing methods must agree; v3 + `fileCategory == RAW` recovered
188/191 and 24.4 TB. **Admission:** `scripts/audit_pool_size_license.py` — the
detector — has the same defect (a `None` fetch is indistinguishable from an empty
listing). **A1: fix + commit control table and archived responses.**

### 4.4 The review instance (why §1's last clause is narrow)

Reviewing this plan, an independent model read the relevant files, cited line
numbers, and concluded the enzyme override reached MSFragger. It does not (§4.1b).
The author had concluded the opposite error from an inaccurate code comment.
Neither reading resolved it; the rendered artifact did, in one command. A single
demonstration that the artifact settles what code-reading did not — **not** a
claim that human/model review generally fails.

---

## 5. Two candidate metrics

1. **Configuration conformance** — fraction of runs whose *rendered* config
   matches the resolved profile (§4.1, both engines, paired). Figure 2.
2. **Execution provenance completeness** — fraction of runs carrying immutable
   code/config/input hashes. `_process_runner.sh` records a git SHA, a
   `working_tree_dirty` boolean, and a config SHA-256. **Measured (A9, n=150):
   100% of runs are `working_tree_dirty: true`** — by the immutable-code-hash
   criterion, completeness is **0%** — and only **12/150** recorded config SHAs
   match the current config. The dirty flag proves ambiguity without identifying
   the diff. A live, quantified instance of the metric, not a hypothetical.

**proteon as contrast, not proof:** it ships validation-first (`ORACLE_SETUP.md`,
stable/experimental tiers, an EVIDENT release-tier claim on CHARMM19+BALL). The
claudius-vs-proteon contrast (drift found post-hoc vs tolerances declared up
front) is illustrative — it is **not** a second audited data point.

---

## 6. Novelty, positioned honestly

Prior art is real and cited, not implied absent:
- **Provenance:** W3C PROV (provenance used to assess trust); **Workflow Run
  RO-Crate** already models workflow-run inputs/outputs/code/execution provenance.
- **Lineage tooling:** MLflow, DVC.
- **Oracle-poor testing:** **differential testing** (cross-implementation
  agreement) and **metamorphic testing** already target scientific software
  without ground truth — the §4.1 cross-engine comparison *is* differential
  testing applied to config, and must be named as such.

So the contribution is **not** a new provenance model or a new testing paradigm.
It is narrower and honest: **a small, executable, domain-grounded conformance
gate** (rendered engine config vs resolved profile) **plus explicit
evidence-status rendering** (Verified/Judged/Absent, deterministic synthesis).
Enough for a poster if presented at that size.

**O6:** confirm the differential/metamorphic framing survives a real lit check
before printing.

---

## 7. Figures (five panels)

1. **Portfolio + problem** — 3 systems × 3 domains, one author (strip); claudius
   at scale as the detailed case. (Portfolio = context, labeled as such.)
2. **Configuration conformance** ⭐ — the §4.1 **paired** result: 128 datasets,
   FragPipe 77.3% vs Sage 96.9%, **one-sided discordance 25 vs 0** (retain both
   marginal rates). The money figure *because it is measured and paired*.
3. **Failure classes** — silent no-op override (§4.1b), silent API null (§4.3),
   deployment drift (§4.2, labeled illustration).
4. **The oracle question** — cross-engine agreement (correlated evidence, not an
   oracle); **cu-ims inset** (simulator truth + CPU/CUDA bit-parity); proteon
   (declared tolerances); EVIDENT (manifest as oracle).
5. **The conformance gate, run prospectively** ⭐ — **not** a label on our own
   finding. A deterministic check: *resolved profile → rendered config → PASS/FAIL,
   the failing config caught before the search runs*, reported as a TrustReport
   (Verified = criterion results read from the artifact; Judged = the mechanism;
   Absent = the O5 invocation path). **✅ Built + demonstrated**
   (`scripts/analysis/conformance_gate.py`): **exits 2 (blocks) on
   non-conformance** — shown failing FragPipe/PXD059168 (Lactyl absent) and
   FragPipe/PXD058376 (HLA rendered tryptic) while passing the conforming Sage
   config; emits `A6_trustreport_example.json` (typed Verified/Judged/Absent,
   engine version + rendered-config sha256, deterministic). This makes EVIDENT
   operational rather than circular.

**Device:** tag every poster claim Verified / Judged / Absent.

---

## 8. Explicitly exploratory (say these out loud)

- **n = 1, no control, author = subject.** This is a case study.
- **No LLM attribution** in commits — "AI-assisted" is setting, not variable.
- Only **claudius** is audited; proteon/cu-ims are context, not generalization.
- The A2 result links rendered config to runs (A9 ✅); the residual is
  **code-pinning — 100% of runs had a dirty working tree**, so the exact code is
  not bit-reconstructable (itself the §5 metric-#2 finding).
- The fig-5 conformance gate is **built and demonstrated** (A6/A10 ✅,
  `conformance_gate.py`, exit-2-on-non-conformance). It emits EVIDENT-schema
  TrustReports directly; the **full EVIDENT engine integration** (manifest → the
  typed-trust Rust crate → orchestrator) **remains future work** — the crate is
  scaffolded, the orchestrator sketched.
- cu-ims is behind SOTA; its parity matrix is partly prospective.

---

## 9. Hook

Lead with §4.1: *one modification profile, two search engines — Sage searched the
PTM, FragPipe searched a generic default, on 25 of 128 datasets, and the
development-time checks in use did not flag it.* Then §4.4: the rendered artifact settled in one command what code
reading did not. Then offer the conformance gate + typed-trust reporting as a
small, concrete answer — and the portfolio as context that this is one person's
operating infrastructure, not a toy.

---

## 10. Claims removed or softened (audit trail)

v1→v3 removals preserved. New in **v4** (from Codex review):
- "amplifies verification gaps" → case-study framing; no causal/baseline claim.
- "invisible to code review by humans and models" → "not resolved by the reviews
  performed here" (one instance, not a generalization).
- A2 "77.3% vs 97.0%" across **different** sets → **paired** 128-dataset
  comparison (25/128 divergent, asymmetric).
- "a real defect, not a design choice" → "non-conformance vs the reconstructed
  `mod_profile` contract; intent not proven."
- "different modification spaces" (unqualified) → scoped to the target PTM/enzyme,
  with rendered≠executed and search-space-breadth caveats added.
- EVIDENT "systematic validation strategy / remedy" → "conformance gate +
  evidence-status rendering," demonstrated prospectively.
- "~250k LOC / team-scale / previously required dedicated teams" → LOC demoted to
  context with a counting policy; sociological claims dropped.
- Novelty: differential + metamorphic testing now cited; contribution narrowed.

---

## 11. Actions and open questions

- **A1** — fix the fetch-failure/empty-listing conflation in
  `audit_pool_size_license.py`; commit control table + archived responses.
- **A2** — ✅ full conformance sweep, both engines (committed `065cc82`).
- **A8** — ✅ paired analysis (128 datasets; 25/128 divergent, asymmetric). Fold
  into `A2_CONFORMANCE_RESULTS.md` + `a2_conformance.py` (`--paired`).
- **A6/A10** — ✅ done: `scripts/analysis/conformance_gate.py` — deterministic
  pre-search gate, exit 2 on non-conformance, emits an EVIDENT TrustReport
  (Verified/Judged/Absent). Demonstrated failing FragPipe/PXD059168 (Lactyl
  absent) + FragPipe/PXD058376 (HLA→tryptic), passing Sage/PXD059168. Example
  report: `A6_trustreport_example.json`. ✅ **Wired into `_process_runner.sh`**
  (soft: writes a TrustReport per rendered config + warns; `CONFORMANCE_GATE_STRICT=1`
  fails the job once the A4 render bug is fixed — hard-blocking now would halt every
  PTM run). Fuller *pre-engine* gating belongs inside `step2_search.py`; full
  typed-trust crate integration is future.
- **A9** — ✅ done: `scripts/analysis/a9_provenance_supplement.py` +
  `A9_provenance_supplement.tsv` (150 rows) — per-(accession,group) hashed chain
  profile→rendered config→engine versions→run→commit, with conformance verdicts.
  Recovered FragPipe 24.0 / MSFragger 4.4 / Sage 0.15; found 100% working_tree_dirty
  and 12/150 config-match-current. Fig 2 is run-level; code-pinning is the stated gap.
  (Regenerate on-cluster for canonical hashes; local hashes are byte-identical.)
- **A3** — fix the inaccurate comment in `fragpipe_job.py:90-95`.
- **A4** — fix/remove the inert `--enzyme` override (§4.1b), with a check.
- **A5** — record the FragPipe version in `DATASET_CARD.md`/`SCHEMA.md`.
- **O1** — abstract editable; **title-edit permission unconfirmed**.
- **O5** — which invocation path selected `Nonspecific-HLA.workflow` for the
  conforming HLA datasets? (§4.1 Absent line).
- **O6** — confirm the differential/metamorphic-testing lit framing.
- **O7** — publishing §4.1 admits a shipped corpus has non-conforming runs.
  Credibility win, or reprocess first and show the fix?
