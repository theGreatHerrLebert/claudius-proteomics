# ECCB 2026 poster — PLAN (v5, reframed)

**Status:** v5, 2026-08-25. **Reframed spine** (author's call): from a
proteomics-infrastructure case study to a general question — *how do you trust
code, configs, and results an AI wrote that you didn't author and don't fully
understand?* — with proteomics as the **sole evaluated case** and proteon/cu-ims
as **design probes**. Then hardened against a third Codex review
(`TO_ECCB_POSTER.codex-review-v3.md`). Supersedes v4 (proteomics-spine). Prior
reviews: `TO_ECCB_POSTER.codex-review{,-v2,-v3}.md`. For ECCB 2026 submission
#1015 (Posters track).

**Why reframe.** ECCB is broad computational biology; MS proteomics is niche
there. The general thesis draws the room; the concrete deliverables + one measured
result reward the specialists who stop by. A poster is a conversation-starter, so
an explorative thesis is fine — **but only if the epistemic status is worn openly**
(the V/J/A tagging device), and only if the scope boundary is unmistakable.

**The scope boundary (the single most important thing — Codex v3 (b)).**
> The conformance gate is **evaluated on ONE proteomics configuration-conformance
> problem.** General AI-assisted-software trust is the **motivating question, not a
> demonstrated conclusion.** proteon/cu-ims are design probes; there is no
> comparative study and no claim of improved reliability. n=1.

**What changed from v4 (all from Codex v3):**
1. Thesis is now a *bounded position + one case study*, not a general solution.
2. proteon/cu-ims relabelled **design probes**, not "the pattern generalizes."
3. Broad factual claims softened (§1) — tags constrain rhetoric, not decorate it.
4. The **one contribution** named (§6): executable claim-level evidence-status
   rendering + the deterministic gate. The 4 "patterns" are framing/design, not
   novelty claims.
5. **Oracle corrected:** the finding's oracle is the *rendered-config-vs-profile
   contract*, NOT cross-engine agreement (which is correlated evidence only).
6. Deliverables must be evidence-bearing or framed as "built with this practice"
   (§7), not proof.
7. New **"What is NOT shown + proposed evaluation"** panel (§8) — answers the
   "where's the evaluation?" puncture before it's asked.

---

## A. Abstract (~165 words, broad-first, softened)

**Title:** *Typed Trust: verifying scientific software you didn't write*
(subtitle: a discipline for AI-assisted research software — proteomics, structural
biology & GPU signal-processing as use-cases). Title editable per #1015.

> Researchers increasingly build scientific software with AI assistance, and a
> recurring difficulty is *trusting* code, configurations, and results one did not
> author and does not fully understand. Failures at configuration, deployment, and
> external-service boundaries are not exposed by unit tests or successful-looking
> runs, and — in the case we document — were not resolved by code review, human or
> model.
>
> We present EVIDENT, a framework for typed trust: every assertion records how it
> was established — Verified (a reproducible procedure ran), Judged (an
> interpretation), or Absent (sought, not found) — with a deterministic final
> report that calls no model. It couples this evidence-status rendering to an
> executable conformance gate that checks a rendered configuration against its
> intended contract.
>
> We evaluate the gate on one proteomics case: a single profile silently drove two
> search engines to search different modification spaces on 25 of 128 datasets,
> caught only by the rendered artifact. Two further AI-assisted systems
> (structural biology; GPU signal-processing) are shown as design probes.

Softened vs the punchier draft: dropped "AI writes faster than anyone can verify"
(unevidenced) and "invisible to humans *and* models alike" (→ "in the case we
document … not resolved by code review"). Keep the punchier line only as the
spoken hook, tagged Judged.

---

## 1. Thesis (a bounded position)

> In AI-assisted scientific software, trusting artifacts you did not author is a
> distinct and under-tooled problem. In one documented proteomics case,
> validation against **rendered run artifacts** exposed a configuration failure
> that unit tests, successful runs, and code review (human and model) did not.

Explicitly a **position + one case study**, not a measured general effect. "AI
writes code faster than we can verify it" is the *motivation* (spoken, tagged
Judged), not an asserted finding. n=1, no control, author = developer = analyst =
subject; "AI-assisted" is the setting, not a variable.

---

## 2. The framework — EVIDENT (design, not novelty claims)

Four organizing ideas; presented as the framework's shape, with the actual
contribution isolated in §6.
1. **Typed trust** — Verified / Judged / Absent; deterministic synthesis (verdict
   calls no model).
2. **Claim ⊥ test** — a manifest binds `claim → oracle → tolerance → command →
   artifact`; a test attaches to / is challenged against a claim without rewriting
   it. (Familiar test architecture unless the manifest shows a capability
   conventional test metadata lacks — do not overclaim.)
3. **The oracle pattern** — what to check against with no ground truth:
   reproducible check · cross-tool agreement (**correlated evidence, not an
   oracle**) · simulator truth · declared tolerance · **the intended-config
   contract** (the oracle in §3).
4. **Verification as learning** — unknown agent output → extract claim → verify
   against an oracle → fold the verified fact into the knowledge base. Presented
   as a *workflow*, not a technical contribution.

---

## 3. The evaluated case (proteomics — the one measured anchor)

**Oracle: the rendered configuration vs the resolved `mod_profile` contract**
(NOT cross-engine agreement). On **128 paired datasets** searched by both engines,
one profile produced two different rendered search configurations: **25/128
divergent, one-sided (25 vs 0)** — Sage carried the profile's PTM, FragPipe a
generic default; 4/128 both-fail (Sage not universally conformant). Verified from
rendered configs; engine versions recovered (FragPipe 24.0 / MSFragger 4.4 / Sage
0.15); run-linked (A9 provenance supplement; caveat: 100% `working_tree_dirty`,
config drift). A **deterministic conformance gate** is built, deployed, exits 2
pre-search, emits a V/J/A TrustReport (A2/A6/A9/A10 done, committed). Two further
failure classes (silent no-op override; silent API null zeroing 191/544 datasets)
support the boundary-error thesis. All Verified; the mechanism is Judged; the
invocation path is Absent.

---

## 4. Design probes (NOT generalization — Judged, n=1)

Two more single-author AI-assisted systems illustrate the same *practice*; they
are **not evaluated** and do not demonstrate that the approach generalizes.
- **proteon** — structural biology, Rust, pip-published; declared tolerances,
  oracle setup, a release-tier typed-trust claim on CHARMM19+BALL electrostatics.
- **cu-ims** — GPU diaPASEF; oracle = simulator truth + CPU/CUDA bit-parity;
  honestly behind SOTA on identification depth.
Report as *"systems built with this practice,"* one author, ~250k LOC — scale
context, not evidence.

---

## 5. Concrete deliverables (real; framed as built-with-this-practice)

- **HF corpus** — public 100% CC0 timsTOF reference layer (v0.2: 58 datasets,
  33.7M precursors, 504M b/y fragments; v0.3 acyl-PTM increment building now).
- **timsTOF property predictors** — RT / ion-mobility / intensity / CCS on the
  imspy_predictors hub (models-v0.6.0; best fine-tune SA 0.798, CCS MAE 0.0148,
  RT r 0.918).
To avoid "bolted-on advertisement" (Codex v3 Q4): **render one corpus claim as a
TrustReport** (claim: "corpus rows carry documented extraction provenance"; oracle:
the manifest; artifact: the parquet) so the deliverables sit inside the trust
frame — else demote to a one-line footer. "Real/usable" is not evidence the
framework works.

---

## 6. Contribution & novelty (modest, honest)

**The contribution is narrow and concrete:** *executable, claim-level
evidence-status rendering (Verified/Judged/Absent, deterministic) coupled to a
deterministic conformance gate over rendered configuration.* Demoable; not a new
provenance model, not a new testing paradigm.

Prior art cited, not implied absent: **differential + metamorphic testing**
(no-ground-truth checks — the §3 comparison *is* differential testing on config);
**W3C PROV**, **Workflow Run RO-Crate** (run provenance); **MLflow/DVC** (lineage).
"Verification as learning" and "claim ⊥ test" are framing, not claimed as novel.

---

## 7. Figures (five panels)

1. **Hook + the framework** (EVIDENT: the 4 ideas; V/J/A device defined).
2. **The evaluated case** ⭐ — the paired conformance result (25 vs 0), oracle =
   config-vs-contract. The one measured panel.
3. **Failure classes + provenance** (no-op override; silent API null; 100% dirty).
4. **Design probes** — proteon + cu-ims, labelled not-evaluated.
5. **Verification-as-learning, worked** — the gate run prospectively (exit 2) as a
   TrustReport (Verified / Judged / Absent). Plus the deliverables band.
Device: tag every claim V/J/A; the motivational/general claims are Judged/Absent.

---

## 8. What is NOT shown (put it on the poster — Codex v3 Q5)

- **No evaluation** that EVIDENT reduces undetected failures or improves review.
- **No comparative study**; **no transfer** beyond one author's proteomics pipeline.
- One measured incident; proteon/cu-ims unevaluated; orchestrator sketched.
- **Proposed evaluation protocol** (future): inject known config non-conformances
  across a dataset panel; measure gate detection rate + false positives; compare
  against reviewer (human/model) detection. State this as the next step.

---

## 9. Actions / open questions

- **A2/A6/A9/A10** ✅ done (measured finding, gate, provenance supplement — committed).
- **A-new1** — render one corpus claim as a TrustReport (§5) or footer-demote.
- **A-new2** — update the mockup: relabel proteon/cu-ims as design probes; fix the
  oracle label (config-contract, not cross-engine agreement); soften the two broad
  claims; add the §8 "not shown" box.
- **A4** — deploy the FragPipe render fixes so the gate can go strict (separate track).
- **O1** — title editability for #1015. **O6** — differential/metamorphic lit check.
- **O7** — publishing a shipped corpus has non-conforming runs: credibility vs reprocess-first.

---

## 10. Audit trail (claims removed/softened)

v1→v4 removals preserved. New in **v5**:
- "the pattern generalizes (3 systems)" → **design probes, not generalization** (Q1).
- "AI writes faster than anyone can verify" / "invisible to humans and models" →
  softened to the documented case (Q2).
- 4 patterns as novelty → **framing; one narrow contribution named** (Q3).
- proteomics oracle "cross-engine agreement" → **config-vs-profile contract** (Q3/c);
  cross-engine agreement is correlated evidence, not an oracle.
- deliverables as standalone value → evidence-bearing or footer (Q4).
- added the explicit **"not shown + proposed evaluation"** (Q5).
