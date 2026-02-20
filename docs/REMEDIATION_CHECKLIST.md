# San Jose Remediation Checklist

Date: 2026-02-20
Scope: Runner + merge logic + scientific quality gates
Goal: Maximize scientific risk reduction per engineering day

## Ranking Method

Each item is ranked by `risk_reduction_per_day` (High, Medium, Low), using:

- scientific impact if unfixed
- probability of silent failure
- implementation effort

## Prioritized Backlog

| Rank | Item | Risk Reduced | Effort | Risk Reduction / Day |
|---|---|---|---:|---|
| 1 | Fail hard when any enabled engine fails in Step 2 | Prevents silent partial-search datasets being treated as valid | 0.5 day | High |
| 2 | Make Step 4 extraction failures fatal (no silent empty fallback) | Prevents empty/partial raw feature stores from propagating | 0.5 day | High |
| 3 | Make Step 3 per-group errors fail accession run (or explicit degraded mode) | Prevents hidden group-level corruption behind success status | 0.5 day | High |
| 4 | Fix local-data config key mismatch (`local_data` vs `local_datasets`) | Prevents wrong data source and accidental re-downloads | 0.25 day | High |
| 5 | Align pipeline order with design (raw-anchored Step 3) or redesign docs + code consistently | Removes core logic contradiction in matching strategy | 1 day | High |
| 6 | Fix engine-only fallback merge to preserve run-specific identity | Prevents inflated false consensus from sequence-level broadcasting | 1.5 days | High |
| 7 | Treat partial step runs as partial state (not `completed=true`) | Prevents invalid checkpoint semantics/resume behavior | 0.5 day | Medium |
| 8 | Add organism resolver support for sheep aliases/taxon | Prevents sample-group misassignment in configured datasets | 0.25 day | Medium |
| 9 | Ensure Step 4 Sage fragment path resolution is group-correct | Restores expected fragment annotation quality in per-group mode | 0.5 day | Medium |
| 10 | Enforce download file filtering in primary path (`pridepy`) | Avoids test-mode/selection drift and unexpected dataset scope | 0.75 day | Medium |
| 11 | Add minimal training gate for “report all” mode (engine-aware decoy/confidence policy) | Reduces contamination of model-training snapshots | 2 days | Medium |
| 12 | Add integration tests for runner failure semantics and merge invariants | Prevents regressions in core scientific data integrity | 2 days | Medium |

## Action Plan (Execution Order)

1. Implement ranks 1-4 in one patch set.
2. Implement ranks 5-7 in one patch set.
3. Implement ranks 8-10 in one patch set.
4. Implement ranks 11-12 with tests and snapshot-gate docs.

## Detailed Acceptance Criteria

### 1) Step 2 hard-fail policy

- Step 2 returns non-success if any enabled engine reports `status != success`.
- `runner/run_dataset.py` must stop and mark pipeline failed.
- Step summary must include per-engine failure reason.

### 2) Step 4 fatal extraction policy

- Any per-file extraction exception fails Step 4 by default.
- Optional degraded mode must be explicit (`--allow-partial-extraction`) and recorded in summary.

### 3) Step 3 group failure policy

- Any required group failure causes Step 3 failure, unless `--allow-partial-groups` is explicitly set.
- QC manifest must record skipped/failed groups with clear status.

### 4) Local-data config key fix

- Step 1 resolves local paths from `config.local_data[accession]`.
- Backward-compatible alias for `local_datasets` can be kept with deprecation warning.

### 5) Pipeline-order consistency

- Either:
- reorder to `1,2,4,3,5` for true raw anchoring,
- or remove raw-anchored claims and make Step 3 purely search-level.
- Docs and code must match exactly.

### 6) Engine-only fallback integrity

- No fallback merge may join DIA-NN/Sage to FragPipe rows solely by `(sequence,charge)` across runs.
- Preserve run identity (`raw_file`) in all fallback joins.
- Add invariant checks for impossible consensus inflation.

### 7) Partial-run checkpoint semantics

- `state.completed=true` only when full requested workflow is complete.
- Partial execution stores `completed=false` with explicit `steps_requested`.

### 8) Sheep organism support

- Add sheep aliases and taxon mapping in `sample_group_resolver`.
- Add test case using representative ovine-style filenames.

### 9) Step 4 group fragment path fix

- Sage fragment lookup in per-group mode resolves to the current group processed directory.
- Add a regression test with synthetic group directory layout.

### 10) Download filter enforcement

- `max_files` and file patterns are honored identically in both primary and fallback download paths.
- Add deterministic logging: selected files list before download starts.

### 11) Training gate for report-all mode

- Define snapshot gate: required confidence constraints per engine and consensus strata.
- Block training snapshot creation if gate metrics fail.

### 12) Integration tests

- Add tests for:
- engine failure propagation,
- extraction failure propagation,
- fallback merge invariants,
- checkpoint correctness after partial runs.

## Suggested Milestones

| Milestone | Includes | ETA |
|---|---|---:|
| M1 | Items 1-4 | 1-2 days |
| M2 | Items 5-7 | 2-3 days |
| M3 | Items 8-10 | 1-2 days |
| M4 | Items 11-12 | 2-3 days |

## Definition of Done

- No silent-success path for failed/partial scientific outputs.
- Merge logic preserves run-level identity in all modes.
- Docs reflect real execution semantics.
- Core runner invariants covered by automated tests.
