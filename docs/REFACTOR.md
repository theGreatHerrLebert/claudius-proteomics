# Refactor Backlog

Date: 2026-02-20
Scope: pipeline + backend + shared libraries
Objective: reduce complexity and drift before large-scale dataset processing.

## Refactor Now (Highest ROI)

1. Canonicalize duplicated matching/normalization modules.
- Problem: diverged copies in `lib/sequence_utils.py` vs `scripts/sequence_utils.py`, and `lib/precursor_matching.py` vs `scripts/precursor_matching.py`.
- Current usage split: runner uses `lib/*` (`runner/steps/step3_stratify.py`), legacy scripts use `scripts/*` (`scripts/build_precursor_index_legacy.py`).
- Action: choose one canonical location and delete/redirect the other.

2. Split Step 3 monolith into focused modules.
- Problem: `runner/steps/step3_stratify.py` is oversized (1665 lines), with `_generate_html_report` alone at 462 lines.
- Action: extract to `matching.py`, `qc.py`, `reporting.py`, and orchestration module.

3. Remove parallel merge/indexing implementations.
- Problem: runner Step3/Step5 path coexists with script stack (`scripts/build_precursor_index.py`, `scripts/build_precursor_index_legacy.py`, `scripts/precursor_merging/*`, `scripts/merge_psm_raw.py` via `rules/merge.smk`).
- Action: pick one production path, mark/remove legacy path.

4. Refactor dashboard backend query hot path.
- Problem: `dashboard/backend/main.py` is large; `list_precursors` is a large endpoint and materializes broad tables before pagination.
- Action: move query logic into service layer and add query-optimized summary pathway.

5. De-duplicate blob parsing logic.
- Problem: `dashboard/backend/main.py` has a local `read_blob` implementation duplicating `sanjose/blob.py`.
- Action: use `sanjose.blob.BlobReader` in API layer.

6. Stream Step 4 output instead of in-memory aggregation.
- Problem: `runner/steps/step4_extract.py` collects all file DataFrames then concatenates.
- Action: incremental Parquet writing (append/chunked writer) to reduce memory pressure.

## Refactor Soon (Medium ROI)

1. Make packaging file collection declarative.
- Problem: `runner/steps/step6_package.py` hard-codes many file patterns and directory rules.
- Action: move include map to config/manifest schema.

2. Remove unused summary helper APIs.
- Problem: `runner/summary.py` includes `create_step1_summary`..`create_step5_summary` and `create_manifest` not used by runner step implementations.
- Action: remove or relocate to tests/examples.

3. Remove dead parser helper.
- Problem: `scripts/engine_parsers/base.py` `_prefix_columns` has no call sites.
- Action: delete helper.

4. Clean stale imports and minor dead code.
- `runner/steps/step4_extract.py` imports unused quality-metric functions.
- `runner/steps/step5_merge.py` imports unused `create_manifest`.
- Action: clean imports and keep modules minimal.

## Stale/Scope-Drift Signals

1. Orchestration split-brain.
- `Snakefile` and `runner/run_dataset.py` both represent high-level execution paths with differing semantics.
- Decision needed: canonical production orchestrator.

2. Legacy artifacts still first-class.
- `scripts/build_precursor_index_legacy.py` and older script flows remain active-looking despite newer runner flow.
- Action: archive or remove explicitly.

3. Checkpoint marker redundancy.
- `.done` marker files in `runner/state.py` are written but not central to scheduling behavior.
- Action: simplify checkpoint model or fully integrate marker checks.

## Suggested Execution Sequence

1. Canonicalize shared modules + remove duplicated code paths.
2. Split Step 3 and dashboard backend into focused modules.
3. Stream extraction output and declarativize packaging.
4. Cleanup pass: dead helpers/imports + legacy archival.
