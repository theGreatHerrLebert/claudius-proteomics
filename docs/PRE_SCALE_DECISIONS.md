# Pre-Scale Decisions (Before Large Dataset Runs)

Date: 2026-02-20
Purpose: lock critical scope and implementation decisions before spending substantial compute and analyst time.

## Why This Exists

Running many public datasets without hard contracts creates expensive ambiguity:

- compute cost is burned on data that may be unusable for training
- silent partial failures poison downstream confidence
- schema and semantics drift make results hard to reproduce

This memo defines decisions to make now, with recommended defaults.

## Decision Summary (Recommended Defaults)

| # | Decision | Recommended Default | Must Decide Today |
|---|---|---|---|
| 1 | Product boundary | Dual artifacts: `archive-all` + `train-ready snapshot` | Yes |
| 2 | Failure policy | Fail-fast by default; degraded mode opt-in only | Yes |
| 3 | Pipeline semantics | Single consistent contract (raw-anchored or search-only, not mixed) | Yes |
| 4 | Snapshot gates | Hard pass/fail gates before training eligibility | Yes |
| 5 | Label policy | Engine-aware confidence policy, not `n_engines` alone | Yes |
| 6 | Initial scope | timsTOF + DDA + standard mods only for wave 1 | Yes |
| 7 | Regression set | 3-5 golden datasets covering edge cases | Yes |
| 8 | Provenance contract | Mandatory checksums + versions + config hash + commit hash | Yes |
| 9 | Schema versioning | Semantic schema version + compatibility tests + migrations | Yes |
| 10 | Storage strategy | Hot query table + cold archive for heavy blobs/intermediates | Yes |
| 11 | Throughput guardrails | Pilot budget and failure-rate thresholds before full queueing | Yes |
| 12 | Publishability rule | No silent degradation + all gates pass + regression pass | Yes |

## Detailed Decisions

## 1) Product Boundary

Decision:

- keep two explicit outputs:
- `archive-all`: broad, minimally filtered scientific archive
- `train-ready snapshot`: strict, gated dataset for model training

Why:

- preserves scientific disagreement information
- prevents low-confidence archival data from contaminating model labels

## 2) Failure Policy

Decision:

- default behavior is fail-fast for production
- degraded operation requires explicit flags and is clearly marked in manifests

Why:

- avoids silent partial successes
- keeps dataset status machine-checkable

## 3) Pipeline Semantics Freeze

Decision:

- choose one and enforce it everywhere:
- Option A: raw-anchored matching (requires ordering and code consistency)
- Option B: search-level merge first, raw integration later

Rule:

- docs, code, and metrics must reflect the same contract

## 4) Snapshot Quality Gates (Hard)

Decision:

- define pass/fail gates required for `train-ready snapshot`

Minimum gate categories:

- engine execution completeness
- contamination/confidence policy per engine
- overlap and consensus sanity
- extraction completeness and blob readability
- schema/provenance completeness

## 5) Label Policy for Training

Decision:

- do not treat `n_engines` as sole confidence proxy

Require:

- engine-specific quality fields
- agreement tier
- coordinate consistency checks
- explicit handling of disagreement strata

## 6) Scope for Wave 1

Decision:

- first scale wave limited to:
- timsTOF
- DDA mode
- standard enzymatic/modification profiles

Defer:

- special digestion modes
- unusual PTM-heavy workflows
- borderline metadata quality datasets

## 7) Golden Dataset Suite

Decision:

- establish fixed regression set before scaling

Recommended composition:

- single-organism canonical dataset
- multi-organism dataset
- multi-instrument or heterogeneous metadata dataset
- one intentionally messy public dataset

Use:

- every major pipeline change must pass this suite

## 8) Provenance Contract

Decision:

- require the following in each final manifest:
- engine versions and executable/container digests
- FASTA identifiers and checksums
- config hash
- code commit hash
- timestamps and run mode flags

Why:

- enables true reproducibility and auditability

## 9) Schema and Compatibility Policy

Decision:

- semantic schema versioning for core parquet outputs
- explicit migration scripts for breaking changes
- compatibility tests for dashboard + training loaders

Why:

- prevents downstream breakage from silent schema drift

## 10) Storage and Data Lifecycle

Decision:

- split storage into:
- hot tier: compact query-focused tables
- cold tier: blobs, intermediates, heavyweight artifacts

Why:

- controls costs and improves query performance at scale

## 11) Throughput and Budget Guardrails

Decision:

- require pilot benchmarking before broad queue submission

Pilot outputs to record:

- CPU-hours per accession
- storage growth per accession
- failure/retry rates
- median wall-clock by step

Gate:

- no full-scale rollout until pilot meets targets

## 12) Publishability Rule

Decision:

- accession/snapshot is publishable only if:
- no silent degradation flags
- all hard gates pass
- golden suite invariants pass

## Concrete Go/No-Go Criteria

## Go

- high-priority remediation items (from `docs/REMEDIATION_CHECKLIST.md`) are completed
- hard snapshot gates implemented and enforced
- provenance manifest contract implemented
- pilot run of golden datasets passes with acceptable cost/failure metrics

## No-Go

- any silent partial-success path remains
- pipeline semantics still contradictory between docs and code
- no reproducible provenance/checksum chain
- golden dataset regressions unresolved

## Immediate Next Steps

1. Ratify the 12 decisions above (owner + date).
2. Implement high-priority remediation items first.
3. Freeze schema/provenance contract.
4. Run pilot on golden suite and record metrics.
5. Reassess scale decision with measured pilot data.
