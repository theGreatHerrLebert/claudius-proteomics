# TO_HF_CORPUS.md — Codex review (gpt-5.5, 2026-06-10)

Independent review of the v3 plan. Verdicts + actions tracked in
TO_HF_CORPUS.md v4.

## Highest-leverage issues

### 1. "1% FDR floor" is not statistically guaranteed
The inclusion rule is the **union** of two engine-specific sets
(`sage_q ≤ 0.01 OR fragpipe_q ≤ 0.01`). Each engine controls FDR within
its own search; their union is not automatically a 1%-FDR corpus
(different search spaces, scoring, decoys, calibration strata). Describe
as "accepted by ≥1 engine at that engine's reported q≤0.01 threshold."
Filtering to `q ≤ 0.001` does not "recompute FDR" — it applies a stricter
reported-q threshold; post-hoc combinations with agreement/PTM/charge/
quality filters do not retain nominal FDR. Verify every exported q is
PSM-level, not peptide/protein-level.

### 2. The row unit is underspecified
"One row per (raw_file, precursor_id)" conflicts with calling it a PSM
corpus. A precursor may yield multiple PASEF events / rank-1 PSMs.
Collapsing needs explicit rules (which PSM supplies q/PEP/RT/mobility/
fragments; aggregate vs discard repeats; conflicting assignments; whether
`collision_energies` covers all or selected events). Ship stable IDs
(`frame_id`, `precursor_number`, `scan_begin/end`, native id, immutable
`psm_id`) + observation counts + aggregation provenance.

### 3. `n_engines` and null semantics need separation
Null conflates: engine not run / artifact unavailable / no candidate /
not rank-1 / failed q / different peptidoform. Add per-engine status
columns (`*_searched`, `*_candidate_present`, `*_rank`, `*_passes_floor`).
Define `n_engines` = number assigning the **same normalized peptidoform +
charge while passing the threshold**. PEPs are complementary but
**correlated** (same spectrum, possibly shared DB/features) — document
that they must not be multiplied/averaged into a consensus.

### 4. Fragment "both" representation risks triple-counting
Keeping sage + fragpipe + `both` rows duplicates agreed fragments; naive
users inflate counts/intensities. Prefer one canonical fragment row
(`sage_matched`, `fragpipe_matched`, per-engine mz_exp + intensity,
`agreement_ppm`). Inputs are asymmetric (Sage native matches vs FragPipe
rematched by imspy) → measures rematcher agreement; apply same matcher to
both or expose `match_method`. Define tolerance, mass convention, NL/
isotope handling, deterministic one-to-one resolution (ppm-range join can
go many-to-many).

### 5. Width reliability is not release-ready
Split `rt_width_reliable` / `im_width_reliable` (a single flag hides
which passed). The SNR definition is an open question yet `SNR≥20` /
"≤6% σ error" are stated as established — freeze + version the baseline
window, MAD scaling, smoothing, fit model, min points, edge handling,
zero-noise treatment; state provenance (simulation vs real, instruments,
intensity range). Include units + fit status/reason; null must
distinguish missing-extraction from failed-fit. Do **not** impute
`gradient_length` with cohort median for calibration — keep missing;
imputation is a separate named column.

### 6. Split construction is currently impossible as written
Sequence-hash assignment and forcing whole `group_id`s into one split
conflict. Build connected components of the bipartite graph linking
`sequence_normalized` and `group_id`, assign components, apply canonical
holdouts first, publish the algorithm (not just a seed). Consider
modified-sequence and spectrum/raw-file leakage; exact fingerprinting
misses converted copies / renamed submissions / repeated acquisitions.

### 7. Manifest does not reproduce the release
Git SHA + mutable MOGON paths is insufficient. Record: SHA-256 + size per
input artifact and output shard; raw-file fingerprints + source URLs/
access dates; engine + MSBooster/PeptideProphet + timsrust/imspy versions;
full search configs + FASTA checksum + contaminant DB + decoy method;
extraction params + container/env digest; schema with units/null
semantics + per-shard row counts; a PRIDE metadata snapshot (not live API).
The "CE [0,1] assertion" looks erroneous — CE is not a unit-interval
quantity (note: it's our /100-normalized convention; ship raw volts).

## Peer-review risks
- "First public" needs a documented search or softer wording.
- "PRIDE-verified CC0 / zero license vetting" — keep per-dataset evidence
  + retrieval dates; the "zero vetting" line is an avoidable provocation.
- Repo names containing `v0.1` complicate semantic evolution — use one
  repo with immutable revision tags + shard checksums.
- Reconcile "v0.1 first release but private" vs "v1.0 first public".
