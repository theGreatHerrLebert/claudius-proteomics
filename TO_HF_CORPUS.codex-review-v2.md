# TO_HF_CORPUS.md — Codex build-readiness review (gpt-5.5, 2026-06-10, v4)

Code-grounded (Codex inspected `BlobReader`, `extract_fragment_peaks.py`,
canonical parquets). Verdicts folded into TO_HF_CORPUS.md v5.

### 1. Blob pass — `BlobReader.read_batch()` is NOT production-usable
Full 6.4 TB compressed I/O (not "a fraction" — each precursor is separately
outer-compressed; any access decompresses the whole blob). Current reader
materializes the full decompressed blob, eagerly loads fragment + MS1 + full
raw 4D arrays, returns a list retaining every `RawSignal`, and silently nulls
decode exceptions → severe memory amplification.
**Fix — projection reader:** one blob at a time; extract only metadata +
`ms1_rt_*`/`ms1_im_*`; compute immediately; release buffers. Sort rows by
`(raw_file, blob_offset)` for sequential reads; open each raw file once; never
retain decoded signals. Validate `offset>=0, size>0, offset+size<=file_size`;
verify metadata `precursor_id` == parquet row; count failures by reason and
quarantine above a frozen threshold (never silent-null); write incremental
Parquet row groups (not an accession-wide pandas frame); checkpoint per raw
file with atomic rename + input checksums; size Slurm concurrency by FS
throughput and memory by **max** uncompressed blob. CE unit is **volts**, not
eV (existing extraction logs mislabel it).

### 2. Tier-1 join + grain coherence
1:1 gate must assert: unique keys in BOTH inputs; equal key **sets** (not just
counts); identical multiplicities; non-null keys; blob refs resolve to the
expected store; joined count == both sources (pre-floor). Do NOT normalize raw
filenames inside the assertion — build an explicit alias table first
(`%20`/underscore/basename rewriting can collide).
**Grain incoherence (my error):** a per-precursor row cannot hold two
peptidoforms as separate accepted PSMs, and `n_engines` is a scalar — "counts
them separately" is wrong. For a conflict row: pick ONE canonical engine
assignment; `n_engines` reflects support for THAT assignment; canonical
`charge`+sequence from the **same** engine (never mixed); add
`sage_assignment_passes_floor`/`fragpipe_assignment_passes_floor` +
`peptidoform_conflict=true`. Union null-fill must be driven by an explicit
Arrow schema; distinguish missing-column vs all-null via `source_column_present`
in the manifest.

### 3. Tier-3 — existing extractor unsuitable as canonical input
`extract_fragment_peaks.py`: emits one row/precursor with list columns; filters
on `is_high_quality` + `n_engines` (doesn't cover the permissive cohort);
`get_sequence_for_psm()` falls back to `consensus_peptide` (can falsely report
an engine match); `_find_peak()` matches each theoretical ion independently
(one peak can satisfy multiple ions — NOT one-to-one); unmatched = `mz=0` not
null.
**Fix:** rerun both passes WITHOUT consensus fallback and WITHOUT the quality
gate; deterministic greedy one-to-one matching, candidates sorted by
`absolute_ppm, desc_peak_intensity, ion_type, ordinal, charge, peak_index`,
consuming each ion and peak once. **Only merge Sage/FragPipe rows under the
canonical fragment key when normalized peptidoform AND precursor charge are
identical** — for conflicting peptidoforms `b3(A) != b3(B)`, so include
assignment identity in the key or emit separate rows, else `agreement_ppm` is
chemically meaningless. Freeze:
`agreement_ppm = 1e6 * abs(sage_mz_exp - fragpipe_mz_exp) / fragment_mz_calculated`,
null unless both matched the **same** theoretical ion.
⚠️ Consequence: the existing 62 `fragment_peaks{,_fragpipe}.parquet` were made
by the flawed extractor → **must be regenerated** for the permissive corpus.
The "62 satisfy matcher-(a) out of the box" claim was wrong.

### 4. Readiness vs inclusion — content gates, not file existence
Add: every included raw file has a blob store + valid offsets; both rematch
files cover the post-floor Tier-1 keys (documented failure counts); identical
imspy version + config between passes; per-raw compression readability + acq
mode + instrument; no mixed excluded acq types within an accession; CC0
evidence attached **before build scheduling**; anomaly quarantine before
advertising totals. Report the ready count AFTER these gates, not "62 minus
exclusions."

### 5. Freeze the SNR estimator — `trace_snr_v1`
1. Reject non-finite coords/intensities; clamp negatives to 0.
2. Require ≥7 trace points and ≥3 each side of the fitted apex.
3. Peak exclusion = apex ± `2.5*sigma`; reject truncated traces where this
   reaches a boundary.
4. Baseline samples = points outside that interval.
5. `baseline = median(baseline_samples)`.
6. `noise_sigma = 1.4826 * median(abs(samples - baseline))`.
7. `snr = max(0, apex_intensity - baseline) / noise_sigma`.
8. If <6 baseline points / `noise_sigma<=0` / fit failed / truncated → SNR
   null + specific status, never infinity.
9. No smoothing before baseline/noise estimation.
Keep `SNR>=20` provisional until this exact estimator is revalidated.
