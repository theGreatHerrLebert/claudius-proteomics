# Bug Report: Precursor Index & Extraction Issues

Date: 2026-01-29

## Summary

Several bugs were identified and fixed in the precursor extraction and indexing pipeline.

---

## Bug 1: MS1 Extraction Returning Empty Data (FIXED)

### Symptoms
- Only ~0.1% of precursors had MS1 raw data
- XIC and mobilogram plots were empty for most precursors

### Root Causes

**1a. NaN "or" Logic Bug**
```python
# BROKEN: np.nan is truthy in Python, so this returns nan
mz = row.get('fragpipe_mz') or row.get('raw_mz') or 0

# FIXED: Use explicit pd.notna() checks
mz_fp = row.get('fragpipe_mz')
mz = float(mz_fp) if pd.notna(mz_fp) else (float(mz_raw) if pd.notna(mz_raw) else 0.0)
```

**1b. RT Unit Conversion Bug**
```python
# BROKEN: Fragment RT is already in seconds, but code multiplied by 60
rt_sec = frag_data.get('rt', 0) * 60.0  # RT ends up as 7200 sec instead of 120 sec

# FIXED: RT is already in seconds, don't multiply
rt_sec = float(frag_data.get('rt', 0))
```

### Impact
With these fixes, MS1 extraction now produces 1.2 billion raw 4D points (vs nearly 0 before).

### Files Changed
- `scripts/precursor_store_parquet.py`

---

## Bug 2: Fragment Extraction OOM (FIXED)

### Symptoms
- Extraction process killed due to 25GB+ memory usage
- `get_pasef_fragments()` loads all ~500k fragments at once

### Root Cause
The Rust API `get_pasef_fragments()` returned all fragments at once, with no batched alternative.

### Fix
Added new Rust API `get_pasef_fragments_for_precursors(precursor_ids, num_threads)` that:
- Extracts fragments only for specified precursor IDs
- Called per batch during extraction
- Properly merges re-targeted precursors (same precursor fragmented multiple times)

### Files Changed
- `rustdf/src/data/dda.rs` - Added `get_pasef_fragments_for_precursors()`
- `imspy_connector/src/py_dda.rs` - Added Python binding
- `scripts/precursor_store_parquet.py` - Use batched API

---

## Bug 3: DIA-NN Not Joined to Precursor Index (PARTIALLY FIXED)

### Symptoms
- `n_engines` max was 2 (only FragPipe + Sage)
- `diann_peptide` column missing from precursor_index.parquet

### Root Cause
DIA-NN matching (`match_diann_to_precursors`) only worked for FragPipe-identified precursors, using sequence-based matching. For unidentified precursors, there was no matching at all.

### Partial Fix
Updated `match_diann_to_precursors` to:
1. First match by sequence + charge + raw_file (for identified precursors)
2. Then match by m/z + charge + RT + IM (for unidentified precursors)

### Remaining Issue
**Unidentified precursors have no RT value in the index**, so DIA-NN is matched by m/z + charge only (no RT/IM filter). This causes:
- DIA-NN matches 170k precursors, but all are different from FragPipe/Sage
- Only 8 precursors have all 3 engines (and only 1 agrees on sequence)

### Required Fix
Need to derive RT for ALL raw precursors from frame metadata:
```python
# In load_raw_precursors(), add:
# 1. Load frame metadata to get frame_id -> RT mapping
# 2. Join raw precursors with RT based on frame_id
```

### Files Changed
- `scripts/build_precursor_index.py` - Improved DIA-NN matching

---

## Bug 4: Fragment IM vs m/z Plot Shows as 1D Line (ADDRESSED)

### Symptoms
- Fragment IM vs m/z scatter plot appeared as a horizontal line
- No 2D distribution visible

### Root Cause
PASEF isolation window is very narrow (~0.02 1/K0), so all fragments from one precursor have nearly identical mobility.

### Solution
Changed visualization from scatter plot to heatmap with scan index on Y-axis:
- X-axis: m/z (binned)
- Y-axis: scan index (discrete integer values)
- Color: accumulated intensity per cell

### Files Changed
- `dashboard/frontend/src/components/PrecursorViz.tsx` - Added HeatmapPlot component
- `scripts/precursor_store_parquet.py` - Added `fragment_scan` column
- `dashboard/backend/main.py` - Added `fragment_scan` to API

---

## Statistics After Fixes

### Extraction Output (frac01)
- Precursors: 214,824
- Fragment peaks: 218 million
- XIC points: 5.3 million
- Raw 4D points: 1.2 billion
- File size: 5.4 GB

### Engine Coverage (3-file dataset)
- FragPipe only: 5,297
- DIA-NN only: 170,682 (no overlap due to RT issue)
- Sage only: 0
- FP + SG: 172,408
- All three: 8 (only 1 with sequence agreement)

---

## Recommendations

1. **Fix RT for raw precursors** - Critical for proper DIA-NN matching
2. **Sequence normalization** - Check I/L and modification format consistency between engines
3. **Add sequence-based validation** - Warn when engines disagree on sequence for same precursor
