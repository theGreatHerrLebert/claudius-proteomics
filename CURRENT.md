# Current Task: PXD046675 Full Pipeline Run on monster2

## What was done

### 1. Instrument model feature (`61db4f4`, `ac4a97f`)
Extended sample group resolution so `group_id` becomes `{organism}_{enzyme}_{instrument}`
when multiple timsTOF instrument models are detected within a PXD. Changes:

- **`scripts/sample_group_resolver.py`**: Added `instrument_model` field to `SampleGroup`,
  instrument normalization via `INSTRUMENT_ALIASES`, TDF `analysis.tdf` reading
  (`_read_instrument_from_d`), `_scan_instruments()`, and `_refine_by_instrument()` which
  splits groups containing runs from different instruments. Single-instrument groups keep
  their original group_id (backward compatible).
- **`runner/steps/step1_download.py`**: Logs multi-instrument info, adds `is_multi_instrument`
  and `group_instruments` to step summary.
- **`runner/steps/step2_search.py`**: Prints instrument model in per-group header.

### 2. Bug fix: instrument_model null for fallback datasets (`c9818a0`)
**Problem**: When `.d` folder names are cryptic (e.g. `L120-1_Slot2-13_1_11628.d`),
`populate_runs` couldn't match them to groups. They ended up as unassigned, then got
rescued into the sole group in step 7. But instrument refinement (step 6) ran before
the rescue, so the group had no runs during TDF scanning and `instrument_model` stayed null.

**Fix**: Swapped steps 6 and 7 in `resolve_sample_groups()` — unassigned run rescue now
happens before instrument refinement, so runs are on the group when TDF files are read.

### 3. Pig (Sus scrofa) organism support (`61db4f4`)
- Added pig aliases, taxon IDs, organism config
- Added `PXD046675` to datasets list and `dataset_metadata` in config

### 4. PXD046675 data preparation
- Downloaded `original_data.zip` (~15 GB) from PRIDE FTP
- Extracted 6 `.d` folders: L120-1/2/3, L240-1/2/3 (pig muscle, 120d vs 240d)
- All single instrument: timsTOF Pro
- Downloaded pig proteome from UniProt (UP000008227, 46,010 proteins)
- Generated decoy FASTA (92,020 sequences) at `resources/fasta/search_db/pig_decoys.fasta`

### 5. Validated locally
- Sample group resolver: `pig_trypsin` group, 6 runs, `instrument_model: timsTOF Pro`
- Step 1 (test mode): passed, instrument info logged correctly
- Step 2: confirmed FASTA found, engines ready to run

### 6. Moved to monster2
- Pushed bug fix commit (`c9818a0`) and pulled on monster2
- Rsynced 6 `.d` folders (~16 GB) to `/globalscratch/dateschn/claudius-data/raw/PXD046675/`
- Rsynced pig FASTA + decoys to monster2's `resources/fasta/search_db/`
- Verified all three search engines (FragPipe, DIA-NN, Sage) present and configured

## What is in progress

### Full pipeline run on monster2 (PID 1860458)
```bash
# Running as nohup on monster2 (128 cores, 503 GB RAM)
.venv/bin/python -u runner/run_dataset.py PXD046675 \
  --local-data /globalscratch/dateschn/claudius-data/raw/PXD046675
```
- Log: `/globalscratch/dateschn/claudius-data/run_PXD046675.log`
- Step 1 completed (0.2s)
- Step 2 in progress: FragPipe running, then DIA-NN, then Sage (all 6 files, no test mode)
- Steps 3-5 will follow automatically

Monitor with:
```bash
ssh monster2 "tail -f /globalscratch/dateschn/claudius-data/run_PXD046675.log"
```

**Note**: There is also a PXD019086 (yeast) run happening on monster2 concurrently.

## What to do next

### 1. Check pipeline results
```bash
ssh monster2 "tail -100 /globalscratch/dateschn/claudius-data/run_PXD046675.log"
```
Look for step 2 (search) completion, step 3 (stratify) results, and step 4/5 output.

### 2. If step 4 (extract) fails with imspy import error
The venv was updated (user confirmed), but if `TimsDatasetDDA` is still missing,
check `runner/steps/step4_extract.py:273` — may need to use `TimsDataset` instead.

### 3. Review outputs
- `data/processed/PXD046675/pig_trypsin/` — per-engine search results
- `data/processed/PXD046675/precursor_index.parquet` — unified index
- `data/merged/PXD046675/precursor_store.parquet` — final merged output

### 4. Consider removing local zip
The 15 GB `data/raw/PXD046675/original_data.zip` on this machine can be deleted
to save space — the extracted `.d` folders and monster2 copy are the sources of truth.
