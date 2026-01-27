# RAW DATA Extraction

**Owner:** TBD
**Priority:** High
**Status:** Not Started

## Overview

Expand the timsTOF feature extraction capabilities beyond current XICs, mobilograms, and isotope patterns. Enable richer raw signal features for model training and quality assessment.

## Current State

- imspy extracts XICs, mobilograms, isotope patterns from .d files
- Implementation in `scripts/extract_raw_features.py` (~463 LOC)
- Parameters configured in `config/config.yaml` under `extraction:`

## TODOs

### 1. Spectral Features

- [ ] Extract MS2 fragment ion intensities (b/y ions)
- [ ] Extract precursor isolation window purity
- [ ] Extract apex intensity vs integrated intensity
- [ ] Extract peak width (FWHM) for RT and mobility

### 2. Raw Signal Quality

- [ ] Extract signal-to-noise ratio per precursor
- [ ] Extract number of MS2 scans per precursor
- [ ] Extract frame-level statistics (TIC, base peak)
- [ ] Extract collision energy values

### 3. timsTOF-Specific Features

- [ ] Extract PASEF efficiency metrics
- [ ] Extract accumulation time
- [ ] Extract mobility calibration coefficients
- [ ] Extract temperature/pressure metadata if available

### 4. Multi-Instrument Support (Future)

- [ ] Thermo .raw via RawFileReader
- [ ] Sciex .wiff via ProteoWizard
- [ ] Design abstraction layer for instrument-agnostic features

### 5. Batch Processing Improvements

- [ ] Parallel extraction across multiple .d files
- [ ] Resume capability for interrupted extractions
- [ ] Memory-efficient streaming for large files

## Files to Modify

- `scripts/extract_raw_features.py` - Main extraction logic
- `rules/extract.smk` - Snakemake rule
- `config/config.yaml` - Add new extraction parameters
- `database/schema.json` - Add new column definitions

## Dependencies

- imspy / rustims stack
- Understanding of timsTOF .d file structure

## Acceptance Criteria

- [ ] All new features extracted and stored in Parquet
- [ ] Unit tests for each extraction function
- [ ] Documentation of new fields in schema.json
- [ ] Benchmark extraction speed (files/hour)
