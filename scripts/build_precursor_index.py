#!/usr/bin/env python3
"""
Build Unified Precursor Index

Creates a single parquet file that maps:
  precursor_id (raw data) → search engine identifications + quality metrics

This is the bridge between raw timsTOF data and search engine results,
enabling easy lookup: "For this precursor ID, what peptide was identified
and by which engines?"

Output columns:
  - precursor_id: timsTOF precursor ID (link to raw data)
  - raw_file: Source .d file name
  - mz, charge, rt, mobility: Precursor properties from raw data

  Search engine IDs (NULL if not identified by that engine):
  - fragpipe_peptide, fragpipe_protein, fragpipe_probability, fragpipe_hyperscore
  - diann_peptide, diann_protein, diann_qvalue, diann_ccs
  - sage_peptide, sage_protein, sage_qvalue, sage_hyperscore

  Consensus:
  - sequence_normalized: I/L normalized sequence for matching
  - n_engines: Number of engines that identified this precursor
  - consensus_peptide: Best peptide (from most engines)
  - confidence_weight: 1.0 if all engines agree, lower if discrepant

Usage:
    python scripts/build_precursor_index.py --accession PXD019086
    python scripts/build_precursor_index.py --accession PXD019086 --raw-file "frac01.d"
"""

import argparse
import re
from pathlib import Path
from typing import Optional, Dict, List, Any

import pandas as pd
import numpy as np

# Import sequence standardization utilities
from sequence_utils import (
    standardize_diann_sequence,
    standardize_sage_sequence,
    standardize_fragpipe_sequence,
    standardize_fragpipe_modified_peptide,
    remove_modifications,
    normalize_sequence_il,
)

# Import precursor matching utilities
from precursor_matching import (
    PrecursorMatcher,
    MatchConfig,
    MatchTier,
    normalize_sequence_for_matching,
)


def normalize_sequence(sequence: str) -> str:
    """Normalize sequence for matching (I→L, uppercase)."""
    if pd.isna(sequence):
        return ""
    # Replace I with L outside of modification brackets
    def replace_match(match):
        if match.group(1):
            return match.group(1)
        return match.group(2).replace("I", "L")
    pattern = r'(\[.*?\])|([^\[\]]+)'
    return re.sub(pattern, replace_match, str(sequence).upper())


def remove_mods(sequence: str) -> str:
    """Remove modification annotations from sequence."""
    if pd.isna(sequence):
        return ""
    clean = re.sub(r'\[UNIMOD:\d+\]', '', str(sequence))
    clean = re.sub(r'\[[+-]?\d+\.?\d*\]', '', clean)
    clean = re.sub(r'\(UniMod:\d+\)', '', clean)
    return clean


def parse_fragpipe_spectrum_id(spectrum: str) -> tuple:
    """Parse FragPipe spectrum ID to get raw file and precursor ID."""
    # Format: rawfile.scannum.scannum.charge
    parts = spectrum.rsplit(".", 3)
    if len(parts) >= 4:
        raw_file = parts[0]
        precursor_id = int(parts[1])
        return raw_file, precursor_id
    return None, None


def load_fragpipe_psms(
    accession: str,
    base_dir: Path,
    raw_file_filter: Optional[str] = None,
) -> pd.DataFrame:
    """Load FragPipe PSMs with precursor_id mapping."""

    # Find all PSM files
    processed_dir = base_dir / accession
    psm_files = list(processed_dir.rglob("psm.tsv"))

    if not psm_files:
        print(f"  No FragPipe PSM files found in {processed_dir}")
        return pd.DataFrame()

    dfs = []
    for psm_file in psm_files:
        df = pd.read_csv(psm_file, sep="\t")

        # Parse spectrum ID to get precursor_id
        parsed = df["Spectrum"].apply(parse_fragpipe_spectrum_id)
        df["raw_file"] = parsed.apply(lambda x: x[0])
        df["precursor_id"] = parsed.apply(lambda x: x[1])

        # Filter by raw file if specified
        if raw_file_filter:
            df = df[df["raw_file"].str.contains(raw_file_filter, case=False, na=False)]

        dfs.append(df)

    if not dfs:
        return pd.DataFrame()

    all_psms = pd.concat(dfs, ignore_index=True)

    # Keep best PSM per precursor (highest probability)
    all_psms = all_psms.sort_values("Probability", ascending=False)
    all_psms = all_psms.drop_duplicates(subset=["raw_file", "precursor_id"], keep="first")

    # Extract charge from Spectrum ID (last part after final dot)
    def extract_charge(spectrum):
        parts = spectrum.rsplit(".", 1)
        if len(parts) >= 2:
            try:
                return int(parts[-1])
            except ValueError:
                return None
        return None

    all_psms["extracted_charge"] = all_psms["Spectrum"].apply(extract_charge)

    # Standardize modifications to UNIMOD format
    # Priority: Use Peptide + Assigned Modifications (includes fixed mods like carbamidomethyl)
    # "Modified Peptide" column often excludes fixed modifications
    if "Assigned Modifications" in all_psms.columns:
        all_psms["fragpipe_modified_std"] = all_psms.apply(
            lambda r: standardize_fragpipe_sequence(
                r["Peptide"],
                r.get("Assigned Modifications", "") if pd.notna(r.get("Assigned Modifications")) else ""
            ),
            axis=1
        )
    elif "Modified Peptide" in all_psms.columns:
        # Fallback: Use "Modified Peptide" column (may miss fixed mods)
        all_psms["fragpipe_modified_std"] = all_psms["Modified Peptide"].apply(
            standardize_fragpipe_modified_peptide
        )
    else:
        all_psms["fragpipe_modified_std"] = all_psms["Peptide"]

    # Rename/select columns
    result = pd.DataFrame({
        "raw_file": all_psms["raw_file"],
        "precursor_id": pd.to_numeric(all_psms["precursor_id"], errors="coerce").astype("Int64"),
        "fragpipe_peptide": all_psms["Peptide"],
        "fragpipe_modified": all_psms["fragpipe_modified_std"],  # Standardized UNIMOD format
        "fragpipe_protein": all_psms["Protein"],
        "fragpipe_probability": all_psms["Probability"],
        "fragpipe_pep": 1.0 - all_psms["Probability"],  # PEP = 1 - probability
        "fragpipe_hyperscore": all_psms.get("Hyperscore", np.nan),
        "fragpipe_qvalue": all_psms.get("Qvalue", np.nan),
        "fragpipe_mz": all_psms.get("Calibrated Observed M/Z", all_psms.get("Observed M/Z")),
        "fragpipe_rt": all_psms.get("Retention"),
        "fragpipe_mobility": all_psms.get("Ion Mobility"),
        "fragpipe_charge": all_psms.get("Charge", all_psms["extracted_charge"]),
    })

    return result


def load_diann_report(
    accession: str,
    base_dir: Path,
    raw_file_filter: Optional[str] = None,
) -> pd.DataFrame:
    """Load DIA-NN report with m/z-based precursor matching."""

    report_path = base_dir / accession / "diann" / "report.parquet"
    if not report_path.exists():
        print(f"  DIA-NN report not found: {report_path}")
        return pd.DataFrame()

    df = pd.read_parquet(report_path)

    # Filter by raw file if specified
    if raw_file_filter and "Run" in df.columns:
        df = df[df["Run"].str.contains(raw_file_filter, case=False, na=False)]

    # Standardize DIA-NN modified sequence to UNIMOD format
    modified_col = df.get("Modified.Sequence", pd.Series(dtype=str))
    df["diann_modified_std"] = modified_col.apply(standardize_diann_sequence)

    # DIA-NN doesn't have direct precursor_id mapping - we'll match by m/z later
    result = pd.DataFrame({
        "diann_run": df.get("Run", ""),
        "diann_peptide": df.get("Stripped.Sequence", df.get("Sequence", "")),
        "diann_modified": df["diann_modified_std"],  # Standardized UNIMOD format
        "diann_protein": df.get("Protein.Ids", ""),
        "diann_charge": df.get("Precursor.Charge"),
        "diann_mz": df.get("Precursor.Mz"),
        "diann_qvalue": df.get("Q.Value"),
        "diann_pep": df.get("PEP"),
        "diann_global_qvalue": df.get("Global.Q.Value"),
        "diann_pg_qvalue": df.get("PG.Q.Value"),
        "diann_rt": df.get("RT"),
        "diann_mobility": df.get("IM"),
        "diann_ccs": df.get("CCS"),
    })

    return result


def load_sage_results(
    accession: str,
    base_dir: Path,
    raw_file_filter: Optional[str] = None,
) -> pd.DataFrame:
    """Load Sage results with precursor_id mapping."""

    sage_path = base_dir / accession / "sage" / "results.sage.parquet"
    if not sage_path.exists():
        print(f"  Sage results not found: {sage_path}")
        return pd.DataFrame()

    df = pd.read_parquet(sage_path)

    # Filter decoys
    if "is_decoy" in df.columns:
        df = df[~df["is_decoy"]]

    # Filter by raw file if specified
    if raw_file_filter and "filename" in df.columns:
        df = df[df["filename"].str.contains(raw_file_filter, case=False, na=False)]

    # Sage uses scannr from mzML conversion - not the timsTOF precursor_id!
    # We'll match by m/z + charge instead
    # Calculate precursor m/z from mass and charge
    df_copy = df.copy()
    df_copy["sage_mz"] = (df_copy["expmass"] + df_copy["charge"] * 1.007276) / df_copy["charge"]

    # Standardize Sage modified sequence to UNIMOD format
    # Sage "peptide" column contains [+mass] format, e.g., "PEPTIDE[+15.9949]K"
    peptide_col = df_copy.get("peptide", pd.Series(dtype=str))
    df_copy["sage_modified_std"] = peptide_col.apply(standardize_sage_sequence)

    result = pd.DataFrame({
        "sage_raw_file": df_copy.get("filename", "").apply(lambda x: Path(x).stem.replace(".d", "") if pd.notna(x) else ""),
        "sage_scannr": df_copy.get("scannr"),  # Keep for reference but don't use for matching
        "sage_peptide": df_copy.get("stripped_peptide"),
        "sage_modified": df_copy["sage_modified_std"],  # Standardized UNIMOD format
        "sage_protein": df_copy.get("proteins"),
        "sage_charge": df_copy.get("charge"),
        "sage_mz": df_copy["sage_mz"],  # Calculated precursor m/z
        "sage_qvalue": df_copy.get("spectrum_q"),
        "sage_peptide_qvalue": df_copy.get("peptide_q"),
        "sage_protein_qvalue": df_copy.get("protein_q"),
        "sage_pep": np.exp(df_copy["posterior_error"]) if "posterior_error" in df_copy.columns else np.nan,
        "sage_hyperscore": df_copy.get("hyperscore"),
        "sage_rt": df_copy.get("rt"),
        "sage_mobility": df_copy.get("ion_mobility"),
    })

    return result


def load_raw_precursors(
    accession: str,
    local_data_dir: Path,
    raw_file_filter: Optional[str] = None,
    use_calibration: bool = True,
) -> pd.DataFrame:
    """Load precursor metadata from raw data (.d files).

    Extracts complete precursor info by joining:
    - fragmented_precursors: precursor_id, mz, charge, intensity
    - pasef_meta_data: frame_id, scan range, isolation window
    - meta_data (frames): retention time

    Args:
        accession: PRIDE accession
        local_data_dir: Path to directory containing .d files
        raw_file_filter: Optional filter for specific raw file
        use_calibration: Use pre-computed IM calibration for accurate values

    Returns one row per unique precursor (aggregated across PASEF events).
    """
    from imspy_core.timstof import TimsDatasetDDA

    # Find .d files
    d_files = list(local_data_dir.glob("*.d"))
    if raw_file_filter:
        d_files = [f for f in d_files if raw_file_filter in f.name]

    if not d_files:
        print(f"  No .d files found in {local_data_dir}")
        return pd.DataFrame()

    # Import calibration utilities if needed
    if use_calibration:
        from extract_calibration import ensure_calibration, get_calibration_path

    all_precursors = []
    for d_file in d_files:
        print(f"  Loading precursors from: {d_file.name}")
        try:
            # Load dataset with optional calibration
            if use_calibration:
                cal_path = get_calibration_path(str(d_file))
                if cal_path.exists():
                    im_lookup = np.load(cal_path)
                else:
                    print(f"    Extracting IM calibration (one-time)...")
                    im_lookup = ensure_calibration(str(d_file), verbose=False)

                # Create calibrated dataset
                from imspy_connector.py_dda import PyTimsDatasetDDA as PyTimsDatasetDDARust
                rust_dataset = PyTimsDatasetDDARust.with_calibration(
                    str(d_file), False, im_lookup.tolist()
                )
                # Also load Python wrapper for metadata access
                dataset = TimsDatasetDDA(str(d_file), in_memory=False, use_bruker_sdk=False)
            else:
                dataset = TimsDatasetDDA(str(d_file), in_memory=False, use_bruker_sdk=False)
                rust_dataset = None

            # Get base precursor info
            precursors = dataset.fragmented_precursors.copy()

            # Get PASEF metadata (links precursor_id to frame_id, scan range, isolation)
            pasef = dataset.pasef_meta_data.copy()

            # Get frame metadata (has Time/RT)
            frames = dataset.meta_data[['frame_id', 'Time']].copy()

            # Join PASEF to frames to get RT
            pasef = pasef.merge(frames, on='frame_id', how='left')

            # Calculate mobility from scan range
            # Use center of scan range for mobility calculation
            scan_center = ((pasef['scan_begin'] + pasef['scan_end']) / 2).astype(np.int32)

            # Convert scans to mobility
            # When calibration is used, rust_dataset has accurate IM values
            # When not, use the dataset (which has linear interpolation)
            mobilities = []
            for frame_id in pasef['frame_id'].unique():
                mask = pasef['frame_id'] == frame_id
                scans = scan_center[mask].values
                try:
                    if use_calibration and rust_dataset is not None:
                        # Use calibrated Rust dataset (accurate + fast)
                        im_values = rust_dataset.scan_to_inverse_mobility(
                            int(frame_id), [int(s) for s in scans]
                        )
                    else:
                        # Use Python dataset (linear interpolation)
                        im_values = dataset.scan_to_inverse_mobility(int(frame_id), scans)
                    mobilities.extend(im_values)
                except Exception:
                    # Fallback: use approximate conversion if calibration fails
                    mobilities.extend([np.nan] * len(scans))

            pasef['mobility'] = mobilities
            pasef['rt_seconds'] = pasef['Time']  # Time is already in seconds

            # Aggregate PASEF events per precursor (same precursor fragmented multiple times)
            # Take first frame_id, mean RT, mean mobility, first isolation params
            pasef_agg = pasef.groupby('precursor_id').agg({
                'frame_id': 'first',
                'rt_seconds': 'mean',
                'mobility': 'mean',
                'isolation_mz': 'first',
                'isolation_width': 'first',
                'scan_begin': 'first',
                'scan_end': 'first',
            }).reset_index()

            # Join precursors with aggregated PASEF data
            merged = precursors.merge(pasef_agg, on='precursor_id', how='left')
            merged["raw_file"] = d_file.stem

            all_precursors.append(merged)

        except Exception as e:
            print(f"    Error loading {d_file}: {e}")
            import traceback
            traceback.print_exc()

    if not all_precursors:
        return pd.DataFrame()

    df = pd.concat(all_precursors, ignore_index=True)

    # Select columns - map timsTOF column names to our schema
    result = pd.DataFrame({
        "raw_file": df["raw_file"],
        "precursor_id": df["precursor_id"].astype(int),
        "raw_mz": df.get("monoisotopic_mz", df.get("mono_mz", df.get("largest_peak_mz"))),
        "raw_charge": df.get("charge"),
        "raw_mobility": df.get("mobility"),  # Now from PASEF scan conversion
        "raw_rt_seconds": df.get("rt_seconds"),  # RT in seconds from frame Time
        "raw_intensity": df.get("intensity", df.get("precuror_total_intensity")),
        "frame_id": df.get("frame_id"),
        "isolation_mz": df.get("isolation_mz"),
        "isolation_width": df.get("isolation_width"),
        "parent_id": df.get("parent_id"),
    })

    return result


def match_by_mz_charge_hash(
    source_df: pd.DataFrame,
    target_df: pd.DataFrame,
    source_mz_col: str,
    source_charge_col: str,
    target_mz_col: str,
    target_charge_col: str,
    target_qvalue_col: str,
    mz_tol_ppm: float = 20.0,
    source_rt_col: str = None,
    target_rt_col: str = None,
    rt_tol_min: float = 1.0,
) -> pd.DataFrame:
    """
    Fast m/z + charge matching using hash bins with optional RT filtering.

    Bins m/z values and uses hash-based lookup. O(n + m).

    For 20 ppm at m/z 1000, tolerance = 0.02 Da
    Use bin size of ~0.01 Da (half the smallest tolerance)
    """
    if target_df.empty or source_df.empty:
        return pd.DataFrame()

    # Work with copies
    source = source_df.copy()
    target = target_df.copy()

    # Ensure numeric types and drop NA
    source[source_mz_col] = pd.to_numeric(source[source_mz_col], errors="coerce")
    source[source_charge_col] = pd.to_numeric(source[source_charge_col], errors="coerce").astype("Int64")
    target[target_mz_col] = pd.to_numeric(target[target_mz_col], errors="coerce")
    target[target_charge_col] = pd.to_numeric(target[target_charge_col], errors="coerce").astype("Int64")

    source = source.dropna(subset=[source_mz_col, source_charge_col])
    target = target.dropna(subset=[target_mz_col, target_charge_col])

    if source.empty or target.empty:
        return pd.DataFrame()

    # Check if RT filtering is available
    use_rt = (source_rt_col and target_rt_col and
              source_rt_col in source.columns and target_rt_col in target.columns)

    # Bin size: At m/z=500 with 20ppm, tolerance=0.01. Use bin_size=0.005
    bin_size = 0.005  # Da

    def mz_to_bin(mz):
        return int(mz / bin_size)

    # Create hash key: (charge, mz_bin)
    source["_hash_key"] = list(zip(source[source_charge_col], source[source_mz_col].apply(mz_to_bin)))
    target["_hash_key"] = list(zip(target[target_charge_col], target[target_mz_col].apply(mz_to_bin)))

    # Build target hash table: key -> list of (idx, mz, qvalue, rt)
    target_hash = {}
    for idx, row in target.iterrows():
        key = row["_hash_key"]
        mz = row[target_mz_col]
        qvalue = row.get(target_qvalue_col, 1.0)
        charge = row[target_charge_col]
        rt = row.get(target_rt_col) if use_rt else None

        # Add to current bin and neighbor bins (for edge cases)
        for offset in [-1, 0, 1]:
            neighbor_key = (charge, key[1] + offset)
            if neighbor_key not in target_hash:
                target_hash[neighbor_key] = []
            target_hash[neighbor_key].append((idx, mz, qvalue, rt))

    # Match sources to targets
    matches = []
    for source_idx, source_row in source.iterrows():
        key = source_row["_hash_key"]
        source_mz = source_row[source_mz_col]
        source_rt = source_row.get(source_rt_col) if use_rt else None
        mz_tol = source_mz * mz_tol_ppm / 1e6

        candidates = target_hash.get(key, [])
        if not candidates:
            continue

        # Find best match within tolerance
        best_idx = None
        best_score = float('inf')  # Lower is better (combines mz diff and qvalue)

        for target_idx, target_mz, qvalue, target_rt in candidates:
            mz_diff = abs(source_mz - target_mz)
            if mz_diff > mz_tol:
                continue

            # Check RT if available
            if use_rt and source_rt is not None and target_rt is not None:
                rt_diff = abs(source_rt - target_rt)
                if rt_diff > rt_tol_min:
                    continue
                # Score: prefer close RT and low qvalue
                score = rt_diff / rt_tol_min + qvalue
            else:
                score = qvalue

            if score < best_score:
                best_score = score
                best_idx = target_idx

        if best_idx is not None:
            match_dict = source_row.to_dict()
            match_dict.update(target.loc[best_idx].to_dict())
            matches.append(match_dict)

    if not matches:
        return pd.DataFrame()

    result = pd.DataFrame(matches)
    # Drop helper columns
    result = result.drop(columns=["_hash_key"], errors="ignore")

    return result


def match_diann_to_precursors(
    raw_df: pd.DataFrame,
    diann_df: pd.DataFrame,
    config: Optional[MatchConfig] = None,
) -> pd.DataFrame:
    """Match DIA-NN results to raw precursors using tiered matching.

    Strategy (via PrecursorMatcher):
    1. Tier 1/2: Match by sequence + charge (highest confidence)
    2. Tier 3: Match by m/z + charge + RT + IM (full coordinate)
    3. Tier 4: Match by m/z + charge only (partial coordinate)

    This allows DIA-NN to identify precursors that FragPipe/Sage missed.
    """
    config = config or MatchConfig()
    matcher = PrecursorMatcher(config)

    raw = raw_df.copy()
    diann = diann_df.copy()

    # DIA-NN has "diann_run" column, convert to raw_file format
    diann["raw_file"] = diann["diann_run"].apply(
        lambda x: Path(x).stem.replace(".d", "") if pd.notna(x) else ""
    )

    matches = []
    tier_counts = {tier.name: 0 for tier in MatchTier}

    # Process per raw file for efficiency
    for raw_file in raw["raw_file"].unique():
        raw_file_df = raw[raw["raw_file"] == raw_file].copy()
        diann_file_df = diann[diann["raw_file"] == raw_file].copy()

        if diann_file_df.empty:
            continue

        # Build sequence index on DIA-NN (using standardized modified sequences)
        diann_seq_index = matcher.create_sequence_index(
            diann_file_df,
            sequence_col='diann_modified',
            charge_col='diann_charge',
            raw_file_col='raw_file',
        )

        # Build coordinate index on DIA-NN
        diann_coord_index = matcher.create_coordinate_index(
            diann_file_df,
            mz_col='diann_mz',
            charge_col='diann_charge',
            raw_file_col='raw_file',
        )

        matched_diann_indices = set()

        for raw_idx, raw_row in raw_file_df.iterrows():
            best_match = None
            best_diann_idx = None
            best_tier = MatchTier.NO_MATCH
            best_score = float('inf')

            # Get source sequence (prefer fragpipe_modified, fallback to fragpipe_peptide)
            raw_seq = raw_row.get('fragpipe_modified', raw_row.get('fragpipe_peptide'))
            raw_charge = raw_row.get('fragpipe_charge', raw_row.get('raw_charge'))

            # Try sequence match first (for identified precursors)
            if not pd.isna(raw_seq) and not pd.isna(raw_charge):
                norm_seq = normalize_sequence_for_matching(raw_seq, config.normalize_il)
                key = (norm_seq, int(raw_charge), raw_file)

                if key in diann_seq_index:
                    for diann_idx in diann_seq_index[key]:
                        if diann_idx in matched_diann_indices:
                            continue
                        diann_row = diann_file_df.loc[diann_idx]
                        best_match = diann_row
                        best_diann_idx = diann_idx
                        best_tier = MatchTier.SEQUENCE_IL_NORM if config.normalize_il else MatchTier.SEQUENCE_EXACT
                        best_score = 0.0
                        break  # Take first sequence match

            # If no sequence match, try coordinate-based matching
            if best_match is None:
                raw_mz = raw_row.get('raw_mz', raw_row.get('fragpipe_mz'))
                # Prefer raw RT (seconds), convert to minutes for DIA-NN comparison
                raw_rt_sec = raw_row.get('raw_rt_seconds')
                if pd.isna(raw_rt_sec):
                    # Fallback to FragPipe RT (already in minutes)
                    raw_rt_min = raw_row.get('fragpipe_rt')
                else:
                    raw_rt_min = float(raw_rt_sec) / 60.0
                # Prefer raw mobility from PASEF, fallback to FragPipe
                raw_im = raw_row.get('raw_mobility')
                if pd.isna(raw_im):
                    raw_im = raw_row.get('fragpipe_mobility')

                if not pd.isna(raw_mz) and not pd.isna(raw_charge):
                    mz_bin = int(float(raw_mz) / 0.01)
                    key = (mz_bin, int(raw_charge), raw_file)

                    if key in diann_coord_index:
                        mz_tol = float(raw_mz) * config.mz_tol_ppm / 1e6
                        rt_tol_min = config.rt_tol_sec / 60.0  # DIA-NN RT is in minutes

                        for diann_idx, diann_mz in diann_coord_index[key]:
                            if diann_idx in matched_diann_indices:
                                continue

                            mz_diff = abs(float(raw_mz) - diann_mz)
                            if mz_diff > mz_tol:
                                continue

                            diann_row = diann_file_df.loc[diann_idx]
                            score = mz_diff / mz_tol

                            # Check RT/IM for tier determination
                            has_rt_match = False
                            has_im_match = False

                            diann_rt = diann_row.get('diann_rt')  # DIA-NN RT in minutes
                            if not pd.isna(raw_rt_min) and not pd.isna(diann_rt):
                                rt_diff = abs(float(raw_rt_min) - float(diann_rt))
                                if rt_diff <= rt_tol_min:
                                    has_rt_match = True
                                    score += rt_diff / rt_tol_min
                                else:
                                    continue  # RT mismatch, skip

                            diann_im = diann_row.get('diann_mobility')
                            if not pd.isna(raw_im) and not pd.isna(diann_im):
                                im_diff = abs(float(raw_im) - float(diann_im))
                                if im_diff <= config.im_tol:
                                    has_im_match = True
                                    score += im_diff / config.im_tol
                                else:
                                    continue  # IM mismatch, skip

                            # Add q-value to score (prefer better matches)
                            qval = diann_row.get('diann_qvalue', 1.0)
                            if pd.notna(qval):
                                score += qval

                            if score < best_score:
                                best_score = score
                                best_match = diann_row
                                best_diann_idx = diann_idx
                                if has_rt_match and has_im_match:
                                    best_tier = MatchTier.COORDINATE_FULL
                                else:
                                    best_tier = MatchTier.COORDINATE_PARTIAL

            # Record match
            if best_match is not None:
                match_dict = {
                    "raw_file": raw_row["raw_file"],
                    "precursor_id": raw_row["precursor_id"],
                    "diann_match_tier": best_tier.name,
                    "diann_match_score": best_score,
                }
                for col in best_match.index:
                    if col.startswith("diann_"):
                        match_dict[col] = best_match[col]
                matches.append(match_dict)
                matched_diann_indices.add(best_diann_idx)
                tier_counts[best_tier.name] += 1

    # Report match tiers
    print(f"    Match tiers:")
    for tier, count in tier_counts.items():
        if count > 0:
            print(f"      {tier}: {count}")

    if not matches:
        return pd.DataFrame()

    return pd.DataFrame(matches)


def match_sage_to_precursors(
    fp_df: pd.DataFrame,
    sage_df: pd.DataFrame,
    config: Optional[MatchConfig] = None,
) -> pd.DataFrame:
    """Match Sage results to FragPipe precursors using tiered matching.

    Strategy (via PrecursorMatcher):
    1. Tier 1/2: Match by sequence + charge (highest confidence)
    2. Tier 3: Match by m/z + charge + RT + IM (full coordinate)
    3. Tier 4: Match by m/z + charge only (partial coordinate)

    Matches are done WITHIN each raw file.

    Returns DataFrame with Sage columns added and match quality info.
    """
    config = config or MatchConfig()
    matcher = PrecursorMatcher(config)

    fp = fp_df.copy()
    sage = sage_df.copy()

    # Rename sage raw_file column for matching
    sage["raw_file"] = sage["sage_raw_file"]

    # Define column mappings for matcher
    fp_cols = {
        'sequence': 'fragpipe_peptide',
        'modified': 'fragpipe_modified',
        'charge': 'fragpipe_charge',
        'mz': 'fragpipe_mz',
        'rt': 'fragpipe_rt',
        'mobility': 'fragpipe_mobility',
    }
    sage_cols = {
        'sequence': 'sage_peptide',
        'modified': 'sage_modified',
        'charge': 'sage_charge',
        'mz': 'sage_mz',
        'rt': 'sage_rt',
        'mobility': 'sage_mobility',
    }

    # Match per raw file for efficiency
    all_matches = []
    tier_counts = {tier.name: 0 for tier in MatchTier}

    for raw_file in fp["raw_file"].unique():
        fp_file = fp[fp["raw_file"] == raw_file].copy()
        sage_file = sage[sage["raw_file"] == raw_file].copy()

        if sage_file.empty or fp_file.empty:
            continue

        # Build sequence index on Sage
        sage_seq_index = matcher.create_sequence_index(
            sage_file,
            sequence_col='sage_modified',
            charge_col='sage_charge',
            raw_file_col='raw_file',
        )

        # Build coordinate index on Sage
        sage_coord_index = matcher.create_coordinate_index(
            sage_file,
            mz_col='sage_mz',
            charge_col='sage_charge',
            raw_file_col='raw_file',
        )

        matched_sage_indices = set()

        for fp_idx, fp_row in fp_file.iterrows():
            # Try sequence match first
            fp_seq = fp_row.get('fragpipe_modified', fp_row.get('fragpipe_peptide'))
            fp_charge = fp_row.get('fragpipe_charge')

            best_match = None
            best_sage_idx = None
            best_tier = MatchTier.NO_MATCH

            if not pd.isna(fp_seq) and not pd.isna(fp_charge):
                norm_seq = normalize_sequence_for_matching(fp_seq, config.normalize_il)
                key = (norm_seq, int(fp_charge), raw_file)

                if key in sage_seq_index:
                    for sage_idx in sage_seq_index[key]:
                        if sage_idx in matched_sage_indices:
                            continue
                        sage_row = sage_file.loc[sage_idx]
                        best_match = sage_row
                        best_sage_idx = sage_idx
                        best_tier = MatchTier.SEQUENCE_IL_NORM if config.normalize_il else MatchTier.SEQUENCE_EXACT
                        break  # Take first sequence match

            # If no sequence match, try coordinate
            if best_match is None:
                fp_mz = fp_row.get('fragpipe_mz')

                if not pd.isna(fp_mz) and not pd.isna(fp_charge):
                    mz_bin = int(float(fp_mz) / 0.01)
                    key = (mz_bin, int(fp_charge), raw_file)

                    if key in sage_coord_index:
                        best_score = float('inf')
                        mz_tol = float(fp_mz) * config.mz_tol_ppm / 1e6

                        for sage_idx, sage_mz in sage_coord_index[key]:
                            if sage_idx in matched_sage_indices:
                                continue

                            mz_diff = abs(float(fp_mz) - sage_mz)
                            if mz_diff > mz_tol:
                                continue

                            sage_row = sage_file.loc[sage_idx]
                            score = mz_diff / mz_tol

                            # Check RT/IM for tier determination
                            has_rt_match = False
                            has_im_match = False

                            fp_rt = fp_row.get('fragpipe_rt')
                            sage_rt = sage_row.get('sage_rt')
                            if not pd.isna(fp_rt) and not pd.isna(sage_rt):
                                rt_diff = abs(float(fp_rt) - float(sage_rt))
                                if rt_diff <= config.rt_tol_sec:
                                    has_rt_match = True
                                    score += rt_diff / config.rt_tol_sec

                            fp_im = fp_row.get('fragpipe_mobility')
                            sage_im = sage_row.get('sage_mobility')
                            if not pd.isna(fp_im) and not pd.isna(sage_im):
                                im_diff = abs(float(fp_im) - float(sage_im))
                                if im_diff <= config.im_tol:
                                    has_im_match = True
                                    score += im_diff / config.im_tol

                            if score < best_score:
                                best_score = score
                                best_match = sage_row
                                best_sage_idx = sage_idx
                                if has_rt_match and has_im_match:
                                    best_tier = MatchTier.COORDINATE_FULL
                                else:
                                    best_tier = MatchTier.COORDINATE_PARTIAL

            if best_match is not None:
                match_dict = fp_row.to_dict()
                for col in best_match.index:
                    if col.startswith("sage_"):
                        match_dict[col] = best_match[col]
                match_dict["sage_match_tier"] = best_tier.name
                all_matches.append(match_dict)
                matched_sage_indices.add(best_sage_idx)
                tier_counts[best_tier.name] += 1

    # Report match tiers
    print(f"    Match tiers:")
    for tier, count in tier_counts.items():
        if count > 0:
            print(f"      {tier}: {count}")

    if not all_matches:
        return pd.DataFrame()

    result = pd.DataFrame(all_matches)
    return result


def build_unified_index(
    accession: str,
    base_dir: Path,
    local_data_dir: Optional[Path] = None,
    raw_file_filter: Optional[str] = None,
    load_raw: bool = True,
    include_unidentified: bool = True,
    match_config: Optional[MatchConfig] = None,
    use_calibration: bool = True,
) -> pd.DataFrame:
    """Build unified precursor index combining all sources.

    Args:
        include_unidentified: If True, include ALL fragmented precursors from raw data,
                              not just those identified by search engines.
        match_config: Configuration for precursor matching tolerances.
                      Defaults: mz_tol_ppm=20, rt_tol_sec=30, im_tol=0.05
        use_calibration: Use pre-computed IM calibration for accurate mobility values.
    """
    # Initialize match config with defaults if not provided
    match_config = match_config or MatchConfig()

    print(f"\nBuilding unified precursor index for {accession}")
    print(f"  Match config: mz_tol={match_config.mz_tol_ppm}ppm, rt_tol={match_config.rt_tol_sec}s, im_tol={match_config.im_tol}")
    print("=" * 60)

    # Load search engine results
    print("\nLoading FragPipe PSMs...")
    fp_df = load_fragpipe_psms(accession, base_dir, raw_file_filter)
    print(f"  Loaded {len(fp_df)} FragPipe PSMs")

    print("\nLoading DIA-NN report...")
    diann_df = load_diann_report(accession, base_dir, raw_file_filter)
    print(f"  Loaded {len(diann_df)} DIA-NN precursors")

    print("\nLoading Sage results...")
    sage_df = load_sage_results(accession, base_dir, raw_file_filter)
    print(f"  Loaded {len(sage_df)} Sage PSMs")

    # Load raw precursors if including unidentified
    raw_df = pd.DataFrame()
    if include_unidentified and local_data_dir and local_data_dir.exists():
        print("\nLoading ALL raw precursors...")
        raw_df = load_raw_precursors(
            accession, local_data_dir, raw_file_filter,
            use_calibration=use_calibration
        )
        print(f"  Loaded {len(raw_df)} raw precursors (ALL fragmented)")

    # Start with FragPipe or raw precursors
    if not raw_df.empty and include_unidentified:
        # Start from ALL raw precursors, then left join search results
        index_df = raw_df.copy()
        print(f"\nStarting from {len(index_df)} raw precursors")

        # Merge FragPipe results
        if not fp_df.empty:
            index_df = index_df.merge(
                fp_df,
                on=["raw_file", "precursor_id"],
                how="left"
            )
            print(f"  After FragPipe merge: {(index_df['fragpipe_peptide'].notna()).sum()} identified")
    elif not fp_df.empty:
        # Fall back to FragPipe-only (old behavior)
        index_df = fp_df.copy()
        print(f"\nStarting from {len(index_df)} FragPipe PSMs")
    else:
        print("\nNo data - cannot build index")
        return pd.DataFrame()

    # Use raw charge if available, otherwise FragPipe
    if "raw_charge" in index_df.columns and "fragpipe_charge" not in index_df.columns:
        index_df["fragpipe_charge"] = index_df["raw_charge"]
    elif "fragpipe_charge" not in index_df.columns:
        index_df["fragpipe_charge"] = 2  # Default

    # Match DIA-NN using tiered strategy
    if not diann_df.empty and "fragpipe_peptide" in index_df.columns:
        print("\nMatching DIA-NN to precursors (tiered strategy)...")
        diann_matches = match_diann_to_precursors(index_df, diann_df, config=match_config)
        if not diann_matches.empty:
            print(f"  Matched {len(diann_matches)} DIA-NN precursors")
            diann_cols = [c for c in diann_matches.columns if c.startswith("diann_")]
            index_df = index_df.merge(
                diann_matches[["raw_file", "precursor_id"] + diann_cols],
                on=["raw_file", "precursor_id"],
                how="left",
                suffixes=("", "_dup")
            )
            # Remove duplicate columns if any
            index_df = index_df.loc[:, ~index_df.columns.str.endswith("_dup")]

    # Merge additional raw metadata if not already done
    if load_raw and local_data_dir and local_data_dir.exists() and "raw_intensity" not in index_df.columns:
        print("\nLoading additional raw precursor metadata...")
        extra_raw_df = load_raw_precursors(accession, local_data_dir, raw_file_filter)
        if not extra_raw_df.empty:
            # Only add columns not already present
            new_cols = [c for c in extra_raw_df.columns if c not in index_df.columns and c not in ["raw_file", "precursor_id"]]
            if new_cols:
                index_df = index_df.merge(
                    extra_raw_df[["raw_file", "precursor_id"] + new_cols],
                    on=["raw_file", "precursor_id"],
                    how="left"
                )
                print(f"  Added columns: {new_cols}")

    # Match Sage using tiered strategy
    if not sage_df.empty and "fragpipe_peptide" in index_df.columns:
        print("\nMatching Sage to precursors (tiered strategy)...")
        sage_matches = match_sage_to_precursors(index_df, sage_df, config=match_config)
        if not sage_matches.empty:
            print(f"  Matched {len(sage_matches)} Sage precursors")
            # sage_matches already has all FP + Sage columns merged
            # Keep only Sage columns from matches, merge back to index_df
            sage_cols = [c for c in sage_matches.columns if c.startswith("sage_") or c == "_match_type"]
            sage_subset = sage_matches[["raw_file", "precursor_id"] + sage_cols].drop_duplicates(
                subset=["raw_file", "precursor_id"]
            )
            index_df = index_df.merge(sage_subset, on=["raw_file", "precursor_id"], how="left")

    # Calculate consensus metrics
    print("\nCalculating consensus metrics...")

    # Normalize peptide sequences for matching (I/L normalized, plain sequence)
    for col in ["fragpipe_peptide", "sage_peptide", "diann_peptide"]:
        if col in index_df.columns:
            index_df[f"{col}_norm"] = index_df[col].apply(
                lambda x: normalize_sequence_il(str(x)) if pd.notna(x) else ""
            )

    # Also normalize modified sequences for consensus
    for col in ["fragpipe_modified", "sage_modified", "diann_modified"]:
        if col in index_df.columns:
            index_df[f"{col}_norm"] = index_df[col].apply(
                lambda x: normalize_sequence_il(str(x)) if pd.notna(x) else ""
            )

    # Count engines
    def count_engines(row):
        n = 0
        if pd.notna(row.get("fragpipe_peptide")) and row.get("fragpipe_peptide"):
            n += 1
        if pd.notna(row.get("diann_peptide")) and row.get("diann_peptide"):
            n += 1
        if pd.notna(row.get("sage_peptide")) and row.get("sage_peptide"):
            n += 1
        return n

    index_df["n_engines"] = index_df.apply(count_engines, axis=1)

    # Get consensus peptide (from most engines, prefer FragPipe if tie)
    def get_consensus_peptide(row):
        peptides = []
        if pd.notna(row.get("fragpipe_peptide_norm")) and row.get("fragpipe_peptide_norm"):
            peptides.append(("fragpipe", row["fragpipe_peptide_norm"], row.get("fragpipe_probability", 0)))
        if pd.notna(row.get("diann_peptide_norm")) and row.get("diann_peptide_norm"):
            peptides.append(("diann", row["diann_peptide_norm"], 1 - row.get("diann_qvalue", 1)))
        if pd.notna(row.get("sage_peptide_norm")) and row.get("sage_peptide_norm"):
            peptides.append(("sage", row["sage_peptide_norm"], 1 - row.get("sage_qvalue", 1)))

        if not peptides:
            return "", 0.0

        # Check agreement
        unique_seqs = set(p[1] for p in peptides)
        if len(unique_seqs) == 1:
            # All agree
            return peptides[0][1], 1.0
        else:
            # Disagreement - return most confident
            best = max(peptides, key=lambda x: x[2])
            return best[1], len(peptides) / (len(unique_seqs) + 1)

    consensus = index_df.apply(get_consensus_peptide, axis=1)
    index_df["consensus_peptide"] = consensus.apply(lambda x: x[0])
    index_df["confidence_weight"] = consensus.apply(lambda x: x[1])

    # Select and order final columns
    final_columns = [
        "raw_file", "precursor_id",
        # Raw properties (from timsTOF metadata + PASEF)
        "raw_mz", "raw_charge", "raw_rt_seconds", "raw_mobility", "raw_intensity",
        "frame_id", "isolation_mz", "isolation_width", "parent_id",
        # Consensus
        "consensus_peptide", "n_engines", "confidence_weight",
        # FragPipe (standardized UNIMOD format)
        "fragpipe_peptide", "fragpipe_modified", "fragpipe_protein",
        "fragpipe_probability", "fragpipe_pep", "fragpipe_hyperscore", "fragpipe_qvalue",
        "fragpipe_mz", "fragpipe_charge", "fragpipe_rt", "fragpipe_mobility",
        # DIA-NN (standardized UNIMOD format) + match quality
        "diann_peptide", "diann_modified", "diann_protein",
        "diann_qvalue", "diann_global_qvalue", "diann_pg_qvalue", "diann_pep", "diann_ccs",
        "diann_mz", "diann_charge", "diann_rt", "diann_mobility",
        "diann_match_tier", "diann_match_score",
        # Sage (standardized UNIMOD format) + match quality
        "sage_peptide", "sage_modified", "sage_protein",
        "sage_qvalue", "sage_peptide_qvalue", "sage_protein_qvalue", "sage_pep", "sage_hyperscore",
        "sage_mz", "sage_charge", "sage_rt", "sage_mobility",
        "sage_match_tier",
    ]

    # Only include columns that exist
    final_columns = [c for c in final_columns if c in index_df.columns]
    index_df = index_df[final_columns]

    # Sort by raw_file and precursor_id
    index_df = index_df.sort_values(["raw_file", "precursor_id"])

    return index_df


def main():
    parser = argparse.ArgumentParser(
        description="Build unified precursor index mapping raw data to search engine IDs"
    )
    parser.add_argument(
        "--accession", "-a",
        required=True,
        help="PRIDE accession"
    )
    parser.add_argument(
        "--base-dir", "-d",
        type=Path,
        default=Path("data/processed"),
        help="Base directory for processed data"
    )
    parser.add_argument(
        "--local-data",
        type=Path,
        help="Path to local raw data directory (with .d files)"
    )
    parser.add_argument(
        "--raw-file",
        help="Filter to specific raw file (partial match)"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        help="Output parquet file path"
    )
    parser.add_argument(
        "--no-raw",
        action="store_true",
        help="Skip loading raw precursor metadata"
    )
    parser.add_argument(
        "--include-unidentified",
        action="store_true",
        help="Include ALL fragmented precursors, not just identified ones"
    )
    parser.add_argument(
        "--no-calibration",
        action="store_true",
        help="Skip IM calibration (use linear interpolation - less accurate)"
    )

    args = parser.parse_args()

    # Build index
    index_df = build_unified_index(
        accession=args.accession,
        base_dir=args.base_dir,
        local_data_dir=args.local_data,
        raw_file_filter=args.raw_file,
        load_raw=not args.no_raw,
        include_unidentified=args.include_unidentified,
        use_calibration=not args.no_calibration,
    )

    if index_df.empty:
        print("\nNo data found - index is empty")
        return

    # Summary
    print(f"\n{'='*60}")
    print("Summary:")
    print(f"  Total precursors: {len(index_df):,}")
    if "fragpipe_peptide" in index_df.columns:
        identified = (index_df['fragpipe_peptide'].notna()).sum()
        unidentified = len(index_df) - identified
        print(f"  With FragPipe ID: {identified:,}")
        print(f"  Unidentified:     {unidentified:,} ({100*unidentified/len(index_df):.1f}%)")
    if "diann_peptide" in index_df.columns:
        print(f"  With DIA-NN ID:   {(index_df['diann_peptide'].notna()).sum():,}")
    if "sage_peptide" in index_df.columns:
        print(f"  With Sage ID:     {(index_df['sage_peptide'].notna()).sum():,}")
    if "raw_intensity" in index_df.columns:
        print(f"\n  Intensity range: {index_df['raw_intensity'].min():.0f} - {index_df['raw_intensity'].max():.0f}")
    print(f"\n  By number of engines:")
    for n in sorted(index_df["n_engines"].unique()):
        count = (index_df["n_engines"] == n).sum()
        print(f"    {n} engine(s): {count:,}")

    # Report match tier statistics
    if "diann_match_tier" in index_df.columns:
        print(f"\n  DIA-NN match tiers:")
        for tier, count in index_df["diann_match_tier"].value_counts().items():
            print(f"    {tier}: {count:,}")

    if "sage_match_tier" in index_df.columns:
        print(f"\n  Sage match tiers:")
        for tier, count in index_df["sage_match_tier"].value_counts().items():
            print(f"    {tier}: {count:,}")

    # Save
    if args.output:
        output_path = args.output
    else:
        output_path = args.base_dir / args.accession / "precursor_index.parquet"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    index_df.to_parquet(output_path, index=False)
    print(f"\nSaved to: {output_path}")


if __name__ == "__main__":
    main()
