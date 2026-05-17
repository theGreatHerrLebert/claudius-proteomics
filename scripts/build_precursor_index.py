#!/usr/bin/env python3
"""
Build Unified Precursor Index

Creates a single parquet file that maps:
  precursor_id (raw data) -> search engine identifications + quality metrics

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
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np

# Import modular components
from engine_parsers import FragPipeParser, DiannParser, SageParser
from precursor_merging import MatchConfig, FragPipeAnchoredMerger


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
                from imspy_connector import py_dda
                rust_dataset = py_dda.PyTimsDatasetDDA.with_calibration(
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

            # Geometric MIDPOINT of the PASEF isolation-scan window — NOT an
            # intensity apex. This yields the coarse 'raw_mobility' column, which
            # is biased toward the window centre and must NOT be used as a
            # measured 1/K0. The real ion-mobility measurement is the 'mobility'
            # column from extract_precursors.py (scan-marginal intensity apex).
            scan_center = ((pasef['scan_begin'] + pasef['scan_end']) / 2).astype(np.int32)

            # Convert scans to mobility
            mobilities = []
            for frame_id in pasef['frame_id'].unique():
                mask = pasef['frame_id'] == frame_id
                scans = scan_center[mask].values
                try:
                    if use_calibration and rust_dataset is not None:
                        im_values = rust_dataset.scan_to_inverse_mobility(
                            int(frame_id), [int(s) for s in scans]
                        )
                    else:
                        im_values = dataset.scan_to_inverse_mobility(int(frame_id), scans)
                    mobilities.extend(im_values)
                except Exception:
                    mobilities.extend([np.nan] * len(scans))

            pasef['mobility'] = mobilities
            pasef['rt_seconds'] = pasef['Time']

            # Aggregate PASEF events per precursor
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
        "raw_mobility": df.get("mobility"),
        "raw_rt_seconds": df.get("rt_seconds"),
        "raw_intensity": df.get("intensity", df.get("precuror_total_intensity")),
        "frame_id": df.get("frame_id"),
        "isolation_mz": df.get("isolation_mz"),
        "isolation_width": df.get("isolation_width"),
        "parent_id": df.get("parent_id"),
    })

    return result


def build_unified_index(
    accession: str,
    base_dir: Path,
    local_data_dir: Optional[Path] = None,
    raw_file_filter: Optional[str] = None,
    include_unidentified: bool = True,
    match_config: Optional[MatchConfig] = None,
    use_calibration: bool = True,
) -> pd.DataFrame:
    """Build unified precursor index combining all sources.

    Args:
        accession: PRIDE accession
        base_dir: Base directory for processed data
        local_data_dir: Path to local raw data directory (with .d files)
        raw_file_filter: Filter to specific raw file (partial match)
        include_unidentified: Include ALL fragmented precursors from raw data
        match_config: Configuration for precursor matching tolerances
        use_calibration: Use pre-computed IM calibration for accurate mobility

    Returns:
        Unified precursor index DataFrame
    """
    match_config = match_config or MatchConfig()

    print(f"\nBuilding unified precursor index for {accession}")
    print("=" * 60)

    # Initialize parsers
    fragpipe_parser = FragPipeParser()
    diann_parser = DiannParser()
    sage_parser = SageParser()

    # Load search engine results
    print("\nLoading FragPipe PSMs...")
    fp_df = fragpipe_parser.parse(base_dir, accession, raw_file_filter)
    print(f"  Loaded {len(fp_df)} FragPipe PSMs")

    print("\nLoading DIA-NN report...")
    diann_df = diann_parser.parse(base_dir, accession, raw_file_filter)
    print(f"  Loaded {len(diann_df)} DIA-NN precursors")

    print("\nLoading Sage results...")
    sage_df = sage_parser.parse(base_dir, accession, raw_file_filter)
    print(f"  Loaded {len(sage_df)} Sage PSMs")

    # Load raw precursors
    raw_df = pd.DataFrame()
    if include_unidentified and local_data_dir and local_data_dir.exists():
        print("\nLoading ALL raw precursors...")
        raw_df = load_raw_precursors(
            accession, local_data_dir, raw_file_filter,
            use_calibration=use_calibration
        )
        print(f"  Loaded {len(raw_df)} raw precursors (ALL fragmented)")

    # Initialize merger
    merger = FragPipeAnchoredMerger(match_config)

    # Handle case where no raw data
    if raw_df.empty and not fp_df.empty:
        print(f"\nNo raw data available, starting from FragPipe PSMs")
        # Add required raw columns from FragPipe
        raw_df = pd.DataFrame({
            "raw_file": fp_df["raw_file"],
            "precursor_id": fp_df["precursor_id"],
            "raw_mz": fp_df.get("fragpipe_mz"),
            "raw_charge": fp_df.get("fragpipe_charge"),
            "raw_rt_seconds": fp_df.get("fragpipe_rt"),
            "raw_mobility": fp_df.get("fragpipe_mobility"),
        })
        # Remove FragPipe columns that will be re-merged
        fp_df_clean = fp_df.drop(columns=["raw_file", "precursor_id"], errors="ignore")
        raw_df = pd.concat([raw_df, fp_df_clean], axis=1)
        # Clear fp_df so merger doesn't re-merge
        fp_df = pd.DataFrame()

    if raw_df.empty:
        print("\nNo data - cannot build index")
        return pd.DataFrame()

    # Merge all sources
    print("\nMerging search engine results...")
    index_df = merger.merge(raw_df, fp_df, diann_df, sage_df)

    # Select and order final columns
    final_columns = [
        "raw_file", "precursor_id",
        # Raw properties
        "raw_mz", "raw_charge", "raw_rt_seconds", "raw_mobility", "raw_intensity",
        "frame_id", "isolation_mz", "isolation_width", "parent_id",
        # Consensus
        "consensus_peptide", "n_engines", "confidence_weight",
        # FragPipe
        "fragpipe_peptide", "fragpipe_modified", "fragpipe_protein",
        "fragpipe_probability", "fragpipe_pep", "fragpipe_hyperscore", "fragpipe_qvalue",
        "fragpipe_mz", "fragpipe_charge", "fragpipe_rt", "fragpipe_mobility",
        # DIA-NN
        "diann_peptide", "diann_modified", "diann_protein",
        "diann_qvalue", "diann_global_qvalue", "diann_pg_qvalue", "diann_pep", "diann_ccs",
        "diann_mz", "diann_charge", "diann_rt", "diann_mobility",
        "diann_match_tier",
        # Sage
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


def print_summary(index_df: pd.DataFrame):
    """Print summary statistics for the index."""
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

    if "diann_match_tier" in index_df.columns:
        print(f"\n  DIA-NN match tiers:")
        for tier, count in index_df["diann_match_tier"].value_counts().items():
            print(f"    {tier}: {count:,}")

    if "sage_match_tier" in index_df.columns:
        print(f"\n  Sage match tiers:")
        for tier, count in index_df["sage_match_tier"].value_counts().items():
            print(f"    {tier}: {count:,}")


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
        "--include-unidentified",
        action="store_true",
        help="Include ALL fragmented precursors, not just identified ones"
    )
    parser.add_argument(
        "--no-calibration",
        action="store_true",
        help="Skip IM calibration (use linear interpolation - less accurate)"
    )
    parser.add_argument(
        "--mz-tol",
        type=float,
        default=20.0,
        help="m/z tolerance in ppm (default: 20)"
    )
    parser.add_argument(
        "--rt-tol",
        type=float,
        default=30.0,
        help="RT tolerance in seconds (default: 30)"
    )
    parser.add_argument(
        "--im-tol",
        type=float,
        default=0.05,
        help="Ion mobility tolerance in 1/K0 (default: 0.05)"
    )

    args = parser.parse_args()

    # Build match config from arguments
    match_config = MatchConfig(
        mz_tol_ppm=args.mz_tol,
        rt_tol_sec=args.rt_tol,
        im_tol=args.im_tol,
    )

    # Build index
    index_df = build_unified_index(
        accession=args.accession,
        base_dir=args.base_dir,
        local_data_dir=args.local_data,
        raw_file_filter=args.raw_file,
        include_unidentified=args.include_unidentified,
        match_config=match_config,
        use_calibration=not args.no_calibration,
    )

    if index_df.empty:
        print("\nNo data found - index is empty")
        return

    # Summary
    print_summary(index_df)

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
