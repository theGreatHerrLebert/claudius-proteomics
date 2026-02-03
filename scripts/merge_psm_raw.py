#!/usr/bin/env python3
"""
Merge search engine identifications with raw extracted features.

Supports three search engines:
- FragPipe: Join on (raw_file, precursor_id parsed from Spectrum column)
- Sage: Join on (filename, scannr) where scannr = precursor_id
- DIA-NN: Match by (raw_file, RT, m/z, charge) within tolerances

Creates a unified dataset with:
- Peptide sequences (UniMod normalized)
- Search engine scores
- Raw signal features (MS1 chromatogram, mobilogram, isotopes)
- Fragment spectrum reference (blob_offset, blob_size)
"""

import sys
import re
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np

# Import sequence standardization utilities
from sequence_utils import (
    standardize_fragpipe_modified_peptide,
    standardize_sage_sequence,
    standardize_diann_sequence,
)


def setup_logging(log_path: Path):
    """Set up logging to file and stdout."""
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


def load_fragpipe_psms(psm_dir: Path, logger) -> pd.DataFrame:
    """
    Load FragPipe PSMs from per-file psm.tsv files.

    Parses precursor_id from Spectrum column format:
    {raw_file}.{precursor_id}.{precursor_id}.{charge}
    """
    psm_files = list(psm_dir.glob("*/psm.tsv"))
    if not psm_files:
        logger.warning(f"No FragPipe PSM files found in {psm_dir}")
        return pd.DataFrame()

    dfs = []
    for psm_file in psm_files:
        df = pd.read_csv(psm_file, sep='\t')
        dfs.append(df)

    psms = pd.concat(dfs, ignore_index=True)
    logger.info(f"Loaded {len(psms):,} FragPipe PSMs from {len(psm_files)} files")

    # Parse precursor_id and raw_file from Spectrum column
    # Format: "20190504_TIMS1_FlMe_SA_HeLa_frac01_A10_1_93.00014.00014.2"
    spectrum_pattern = r'^(.+?)\.(\d+)\.\d+\.\d+$'
    matches = psms['Spectrum'].str.extract(spectrum_pattern)
    psms['raw_file'] = matches[0] + '.d'
    psms['precursor_id'] = matches[1].astype(int)

    # Normalize column names
    psms = psms.rename(columns={
        'Peptide': 'sequence',
        'Modified Peptide': 'modified_sequence',
        'Charge': 'charge',
        'Retention': 'rt_minutes',
        'Ion Mobility': 'mobility',
        'Observed M/Z': 'mz',
        'Hyperscore': 'fragpipe_hyperscore',
        'Expectation': 'fragpipe_expectation',
        'Probability': 'fragpipe_probability',
    })

    psms['rt_seconds'] = psms['rt_minutes'] * 60
    psms['engine'] = 'fragpipe'

    # Standardize modified_sequence to UNIMOD format
    psms['modified_sequence_unimod'] = psms['modified_sequence'].apply(
        lambda x: standardize_fragpipe_modified_peptide(x) if pd.notna(x) else ''
    )

    return psms[['raw_file', 'precursor_id', 'sequence', 'modified_sequence',
                 'modified_sequence_unimod', 'charge',
                 'mz', 'rt_seconds', 'mobility', 'fragpipe_hyperscore',
                 'fragpipe_expectation', 'fragpipe_probability', 'engine']]


def load_sage_psms(sage_dir: Path, logger) -> pd.DataFrame:
    """
    Load Sage PSMs from results.sage.parquet.

    scannr column = precursor_id
    """
    sage_file = sage_dir / "results.sage.parquet"
    if not sage_file.exists():
        logger.warning(f"Sage results not found: {sage_file}")
        return pd.DataFrame()

    psms = pd.read_parquet(sage_file)
    logger.info(f"Loaded {len(psms):,} Sage PSMs")

    # Rename columns
    psms = psms.rename(columns={
        'filename': 'raw_file',
        'scannr': 'precursor_id',
        'peptide': 'modified_sequence',
        'stripped_peptide': 'sequence',
        'charge': 'charge',
        'expmass': 'exp_mass',
        'hyperscore': 'sage_hyperscore',
    })

    psms['precursor_id'] = psms['precursor_id'].astype(int)
    psms['engine'] = 'sage'

    # Standardize modified_sequence to UNIMOD format
    psms['modified_sequence_unimod'] = psms['modified_sequence'].apply(
        lambda x: standardize_sage_sequence(x) if pd.notna(x) else ''
    )

    # Sage doesn't have RT/mz directly in output, will get from raw features
    return psms[['raw_file', 'precursor_id', 'sequence', 'modified_sequence',
                 'modified_sequence_unimod', 'charge', 'sage_hyperscore', 'engine']]


def load_diann_psms(diann_dir: Path, logger) -> pd.DataFrame:
    """
    Load DIA-NN results from report.parquet.

    DIA-NN doesn't have precursor_id, so we return RT/mz/IM for tolerance matching.
    """
    diann_file = diann_dir / "report.parquet"
    if not diann_file.exists():
        logger.warning(f"DIA-NN results not found: {diann_file}")
        return pd.DataFrame()

    psms = pd.read_parquet(diann_file)
    logger.info(f"Loaded {len(psms):,} DIA-NN precursors")

    # Rename columns
    psms = psms.rename(columns={
        'Run': 'raw_file_base',
        'Modified.Sequence': 'modified_sequence',
        'Stripped.Sequence': 'sequence',
        'Precursor.Charge': 'charge',
        'Precursor.Mz': 'mz',
        'RT': 'rt_minutes',
        'IM': 'mobility',
        'Lib.Q.Value': 'diann_qvalue',
    })

    # Add .d extension to raw file
    psms['raw_file'] = psms['raw_file_base'] + '.d'
    psms['rt_seconds'] = psms['rt_minutes'] * 60
    psms['engine'] = 'diann'

    # Standardize modified_sequence to UNIMOD format
    psms['modified_sequence_unimod'] = psms['modified_sequence'].apply(
        lambda x: standardize_diann_sequence(x) if pd.notna(x) else ''
    )

    # DIA-NN needs tolerance-based matching
    psms['precursor_id'] = None  # Will be filled during merge

    return psms[['raw_file', 'precursor_id', 'sequence', 'modified_sequence',
                 'modified_sequence_unimod', 'charge', 'mz', 'rt_seconds', 'mobility', 'engine']]


def load_raw_features(features_path: Path, logger) -> pd.DataFrame:
    """Load raw extracted features."""
    features = pd.read_parquet(features_path)
    logger.info(f"Loaded {len(features):,} raw feature records")
    return features


def merge_by_precursor_id(psms: pd.DataFrame, features: pd.DataFrame, logger) -> pd.DataFrame:
    """
    Merge PSMs with features by (raw_file, precursor_id).

    Works for FragPipe and Sage which have direct precursor_id.
    """
    logger.info(f"Merging {len(psms):,} PSMs by precursor_id...")

    merged = psms.merge(
        features,
        on=['raw_file', 'precursor_id'],
        how='inner',
        suffixes=('', '_raw')
    )

    logger.info(f"Matched {len(merged):,} PSMs ({len(merged)/len(psms)*100:.1f}%)")
    return merged


def merge_by_tolerance(psms: pd.DataFrame, features: pd.DataFrame,
                       rt_tol_sec: float = 30.0, mz_tol_ppm: float = 20.0,
                       im_tol: float = 0.05, logger=None) -> pd.DataFrame:
    """
    Merge PSMs with features by RT/mz/IM tolerance matching.

    Used for DIA-NN which doesn't have precursor_id.
    """
    if logger:
        logger.info(f"Merging {len(psms):,} PSMs by tolerance matching...")
        logger.info(f"  RT tolerance: {rt_tol_sec}s, m/z tolerance: {mz_tol_ppm}ppm, IM tolerance: {im_tol}")

    matched = []

    for raw_file in psms['raw_file'].unique():
        psms_file = psms[psms['raw_file'] == raw_file]
        features_file = features[features['raw_file'] == raw_file]

        if len(features_file) == 0:
            continue

        for _, psm in psms_file.iterrows():
            # Find matching features within tolerances
            rt_match = np.abs(features_file['rt_seconds'] - psm['rt_seconds']) <= rt_tol_sec
            mz_match = np.abs((features_file['mz'] - psm['mz']) / psm['mz'] * 1e6) <= mz_tol_ppm
            im_match = np.abs(features_file['mobility'] - psm['mobility']) <= im_tol
            charge_match = features_file['charge'] == psm['charge']

            matches = features_file[rt_match & mz_match & im_match & charge_match]

            if len(matches) == 1:
                # Unique match
                match = matches.iloc[0].to_dict()
                match.update(psm.to_dict())
                matched.append(match)
            elif len(matches) > 1:
                # Multiple matches - take closest by RT
                rt_diffs = np.abs(matches['rt_seconds'] - psm['rt_seconds'])
                best_idx = rt_diffs.idxmin()
                match = matches.loc[best_idx].to_dict()
                match.update(psm.to_dict())
                matched.append(match)

    result = pd.DataFrame(matched)
    if logger:
        logger.info(f"Matched {len(result):,} PSMs ({len(result)/len(psms)*100:.1f}%)")

    return result


def merge_diann_by_sequence(diann_psms: pd.DataFrame, existing_merged: pd.DataFrame,
                            features: pd.DataFrame, logger) -> pd.DataFrame:
    """
    Merge DIA-NN results using a two-step approach:

    1. First: Match by (raw_file, sequence, charge) to existing PSMs → inherit precursor_id
    2. Fallback: For unmatched, use RT/mz/IM tolerance matching to raw features

    This is more reliable than pure tolerance matching.
    """
    logger.info(f"Merging {len(diann_psms):,} DIA-NN precursors...")

    # Step 1: Match by modified_sequence_unimod+charge to existing merged PSMs
    # Build lookup from existing merged data (use UNIMOD-standardized sequences)
    if len(existing_merged) > 0:
        seq_lookup = existing_merged[['raw_file', 'modified_sequence_unimod', 'charge', 'precursor_id']].drop_duplicates()
        seq_lookup = seq_lookup.rename(columns={'precursor_id': 'matched_precursor_id'})

        diann_with_seq = diann_psms.merge(
            seq_lookup,
            on=['raw_file', 'modified_sequence_unimod', 'charge'],
            how='left'
        )

        # Split into matched and unmatched
        matched_by_seq = diann_with_seq[diann_with_seq['matched_precursor_id'].notna()].copy()
        unmatched = diann_with_seq[diann_with_seq['matched_precursor_id'].isna()].copy()

        logger.info(f"  Matched by sequence+charge: {len(matched_by_seq):,}")
        logger.info(f"  Unmatched (DIA-NN unique): {len(unmatched):,}")

        # For matched, use the precursor_id from existing PSMs
        if len(matched_by_seq) > 0:
            matched_by_seq['precursor_id'] = matched_by_seq['matched_precursor_id'].astype(int)
            matched_by_seq = matched_by_seq.drop(columns=['matched_precursor_id'])
            matched_merged = merge_by_precursor_id(matched_by_seq, features, logger)
        else:
            matched_merged = pd.DataFrame()
    else:
        unmatched = diann_psms.copy()
        matched_merged = pd.DataFrame()
        logger.info("  No existing PSMs to match against")

    # Step 2: For unmatched, use tolerance matching
    if len(unmatched) > 0:
        unmatched = unmatched.drop(columns=['matched_precursor_id'], errors='ignore')
        logger.info(f"  Tolerance matching for {len(unmatched):,} DIA-NN unique precursors...")
        unmatched_merged = merge_by_tolerance(unmatched, features, logger=logger)
    else:
        unmatched_merged = pd.DataFrame()

    # Combine
    result_parts = [df for df in [matched_merged, unmatched_merged] if len(df) > 0]
    if result_parts:
        result = pd.concat(result_parts, ignore_index=True)
        logger.info(f"  Total DIA-NN merged: {len(result):,}")
        return result
    else:
        return pd.DataFrame()


def merge_all_engines(processed_dir: Path, features_path: Path, output_path: Path, logger):
    """
    Merge all search engine results with raw features.

    Strategy:
    1. FragPipe + Sage: Direct join on (raw_file, precursor_id)
    2. DIA-NN: First match by (sequence, charge) to inherit precursor_id,
               then fallback to tolerance matching for unique sequences
    """
    # Load raw features
    features = load_raw_features(features_path, logger)

    all_merged = []

    # FragPipe (direct precursor_id match)
    fragpipe_psms = load_fragpipe_psms(processed_dir, logger)
    if len(fragpipe_psms) > 0:
        fragpipe_merged = merge_by_precursor_id(fragpipe_psms, features, logger)
        all_merged.append(fragpipe_merged)

    # Sage (direct precursor_id match via scannr)
    sage_dir = processed_dir / "sage"
    sage_psms = load_sage_psms(sage_dir, logger)
    if len(sage_psms) > 0:
        sage_merged = merge_by_precursor_id(sage_psms, features, logger)
        all_merged.append(sage_merged)

    # DIA-NN (sequence+charge match first, then tolerance fallback)
    diann_dir = processed_dir / "diann"
    diann_psms = load_diann_psms(diann_dir, logger)
    if len(diann_psms) > 0:
        # Combine existing merged for sequence lookup
        existing_merged = pd.concat(all_merged, ignore_index=True) if all_merged else pd.DataFrame()
        diann_merged = merge_diann_by_sequence(diann_psms, existing_merged, features, logger)
        if len(diann_merged) > 0:
            all_merged.append(diann_merged)

    if not all_merged:
        logger.error("No PSMs found from any search engine!")
        return

    # Combine all
    combined = pd.concat(all_merged, ignore_index=True)
    logger.info(f"\nTotal merged: {len(combined):,} PSMs")
    logger.info(f"By engine: {combined['engine'].value_counts().to_dict()}")

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(output_path, index=False)
    logger.info(f"Saved to {output_path}")


# Snakemake entry point
if __name__ == "__main__" and 'snakemake' in dir():
    logger = setup_logging(Path(snakemake.log[0]))

    processed_dir = Path(snakemake.input.peptides).parent
    features_path = Path(snakemake.input.raw_features)
    output_path = Path(snakemake.output.merged)

    merge_all_engines(processed_dir, features_path, output_path, logger)


# CLI entry point for testing
if __name__ == "__main__" and 'snakemake' not in dir():
    import argparse
    parser = argparse.ArgumentParser(description="Merge search engine results with raw features")
    parser.add_argument("--processed", required=True, help="Processed data directory")
    parser.add_argument("--features", required=True, help="Raw features parquet")
    parser.add_argument("--output", required=True, help="Output parquet")
    parser.add_argument("--log", default="merge.log", help="Log file")
    args = parser.parse_args()

    logger = setup_logging(Path(args.log))
    merge_all_engines(Path(args.processed), Path(args.features), Path(args.output), logger)
