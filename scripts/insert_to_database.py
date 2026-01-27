#!/usr/bin/env python3
"""
Insert processed data into the peptide property database.

Merges FragPipe PSMs with raw extracted features and stores
in the database structure with full provenance tracking.
"""

import json
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import numpy as np

# Snakemake provides these variables
psm_file = Path(snakemake.input.psm)
raw_features_file = Path(snakemake.input.raw_features)
output_peptides = Path(snakemake.output.peptides)
output_manifest = Path(snakemake.output.manifest)
log_file = Path(snakemake.log[0])


def setup_logging(log_path: Path):
    """Set up logging."""
    import logging

    log_path.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    fh = logging.FileHandler(log_path)
    ch = logging.StreamHandler(sys.stdout)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger


def load_psms(psm_path: Path) -> pd.DataFrame:
    """Load and normalize PSM data."""
    logger.info(f"Loading PSMs from {psm_path}")
    psms = pd.read_csv(psm_path, sep='\t')

    # Normalize column names
    column_mapping = {
        'Peptide': 'sequence',
        'peptide': 'sequence',
        'Modified Sequence': 'modified_sequence',
        'Charge': 'charge',
        'charge': 'charge',
        'Retention': 'retention_time',
        'RT': 'retention_time',
        'Ion Mobility': 'mobility',
        'IonMobility': 'mobility',
        '1/K0': 'mobility',
        'CCS': 'ccs',
        'Calculated M/Z': 'mz',
        'm/z': 'mz',
        'Calculated Peptide Mass': 'mass',
        'Spectrum File': 'raw_file',
    }

    psms = psms.rename(columns={k: v for k, v in column_mapping.items() if k in psms.columns})
    logger.info(f"Loaded {len(psms)} PSMs")

    return psms


def load_raw_features(features_path: Path) -> pd.DataFrame:
    """Load raw extracted features."""
    logger.info(f"Loading raw features from {features_path}")
    features = pd.read_parquet(features_path)
    logger.info(f"Loaded {len(features)} feature records")
    return features


def merge_data(psms: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    """Merge PSMs with raw features."""
    logger.info("Merging PSMs with raw features...")

    # Add index to PSMs for merging
    psms = psms.reset_index().rename(columns={'index': 'psm_id'})

    # Merge on psm_id
    merged = psms.merge(features, on='psm_id', how='left', suffixes=('', '_raw'))

    # Handle duplicate columns from merge
    for col in merged.columns:
        if col.endswith('_raw'):
            base_col = col[:-4]
            if base_col in merged.columns:
                # Fill missing values from raw
                merged[base_col] = merged[base_col].fillna(merged[col])
            merged = merged.drop(columns=[col])

    logger.info(f"Merged dataset: {len(merged)} records")
    return merged


def clean_and_validate(df: pd.DataFrame, accession: str) -> pd.DataFrame:
    """Clean data and add metadata."""
    logger.info("Cleaning and validating data...")

    initial_count = len(df)

    # Add accession
    df['accession'] = accession

    # Ensure required columns exist
    required_cols = ['sequence', 'charge']
    for col in required_cols:
        if col not in df.columns:
            logger.error(f"Missing required column: {col}")
            raise ValueError(f"Missing required column: {col}")

    # Filter valid sequences
    df = df[df['sequence'].notna() & (df['sequence'] != '')]
    df = df[df['sequence'].str.match(r'^[A-Z]+$', na=False)]  # Only standard AAs

    # Filter valid charges
    df = df[df['charge'].notna() & (df['charge'] > 0) & (df['charge'] <= 10)]

    # Filter CCS if present
    if 'ccs' in df.columns:
        # Keep entries without CCS (they may have other useful features)
        # But filter out clearly invalid CCS values
        valid_ccs = df['ccs'].isna() | ((df['ccs'] >= 100) & (df['ccs'] <= 2000))
        df = df[valid_ccs]

    logger.info(f"After cleaning: {len(df)} records (removed {initial_count - len(df)})")

    return df


def compute_statistics(df: pd.DataFrame) -> dict:
    """Compute dataset statistics for manifest."""
    stats = {
        'n_psms': len(df),
        'n_unique_peptides': df['sequence'].nunique(),
        'n_unique_modified': df['modified_sequence'].nunique() if 'modified_sequence' in df.columns else 0,
        'charge_distribution': df['charge'].value_counts().to_dict(),
    }

    if 'ccs' in df.columns:
        ccs_valid = df['ccs'].dropna()
        if len(ccs_valid) > 0:
            stats['ccs_stats'] = {
                'n_with_ccs': len(ccs_valid),
                'mean': float(ccs_valid.mean()),
                'std': float(ccs_valid.std()),
                'min': float(ccs_valid.min()),
                'max': float(ccs_valid.max()),
            }

    if 'retention_time' in df.columns:
        rt_valid = df['retention_time'].dropna()
        if len(rt_valid) > 0:
            stats['rt_stats'] = {
                'mean': float(rt_valid.mean()),
                'min': float(rt_valid.min()),
                'max': float(rt_valid.max()),
            }

    return stats


def main():
    global logger
    logger = setup_logging(log_file)

    logger.info("=" * 60)
    logger.info("Inserting data into peptide database")
    logger.info("=" * 60)

    # Get accession from path
    accession = output_peptides.parent.name
    logger.info(f"Accession: {accession}")

    # Load data
    psms = load_psms(psm_file)
    features = load_raw_features(raw_features_file)

    # Merge
    merged = merge_data(psms, features)

    # Clean and validate
    cleaned = clean_and_validate(merged, accession)

    # Compute statistics
    stats = compute_statistics(cleaned)

    # Create manifest
    manifest = {
        'accession': accession,
        'created': datetime.now().isoformat(),
        'source_files': {
            'psm': str(psm_file),
            'raw_features': str(raw_features_file),
        },
        'processing': {
            'pipeline_version': '1.0',
        },
        **stats,
    }

    # Save
    output_peptides.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Saving peptides to {output_peptides}")
    cleaned.to_parquet(output_peptides, index=False)

    logger.info(f"Saving manifest to {output_manifest}")
    with open(output_manifest, 'w') as f:
        json.dump(manifest, f, indent=2)

    logger.info("=" * 60)
    logger.info(f"Database insertion complete for {accession}")
    logger.info(f"  PSMs: {stats['n_psms']:,}")
    logger.info(f"  Unique peptides: {stats['n_unique_peptides']:,}")
    if 'ccs_stats' in stats:
        logger.info(f"  With CCS: {stats['ccs_stats']['n_with_ccs']:,}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
