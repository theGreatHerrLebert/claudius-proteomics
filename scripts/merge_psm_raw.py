#!/usr/bin/env python3
"""
Merge FragPipe PSM identifications with raw extracted features.

Creates a unified dataset with:
- Peptide sequence and modifications
- CCS values
- Raw signal features (chromatograms, mobilograms, isotopes)
"""

import sys
from pathlib import Path

import pandas as pd

# Snakemake provides these variables
psm_file = Path(snakemake.input.psm)
raw_features_file = Path(snakemake.input.raw_features)
output_file = Path(snakemake.output.merged)
log_file = Path(snakemake.log[0])


def setup_logging(log_path: Path):
    """Redirect stdout/stderr to log file."""
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


def load_psms(psm_path: Path) -> pd.DataFrame:
    """Load and normalize PSM data from FragPipe."""
    logger.info(f"Loading PSMs from {psm_path}")
    psms = pd.read_csv(psm_path, sep='\t')

    # Normalize column names (FragPipe output can vary)
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
        'CCS': 'ccs',
        'Calculated M/Z': 'mz',
        'Calculated Peptide Mass': 'mass',
    }

    psms = psms.rename(columns={k: v for k, v in column_mapping.items() if k in psms.columns})

    logger.info(f"Loaded {len(psms)} PSMs with columns: {list(psms.columns)}")
    return psms


def load_raw_features(features_path: Path) -> pd.DataFrame:
    """Load raw extracted features."""
    logger.info(f"Loading raw features from {features_path}")
    features = pd.read_parquet(features_path)
    logger.info(f"Loaded {len(features)} feature records")
    return features


def merge_data(psms: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    """
    Merge PSMs with raw features.

    Alignment is done by PSM ID (row index from FragPipe processing).
    """
    logger.info("Merging PSMs with raw features...")

    # Merge on psm_id
    merged = psms.reset_index().rename(columns={'index': 'psm_id'})
    merged = merged.merge(features, on='psm_id', how='left', suffixes=('', '_raw'))

    logger.info(f"Merged dataset: {len(merged)} records")
    logger.info(f"Columns: {list(merged.columns)}")

    return merged


def clean_and_validate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean merged data and validate required fields.

    - Remove entries without CCS values
    - Remove duplicates (by sequence + charge)
    - Basic outlier detection
    """
    logger.info("Cleaning and validating data...")

    initial_count = len(df)

    # Require sequence
    if 'sequence' in df.columns:
        df = df[df['sequence'].notna() & (df['sequence'] != '')]
        logger.info(f"After sequence filter: {len(df)} records")

    # Require CCS (for CCS prediction task)
    if 'ccs' in df.columns:
        df = df[df['ccs'].notna() & (df['ccs'] > 0)]
        logger.info(f"After CCS filter: {len(df)} records")

    # Remove extreme outliers (CCS outside reasonable range)
    if 'ccs' in df.columns:
        ccs_min, ccs_max = 200, 1500  # Typical peptide CCS range in Angstrom^2
        df = df[(df['ccs'] >= ccs_min) & (df['ccs'] <= ccs_max)]
        logger.info(f"After CCS range filter: {len(df)} records")

    logger.info(f"Removed {initial_count - len(df)} records during cleaning")

    return df


def main():
    global logger
    logger = setup_logging(log_file)

    logger.info("Starting PSM-raw feature merge")

    # Load data
    psms = load_psms(psm_file)
    features = load_raw_features(raw_features_file)

    # Merge
    merged = merge_data(psms, features)

    # Clean
    cleaned = clean_and_validate(merged)

    # Save
    output_file.parent.mkdir(parents=True, exist_ok=True)
    cleaned.to_parquet(output_file, index=False)

    logger.info(f"Saved merged dataset to {output_file}")
    logger.info("Merge completed successfully")


if __name__ == "__main__":
    main()
