#!/usr/bin/env python3
"""
Combine multiple merged datasets into a single training set.

Handles:
- Deduplication across studies
- Alignment of CCS values (if from different instruments/calibrations)
- Train/test split preparation
"""

import sys
from pathlib import Path
from typing import List

import pandas as pd

# Snakemake provides these variables
input_files = [Path(f) for f in snakemake.input.datasets]
output_file = Path(snakemake.output.combined)
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


def load_datasets(paths: List[Path]) -> pd.DataFrame:
    """Load and concatenate all dataset files."""
    dfs = []

    for path in paths:
        logger.info(f"Loading {path}")
        df = pd.read_parquet(path)
        # Add source column
        df['source'] = path.parent.name  # Use accession as source
        dfs.append(df)
        logger.info(f"  -> {len(df)} records")

    combined = pd.concat(dfs, ignore_index=True)
    logger.info(f"Combined: {len(combined)} total records")

    return combined


def deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Deduplicate peptides across studies.

    For peptides seen in multiple studies, we can:
    1. Keep all observations (more training data)
    2. Average CCS values (reduce noise)
    3. Keep highest-confidence observation

    For now, we keep all observations but flag duplicates.
    """
    logger.info("Checking for duplicates...")

    # Define unique peptide as sequence + charge + modifications
    if 'modified_sequence' in df.columns:
        df['peptide_key'] = df['modified_sequence'].astype(str) + '_' + df['charge'].astype(str)
    else:
        df['peptide_key'] = df['sequence'].astype(str) + '_' + df['charge'].astype(str)

    # Count occurrences
    peptide_counts = df['peptide_key'].value_counts()
    duplicates = peptide_counts[peptide_counts > 1]

    logger.info(f"Unique peptides: {len(peptide_counts)}")
    logger.info(f"Peptides with multiple observations: {len(duplicates)}")

    # Flag duplicates but keep all
    df['is_duplicate'] = df['peptide_key'].isin(duplicates.index)

    return df


def add_split_column(df: pd.DataFrame, test_ratio: float = 0.1) -> pd.DataFrame:
    """
    Add train/test split column based on peptide (not PSM).

    Ensures same peptide doesn't appear in both train and test.
    """
    logger.info("Adding train/test split...")

    unique_peptides = df['peptide_key'].unique()
    n_test = int(len(unique_peptides) * test_ratio)

    # Random split
    import numpy as np
    np.random.seed(42)  # Reproducibility
    test_peptides = set(np.random.choice(unique_peptides, n_test, replace=False))

    df['split'] = df['peptide_key'].apply(lambda x: 'test' if x in test_peptides else 'train')

    train_count = (df['split'] == 'train').sum()
    test_count = (df['split'] == 'test').sum()
    logger.info(f"Train: {train_count} records, Test: {test_count} records")

    return df


def main():
    global logger
    logger = setup_logging(log_file)

    logger.info("Starting dataset combination")
    logger.info(f"Input files: {input_files}")

    # Load all datasets
    combined = load_datasets(input_files)

    # Deduplicate
    deduped = deduplicate(combined)

    # Add split column
    final = add_split_column(deduped)

    # Summary statistics
    logger.info("\n=== Dataset Summary ===")
    logger.info(f"Total records: {len(final)}")
    logger.info(f"Unique peptides: {final['peptide_key'].nunique()}")
    if 'ccs' in final.columns:
        logger.info(f"CCS range: {final['ccs'].min():.1f} - {final['ccs'].max():.1f}")
        logger.info(f"CCS mean: {final['ccs'].mean():.1f}")
    logger.info(f"Sources: {final['source'].unique().tolist()}")

    # Save
    output_file.parent.mkdir(parents=True, exist_ok=True)
    final.to_parquet(output_file, index=False)

    logger.info(f"Saved combined dataset to {output_file}")


if __name__ == "__main__":
    main()
