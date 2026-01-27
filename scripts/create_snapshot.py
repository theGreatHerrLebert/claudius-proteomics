#!/usr/bin/env python3
"""
Create a versioned snapshot from the peptide database for model training.

Snapshots are immutable exports that combine data from multiple datasets
and include train/val/test splits for reproducible experiments.
"""

import json
import sys
from datetime import datetime
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

# Snakemake provides these variables
metadata_file = Path(snakemake.input.metadata)
peptide_files = [Path(f) for f in snakemake.input.peptides]
output_training_set = Path(snakemake.output.training_set)
output_manifest = Path(snakemake.output.manifest)
log_file = Path(snakemake.log[0])

# Parameters
test_split = snakemake.params.test_split
min_ccs = snakemake.params.min_ccs
max_ccs = snakemake.params.max_ccs


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


def load_all_datasets(peptide_files: List[Path]) -> pd.DataFrame:
    """Load and concatenate all dataset files."""
    logger.info(f"Loading {len(peptide_files)} datasets...")

    dfs = []
    for path in peptide_files:
        logger.info(f"  Loading {path.parent.name}...")
        df = pd.read_parquet(path)
        dfs.append(df)
        logger.info(f"    {len(df):,} records")

    combined = pd.concat(dfs, ignore_index=True)
    logger.info(f"Combined: {len(combined):,} total records")

    return combined


def filter_for_training(df: pd.DataFrame, min_ccs: float, max_ccs: float) -> pd.DataFrame:
    """Apply filters for training data quality."""
    logger.info("Applying training data filters...")

    initial_count = len(df)

    # Require CCS for CCS prediction training
    df = df[df['ccs'].notna()]
    logger.info(f"  After CCS required: {len(df):,}")

    # Filter CCS range
    df = df[(df['ccs'] >= min_ccs) & (df['ccs'] <= max_ccs)]
    logger.info(f"  After CCS range [{min_ccs}, {max_ccs}]: {len(df):,}")

    # Require sequence
    df = df[df['sequence'].notna() & (df['sequence'] != '')]
    logger.info(f"  After sequence required: {len(df):,}")

    # Require charge
    df = df[df['charge'].notna() & (df['charge'] > 0)]
    logger.info(f"  After charge required: {len(df):,}")

    logger.info(f"Filtered: {initial_count - len(df):,} records removed")

    return df


def deduplicate_peptides(df: pd.DataFrame) -> pd.DataFrame:
    """
    Handle duplicate peptides across datasets.

    For peptides seen multiple times, we aggregate CCS values
    and keep the most common/median value.
    """
    logger.info("Deduplicating peptides...")

    # Create peptide key (sequence + charge)
    df['peptide_key'] = df['sequence'] + '_' + df['charge'].astype(str)

    initial_records = len(df)
    initial_peptides = df['peptide_key'].nunique()

    # Group by peptide and aggregate
    agg_funcs = {
        'sequence': 'first',
        'modified_sequence': 'first',
        'charge': 'first',
        'ccs': 'median',  # Use median CCS across observations
        'mz': 'first',
        'mass': 'first',
        'retention_time': 'median',
        'mobility': 'median',
        'accession': lambda x: ','.join(x.unique()),  # Track source datasets
    }

    # Only aggregate columns that exist
    agg_funcs = {k: v for k, v in agg_funcs.items() if k in df.columns}

    # Add count column
    df['n_observations'] = 1

    deduped = df.groupby('peptide_key').agg({
        **agg_funcs,
        'n_observations': 'sum'
    }).reset_index()

    logger.info(f"Deduplicated: {initial_records:,} records -> {len(deduped):,} unique peptides")
    logger.info(f"  Multi-observation peptides: {(deduped['n_observations'] > 1).sum():,}")

    return deduped


def create_splits(df: pd.DataFrame, test_ratio: float, val_ratio: float = 0.1) -> pd.DataFrame:
    """
    Create train/val/test splits by peptide.

    Ensures same peptide doesn't appear in multiple splits.
    """
    logger.info(f"Creating splits (test={test_ratio}, val={val_ratio})...")

    np.random.seed(42)  # Reproducibility

    peptides = df['peptide_key'].unique()
    n_peptides = len(peptides)

    n_test = int(n_peptides * test_ratio)
    n_val = int(n_peptides * val_ratio)
    n_train = n_peptides - n_test - n_val

    # Shuffle and split
    np.random.shuffle(peptides)

    test_peptides = set(peptides[:n_test])
    val_peptides = set(peptides[n_test:n_test + n_val])
    train_peptides = set(peptides[n_test + n_val:])

    # Assign splits
    def assign_split(key):
        if key in test_peptides:
            return 'test'
        elif key in val_peptides:
            return 'val'
        else:
            return 'train'

    df['split'] = df['peptide_key'].apply(assign_split)

    # Log split statistics
    for split in ['train', 'val', 'test']:
        split_df = df[df['split'] == split]
        logger.info(f"  {split}: {len(split_df):,} peptides")

    return df


def compute_snapshot_stats(df: pd.DataFrame) -> dict:
    """Compute statistics for snapshot manifest."""
    stats = {
        'n_psms': len(df),
        'n_unique_peptides': df['peptide_key'].nunique(),
        'source_datasets': df['accession'].str.split(',').explode().unique().tolist(),
        'split_counts': df['split'].value_counts().to_dict(),
    }

    if 'ccs' in df.columns:
        stats['ccs_min'] = float(df['ccs'].min())
        stats['ccs_max'] = float(df['ccs'].max())
        stats['ccs_mean'] = float(df['ccs'].mean())
        stats['ccs_std'] = float(df['ccs'].std())

    stats['charge_distribution'] = df['charge'].value_counts().to_dict()

    return stats


def main():
    global logger
    logger = setup_logging(log_file)

    version = output_training_set.parent.name

    logger.info("=" * 60)
    logger.info(f"Creating database snapshot: {version}")
    logger.info("=" * 60)

    # Load all datasets
    combined = load_all_datasets(peptide_files)

    # Filter for training
    filtered = filter_for_training(combined, min_ccs, max_ccs)

    # Deduplicate
    deduped = deduplicate_peptides(filtered)

    # Create splits
    with_splits = create_splits(deduped, test_split)

    # Compute statistics
    stats = compute_snapshot_stats(with_splits)

    # Create manifest
    manifest = {
        'version': version,
        'created': datetime.now().isoformat(),
        'split_ratio': {
            'test': test_split,
            'val': 0.1,
            'train': 1.0 - test_split - 0.1,
        },
        'filters': {
            'min_ccs': min_ccs,
            'max_ccs': max_ccs,
            'require_ccs': True,
        },
        **stats,
    }

    # Save
    output_training_set.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Saving training set to {output_training_set}")
    with_splits.to_parquet(output_training_set, index=False)

    logger.info(f"Saving manifest to {output_manifest}")
    with open(output_manifest, 'w') as f:
        json.dump(manifest, f, indent=2)

    logger.info("=" * 60)
    logger.info(f"Snapshot {version} created successfully")
    logger.info(f"  Total peptides: {stats['n_unique_peptides']:,}")
    logger.info(f"  Source datasets: {stats['source_datasets']}")
    logger.info(f"  CCS range: {stats['ccs_min']:.1f} - {stats['ccs_max']:.1f} Å²")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
