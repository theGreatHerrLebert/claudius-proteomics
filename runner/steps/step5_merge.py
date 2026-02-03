#!/usr/bin/env python3
"""
Step 5: Merge Search Results with Extracted Raw Data

Joins engine IDs with raw features into final dashboard-ready dataset.

Input:
- data/processed/{accession}/precursor_index.parquet (Step 3)
- data/extracted/{accession}/raw_features.parquet (Step 4)

Outputs:
- data/merged/{accession}/precursor_store.parquet
- data/merged/{accession}/manifest.json
- step5_summary.json
"""

import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional

import pandas as pd
import numpy as np

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from runner.summary import StepSummary, write_step_summary, create_manifest


def run_step5_merge(
    accession: str,
    config: Dict[str, Any],
    output_base_dir: Path,
    precursor_index_path: Optional[Path] = None,
    raw_features_path: Optional[Path] = None,
) -> StepSummary:
    """
    Execute Step 5: Final merge of search + raw data.

    Args:
        accession: PRIDE accession
        config: Pipeline configuration dict
        output_base_dir: Base directory for outputs
        precursor_index_path: Path to precursor index (default: from Step 3)
        raw_features_path: Path to raw features (default: from Step 4)

    Returns:
        StepSummary with results
    """
    summary = StepSummary(
        step_name="step5",
        accession=accession,
    )

    processed_dir = output_base_dir / "processed" / accession
    extracted_dir = output_base_dir / "extracted" / accession
    merged_dir = output_base_dir / "merged" / accession

    try:
        # Resolve input paths
        if precursor_index_path is None:
            precursor_index_path = processed_dir / "precursor_index.parquet"
        if raw_features_path is None:
            raw_features_path = extracted_dir / "raw_features.parquet"

        # Load inputs
        print("  Loading precursor index from Step 3...")
        if not precursor_index_path.exists():
            raise FileNotFoundError(f"Precursor index not found: {precursor_index_path}")
        precursor_index = pd.read_parquet(precursor_index_path)
        print(f"    {len(precursor_index)} precursors")

        print("  Loading raw features from Step 4...")
        if raw_features_path.exists():
            raw_features = pd.read_parquet(raw_features_path)
            print(f"    {len(raw_features)} extracted precursors")
        else:
            print("    No raw features found - creating from precursor index only")
            raw_features = pd.DataFrame()

        # Merge datasets
        print("  Merging search results with raw features...")
        merged_df = _merge_datasets(precursor_index, raw_features)

        # Add consensus columns
        merged_df = _add_consensus_columns(merged_df)

        # Add quality summary columns
        merged_df = _add_quality_columns(merged_df)

        # Save merged dataset
        merged_dir.mkdir(parents=True, exist_ok=True)
        output_path = merged_dir / "precursor_store.parquet"
        merged_df.to_parquet(output_path, index=False)
        print(f"  Saved precursor store: {output_path}")

        # Compute statistics
        n_total = len(merged_df)
        n_per_engine = {
            "fragpipe": int((merged_df.get("fragpipe_peptide", pd.Series()).notna()).sum()),
            "diann": int((merged_df.get("diann_peptide", pd.Series()).notna()).sum()),
            "sage": int((merged_df.get("sage_peptide", pd.Series()).notna()).sum()),
        }
        n_unidentified = int((merged_df["n_engines"] == 0).sum())

        # Quality summary
        quality_summary = _compute_quality_summary(merged_df)

        # Create manifest
        manifest = {
            "accession": accession,
            "pipeline_version": "1.0",
            "generated_at": datetime.now().isoformat(),
            "n_total_precursors": n_total,
            "n_per_engine": n_per_engine,
            "n_unidentified": n_unidentified,
            "quality_summary": quality_summary,
            "output_file": str(output_path),
        }
        manifest_path = merged_dir / "manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2)
        print(f"  Saved manifest: {manifest_path}")

        # Update summary
        summary.data = {
            "n_total_precursors": n_total,
            "n_per_engine": n_per_engine,
            "n_unidentified": n_unidentified,
            "quality_summary": quality_summary,
        }
        summary.outputs = [str(output_path), str(manifest_path)]
        summary.complete(success=True)

    except Exception as e:
        summary.complete(success=False, error_message=str(e))
        raise

    # Write summary file
    write_step_summary(summary, processed_dir)

    return summary


def _merge_datasets(
    precursor_index: pd.DataFrame,
    raw_features: pd.DataFrame,
) -> pd.DataFrame:
    """Merge precursor index with raw features."""

    if raw_features.empty:
        # No raw features - use precursor index only
        merged = precursor_index.copy()
        # Add placeholder columns for raw data
        merged["raw_file"] = None
        merged["precursor_id"] = None
        return merged

    # Create join key from sequence + charge
    precursor_index = precursor_index.copy()
    raw_features = raw_features.copy()

    # Normalize join keys
    if "sequence_normalized" not in raw_features.columns:
        # Need to create normalized sequence from raw data
        # This would require re-matching - for now, join on precursor_id if available
        pass

    # Try joining on (raw_file, precursor_id) if available in both
    if "precursor_id" in raw_features.columns and "precursor_id" in precursor_index.columns:
        merged = raw_features.merge(
            precursor_index,
            on=["raw_file", "precursor_id"],
            how="outer",
            suffixes=("", "_idx"),
        )
    elif "sequence_normalized" in raw_features.columns:
        # Join on sequence + charge
        merged = raw_features.merge(
            precursor_index,
            on=["sequence_normalized", "charge"],
            how="outer",
            suffixes=("", "_idx"),
        )
    else:
        # Cannot merge - return combined
        merged = pd.concat([raw_features, precursor_index], ignore_index=True)

    # Clean up duplicate columns
    for col in merged.columns:
        if col.endswith("_idx"):
            base_col = col[:-4]
            if base_col in merged.columns:
                merged.drop(columns=[col], inplace=True)

    return merged


def _add_consensus_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add consensus peptide and confidence weight columns."""
    df = df.copy()

    # Ensure n_engines column exists
    if "n_engines" not in df.columns:
        def count_engines(row):
            n = 0
            if pd.notna(row.get("fragpipe_peptide")) and row.get("fragpipe_peptide"):
                n += 1
            if pd.notna(row.get("diann_peptide")) and row.get("diann_peptide"):
                n += 1
            if pd.notna(row.get("sage_peptide")) and row.get("sage_peptide"):
                n += 1
            return n

        df["n_engines"] = df.apply(count_engines, axis=1)

    # Get consensus peptide (prefer more engines, then FragPipe)
    def get_consensus(row):
        peptides = []
        if pd.notna(row.get("fragpipe_modified")) and row.get("fragpipe_modified"):
            peptides.append(("fragpipe", row["fragpipe_modified"]))
        if pd.notna(row.get("diann_modified")) and row.get("diann_modified"):
            peptides.append(("diann", row["diann_modified"]))
        if pd.notna(row.get("sage_modified")) and row.get("sage_modified"):
            peptides.append(("sage", row["sage_modified"]))

        if not peptides:
            return ""
        # Return first available (FragPipe preferred due to order)
        return peptides[0][1]

    df["consensus_peptide"] = df.apply(get_consensus, axis=1)

    # Confidence weight based on engine agreement
    def get_confidence(row):
        n = row.get("n_engines", 0)
        if n >= 3:
            return 1.0
        elif n == 2:
            return 0.9
        elif n == 1:
            return 0.7
        else:
            return 0.0

    df["confidence_weight"] = df.apply(get_confidence, axis=1)

    return df


def _add_quality_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add summary quality columns."""
    df = df.copy()

    # Quality flag based on RT R², IM R², isotope cosim
    def is_high_quality(row):
        rt_r2 = row.get("ms1_rt_r2", 0)
        im_r2 = row.get("ms1_im_r2", 0)
        iso_cosim = row.get("isotope_cosim", 0)

        return (
            (pd.notna(rt_r2) and rt_r2 >= 0.8) and
            (pd.notna(im_r2) and im_r2 >= 0.8) and
            (pd.notna(iso_cosim) and iso_cosim >= 0.9)
        )

    if "ms1_rt_r2" in df.columns:
        df["is_high_quality"] = df.apply(is_high_quality, axis=1)
    else:
        df["is_high_quality"] = False

    return df


def _compute_quality_summary(df: pd.DataFrame) -> Dict[str, Any]:
    """Compute quality summary statistics."""
    summary = {}

    # Engine distribution
    n_total = len(df)
    for n_eng in range(4):
        count = (df["n_engines"] == n_eng).sum()
        summary[f"n_{n_eng}_engines"] = int(count)
        summary[f"pct_{n_eng}_engines"] = round(count / max(n_total, 1) * 100, 1)

    # Quality metrics if available
    for col in ["ms1_rt_r2", "ms1_im_r2", "isotope_cosim"]:
        if col in df.columns:
            values = df[col].dropna()
            if len(values) > 0:
                summary[f"{col}_mean"] = round(float(values.mean()), 4)
                summary[f"{col}_median"] = round(float(values.median()), 4)

    # High quality count
    if "is_high_quality" in df.columns:
        n_high = df["is_high_quality"].sum()
        summary["n_high_quality"] = int(n_high)
        summary["pct_high_quality"] = round(n_high / max(n_total, 1) * 100, 1)

    return summary


if __name__ == "__main__":
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description="Step 5: Final merge")
    parser.add_argument("--accession", "-a", required=True, help="PRIDE accession")
    parser.add_argument("--config", "-c", default="config/config.yaml", help="Config file")
    parser.add_argument("--output-dir", "-o", default="data", help="Output base directory")
    parser.add_argument("--precursor-index", type=Path, help="Precursor index path")
    parser.add_argument("--raw-features", type=Path, help="Raw features path")

    args = parser.parse_args()

    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Run step
    summary = run_step5_merge(
        accession=args.accession,
        config=config,
        output_base_dir=Path(args.output_dir),
        precursor_index_path=args.precursor_index,
        raw_features_path=args.raw_features,
    )

    print(f"\nStep 5 completed: {summary.status}")
    print(f"  Total precursors: {summary.data['n_total_precursors']}")
    print(f"  Per engine: {summary.data['n_per_engine']}")
    print(f"  Unidentified: {summary.data['n_unidentified']}")
