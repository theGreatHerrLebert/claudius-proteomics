#!/usr/bin/env python3
"""
Step 4: Raw Data Extraction with Quality Metrics

Extracts 4D raw signals with Gaussian fits and quality metrics.

Input: data/raw/{accession}/*.d

Outputs (per-group when sample_groups.yaml exists):
- data/extracted/{accession}/{group_id}/{raw_file}/
    - index.parquet
    - blobs.bin
- data/extracted/{accession}/{group_id}/raw_features.parquet (merged)
- step4_summary.json

Legacy (no sample_groups.yaml):
- data/extracted/{accession}/{raw_file}/
- data/extracted/{accession}/raw_features.parquet
"""

import sys
from pathlib import Path
from typing import Dict, Any, List, Optional

import numpy as np
import pandas as pd

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from lib.quality_metrics import (
    calculate_moments,
    fit_gaussian,
    compute_isotope_cosine_similarity,
)
from runner.summary import StepSummary, write_step_summary


def run_step4_extract(
    accession: str,
    config: Dict[str, Any],
    output_base_dir: Path,
    raw_dir: Optional[Path] = None,
    num_threads: int = 16,
    batch_size: int = 10000,
    max_files: int = 0,
) -> StepSummary:
    """
    Execute Step 4: Extract raw 4D signals with quality metrics.

    If sample_groups.yaml exists, extracts per-group into
    data/extracted/{accession}/{group_id}/. Otherwise falls back
    to legacy single-group output at data/extracted/{accession}/.

    Args:
        accession: PRIDE accession
        config: Pipeline configuration dict
        output_base_dir: Base directory for outputs
        raw_dir: Path to raw data directory (default: data/raw/{accession})
        num_threads: Number of threads for extraction
        batch_size: Batch size for memory-efficient extraction
        max_files: Maximum files to process (0 = all)

    Returns:
        StepSummary with results
    """
    summary = StepSummary(
        step_name="step4",
        accession=accession,
    )

    if raw_dir is None:
        raw_dir = output_base_dir / "raw" / accession

    extracted_dir = output_base_dir / "extracted" / accession

    try:
        # Check for sample group manifest
        metadata_dir = output_base_dir / "metadata" / accession
        sg_path = metadata_dir / "sample_groups.yaml"

        if sg_path.exists():
            from scripts.sample_group_resolver import SampleGroupManifest

            manifest = SampleGroupManifest.from_yaml(sg_path)
            print(f"  Loaded {len(manifest.groups)} sample groups from {sg_path}")

            group_results = {}
            total_precursors = 0
            total_files = 0
            total_blob_size = 0
            all_quality_metrics = {"rt_r2": [], "im_r2": [], "isotope_cosim": []}
            all_outputs = []

            for group in manifest.groups:
                print(f"\n  === Sample group: {group.group_id} ===")
                print(f"      Organism: {group.organism_name} ({group.organism_key})")
                print(f"      Enzyme:   {group.enzyme}")
                print(f"      Runs:     {group.n_runs}")

                if group.n_runs == 0:
                    print(f"      Skipping (no runs)")
                    group_results[group.group_id] = {"status": "skipped", "reason": "no runs"}
                    continue

                # Filter .d files to this group's runs
                group_d_files = [raw_dir / run for run in group.runs]
                group_d_files = [f for f in group_d_files if f.exists()]

                if max_files > 0:
                    group_d_files = group_d_files[:max_files]

                if not group_d_files:
                    print(f"      Warning: No .d files found for this group")
                    group_results[group.group_id] = {"status": "skipped", "reason": "no .d files found"}
                    continue

                # Per-group output directory
                group_extracted_dir = extracted_dir / group.group_id

                print(f"      Extracting from {len(group_d_files)} .d files")

                group_features, group_stats = _extract_file_list(
                    d_files=group_d_files,
                    extracted_dir=group_extracted_dir,
                    num_threads=num_threads,
                    batch_size=batch_size,
                )

                total_precursors += group_stats["n_precursors"]
                total_files += group_stats["n_files"]
                total_blob_size += group_stats["blob_size_bytes"]

                # Accumulate quality metrics
                for key in all_quality_metrics:
                    all_quality_metrics[key].extend(group_stats["quality_metrics"].get(key, []))

                all_outputs.append(str(group_extracted_dir))

                group_results[group.group_id] = {
                    "status": "success",
                    "n_precursors": group_stats["n_precursors"],
                    "n_files": group_stats["n_files"],
                    "blob_size_bytes": group_stats["blob_size_bytes"],
                }

            # Compute quality statistics
            quality_stats = _compute_quality_stats(all_quality_metrics)

            summary.data = {
                "mode": "per_group",
                "n_precursors_extracted": total_precursors,
                "n_files_processed": total_files,
                "blob_size_gb": round(total_blob_size / (1024 ** 3), 3),
                "quality_stats": quality_stats,
                "group_results": group_results,
            }
            summary.outputs = all_outputs
            summary.complete(success=True)

        else:
            # Legacy single-group behavior
            print(f"  No sample_groups.yaml found — using single-group mode")

            d_files = sorted(raw_dir.glob("*.d"))
            if max_files > 0:
                d_files = d_files[:max_files]

            if not d_files:
                raise FileNotFoundError(f"No .d files found in {raw_dir}")

            print(f"  Extracting from {len(d_files)} .d files")

            group_features, group_stats = _extract_file_list(
                d_files=d_files,
                extracted_dir=extracted_dir,
                num_threads=num_threads,
                batch_size=batch_size,
            )

            quality_stats = _compute_quality_stats(group_stats["quality_metrics"])

            summary.data = {
                "n_precursors_extracted": group_stats["n_precursors"],
                "n_files_processed": group_stats["n_files"],
                "blob_size_gb": round(group_stats["blob_size_bytes"] / (1024 ** 3), 3),
                "quality_stats": quality_stats,
            }
            summary.outputs = [str(extracted_dir)]
            summary.complete(success=True)

    except Exception as e:
        summary.complete(success=False, error_message=str(e))
        raise

    # Write summary file
    write_step_summary(summary, output_base_dir / "processed" / accession)

    return summary


def _extract_file_list(
    d_files: List[Path],
    extracted_dir: Path,
    num_threads: int = 16,
    batch_size: int = 10000,
) -> tuple:
    """
    Extract precursors from a list of .d files and merge features.

    Returns:
        (merged_features DataFrame, stats dict)
    """
    all_features = []
    total_precursors = 0
    total_blob_size = 0
    quality_metrics = {"rt_r2": [], "im_r2": [], "isotope_cosim": []}

    for i, d_file in enumerate(d_files):
        print(f"  [{i + 1}/{len(d_files)}] Processing {d_file.name}...")

        file_features, file_stats = _extract_single_file(
            d_file,
            extracted_dir / d_file.name,
            num_threads=num_threads,
            batch_size=batch_size,
        )

        if file_features is not None and not file_features.empty:
            all_features.append(file_features)
            total_precursors += len(file_features)
            total_blob_size += file_stats.get("blob_size_bytes", 0)

            for col in ["ms1_rt_r2", "ms1_im_r2", "isotope_cosim"]:
                if col in file_features.columns:
                    quality_metrics[col.replace("ms1_", "")].extend(
                        file_features[col].dropna().tolist()
                    )

    # Merge all features
    if all_features:
        merged_features = pd.concat(all_features, ignore_index=True)
        merged_path = extracted_dir / "raw_features.parquet"
        extracted_dir.mkdir(parents=True, exist_ok=True)
        merged_features.to_parquet(merged_path, index=False)
        print(f"  Merged features saved: {merged_path}")
    else:
        merged_features = pd.DataFrame()

    stats = {
        "n_precursors": total_precursors,
        "n_files": len(d_files),
        "blob_size_bytes": total_blob_size,
        "quality_metrics": quality_metrics,
    }

    return merged_features, stats


def _extract_single_file(
    d_file: Path,
    output_dir: Path,
    num_threads: int = 16,
    batch_size: int = 10000,
) -> tuple:
    """Extract precursors from a single .d file."""
    try:
        from scripts.extract_precursors import (
            extract_precursors_batched,
            setup_logging,
        )
        from scripts.extract_calibration import ensure_calibration, get_calibration_path
        from imspy_core.timstof import TimsDatasetDDA

        output_dir.mkdir(parents=True, exist_ok=True)

        # Load dataset
        dataset = TimsDatasetDDA(str(d_file), in_memory=False, use_bruker_sdk=False)

        # Get calibration
        cal_path = get_calibration_path(str(d_file))
        if cal_path.exists():
            calibration = np.load(cal_path)
        else:
            calibration = ensure_calibration(str(d_file), verbose=False)

        # Extract with batching
        stats = extract_precursors_batched(
            dataset=dataset,
            raw_file_name=d_file.name,
            output_dir=output_dir,
            batch_size=batch_size,
            num_threads=num_threads,
            calibration=calibration,
        )

        # Load extracted features
        index_path = output_dir / "index.parquet"
        if index_path.exists():
            features = pd.read_parquet(index_path)
            features["raw_file"] = d_file.name
            return features, stats
        else:
            return pd.DataFrame(), {"n_precursors": 0, "blob_size_bytes": 0}

    except ImportError as e:
        print(f"    Warning: imspy not available, using fallback extraction: {e}")
        return _extract_fallback(d_file, output_dir)

    except Exception as e:
        print(f"    Error extracting {d_file.name}: {e}")
        return pd.DataFrame(), {"n_precursors": 0, "blob_size_bytes": 0}


def _extract_fallback(d_file: Path, output_dir: Path) -> tuple:
    """Fallback extraction when imspy is not available."""
    # This would use a simpler extraction method
    # For now, return empty results
    return pd.DataFrame(), {"n_precursors": 0, "blob_size_bytes": 0}


def _compute_quality_stats(quality_metrics: Dict[str, List[float]]) -> Dict[str, Any]:
    """Compute summary statistics for quality metrics."""
    stats = {}

    for metric_name, values in quality_metrics.items():
        if values:
            arr = np.array(values)
            stats[metric_name] = {
                "mean": round(float(np.mean(arr)), 4),
                "median": round(float(np.median(arr)), 4),
                "q10": round(float(np.percentile(arr, 10)), 4),
                "q90": round(float(np.percentile(arr, 90)), 4),
                "n_high_quality": int(np.sum(arr >= 0.8)),
                "pct_high_quality": round(float(np.sum(arr >= 0.8) / len(arr) * 100), 1),
            }
        else:
            stats[metric_name] = {
                "mean": 0.0, "median": 0.0, "q10": 0.0, "q90": 0.0,
                "n_high_quality": 0, "pct_high_quality": 0.0,
            }

    return stats


if __name__ == "__main__":
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description="Step 4: Extract raw features")
    parser.add_argument("--accession", "-a", required=True, help="PRIDE accession")
    parser.add_argument("--config", "-c", default="config/config.yaml", help="Config file")
    parser.add_argument("--output-dir", "-o", default="data", help="Output base directory")
    parser.add_argument("--raw-dir", type=Path, help="Raw data directory")
    parser.add_argument("--threads", type=int, default=16, help="Number of threads")
    parser.add_argument("--batch-size", type=int, default=10000, help="Batch size")
    parser.add_argument("--max-files", type=int, default=0, help="Max files (0=all)")

    args = parser.parse_args()

    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Run step
    summary = run_step4_extract(
        accession=args.accession,
        config=config,
        output_base_dir=Path(args.output_dir),
        raw_dir=args.raw_dir,
        num_threads=args.threads,
        batch_size=args.batch_size,
        max_files=args.max_files,
    )

    print(f"\nStep 4 completed: {summary.status}")
    print(f"  Precursors extracted: {summary.data['n_precursors_extracted']}")
    print(f"  Blob size: {summary.data['blob_size_gb']:.3f} GB")
