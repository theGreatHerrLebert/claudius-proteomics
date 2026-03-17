#!/usr/bin/env python3
"""
Rebuild precursor_index and precursor_store for /scratch datasets.

Re-runs step 3 (stratify) and step 5 (merge) using existing engine outputs
and extracted per-file index.parquet files. Useful when merging parameters
change (e.g., RT tolerance).

Usage:
    python scripts/rebuild_index.py /scratch/claudius-proteomics/PXD019086
    python scripts/rebuild_index.py --all   # rebuild all datasets in /scratch
"""

import argparse
import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from runner.steps.step3_stratify import (
    _build_precursor_index_anchored,
    _compute_overlap_stats_from_parsers,
)
from runner.steps.step5_merge import _merge_single_group

SCRATCH_BASE = Path("/scratch/claudius-proteomics")


def reconstruct_raw_features(dataset_dir: Path) -> pd.DataFrame:
    """Concatenate per-file index.parquet into raw_features DataFrame."""
    extracted_dir = dataset_dir / "extracted"
    if not extracted_dir.exists():
        print(f"  WARNING: No extracted/ directory in {dataset_dir}")
        return pd.DataFrame()

    index_files = sorted(extracted_dir.glob("*.d/index.parquet"))
    if not index_files:
        print(f"  WARNING: No index.parquet files found in {extracted_dir}")
        return pd.DataFrame()

    print(f"  Loading {len(index_files)} per-file index.parquet files...")
    frames = []
    for idx_file in index_files:
        raw_file = idx_file.parent.name.replace(".d", "")
        df = pd.read_parquet(idx_file)
        if "raw_file" not in df.columns:
            df["raw_file"] = raw_file
        frames.append(df)

    raw_features = pd.concat(frames, ignore_index=True)
    print(f"  Reconstructed raw_features: {len(raw_features)} precursors from {len(index_files)} files")
    return raw_features


def rebuild_dataset(dataset_dir: Path):
    """Rebuild precursor_index and precursor_store for a single dataset."""
    print(f"\n{'='*60}")
    print(f"Rebuilding: {dataset_dir.name}")
    print(f"{'='*60}")

    # Read existing manifest for metadata
    manifest_path = dataset_dir / "manifest.json"
    if manifest_path.exists():
        with open(manifest_path) as f:
            manifest = json.load(f)
        accession = manifest.get("accession", dataset_dir.name)
    else:
        accession = dataset_dir.name

    # Back up existing index
    old_index = dataset_dir / "precursor_index.parquet"
    if old_index.exists():
        backup = dataset_dir / "precursor_index.parquet.bak"
        shutil.copy2(old_index, backup)
        print(f"  Backed up old index to {backup.name}")

    # Step 1: Reconstruct raw_features from per-file index.parquet
    raw_features = reconstruct_raw_features(dataset_dir)

    # Keep only anchor columns
    anchor_cols = ["precursor_id", "raw_file", "mz", "charge", "rt_seconds", "mobility"]
    available = [c for c in anchor_cols if c in raw_features.columns]
    raw_precursors = raw_features[available].copy() if not raw_features.empty else pd.DataFrame()

    # Step 2: Load engine results via parsers adapted to /scratch layout
    print("  Loading engine results...")
    fp_df = _load_engine_from_scratch(dataset_dir, "fragpipe")
    dn_df = _load_engine_from_scratch(dataset_dir, "diann")
    sg_df = _load_engine_from_scratch(dataset_dir, "sage")

    print(f"    FragPipe: {len(fp_df) if fp_df is not None and not fp_df.empty else 0} PSMs")
    print(f"    DIA-NN: {len(dn_df) if dn_df is not None and not dn_df.empty else 0} precursors")
    print(f"    Sage: {len(sg_df) if sg_df is not None and not sg_df.empty else 0} PSMs")

    # Step 3: Build new precursor index with current RT tolerance
    from scripts.precursor_merging.config import MatchConfig
    print(f"  RT tolerance: {MatchConfig().rt_tol_sec}s")

    overlap_stats = _compute_overlap_stats_from_parsers(fp_df, dn_df, sg_df)
    print(f"  Overlap: all_three={overlap_stats['n_all_three']}, "
          f"at_least_two={overlap_stats['n_at_least_two']}, "
          f"union={overlap_stats['n_union']}")

    precursor_index = _build_precursor_index_anchored(
        raw_precursors=raw_precursors,
        fp_df=fp_df,
        dn_df=dn_df,
        sg_df=sg_df,
        config={},
    )

    # Save new index
    index_path = dataset_dir / "precursor_index.parquet"
    precursor_index.to_parquet(index_path, index=False)
    print(f"  Wrote new precursor_index: {len(precursor_index)} rows → {index_path}")

    # Step 4: Rebuild precursor_store (step 5 merge)
    old_store = dataset_dir / "precursor_store.parquet"
    if old_store.exists():
        backup = dataset_dir / "precursor_store.parquet.bak"
        shutil.copy2(old_store, backup)
        print(f"  Backed up old store to {backup.name}")

    # Save raw_features temporarily for step 5
    raw_features_path = dataset_dir / "extracted" / "raw_features.parquet"
    wrote_temp = False
    if not raw_features_path.exists() and not raw_features.empty:
        raw_features.to_parquet(raw_features_path, index=False)
        wrote_temp = True

    result = _merge_single_group(
        accession=accession,
        group_id=None,
        precursor_index_path=index_path,
        raw_features_path=raw_features_path,
        merged_dir=dataset_dir,
    )

    # Clean up temp raw_features if we created it
    if wrote_temp and raw_features_path.exists():
        raw_features_path.unlink()
        print(f"  Cleaned up temporary raw_features.parquet")

    # Update manifest
    if manifest_path.exists():
        with open(manifest_path) as f:
            manifest = json.load(f)
        manifest["generated_at"] = datetime.now().isoformat()
        manifest["n_total_precursors"] = result["n_total_precursors"]
        manifest["n_per_engine"] = result["n_per_engine"]
        manifest["n_unidentified"] = result["n_unidentified"]
        manifest["quality_summary"] = result["quality_summary"]
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2)
        print(f"  Updated manifest.json")

    print(f"\n  DONE: {result['n_total_precursors']} precursors")
    return result


def _load_engine_from_scratch(dataset_dir: Path, engine_name: str):
    """Load engine results from /scratch flat layout using parsers.

    /scratch layout: dataset_dir/engines/{engine}/{files}
    Parser expects: base_dir/accession/{engine}/{files}
    So we call parser.parse(dataset_dir, "engines") to resolve correctly.
    """
    engine_dir = dataset_dir / "engines" / engine_name
    if not engine_dir.exists():
        return None

    # Check for canonical parquet first
    canonical = engine_dir / f"{engine_name}_canonical.parquet"
    if canonical.exists():
        return pd.read_parquet(canonical)

    # Use parser with path mapping: parser looks in base_dir/accession/
    # so parse(dataset_dir, "engines") → dataset_dir/engines/{engine}/
    from scripts.engine_parsers.fragpipe_parser import FragPipeParser
    from scripts.engine_parsers.diann_parser import DiannParser
    from scripts.engine_parsers.sage_parser import SageParser

    parsers = {"fragpipe": FragPipeParser, "diann": DiannParser, "sage": SageParser}
    parser_cls = parsers.get(engine_name)
    if parser_cls is None:
        return None

    try:
        return parser_cls().parse(dataset_dir, "engines")
    except Exception as e:
        print(f"    WARNING: Could not parse {engine_name}: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Rebuild precursor indices with updated merging parameters")
    parser.add_argument("dataset_dir", nargs="?", help="Path to dataset directory")
    parser.add_argument("--all", action="store_true", help="Rebuild all datasets in /scratch/claudius-proteomics")
    args = parser.parse_args()

    if args.all:
        datasets = sorted(SCRATCH_BASE.glob("PXD*"))
        # Also include SIM_SMOKE
        sim = SCRATCH_BASE / "SIM_SMOKE_human_trypsin"
        if sim.exists():
            datasets.append(sim)
    elif args.dataset_dir:
        datasets = [Path(args.dataset_dir)]
    else:
        parser.error("Provide a dataset path or --all")

    print(f"Will rebuild {len(datasets)} dataset(s)")
    for d in datasets:
        if (d / "engines").exists():
            rebuild_dataset(d)
        else:
            print(f"\nSkipping {d.name} (no engines/ directory)")


if __name__ == "__main__":
    main()
