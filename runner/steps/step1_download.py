#!/usr/bin/env python3
"""
Step 1: Download Raw Data

Downloads timsTOF .d files from PRIDE or symlinks local data.

Input: PRIDE accession (e.g., PXD019086)

Outputs:
- data/raw/{accession}/*.d/
- data/metadata/{accession}/metadata.yaml
- step1_summary.json
"""

import sys
from pathlib import Path
from typing import Dict, Any, Optional, List

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from runner.summary import StepSummary, write_step_summary
from scripts.sample_group_resolver import resolve_sample_groups


def run_step1_download(
    accession: str,
    config: Dict[str, Any],
    output_base_dir: Path,
    test_mode: bool = False,
    max_files: int = 0,
    local_data_path: Optional[Path] = None,
) -> StepSummary:
    """
    Execute Step 1: Download raw data from PRIDE.

    Args:
        accession: PRIDE accession
        config: Pipeline configuration dict
        output_base_dir: Base directory for outputs
        test_mode: If True, limit files for testing
        max_files: Maximum number of files (0 = unlimited)
        local_data_path: Optional path to local data (skip download)

    Returns:
        StepSummary with results
    """
    summary = StepSummary(
        step_name="step1",
        accession=accession,
    )

    raw_dir = output_base_dir / "raw" / accession
    metadata_dir = output_base_dir / "metadata" / accession

    try:
        # Check for local data path in config or argument
        local_path = local_data_path
        if local_path is None:
            local_datasets = config.get("local_data", {})
            if not local_datasets:
                local_datasets = config.get("local_datasets", {})
                if local_datasets:
                    print("  Warning: 'local_datasets' config key is deprecated, use 'local_data'")
            if accession in local_datasets:
                local_path = Path(local_datasets[accession])

        if local_path and local_path.exists():
            # Use local data - create symlinks
            raw_files = _link_local_data(local_path, raw_dir)
            metadata = _create_local_metadata(accession, local_path, metadata_dir)
        else:
            # Download from PRIDE
            raw_files = _download_from_pride(
                accession, raw_dir, metadata_dir, max_files, test_mode
            )
            metadata = _load_metadata(metadata_dir)

        # Apply max_files limit if in test mode
        if test_mode and max_files > 0:
            raw_files = raw_files[:max_files]

        # Resolve sample groups (organism + enzyme per file)
        print("  Resolving sample groups...")
        sg_manifest = resolve_sample_groups(accession, raw_dir, metadata_dir, config)
        sg_manifest.to_yaml(metadata_dir / "sample_groups.yaml")
        print(f"  Sample groups: {len(sg_manifest.groups)} groups, "
              f"multi-organism={sg_manifest.is_multi_organism}")
        if sg_manifest.is_multi_instrument:
            instruments = {g.instrument_model for g in sg_manifest.groups if g.instrument_model}
            print(f"  Multi-instrument: {', '.join(sorted(instruments))}")

        # Paper metadata extraction (optional, graceful)
        paper_result = None
        print("  Paper metadata extraction...")
        try:
            from scripts.paper_metadata import run_paper_extraction

            paper_result = run_paper_extraction(accession, metadata_dir, config)
            if paper_result and paper_result.get("status") == "success":
                n_fields = len(paper_result.get("fields", {}))
                print(f"  Paper extraction: {n_fields} fields extracted")
            elif paper_result:
                print(f"  Paper extraction: {paper_result.get('status', 'unknown')}")
            else:
                print("  Paper extraction: skipped (disabled in config)")
        except Exception as e:
            print(f"  Paper extraction: failed ({e}), continuing...")

        # Calculate total size
        total_size_gb = sum(
            sum(f.stat().st_size for f in Path(rf).rglob("*") if f.is_file())
            for rf in raw_files if Path(rf).exists()
        ) / (1024 ** 3)

        # Update summary
        summary.data = {
            "n_raw_files": len(raw_files),
            "total_size_gb": round(total_size_gb, 2),
            "raw_files": [str(f) for f in raw_files],
            "metadata": {
                "organism": metadata.get("organism", "unknown"),
                "instrument": metadata.get("instrument", "unknown"),
                "lab_id": metadata.get("lab_id", "unknown"),
                "completeness": metadata.get("completeness", 0),
            },
            "paper_extraction": {
                "status": paper_result.get("status") if paper_result else "skipped",
                "n_fields": len(paper_result.get("fields", {})) if paper_result else 0,
            },
            "sample_groups": {
                "n_groups": len(sg_manifest.groups),
                "is_multi_organism": sg_manifest.is_multi_organism,
                "is_multi_instrument": sg_manifest.is_multi_instrument,
                "groups": [g.group_id for g in sg_manifest.groups],
                "group_instruments": {
                    g.group_id: g.instrument_model
                    for g in sg_manifest.groups
                },
                "unassigned_runs": len(sg_manifest.unassigned_runs),
            },
        }
        summary.outputs = [str(raw_dir), str(metadata_dir)]
        summary.complete(success=True)

    except Exception as e:
        summary.complete(success=False, error_message=str(e))
        raise

    # Write summary file
    write_step_summary(summary, output_base_dir / "processed" / accession)

    return summary


def _link_local_data(local_path: Path, raw_dir: Path) -> List[Path]:
    """Create symlinks to local .d files."""
    raw_dir.mkdir(parents=True, exist_ok=True)

    d_folders = list(local_path.glob("*.d"))
    raw_files = []

    for d_folder in d_folders:
        link_path = raw_dir / d_folder.name
        if link_path.exists():
            # Already linked or exists
            if link_path.is_symlink():
                link_path.unlink()
            elif link_path.is_dir():
                # Real directory exists, skip
                raw_files.append(link_path)
                continue

        link_path.symlink_to(d_folder.absolute())
        raw_files.append(link_path)

    print(f"  Linked {len(raw_files)} .d folders from {local_path}")
    return raw_files


def _download_from_pride(
    accession: str,
    raw_dir: Path,
    metadata_dir: Path,
    max_files: int,
    test_mode: bool,
) -> List[Path]:
    """Download files from PRIDE Archive."""
    from scripts.download_pride import download_pride

    print(f"  Downloading from PRIDE: {accession}")

    # Apply max_files for test mode
    download_max = max_files if test_mode and max_files > 0 else 0

    success = download_pride(
        accession=accession,
        output_dir=raw_dir,
        metadata_dir=metadata_dir,
        max_files=download_max,
        extract=True,
    )

    if not success:
        raise RuntimeError(f"Failed to download {accession} from PRIDE")

    # Find downloaded .d folders
    raw_files = list(raw_dir.glob("*.d"))
    return raw_files


def _create_local_metadata(
    accession: str,
    local_path: Path,
    metadata_dir: Path,
) -> Dict[str, Any]:
    """Create metadata for local data."""
    import yaml

    metadata_dir.mkdir(parents=True, exist_ok=True)

    # Try to read existing metadata
    existing_metadata_file = local_path / "metadata.yaml"
    if existing_metadata_file.exists():
        with open(existing_metadata_file) as f:
            metadata = yaml.safe_load(f) or {}
    else:
        metadata = {}

    # Fill in defaults
    metadata.setdefault("accession", accession)
    metadata.setdefault("source", "local")
    metadata.setdefault("local_path", str(local_path))

    # Write metadata
    output_file = metadata_dir / "metadata.yaml"
    with open(output_file, "w") as f:
        yaml.dump(metadata, f, default_flow_style=False)

    return metadata


def _load_metadata(metadata_dir: Path) -> Dict[str, Any]:
    """Load metadata from YAML file."""
    import yaml

    metadata_file = metadata_dir / "pride_metadata.yaml"
    if not metadata_file.exists():
        metadata_file = metadata_dir / "metadata.yaml"

    if metadata_file.exists():
        with open(metadata_file) as f:
            data = yaml.safe_load(f) or {}

        # Extract relevant fields from PRIDE metadata format
        fields = data.get("fields", {})
        return {
            "organism": fields.get("organism", {}).get("value", "unknown"),
            "instrument": fields.get("instrument", {}).get("value", "unknown"),
            "lab_id": fields.get("lab_id", {}).get("value", "unknown"),
            "completeness": 1.0 if data.get("validation", {}).get("complete") else 0.5,
        }

    return {}


if __name__ == "__main__":
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description="Step 1: Download raw data")
    parser.add_argument("--accession", "-a", required=True, help="PRIDE accession")
    parser.add_argument("--config", "-c", default="config/config.yaml", help="Config file")
    parser.add_argument("--output-dir", "-o", default="data", help="Output base directory")
    parser.add_argument("--test-mode", action="store_true", help="Test mode")
    parser.add_argument("--max-files", type=int, default=3, help="Max files in test mode")
    parser.add_argument("--local-data", type=Path, help="Local data path (skip download)")

    args = parser.parse_args()

    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Run step
    summary = run_step1_download(
        accession=args.accession,
        config=config,
        output_base_dir=Path(args.output_dir),
        test_mode=args.test_mode,
        max_files=args.max_files if args.test_mode else 0,
        local_data_path=args.local_data,
    )

    print(f"\nStep 1 completed: {summary.status}")
    print(f"  Raw files: {summary.data['n_raw_files']}")
    print(f"  Size: {summary.data['total_size_gb']:.2f} GB")
