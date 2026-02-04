#!/usr/bin/env python3
"""
Step 6: Package - Create Distributable Archive

Creates a self-contained zip archive of runner results for collection upload.

Input:
- data/merged/{accession}/precursor_store.parquet (Step 5)
- data/merged/{accession}/manifest.json (Step 5)
- data/processed/{accession}/precursor_index.parquet (Step 3)
- data/processed/{accession}/consensus/ (Step 3)
- data/processed/{accession}/step*_summary.json (Steps 1-5)

Outputs:
- data/packages/{accession}_v{version}.zip
- step6_summary.json

Archive Structure:
{accession}/
├── manifest.json              # Required - metadata, QC, versions
├── precursor_store.parquet    # Required - main data (IDs + features)
├── precursor_index.parquet    # Engine-level identification details
├── consensus/
│   ├── overlap_stats.json     # 3-way Venn statistics
│   ├── overlap_report.html    # Visual QC report
│   └── stratified/            # PSMs by agreement tier
│       ├── all_three.parquet
│       ├── two_plus.parquet
│       ├── fragpipe_only.parquet
│       ├── diann_only.parquet
│       └── sage_only.parquet
├── extracted/                 # Raw 4D signal data per raw file
│   └── {raw_file}.d/
│       ├── blobs.bin          # Raw 4D data (MS1 + MS2)
│       └── index.parquet      # Blob offsets/sizes
├── engines/                   # Search engine raw outputs
│   ├── fragpipe/
│   │   ├── fragpipe.workflow  # Settings
│   │   ├── fragger.params
│   │   ├── combined_*.tsv     # Combined results
│   │   └── per_file/{raw}/    # Per-file psm/ion/peptide/protein.tsv
│   ├── diann/
│   │   ├── report.parquet     # Main results
│   │   ├── report.log.txt     # Settings/log
│   │   └── report.stats.tsv
│   └── sage/
│       ├── results.sage.parquet
│       ├── matched_fragments.sage.parquet
│       ├── lfq.parquet
│       └── results.json       # Settings
└── summaries/
    ├── step1_summary.json
    ├── step2_summary.json
    ├── step3_summary.json
    ├── step4_summary.json
    └── step5_summary.json
"""

import hashlib
import json
import sys
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from runner.summary import StepSummary, write_step_summary


def run_step6_package(
    accession: str,
    config: Dict[str, Any],
    output_base_dir: Path,
    version: str = "1.0",
) -> StepSummary:
    """
    Execute Step 6: Create distributable zip archive.

    Args:
        accession: PRIDE accession
        config: Pipeline configuration dict
        output_base_dir: Base directory for outputs
        version: Archive version string (default: 1.0)

    Returns:
        StepSummary with results
    """
    summary = StepSummary(
        step_name="step6",
        accession=accession,
    )

    processed_dir = output_base_dir / "processed" / accession
    merged_dir = output_base_dir / "merged" / accession
    packages_dir = output_base_dir / "packages"

    try:
        # Create packages directory
        packages_dir.mkdir(parents=True, exist_ok=True)

        # Collect files to include
        extracted_dir = output_base_dir / "extracted" / accession
        files_to_include = _collect_files(
            accession=accession,
            processed_dir=processed_dir,
            merged_dir=merged_dir,
            extracted_dir=extracted_dir,
        )

        # Validate required files exist
        _validate_required_files(files_to_include)

        # Create archive
        archive_name = f"{accession}_v{version}.zip"
        archive_path = packages_dir / archive_name

        print(f"  Creating archive: {archive_path}")
        archive_info = _create_archive(
            archive_path=archive_path,
            accession=accession,
            files_to_include=files_to_include,
        )

        # Compute checksum
        checksum = _compute_checksum(archive_path)

        # Update manifest with archive info
        _update_manifest_with_archive_info(
            merged_dir=merged_dir,
            version=version,
            checksum=checksum,
            archive_size_bytes=archive_info["total_size"],
        )

        # Update summary
        summary.data = {
            "archive_name": archive_name,
            "archive_path": str(archive_path),
            "version": version,
            "checksum_sha256": checksum,
            "archive_size_bytes": archive_info["total_size"],
            "archive_size_mb": round(archive_info["total_size"] / (1024 * 1024), 2),
            "n_files_included": archive_info["n_files"],
            "files_included": archive_info["files"],
        }
        summary.outputs = [str(archive_path)]
        summary.complete(success=True)

        print(f"  Archive size: {summary.data['archive_size_mb']:.1f} MB")
        print(f"  Files included: {archive_info['n_files']}")
        print(f"  Checksum: {checksum[:16]}...")

    except Exception as e:
        summary.complete(success=False, error_message=str(e))
        raise

    # Write summary file
    write_step_summary(summary, processed_dir)

    return summary


def _collect_files(
    accession: str,
    processed_dir: Path,
    merged_dir: Path,
    extracted_dir: Path,
) -> Dict[str, Path]:
    """
    Collect all files to include in the archive.

    Returns dict mapping archive paths to source paths.
    """
    files = {}

    # Required: manifest.json from merged dir
    manifest_path = merged_dir / "manifest.json"
    if manifest_path.exists():
        files[f"{accession}/manifest.json"] = manifest_path

    # Required: precursor_store.parquet from merged dir
    store_path = merged_dir / "precursor_store.parquet"
    if store_path.exists():
        files[f"{accession}/precursor_store.parquet"] = store_path

    # Optional: precursor_index.parquet from processed dir
    index_path = processed_dir / "precursor_index.parquet"
    if index_path.exists():
        files[f"{accession}/precursor_index.parquet"] = index_path

    # Optional: consensus folder
    consensus_dir = processed_dir / "consensus"
    if consensus_dir.exists():
        # overlap_stats.json
        overlap_stats = consensus_dir / "overlap_stats.json"
        if overlap_stats.exists():
            files[f"{accession}/consensus/overlap_stats.json"] = overlap_stats

        # overlap_report.html
        overlap_report = consensus_dir / "overlap_report.html"
        if overlap_report.exists():
            files[f"{accession}/consensus/overlap_report.html"] = overlap_report

        # stratified folder
        stratified_dir = consensus_dir / "stratified"
        if stratified_dir.exists():
            for parquet_file in stratified_dir.glob("*.parquet"):
                archive_path = f"{accession}/consensus/stratified/{parquet_file.name}"
                files[archive_path] = parquet_file

    # Optional: step summaries
    for step_num in range(1, 6):
        summary_file = processed_dir / f"step{step_num}_summary.json"
        if summary_file.exists():
            files[f"{accession}/summaries/step{step_num}_summary.json"] = summary_file

    # Raw 4D signal data from extracted dir (blobs + index per raw file)
    if extracted_dir.exists():
        for raw_file_dir in extracted_dir.iterdir():
            if raw_file_dir.is_dir() and raw_file_dir.name.endswith('.d'):
                # Include blobs.bin (raw 4D data)
                blobs_file = raw_file_dir / "blobs.bin"
                if blobs_file.exists():
                    archive_path = f"{accession}/extracted/{raw_file_dir.name}/blobs.bin"
                    files[archive_path] = blobs_file

                # Include index.parquet (blob offsets)
                index_file = raw_file_dir / "index.parquet"
                if index_file.exists():
                    archive_path = f"{accession}/extracted/{raw_file_dir.name}/index.parquet"
                    files[archive_path] = index_file

    # ==========================================================================
    # Search Engine Raw Outputs
    # ==========================================================================

    # FragPipe output
    fragpipe_dir = processed_dir / "fragpipe_output"
    if fragpipe_dir.exists():
        # Settings/config files
        for config_file in ["fragpipe.workflow", "workflow.workflow", "fragger.params"]:
            config_path = fragpipe_dir / config_file
            if config_path.exists():
                files[f"{accession}/engines/fragpipe/{config_file}"] = config_path

        # Combined results
        for combined_file in fragpipe_dir.glob("combined_*.tsv"):
            files[f"{accession}/engines/fragpipe/{combined_file.name}"] = combined_file

        # Per-file results (psm.tsv, ion.tsv, peptide.tsv, protein.tsv)
        for raw_file_dir in fragpipe_dir.iterdir():
            if raw_file_dir.is_dir() and not raw_file_dir.name.startswith('.'):
                # Skip MSBooster subdirectory
                if raw_file_dir.name == "MSBooster":
                    continue
                for result_file in ["psm.tsv", "ion.tsv", "peptide.tsv", "protein.tsv"]:
                    result_path = raw_file_dir / result_file
                    if result_path.exists():
                        archive_path = f"{accession}/engines/fragpipe/per_file/{raw_file_dir.name}/{result_file}"
                        files[archive_path] = result_path

    # DIA-NN output
    diann_dir = processed_dir / "diann"
    if diann_dir.exists():
        # Main results
        for diann_file in ["report.parquet", "report.log.txt", "report.stats.tsv"]:
            diann_path = diann_dir / diann_file
            if diann_path.exists():
                files[f"{accession}/engines/diann/{diann_file}"] = diann_path

    # Sage output
    sage_dir = processed_dir / "sage"
    if sage_dir.exists():
        # Results and settings
        for sage_file in ["results.sage.parquet", "matched_fragments.sage.parquet",
                          "lfq.parquet", "results.json"]:
            sage_path = sage_dir / sage_file
            if sage_path.exists():
                files[f"{accession}/engines/sage/{sage_file}"] = sage_path

    return files


def _validate_required_files(files: Dict[str, Path]) -> None:
    """Validate that required files are present."""
    required_patterns = [
        "manifest.json",
        "precursor_store.parquet",
    ]

    missing = []
    for pattern in required_patterns:
        if not any(pattern in path for path in files.keys()):
            missing.append(pattern)

    if missing:
        raise FileNotFoundError(
            f"Required files missing: {missing}. "
            "Ensure steps 1-5 completed successfully."
        )


def _create_archive(
    archive_path: Path,
    accession: str,
    files_to_include: Dict[str, Path],
) -> Dict[str, Any]:
    """
    Create the zip archive.

    Returns info about the created archive.
    """
    total_size = 0
    files_added = []

    with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED) as zf:
        for archive_path_str, source_path in sorted(files_to_include.items()):
            if source_path.exists():
                zf.write(source_path, archive_path_str)
                total_size += source_path.stat().st_size
                files_added.append(archive_path_str)
                print(f"    + {archive_path_str}")

    # Get compressed size
    compressed_size = archive_path.stat().st_size

    return {
        "n_files": len(files_added),
        "files": files_added,
        "total_size": compressed_size,
        "uncompressed_size": total_size,
        "compression_ratio": round(compressed_size / max(total_size, 1), 3),
    }


def _compute_checksum(file_path: Path, algorithm: str = "sha256") -> str:
    """Compute file checksum."""
    hasher = hashlib.new(algorithm)

    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            hasher.update(chunk)

    return hasher.hexdigest()


def _update_manifest_with_archive_info(
    merged_dir: Path,
    version: str,
    checksum: str,
    archive_size_bytes: int,
) -> None:
    """Update manifest.json with archive information."""
    manifest_path = merged_dir / "manifest.json"

    if manifest_path.exists():
        with open(manifest_path, 'r') as f:
            manifest = json.load(f)
    else:
        manifest = {}

    manifest["archive_info"] = {
        "version": version,
        "checksum_sha256": checksum,
        "archive_size_bytes": archive_size_bytes,
        "packaged_at": datetime.now().isoformat(),
    }

    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)


if __name__ == "__main__":
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description="Step 6: Package runner output")
    parser.add_argument("--accession", "-a", required=True, help="PRIDE accession")
    parser.add_argument("--config", "-c", default="config/config.yaml", help="Config file")
    parser.add_argument("--output-dir", "-o", default="data", help="Output base directory")
    parser.add_argument("--version", "-v", default="1.0", help="Archive version")

    args = parser.parse_args()

    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Run step
    summary = run_step6_package(
        accession=args.accession,
        config=config,
        output_base_dir=Path(args.output_dir),
        version=args.version,
    )

    print(f"\nStep 6 completed: {summary.status}")
    print(f"  Archive: {summary.data['archive_name']}")
    print(f"  Size: {summary.data['archive_size_mb']} MB")
    print(f"  Files: {summary.data['n_files_included']}")
