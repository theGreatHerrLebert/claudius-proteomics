#!/usr/bin/env python3
"""
Step 6: Package - Create Distributable Archive

Creates self-contained zip archive(s) of runner results for collection upload.
When sample_groups.yaml exists, creates one archive per group.

Outputs (per-group):
- data/packages/{accession}_{group_id}_v{version}.zip (one per group)
- step6_summary.json

Legacy (no sample_groups.yaml):
- data/packages/{accession}_v{version}.zip

Per-group Archive Structure:
{accession}_{group_id}/
├── manifest.json              # Required - metadata, QC, versions
├── precursor_store.parquet    # Required - main data (IDs + features)
├── precursor_index.parquet    # Engine-level identification details
├── consensus/
│   ├── overlap_stats.json     # 3-way Venn statistics
│   ├── overlap_report.html    # Visual QC report
│   └── stratified/            # PSMs by agreement tier
├── extracted/                 # Raw 4D signal data per raw file
│   └── {raw_file}.d/
│       ├── blobs.bin          # Raw 4D data (MS1 + MS2)
│       └── index.parquet      # Blob offsets/sizes
├── engines/                   # Search engine raw outputs
│   ├── fragpipe/
│   ├── diann/
│   └── sage/
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
    Execute Step 6: Create distributable zip archive(s).

    If sample_groups.yaml exists, creates one archive per group:
    {accession}_{group_id}_v{version}.zip. Otherwise creates a single
    legacy archive: {accession}_v{version}.zip.

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
    extracted_dir = output_base_dir / "extracted" / accession
    packages_dir = output_base_dir / "packages"

    try:
        packages_dir.mkdir(parents=True, exist_ok=True)

        # Check for sample group manifest
        metadata_dir = output_base_dir / "metadata" / accession
        sg_path = metadata_dir / "sample_groups.yaml"

        if sg_path.exists():
            from scripts.sample_group_resolver import SampleGroupManifest

            manifest = SampleGroupManifest.from_yaml(sg_path)
            print(f"  Loaded {len(manifest.groups)} sample groups from {sg_path}")

            group_results = {}
            all_archive_paths = []

            for group in manifest.groups:
                print(f"\n  === Sample group: {group.group_id} ===")

                if group.n_runs == 0:
                    print(f"      Skipping (no runs)")
                    group_results[group.group_id] = {"status": "skipped", "reason": "no runs"}
                    continue

                archive_prefix = f"{accession}_{group.group_id}"

                # Collect files scoped to this group
                files_to_include = _collect_files(
                    archive_prefix=archive_prefix,
                    processed_dir=processed_dir,
                    merged_dir=merged_dir,
                    extracted_dir=extracted_dir,
                    group_id=group.group_id,
                )

                # Validate required files exist
                _validate_required_files(files_to_include)

                # Create archive
                archive_name = f"{archive_prefix}_v{version}.zip"
                archive_path = packages_dir / archive_name

                print(f"  Creating archive: {archive_path}")
                archive_info = _create_archive(
                    archive_path=archive_path,
                    files_to_include=files_to_include,
                )

                checksum = _compute_checksum(archive_path)

                # Update manifest with archive info
                _update_manifest_with_archive_info(
                    merged_dir=merged_dir / group.group_id,
                    version=version,
                    checksum=checksum,
                    archive_size_bytes=archive_info["total_size"],
                )

                all_archive_paths.append(str(archive_path))

                group_results[group.group_id] = {
                    "status": "success",
                    "archive_name": archive_name,
                    "archive_path": str(archive_path),
                    "checksum_sha256": checksum,
                    "archive_size_mb": round(archive_info["total_size"] / (1024 * 1024), 2),
                    "n_files_included": archive_info["n_files"],
                }

                print(f"  Archive size: {group_results[group.group_id]['archive_size_mb']:.1f} MB")
                print(f"  Files included: {archive_info['n_files']}")
                print(f"  Checksum: {checksum[:16]}...")

            summary.data = {
                "mode": "per_group",
                "version": version,
                "n_archives": len(all_archive_paths),
                "group_results": group_results,
            }
            summary.outputs = all_archive_paths
            summary.complete(success=True)

        else:
            # Legacy single-group behavior
            print(f"  No sample_groups.yaml found — using single-group mode")

            archive_prefix = accession

            files_to_include = _collect_files(
                archive_prefix=archive_prefix,
                processed_dir=processed_dir,
                merged_dir=merged_dir,
                extracted_dir=extracted_dir,
            )

            _validate_required_files(files_to_include)

            archive_name = f"{accession}_v{version}.zip"
            archive_path = packages_dir / archive_name

            print(f"  Creating archive: {archive_path}")
            archive_info = _create_archive(
                archive_path=archive_path,
                files_to_include=files_to_include,
            )

            checksum = _compute_checksum(archive_path)

            _update_manifest_with_archive_info(
                merged_dir=merged_dir,
                version=version,
                checksum=checksum,
                archive_size_bytes=archive_info["total_size"],
            )

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
    archive_prefix: str,
    processed_dir: Path,
    merged_dir: Path,
    extracted_dir: Path,
    group_id: Optional[str] = None,
) -> Dict[str, Path]:
    """
    Collect all files to include in the archive.

    When group_id is provided, scopes all paths to the group's subdirectories.
    The archive_prefix is used as the top-level folder name in the zip
    (e.g. "PXD019086" for legacy or "PXD019086_yeast_trypsin" for per-group).

    Returns dict mapping archive paths to source paths.
    """
    files = {}

    # Resolve source directories based on whether this is per-group
    if group_id:
        src_merged = merged_dir / group_id
        src_processed = processed_dir / group_id
        src_extracted = extracted_dir / group_id
    else:
        src_merged = merged_dir
        src_processed = processed_dir
        src_extracted = extracted_dir

    # Required: manifest.json from merged dir
    manifest_path = src_merged / "manifest.json"
    if manifest_path.exists():
        files[f"{archive_prefix}/manifest.json"] = manifest_path

    # Required: precursor_store.parquet from merged dir
    store_path = src_merged / "precursor_store.parquet"
    if store_path.exists():
        files[f"{archive_prefix}/precursor_store.parquet"] = store_path

    # Optional: precursor_index.parquet from processed dir
    index_path = src_processed / "precursor_index.parquet"
    if index_path.exists():
        files[f"{archive_prefix}/precursor_index.parquet"] = index_path

    # Optional: consensus folder
    consensus_dir = src_processed / "consensus"
    if consensus_dir.exists():
        overlap_stats = consensus_dir / "overlap_stats.json"
        if overlap_stats.exists():
            files[f"{archive_prefix}/consensus/overlap_stats.json"] = overlap_stats

        overlap_report = consensus_dir / "overlap_report.html"
        if overlap_report.exists():
            files[f"{archive_prefix}/consensus/overlap_report.html"] = overlap_report

        stratified_dir = consensus_dir / "stratified"
        if stratified_dir.exists():
            for parquet_file in stratified_dir.glob("*.parquet"):
                archive_path = f"{archive_prefix}/consensus/stratified/{parquet_file.name}"
                files[archive_path] = parquet_file

    # Optional: step summaries (accession-level, shared across groups)
    for step_num in range(1, 6):
        summary_file = processed_dir / f"step{step_num}_summary.json"
        if summary_file.exists():
            files[f"{archive_prefix}/summaries/step{step_num}_summary.json"] = summary_file

    # Raw 4D signal data from extracted dir (blobs + index per raw file)
    if src_extracted.exists():
        for raw_file_dir in src_extracted.iterdir():
            if raw_file_dir.is_dir() and raw_file_dir.name.endswith('.d'):
                blobs_file = raw_file_dir / "blobs.bin"
                if blobs_file.exists():
                    archive_path = f"{archive_prefix}/extracted/{raw_file_dir.name}/blobs.bin"
                    files[archive_path] = blobs_file

                index_file = raw_file_dir / "index.parquet"
                if index_file.exists():
                    archive_path = f"{archive_prefix}/extracted/{raw_file_dir.name}/index.parquet"
                    files[archive_path] = index_file

    # ==========================================================================
    # Search Engine Raw Outputs
    # ==========================================================================

    # FragPipe output
    fragpipe_dir = src_processed / "fragpipe_output"
    if fragpipe_dir.exists():
        for config_file in ["fragpipe.workflow", "workflow.workflow", "fragger.params"]:
            config_path = fragpipe_dir / config_file
            if config_path.exists():
                files[f"{archive_prefix}/engines/fragpipe/{config_file}"] = config_path

        for combined_file in fragpipe_dir.glob("combined_*.tsv"):
            files[f"{archive_prefix}/engines/fragpipe/{combined_file.name}"] = combined_file

        for raw_file_dir in fragpipe_dir.iterdir():
            if raw_file_dir.is_dir() and not raw_file_dir.name.startswith('.'):
                if raw_file_dir.name == "MSBooster":
                    continue
                for result_file in ["psm.tsv", "ion.tsv", "peptide.tsv", "protein.tsv"]:
                    result_path = raw_file_dir / result_file
                    if result_path.exists():
                        archive_path = f"{archive_prefix}/engines/fragpipe/per_file/{raw_file_dir.name}/{result_file}"
                        files[archive_path] = result_path

    # DIA-NN output
    diann_dir = src_processed / "diann"
    if diann_dir.exists():
        for diann_file in ["report.parquet", "report.log.txt", "report.stats.tsv"]:
            diann_path = diann_dir / diann_file
            if diann_path.exists():
                files[f"{archive_prefix}/engines/diann/{diann_file}"] = diann_path

    # Sage output
    sage_dir = src_processed / "sage"
    if sage_dir.exists():
        for sage_file in ["results.sage.parquet", "matched_fragments.sage.parquet",
                          "lfq.parquet", "results.json"]:
            sage_path = sage_dir / sage_file
            if sage_path.exists():
                files[f"{archive_prefix}/engines/sage/{sage_file}"] = sage_path

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
    if summary.data.get("mode") == "per_group":
        print(f"  Archives created: {summary.data['n_archives']}")
        for gid, gres in summary.data["group_results"].items():
            if gres.get("status") == "success":
                print(f"    {gres['archive_name']}: {gres['archive_size_mb']} MB, {gres['n_files_included']} files")
    else:
        print(f"  Archive: {summary.data['archive_name']}")
        print(f"  Size: {summary.data['archive_size_mb']} MB")
        print(f"  Files: {summary.data['n_files_included']}")
