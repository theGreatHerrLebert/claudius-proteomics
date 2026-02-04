#!/usr/bin/env python3
"""
Add Dataset Archive to Collection

Extracts a runner archive to a study folder and updates the collection manifest.

Usage:
    python scripts/add_to_collection.py \\
        --archive data/packages/PXD019086_v1.0.zip \\
        --collection /path/to/san_jose_collection/ \\
        --study meier_2021_ccs

Collection Structure:
    /san_jose_collection/
    ├── studies.yaml                   # Study definitions (manual)
    ├── collection_manifest.json       # Auto-generated index
    ├── meier_2021_ccs/               # Study folder
    │   └── PXD019086_v1.0/           # Extracted dataset
    └── archives/                     # Original zips (optional backup)
"""

import argparse
import json
import shutil
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional

import yaml


def add_to_collection(
    archive_path: Path,
    collection_path: Path,
    study_id: str,
    keep_archive: bool = True,
) -> Dict[str, Any]:
    """
    Add a dataset archive to a collection.

    Args:
        archive_path: Path to the zip archive
        collection_path: Path to the collection root
        study_id: ID of the study to add the dataset to
        keep_archive: If True, copy archive to archives/ folder

    Returns:
        Dict with info about the added dataset
    """
    # Validate inputs
    if not archive_path.exists():
        raise FileNotFoundError(f"Archive not found: {archive_path}")

    if not collection_path.exists():
        raise FileNotFoundError(f"Collection path not found: {collection_path}")

    # Load studies.yaml
    studies_config = _load_studies_config(collection_path)
    if study_id not in [s["id"] for s in studies_config.get("studies", [])]:
        available = [s["id"] for s in studies_config.get("studies", [])]
        raise ValueError(
            f"Study '{study_id}' not found in studies.yaml. "
            f"Available studies: {available}"
        )

    # Parse archive name to get accession and version
    archive_name = archive_path.stem  # e.g., "PXD019086_v1.0"
    parts = archive_name.rsplit("_v", 1)
    if len(parts) != 2:
        raise ValueError(
            f"Invalid archive name format: {archive_name}. "
            "Expected format: {accession}_v{version}.zip"
        )
    accession, version = parts

    print(f"Adding {accession} v{version} to study '{study_id}'")

    # Create study folder if needed
    study_dir = collection_path / study_id
    study_dir.mkdir(parents=True, exist_ok=True)

    # Extract archive
    dataset_dir = study_dir / f"{accession}_v{version}"
    if dataset_dir.exists():
        print(f"  Removing existing dataset at {dataset_dir}")
        shutil.rmtree(dataset_dir)

    print(f"  Extracting to {dataset_dir}")
    _extract_archive(archive_path, dataset_dir, accession)

    # Load manifest from extracted dataset
    manifest_path = dataset_dir / "manifest.json"
    if manifest_path.exists():
        with open(manifest_path) as f:
            manifest = json.load(f)
    else:
        manifest = {"accession": accession}

    # Compute dataset summary
    dataset_info = _compute_dataset_info(
        accession=accession,
        version=version,
        study_id=study_id,
        dataset_dir=dataset_dir,
        manifest=manifest,
    )

    # Copy archive to archives folder if requested
    if keep_archive:
        archives_dir = collection_path / "archives"
        archives_dir.mkdir(parents=True, exist_ok=True)
        archive_dest = archives_dir / archive_path.name
        if not archive_dest.exists():
            shutil.copy2(archive_path, archive_dest)
            print(f"  Archived to {archive_dest}")

    # Update collection manifest
    _update_collection_manifest(collection_path, dataset_info, studies_config)

    print(f"  Added {accession} with {dataset_info['n_precursors']} precursors")
    return dataset_info


def _load_studies_config(collection_path: Path) -> Dict[str, Any]:
    """Load studies.yaml configuration."""
    studies_path = collection_path / "studies.yaml"
    if not studies_path.exists():
        raise FileNotFoundError(
            f"studies.yaml not found in {collection_path}. "
            "Create it with study definitions first."
        )

    with open(studies_path) as f:
        return yaml.safe_load(f)


def _extract_archive(
    archive_path: Path,
    dest_dir: Path,
    accession: str,
) -> None:
    """Extract archive contents to destination directory."""
    dest_dir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(archive_path, 'r') as zf:
        # Get all members
        members = zf.namelist()

        # Extract each member, stripping the accession prefix
        for member in members:
            # Member paths are like "PXD019086/manifest.json"
            if member.startswith(f"{accession}/"):
                # Strip prefix
                rel_path = member[len(f"{accession}/"):]
                if not rel_path:
                    continue

                dest_path = dest_dir / rel_path

                # Create parent directories
                dest_path.parent.mkdir(parents=True, exist_ok=True)

                # Extract file
                if not member.endswith('/'):
                    with zf.open(member) as src, open(dest_path, 'wb') as dst:
                        shutil.copyfileobj(src, dst)


def _compute_dataset_info(
    accession: str,
    version: str,
    study_id: str,
    dataset_dir: Path,
    manifest: Dict[str, Any],
) -> Dict[str, Any]:
    """Compute dataset summary information."""
    import pandas as pd

    info = {
        "accession": accession,
        "version": version,
        "study_id": study_id,
        "path": f"{study_id}/{accession}_v{version}",
        "added_at": datetime.now().isoformat(),
    }

    # Try to load precursor store for statistics
    store_path = dataset_dir / "precursor_store.parquet"
    if store_path.exists():
        try:
            df = pd.read_parquet(store_path)
            info["n_precursors"] = len(df)

            # Engine agreement stats
            if "n_engines" in df.columns:
                info["n_all_three"] = int((df["n_engines"] == 3).sum())
                info["n_at_least_two"] = int((df["n_engines"] >= 2).sum())

            # Quality stats
            quality = {}
            if "ms1_rt_r2" in df.columns:
                quality["rt_r2_median"] = round(float(df["ms1_rt_r2"].median()), 4)
            if "ms1_im_r2" in df.columns:
                quality["im_r2_median"] = round(float(df["ms1_im_r2"].median()), 4)
            if "is_high_quality" in df.columns:
                pct = df["is_high_quality"].sum() / max(len(df), 1) * 100
                quality["pct_high_quality"] = round(pct, 1)

            if quality:
                info["quality"] = quality

        except Exception as e:
            print(f"  Warning: Could not read precursor store: {e}")
            info["n_precursors"] = manifest.get("n_total_precursors", 0)
    else:
        info["n_precursors"] = manifest.get("n_total_precursors", 0)

    # Copy quality summary from manifest if available
    if "quality_summary" in manifest and "quality" not in info:
        info["quality"] = manifest["quality_summary"]

    return info


def _update_collection_manifest(
    collection_path: Path,
    dataset_info: Dict[str, Any],
    studies_config: Dict[str, Any],
) -> None:
    """Update the collection manifest with new dataset info."""
    manifest_path = collection_path / "collection_manifest.json"

    # Load existing manifest or create new
    if manifest_path.exists():
        with open(manifest_path) as f:
            manifest = json.load(f)
    else:
        manifest = {
            "version": "1.0",
            "updated_at": None,
            "studies": [],
            "datasets": [],
        }

    # Update timestamp
    manifest["updated_at"] = datetime.now().isoformat()

    # Update datasets list
    datasets = manifest.get("datasets", [])

    # Remove existing entry for this accession (if re-adding)
    datasets = [d for d in datasets if d["accession"] != dataset_info["accession"]]

    # Add new entry
    datasets.append(dataset_info)

    # Sort by accession
    datasets.sort(key=lambda d: d["accession"])
    manifest["datasets"] = datasets

    # Rebuild studies summary
    studies_summary = []
    for study in studies_config.get("studies", []):
        study_datasets = [d for d in datasets if d["study_id"] == study["id"]]
        n_precursors = sum(d.get("n_precursors", 0) for d in study_datasets)

        studies_summary.append({
            "id": study["id"],
            "title": study.get("title", study["id"]),
            "organism": study.get("organism", ""),
            "n_datasets": len(study_datasets),
            "n_total_precursors": n_precursors,
            "datasets": [d["accession"] for d in study_datasets],
        })

    manifest["studies"] = studies_summary

    # Write manifest
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"  Updated collection manifest: {manifest_path}")


def rebuild_collection_manifest(collection_path: Path) -> None:
    """
    Rebuild collection manifest from disk contents.

    Scans all study folders and regenerates the manifest.
    """
    import pandas as pd

    print(f"Rebuilding collection manifest at {collection_path}")

    # Load studies config
    studies_config = _load_studies_config(collection_path)

    datasets = []

    # Scan each study folder
    for study in studies_config.get("studies", []):
        study_id = study["id"]
        study_dir = collection_path / study_id

        if not study_dir.exists():
            continue

        # Find all dataset folders
        for dataset_dir in study_dir.iterdir():
            if not dataset_dir.is_dir():
                continue

            # Parse folder name
            name = dataset_dir.name
            parts = name.rsplit("_v", 1)
            if len(parts) != 2:
                continue

            accession, version = parts

            # Load manifest
            manifest_path = dataset_dir / "manifest.json"
            if manifest_path.exists():
                with open(manifest_path) as f:
                    manifest = json.load(f)
            else:
                manifest = {"accession": accession}

            # Compute info
            info = _compute_dataset_info(
                accession=accession,
                version=version,
                study_id=study_id,
                dataset_dir=dataset_dir,
                manifest=manifest,
            )
            datasets.append(info)

            print(f"  Found {accession} in {study_id}")

    # Build manifest
    manifest = {
        "version": "1.0",
        "updated_at": datetime.now().isoformat(),
        "studies": [],
        "datasets": sorted(datasets, key=lambda d: d["accession"]),
    }

    # Build studies summary
    for study in studies_config.get("studies", []):
        study_datasets = [d for d in datasets if d["study_id"] == study["id"]]
        n_precursors = sum(d.get("n_precursors", 0) for d in study_datasets)

        manifest["studies"].append({
            "id": study["id"],
            "title": study.get("title", study["id"]),
            "organism": study.get("organism", ""),
            "n_datasets": len(study_datasets),
            "n_total_precursors": n_precursors,
            "datasets": [d["accession"] for d in study_datasets],
        })

    # Write manifest
    manifest_path = collection_path / "collection_manifest.json"
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"\nRebuilt manifest with {len(datasets)} datasets")


def main():
    parser = argparse.ArgumentParser(
        description="Add dataset archive to San José collection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Add a single dataset
    python scripts/add_to_collection.py \\
        --archive data/packages/PXD019086_v1.0.zip \\
        --collection /path/to/collection \\
        --study meier_2021_ccs

    # Rebuild manifest from existing files
    python scripts/add_to_collection.py \\
        --collection /path/to/collection \\
        --rebuild
        """,
    )

    parser.add_argument(
        "--archive", "-a",
        type=Path,
        help="Path to dataset archive (.zip)",
    )
    parser.add_argument(
        "--collection", "-c",
        type=Path,
        required=True,
        help="Path to collection root directory",
    )
    parser.add_argument(
        "--study", "-s",
        help="Study ID to add the dataset to",
    )
    parser.add_argument(
        "--no-archive",
        action="store_true",
        help="Don't copy archive to archives/ folder",
    )
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Rebuild collection manifest from disk contents",
    )

    args = parser.parse_args()

    if args.rebuild:
        rebuild_collection_manifest(args.collection)
    else:
        if not args.archive:
            parser.error("--archive is required when not using --rebuild")
        if not args.study:
            parser.error("--study is required when not using --rebuild")

        add_to_collection(
            archive_path=args.archive,
            collection_path=args.collection,
            study_id=args.study,
            keep_archive=not args.no_archive,
        )


if __name__ == "__main__":
    main()
