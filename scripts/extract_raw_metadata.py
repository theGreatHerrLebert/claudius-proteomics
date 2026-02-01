#!/usr/bin/env python3
"""
Merge metadata from PRIDE API and raw timsTOF .d folders.

This script combines metadata from multiple sources with priority handling:
1. PRIDE API (highest priority for identifiers and scientific context)
2. Raw .d folders (best source for LC/MS parameters)
3. Manual entry (fallback for missing fields)
"""

import argparse
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

from pride_metadata import DatasetMetadata, MetadataField, fetch_pride_metadata
from raw_metadata import (
    RawMetadata,
    aggregate_metadata,
    extract_all_raw_metadata,
    extract_raw_metadata,
)


def merge_metadata(
    pride_meta: Optional[DatasetMetadata],
    raw_meta: Dict[str, Any],
    manual_overrides: Optional[Dict[str, Any]] = None,
) -> DatasetMetadata:
    """
    Merge metadata from multiple sources with priority.

    Priority: PRIDE API > raw .d folder > manual entry

    For each field, the highest-priority non-null value is kept.
    Raw data is best for LC/MS parameters (gradient, mode, instrument confirmation).
    PRIDE API is best for identifiers and scientific context.

    Args:
        pride_meta: Metadata from PRIDE API (or None)
        raw_meta: Aggregated metadata from raw .d folders
        manual_overrides: Optional manual values to apply

    Returns:
        Merged DatasetMetadata
    """
    # Start with PRIDE metadata or create empty
    if pride_meta:
        merged = pride_meta
    else:
        merged = DatasetMetadata(accession="unknown")

    # Update from raw data where PRIDE is missing or raw is more accurate
    # Raw data is authoritative for LC/MS parameters

    # Instrument: confirm/update from raw data
    raw_instrument = raw_meta.get("instrument", {}).get("value")
    if raw_instrument:
        if merged.instrument.status == "missing":
            merged.instrument = MetadataField.auto(
                raw_instrument, "raw_data.GlobalMetadata.InstrumentName"
            )
        elif merged.instrument.value != raw_instrument:
            # Raw data confirms/updates instrument
            merged.instrument = MetadataField.auto(
                raw_instrument,
                "raw_data.GlobalMetadata.InstrumentName",
            )

    # Gradient length: raw data is primary source
    raw_gradient = raw_meta.get("gradient_length", {}).get("value")
    if raw_gradient:
        merged.gradient_length = MetadataField.auto(
            raw_gradient, "raw_data.lc_method_name"
        )
        # Add confirmed_by if we have run duration
        run_duration = raw_meta.get("run_duration", {})
        if run_duration.get("mean"):
            # Store as additional context (could add to source string)
            pass

    # Acquisition mode: raw data is primary source
    raw_mode = raw_meta.get("acquisition_mode", {}).get("value")
    if raw_mode:
        merged.acquisition_mode = MetadataField.auto(
            raw_mode, "raw_data.ms_method_name"
        )

    # LC system: raw data only
    raw_lc = raw_meta.get("lc_system", {}).get("value")
    if raw_lc:
        merged.lc_system = MetadataField.auto(
            raw_lc, "raw_data.chromatography.TraceSources"
        )

    # Apply manual overrides
    if manual_overrides:
        for field_name, value in manual_overrides.items():
            if value is not None and hasattr(merged, field_name):
                setattr(merged, field_name, MetadataField.manual(value))

    # Update timestamp
    merged.generated_at = datetime.now().isoformat()

    return merged


def extract_from_all_raw_files(raw_dir: Path) -> Dict[str, Any]:
    """
    Extract and aggregate metadata from all .d folders in a directory.

    Validates consistency across files and returns aggregated values.

    Args:
        raw_dir: Directory containing .d folders

    Returns:
        Aggregated metadata dictionary
    """
    raw_dir = Path(raw_dir)

    if not raw_dir.exists():
        print(f"Error: Raw directory does not exist: {raw_dir}")
        return {}

    # Check if it's a single .d folder or directory of .d folders
    if raw_dir.suffix == ".d" and raw_dir.is_dir():
        # Single .d folder
        metadata = extract_raw_metadata(raw_dir)
        return aggregate_metadata([metadata])

    # Directory with multiple .d folders
    d_folders = list(raw_dir.glob("*.d"))
    if not d_folders:
        print(f"Warning: No .d folders found in {raw_dir}")
        return {}

    print(f"Found {len(d_folders)} .d folders in {raw_dir}")

    metadata_list = extract_all_raw_metadata(raw_dir)

    if not metadata_list:
        print("Warning: Could not extract metadata from any .d folder")
        return {}

    aggregated = aggregate_metadata(metadata_list)

    # Report consistency issues
    for field in ["instrument", "lc_system", "gradient_length", "acquisition_mode"]:
        field_info = aggregated.get(field, {})
        if field_info.get("unique_values"):
            print(
                f"Warning: Inconsistent {field} across files: "
                f"{field_info['unique_values']}"
            )

    return aggregated


def load_pride_metadata(path: Path) -> Optional[DatasetMetadata]:
    """Load PRIDE metadata from YAML file."""
    if not path.exists():
        return None
    return DatasetMetadata.from_yaml(path)


def main():
    parser = argparse.ArgumentParser(
        description="Extract and merge metadata from PRIDE API and raw .d folders"
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        help="Directory containing .d folders (or single .d folder)"
    )
    parser.add_argument(
        "--pride-meta",
        type=Path,
        help="PRIDE metadata YAML file (from download_pride.py --metadata-only)"
    )
    parser.add_argument(
        "--accession",
        help="PRIDE accession to fetch metadata from API (if --pride-meta not provided)"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        help="Output YAML file for merged metadata"
    )
    parser.add_argument(
        "--manual",
        type=Path,
        help="YAML file with manual overrides"
    )

    args = parser.parse_args()

    # Load or fetch PRIDE metadata
    pride_meta = None
    if args.pride_meta:
        print(f"Loading PRIDE metadata from {args.pride_meta}")
        pride_meta = load_pride_metadata(args.pride_meta)
    elif args.accession:
        print(f"Fetching PRIDE metadata for {args.accession}")
        pride_meta = fetch_pride_metadata(args.accession)

    # Extract raw metadata
    raw_meta = {}
    if args.raw_dir:
        print(f"Extracting metadata from raw files in {args.raw_dir}")
        raw_meta = extract_from_all_raw_files(args.raw_dir)

        if raw_meta:
            print(f"  Extracted from {raw_meta.get('num_files', 0)} files")
            if raw_meta.get("instrument", {}).get("value"):
                print(f"  Instrument: {raw_meta['instrument']['value']}")
            if raw_meta.get("gradient_length", {}).get("value"):
                print(f"  Gradient: {raw_meta['gradient_length']['value']} min")
            if raw_meta.get("acquisition_mode", {}).get("value"):
                print(f"  Mode: {raw_meta['acquisition_mode']['value']}")
            if raw_meta.get("lc_system", {}).get("value"):
                print(f"  LC system: {raw_meta['lc_system']['value']}")

    # Load manual overrides
    manual_overrides = None
    if args.manual and args.manual.exists():
        with open(args.manual) as f:
            manual_overrides = yaml.safe_load(f)

    # Check we have at least some input
    if not pride_meta and not raw_meta:
        print("Error: No metadata sources provided. Use --pride-meta, --accession, or --raw-dir")
        sys.exit(1)

    # Merge metadata
    print("\nMerging metadata...")
    merged = merge_metadata(pride_meta, raw_meta, manual_overrides)

    # Write output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_yaml(args.output)
    print(f"\nWrote merged metadata to {args.output}")

    # Print validation summary
    summary = merged.validation_summary()
    print(f"\nValidation summary:")
    print(f"  Complete: {summary['complete']}")
    print(f"  Auto fields: {summary['auto_fields']}")
    print(f"  Inferred fields: {summary['inferred_fields']}")
    print(f"  Manual fields: {summary['manual_fields']}")
    print(f"  Missing fields: {summary['missing_fields']}")

    if summary["missing_field_names"]:
        print(f"\nMissing fields that may need manual entry:")
        for field_name in summary["missing_field_names"]:
            field = getattr(merged, field_name, None)
            if field and field.hint:
                print(f"  - {field_name}: {field.hint}")
            else:
                print(f"  - {field_name}")

    if not summary["complete"]:
        print(f"\nTo complete metadata, edit {args.output} and fill missing fields,")
        print("then run validate_metadata to confirm.")

    return 0 if summary["complete"] else 1


if __name__ == "__main__":
    sys.exit(main())
