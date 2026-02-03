#!/usr/bin/env python3
"""
Extract Ion Mobility Calibration Lookup Tables

Probes the Bruker SDK once per dataset to extract a complete scan→1/K0 lookup table.
This enables accurate ion mobility calibration with fast parallel extraction.

Background:
- The Bruker calibration formula is patented and proprietary
- Using use_bruker_sdk=True gives accurate values but is slow (not thread-safe)
- Using use_bruker_sdk=False uses linear interpolation (fast but inaccurate)
- This script extracts the calibration once, then it can be used for fast parallel extraction

The lookup table is small (~8KB for 1000 scans) and constant across all frames.

Usage:
    python scripts/extract_calibration.py data/raw/PXD019086/sample.d
    python scripts/extract_calibration.py data/raw/PXD019086/sample.d --output calibration.npy
"""

import argparse
from pathlib import Path
from typing import Optional, Tuple
import numpy as np


def extract_im_calibration(d_path: str, verbose: bool = False) -> np.ndarray:
    """Extract ion mobility calibration lookup table using Bruker SDK.

    Opens the dataset with the Bruker SDK to get accurate scan→1/K0 conversion,
    then builds a complete lookup table for all possible scan indices.

    Args:
        d_path: Path to the .d folder
        verbose: Print progress information

    Returns:
        Array of shape (num_scans,) mapping scan_index → 1/K0 value

    Raises:
        RuntimeError: If Bruker SDK is not available or fails
    """
    from imspy_core.timstof import TimsDatasetDDA

    if verbose:
        print(f"Loading dataset with Bruker SDK: {d_path}")

    # Open with Bruker SDK for accurate calibration
    # This is slow and not thread-safe, but we only do it once
    dataset = TimsDatasetDDA(str(d_path), in_memory=False, use_bruker_sdk=True)

    # Get the number of scans from metadata
    # The scan range is typically 0 to scan_max_index (inclusive)
    meta_data = dataset.meta_data
    # Column name may be 'NumScans' or 'num_scans' depending on renaming
    num_scans_col = 'NumScans' if 'NumScans' in meta_data.columns else 'num_scans'
    scan_max_index = meta_data[num_scans_col].max()

    if verbose:
        print(f"  Frame count: {dataset.frame_count}")
        print(f"  Max scan index: {scan_max_index}")

    # Build lookup table: convert ALL scan indices to 1/K0
    # Use frame_id=1 (any MS1 frame works - calibration is typically frame-independent)
    # The calibration parameters are stored once per acquisition in TimsCalibration table
    num_scans = int(scan_max_index) + 1
    scans = np.arange(num_scans, dtype=np.int32)

    if verbose:
        print(f"  Converting {num_scans} scan indices to 1/K0...")

    # Convert using the Bruker SDK
    im_lookup = np.array(dataset.scan_to_inverse_mobility(1, scans.tolist()), dtype=np.float64)

    if verbose:
        print(f"  1/K0 range: {im_lookup.min():.4f} - {im_lookup.max():.4f}")
        print(f"  Lookup table size: {im_lookup.nbytes / 1024:.1f} KB")

    return im_lookup


def get_calibration_path(d_path: str) -> Path:
    """Get the standard calibration file path for a .d folder.

    Convention: {d_folder}.im_calibration.npy placed alongside the .d folder

    Args:
        d_path: Path to the .d folder

    Returns:
        Path to the calibration file
    """
    d_path = Path(d_path)
    return d_path.parent / f"{d_path.stem}.im_calibration.npy"


def extract_and_save_calibration(
    d_path: str,
    output_path: Optional[str] = None,
    verbose: bool = False,
) -> Path:
    """Extract calibration and save to .npy file.

    Args:
        d_path: Path to the .d folder
        output_path: Custom output path (default: alongside .d folder)
        verbose: Print progress information

    Returns:
        Path to the saved calibration file
    """
    if output_path is None:
        output_path = get_calibration_path(d_path)
    else:
        output_path = Path(output_path)

    # Extract calibration
    im_lookup = extract_im_calibration(d_path, verbose=verbose)

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(output_path, im_lookup)

    if verbose:
        print(f"  Saved calibration to: {output_path}")

    return output_path


def load_calibration(d_path: str) -> Optional[np.ndarray]:
    """Load cached calibration if available.

    Args:
        d_path: Path to the .d folder

    Returns:
        Calibration array if cached, None otherwise
    """
    cal_path = get_calibration_path(d_path)
    if cal_path.exists():
        return np.load(cal_path)
    return None


def ensure_calibration(d_path: str, verbose: bool = False) -> np.ndarray:
    """Ensure calibration exists, extracting if necessary.

    Args:
        d_path: Path to the .d folder
        verbose: Print progress information

    Returns:
        Calibration lookup array
    """
    cal_path = get_calibration_path(d_path)

    if cal_path.exists():
        if verbose:
            print(f"Loading cached calibration: {cal_path}")
        return np.load(cal_path)
    else:
        if verbose:
            print(f"Extracting calibration (this is slow, but only done once)...")
        extract_and_save_calibration(d_path, verbose=verbose)
        return np.load(cal_path)


def verify_calibration(d_path: str, im_lookup: np.ndarray) -> Tuple[float, float]:
    """Verify calibration accuracy by comparing with SDK.

    Args:
        d_path: Path to the .d folder
        im_lookup: Pre-computed calibration lookup table

    Returns:
        Tuple of (max_absolute_error, max_relative_error)
    """
    from imspy_core.timstof import TimsDatasetDDA

    # Open with SDK for reference
    dataset = TimsDatasetDDA(str(d_path), in_memory=False, use_bruker_sdk=True)

    # Test random scan values
    test_scans = np.random.choice(len(im_lookup), size=min(100, len(im_lookup)), replace=False)

    sdk_values = np.array(dataset.scan_to_inverse_mobility(1, test_scans.tolist()))
    lookup_values = im_lookup[test_scans]

    abs_errors = np.abs(sdk_values - lookup_values)
    rel_errors = abs_errors / np.maximum(sdk_values, 1e-10)

    return float(abs_errors.max()), float(rel_errors.max())


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Extract ion mobility calibration lookup table from timsTOF data"
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Path to .d folder"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        help="Output .npy file path (default: {input}.im_calibration.npy)"
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Verify calibration accuracy after extraction"
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress progress output"
    )

    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: Input path does not exist: {args.input}")
        return 1

    verbose = not args.quiet

    # Extract calibration
    cal_path = extract_and_save_calibration(
        str(args.input),
        str(args.output) if args.output else None,
        verbose=verbose,
    )

    # Optionally verify
    if args.verify:
        if verbose:
            print("\nVerifying calibration accuracy...")
        im_lookup = np.load(cal_path)
        max_abs, max_rel = verify_calibration(str(args.input), im_lookup)
        if verbose:
            print(f"  Max absolute error: {max_abs:.2e} 1/K0")
            print(f"  Max relative error: {max_rel:.2e}")
            if max_rel < 1e-10:
                print("  Calibration is exact (as expected)")
            else:
                print("  WARNING: Calibration has unexpected errors!")

    if verbose:
        print(f"\nCalibration file: {cal_path}")

    return 0


if __name__ == "__main__":
    exit(main())
