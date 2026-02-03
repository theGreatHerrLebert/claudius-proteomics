#!/usr/bin/env python3
"""
Extract ion mobility calibration lookup table from timsTOF .d files.

The Bruker SDK provides accurate scan→1/K0 conversion but is:
1. Not thread-safe (single-threaded only)
2. Requires Bruker libraries at runtime

This script extracts the calibration once using the SDK, saving it as a numpy
array. The calibration can then be loaded via PyTimsDatasetDDA.with_calibration()
for fast, thread-safe parallel extraction with accurate mobility values.

Usage:
    python extract_calibration.py /path/to/data.d /path/to/output.npy

The output is a 1D numpy array where index=scan_number, value=1/K0.
"""

import sys
import argparse
from pathlib import Path

import numpy as np


def extract_calibration(data_path: str, use_sdk: bool = True) -> np.ndarray:
    """
    Extract ion mobility calibration from a timsTOF .d file.

    Args:
        data_path: Path to the .d folder
        use_sdk: If True, use Bruker SDK for accurate calibration.
                 If False, returns empty array (linear interpolation fallback).

    Returns:
        1D numpy array mapping scan_number → 1/K0 value
    """
    from imspy_core.timstof.dda import TimsDatasetDDA

    # Open dataset with Bruker SDK for accurate calibration
    dataset = TimsDatasetDDA(data_path, in_memory=False, use_bruker_sdk=use_sdk)

    # Get max scan number from metadata
    # The scan range is typically 0 to ~1000 for timsTOF
    max_scan = int(dataset.meta_data['NumScans'].max())

    # Build lookup table: scan → 1/K0
    # Use frame 1 (first MS1 frame) as reference - calibration is constant across frames
    frame_id = 1
    scan_values = list(range(max_scan + 1))

    # Get 1/K0 values from SDK
    try:
        # Access the internal Rust dataset for calibration
        im_values = dataset._TimsDataset__dataset.scan_to_inverse_mobility(frame_id, scan_values)
        calibration = np.array(im_values, dtype=np.float64)
    except Exception as e:
        print(f"Warning: Could not extract calibration via SDK: {e}")
        print("Returning empty calibration (will use linear interpolation)")
        calibration = np.array([])

    return calibration


def extract_and_save_calibration(data_path: str, output_path: str, verbose: bool = True) -> None:
    """
    Extract calibration and save to numpy file.

    Args:
        data_path: Path to the .d folder
        output_path: Path for output .npy file
        verbose: Print progress messages
    """
    if verbose:
        print(f"Extracting IM calibration from: {data_path}")

    calibration = extract_calibration(data_path, use_sdk=True)

    if len(calibration) > 0:
        np.save(output_path, calibration)
        if verbose:
            print(f"Saved calibration ({len(calibration)} scans) to: {output_path}")
            print(f"  1/K0 range: {calibration.min():.4f} - {calibration.max():.4f}")
    else:
        # Save empty array as marker
        np.save(output_path, np.array([]))
        if verbose:
            print("Saved empty calibration (will use linear interpolation)")


def main():
    parser = argparse.ArgumentParser(
        description="Extract ion mobility calibration from timsTOF .d files"
    )
    parser.add_argument("data_path", help="Path to the .d folder")
    parser.add_argument("output_path", help="Path for output .npy file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    if not Path(args.data_path).exists():
        print(f"Error: Data path does not exist: {args.data_path}")
        sys.exit(1)

    extract_and_save_calibration(args.data_path, args.output_path, verbose=args.verbose)


if __name__ == "__main__":
    main()
