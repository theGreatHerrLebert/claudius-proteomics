#!/usr/bin/env python3
"""
Test Ion Mobility Calibration Accuracy

Validates that the pre-computed calibration lookup table produces identical
results to the Bruker SDK for scan→1/K0 conversion.

Usage:
    python scripts/test_calibration_accuracy.py data/raw/PXD019086/sample.d
    python scripts/test_calibration_accuracy.py data/raw/PXD019086/sample.d --verbose
"""

import argparse
from pathlib import Path
from typing import Tuple, Optional
import numpy as np


def test_calibration_accuracy(
    d_path: str,
    calibration_path: Optional[str] = None,
    n_test_frames: int = 10,
    n_test_scans: int = 100,
    verbose: bool = False,
) -> Tuple[bool, dict]:
    """Test calibration accuracy by comparing lookup vs SDK.

    Args:
        d_path: Path to the .d folder
        calibration_path: Path to calibration .npy file (auto-detect if None)
        n_test_frames: Number of frames to test
        n_test_scans: Number of scan values to test per frame
        verbose: Print detailed output

    Returns:
        Tuple of (passed, results_dict)
    """
    from imspy_core.timstof import TimsDatasetDDA

    d_path = Path(d_path)

    # Load or extract calibration
    if calibration_path is None:
        from extract_calibration import get_calibration_path, ensure_calibration
        calibration_path = get_calibration_path(str(d_path))
        if not calibration_path.exists():
            if verbose:
                print(f"Extracting calibration...")
            ensure_calibration(str(d_path), verbose=verbose)

    if verbose:
        print(f"\nTesting calibration accuracy")
        print(f"  Dataset: {d_path}")
        print(f"  Calibration: {calibration_path}")

    # Load calibration lookup table
    im_lookup = np.load(calibration_path)
    if len(im_lookup) == 0:
        print("  ERROR: Empty calibration file (SDK was disabled)")
        return False, {"error": "empty_calibration"}

    if verbose:
        print(f"  Lookup table size: {len(im_lookup)} scans")
        print(f"  1/K0 range: {im_lookup.min():.4f} - {im_lookup.max():.4f}")

    # Load dataset with Bruker SDK for reference
    if verbose:
        print(f"\n  Loading dataset with Bruker SDK for reference...")
    sdk_dataset = TimsDatasetDDA(str(d_path), in_memory=False, use_bruker_sdk=True)

    # Also load with calibration for comparison
    from imspy_connector import py_dda
    cal_dataset = py_dda.PyTimsDatasetDDA.with_calibration(
        str(d_path), False, im_lookup.tolist()
    )

    # Test across multiple frames and scan values
    all_errors = []
    all_rel_errors = []
    test_count = 0

    # Get frame IDs to test
    frame_ids = sdk_dataset.precursor_frames[:n_test_frames]
    if verbose:
        print(f"\n  Testing {len(frame_ids)} frames...")

    for frame_id in frame_ids:
        # Generate random scan values within the valid range
        max_scan = min(len(im_lookup) - 1, 1000)
        test_scans = np.random.choice(max_scan, size=min(n_test_scans, max_scan), replace=False)
        test_scans = np.sort(test_scans)  # Keep as numpy array

        # Get SDK values (ground truth)
        sdk_values = np.array(sdk_dataset.scan_to_inverse_mobility(int(frame_id), test_scans.tolist()))

        # Get lookup values
        lookup_values = np.array(cal_dataset.scan_to_inverse_mobility(int(frame_id), test_scans.tolist()))

        # Calculate errors
        abs_errors = np.abs(sdk_values - lookup_values)
        rel_errors = abs_errors / np.maximum(sdk_values, 1e-10)

        all_errors.extend(abs_errors)
        all_rel_errors.extend(rel_errors)
        test_count += len(test_scans)

        if verbose and len(frame_ids) <= 5:
            print(f"\n    Frame {frame_id}:")
            print(f"      Max absolute error: {abs_errors.max():.2e}")
            print(f"      Max relative error: {rel_errors.max():.2e}")

    # Calculate summary statistics
    all_errors = np.array(all_errors)
    all_rel_errors = np.array(all_rel_errors)

    results = {
        "n_tests": test_count,
        "n_frames": len(frame_ids),
        "max_abs_error": float(all_errors.max()),
        "mean_abs_error": float(all_errors.mean()),
        "max_rel_error": float(all_rel_errors.max()),
        "mean_rel_error": float(all_rel_errors.mean()),
    }

    # Check if within tolerance
    # The lookup should be exact (< 1e-10 relative error)
    passed = results["max_rel_error"] < 1e-8

    if verbose:
        print(f"\n  Results:")
        print(f"    Total tests: {test_count:,}")
        print(f"    Max absolute error: {results['max_abs_error']:.2e} 1/K0")
        print(f"    Mean absolute error: {results['mean_abs_error']:.2e} 1/K0")
        print(f"    Max relative error: {results['max_rel_error']:.2e}")
        print(f"    Mean relative error: {results['mean_rel_error']:.2e}")
        print(f"\n  {'PASSED' if passed else 'FAILED'}: Calibration is {'exact' if passed else 'NOT accurate'}")

    return passed, results


def test_linear_vs_calibration(
    d_path: str,
    verbose: bool = False,
) -> Tuple[float, float]:
    """Compare linear interpolation vs calibrated values.

    This shows the error introduced by using linear interpolation
    instead of the proper Bruker SDK calibration.

    Args:
        d_path: Path to the .d folder
        verbose: Print detailed output

    Returns:
        Tuple of (max_error, mean_error) in 1/K0 units
    """
    from imspy_core.timstof import TimsDatasetDDA

    d_path = Path(d_path)

    if verbose:
        print(f"\nComparing linear interpolation vs SDK calibration")
        print(f"  Dataset: {d_path}")

    # Load with SDK (accurate)
    sdk_dataset = TimsDatasetDDA(str(d_path), in_memory=False, use_bruker_sdk=True)

    # Load with linear interpolation (fast but potentially inaccurate)
    linear_dataset = TimsDatasetDDA(str(d_path), in_memory=False, use_bruker_sdk=False)

    # Test across scan range
    meta = sdk_dataset.meta_data
    num_scans_col = 'NumScans' if 'NumScans' in meta.columns else 'num_scans'
    max_scan = min(int(meta[num_scans_col].max()), 1000)
    test_scans = np.arange(0, max_scan, max(1, max_scan // 100))

    if verbose:
        print(f"  Testing {len(test_scans)} scan values...")

    # Get values from both methods
    frame_id = int(sdk_dataset.precursor_frames[0])
    sdk_values = np.array(sdk_dataset.scan_to_inverse_mobility(frame_id, test_scans.tolist()))
    linear_values = np.array(linear_dataset.scan_to_inverse_mobility(frame_id, test_scans.tolist()))

    # Calculate errors
    abs_errors = np.abs(sdk_values - linear_values)

    max_error = float(abs_errors.max())
    mean_error = float(abs_errors.mean())

    if verbose:
        print(f"\n  Linear interpolation errors:")
        print(f"    Max error:  {max_error:.4f} 1/K0")
        print(f"    Mean error: {mean_error:.4f} 1/K0")
        print(f"\n  For reference, typical 1/K0 range is ~0.6-1.6")
        print(f"  A {max_error:.4f} error at 1/K0=1.0 is {max_error*100:.2f}% relative error")

    return max_error, mean_error


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Test ion mobility calibration accuracy"
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Path to .d folder"
    )
    parser.add_argument(
        "--calibration", "-c",
        type=Path,
        help="Path to calibration .npy file (auto-detect if not specified)"
    )
    parser.add_argument(
        "--n-frames",
        type=int,
        default=10,
        help="Number of frames to test (default: 10)"
    )
    parser.add_argument(
        "--n-scans",
        type=int,
        default=100,
        help="Number of scan values to test per frame (default: 100)"
    )
    parser.add_argument(
        "--compare-linear",
        action="store_true",
        help="Also compare linear interpolation vs SDK"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )

    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: Input path does not exist: {args.input}")
        return 1

    # Run calibration accuracy test
    passed, results = test_calibration_accuracy(
        str(args.input),
        str(args.calibration) if args.calibration else None,
        n_test_frames=args.n_frames,
        n_test_scans=args.n_scans,
        verbose=args.verbose,
    )

    # Optionally compare linear interpolation
    if args.compare_linear:
        test_linear_vs_calibration(str(args.input), verbose=args.verbose)

    return 0 if passed else 1


if __name__ == "__main__":
    exit(main())
