#!/usr/bin/env python3
"""diff_reader_axes.py — compare two probe_reader_axes.py dumps in ppm.

    python diff_reader_axes.py old.npz new.npz

Reports the m/z and 1/K0 axes separately, because the two failure modes have
very different sizes: the ModelType-2 m/z bias is single-digit ppm, while the
linear scan->1/K0 fallback is thousands. Prints per-frame worst cases too, since
a per-frame calibration fix shows up as drift across the run rather than as a
constant offset.
"""
import argparse
import json
import sys

import numpy as np


def _ppm(new: np.ndarray, old: np.ndarray) -> np.ndarray:
    """Relative difference in ppm, guarding the zeros in both grids."""
    denom = np.where(np.abs(old) > 0, np.abs(old), np.nan)
    return (new - old) / denom * 1e6


def _report(name: str, old: np.ndarray, new: np.ndarray, frame_ids: np.ndarray) -> float:
    d = _ppm(new, old)
    finite = d[np.isfinite(d)]
    if finite.size == 0:
        print(f"\n{name}: no comparable values")
        return 0.0
    med, p99 = np.median(np.abs(finite)), np.percentile(np.abs(finite), 99)
    mx = np.max(np.abs(finite))
    print(f"\n{name}  (n={finite.size})")
    print(f"  |delta| ppm : median {med:12.4f}   p99 {p99:12.4f}   max {mx:12.4f}")
    print(f"  signed ppm  : min {np.min(finite):12.4f}   max {np.max(finite):12.4f}")
    per_frame = np.nanmax(np.abs(d), axis=1)
    worst = np.argsort(per_frame)[::-1][:3]
    print("  worst frames: " + ", ".join(
        f"frame {int(frame_ids[i])} {per_frame[i]:.4f} ppm" for i in worst))
    # Drift across the run is the per-frame-calibration signature.
    if len(frame_ids) > 2 and np.isfinite(per_frame).all():
        print(f"  first frame {per_frame[0]:.4f} ppm -> last frame {per_frame[-1]:.4f} ppm")
    return float(mx)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("old", help="npz from the currently deployed build")
    ap.add_argument("new", help="npz from the rebuilt reader")
    ap.add_argument("--threshold-ppm", type=float, default=1.0,
                    help="max |delta| below which the axes are called unchanged")
    args = ap.parse_args()

    a, b = np.load(args.old, allow_pickle=True), np.load(args.new, allow_pickle=True)
    ma, mb = json.loads(str(a["meta"])), json.loads(str(b["meta"]))
    if ma["dataset"] != mb["dataset"]:
        print(f"REFUSING: different datasets\n  old {ma['dataset']}\n  new {mb['dataset']}",
              file=sys.stderr)
        return 2
    if not np.array_equal(a["frame_ids"], b["frame_ids"]):
        print("REFUSING: frame grids differ", file=sys.stderr)
        return 2

    print(f"dataset   : {ma['dataset']}")
    print(f"calibration: {ma['calibration']}")
    print(f"old build : {ma['build']}")
    print(f"new build : {mb['build']}")

    worst_mz = _report("TOF -> m/z", a["mz"], b["mz"], a["frame_ids"])
    worst_im = _report("scan -> 1/K0", a["inv_mobility"], b["inv_mobility"], a["frame_ids"])

    moved = max(worst_mz, worst_im) > args.threshold_ppm
    print("\n" + ("=" * 60))
    print(f"VERDICT: axes {'MOVED' if moved else 'unchanged'} "
          f"(worst m/z {worst_mz:.4f} ppm, worst 1/K0 {worst_im:.4f} ppm, "
          f"threshold {args.threshold_ppm} ppm)")
    if moved:
        print("The merged corpus was extracted with the OLD build — size the")
        print("re-extraction from these numbers before building any more tiers.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
