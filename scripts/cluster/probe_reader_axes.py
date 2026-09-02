#!/usr/bin/env python3
"""probe_reader_axes.py — dump a timsTOF file's TOF->m/z and scan->1/K0 curves.

Used to A/B two rustims builds against each other. The reader's two calibrated
axes are the whole product here: fragment m/z feeds the b/y matching and 1/K0
feeds every CCS label, so a reader upgrade has to be measured, not assumed.

The conversions are pure functions of the file's stored calibration, so we do
not need real peaks — a fixed grid of TOF indices and scan numbers, evaluated on
frames spread across the run, pins the curve exactly. Frames are sampled across
the run because the calibration is resolved per frame (digitizer temperature
drifts, ~20 ppm/degC), so a one-frame probe would hide exactly that.

    python probe_reader_axes.py <dataset.d> -o out.npz [--frames 9]

Writes an .npz with the two curves plus the build fingerprint, for
diff_reader_axes.py to compare.
"""
import argparse
import json
import sqlite3
import sys
from pathlib import Path

import numpy as np

TOF_GRID = np.arange(1, 350_001, 2_500, dtype=np.int32)
SCAN_STEP = 5


def _build_fingerprint() -> dict:
    """Whatever identifies the reader that produced these numbers."""
    info = {}
    try:
        import importlib.metadata as md
        for pkg in ("imspy_connector", "imspy_core"):
            try:
                info[pkg] = md.version(pkg)
            except Exception:
                info[pkg] = "?"
    except Exception:
        pass
    try:
        import imspy_connector
        info["connector_file"] = getattr(imspy_connector, "__file__", "?")
    except Exception:
        pass
    return info


def _calibration_model_types(d_path: Path) -> dict:
    """ModelType of both calibrations, straight from analysis.tdf.

    ModelType 2 on the m/z axis is the case the C3/C4-duplicate fix touches.
    """
    tdf = d_path / "analysis.tdf"
    out = {"mz_model_type": None, "tims_model_type": None}
    if not tdf.exists():
        return out
    try:
        con = sqlite3.connect(f"file:{tdf}?mode=ro", uri=True)
        for key, table in (("mz_model_type", "MzCalibration"),
                           ("tims_model_type", "TimsCalibration")):
            rows = con.execute(f"SELECT DISTINCT ModelType FROM {table}").fetchall()
            out[key] = sorted(r[0] for r in rows)
        con.close()
    except Exception as e:
        out["error"] = str(e)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset", type=Path, help="path to a .d folder")
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--frames", type=int, default=9,
                    help="how many frames to sample across the run (default 9)")
    ap.add_argument("--max-scan", type=int, default=900)
    args = ap.parse_args()

    from imspy_core.timstof import TimsDatasetDDA

    # SDK-free is how the runner opens files on this cluster, so probe that path.
    ds = TimsDatasetDDA(str(args.dataset), use_bruker_sdk=False)

    n_frames = int(getattr(ds, "frame_count", 0) or 0)
    if n_frames <= 0:  # older builds expose it differently
        n_frames = len(getattr(ds, "meta_data", []))
    if n_frames <= 0:
        print("could not determine frame count", file=sys.stderr)
        return 1

    frame_ids = np.unique(
        np.linspace(1, n_frames, num=min(args.frames, n_frames)).astype(np.int32)
    )
    scan_grid = np.arange(0, args.max_scan + 1, SCAN_STEP, dtype=np.int32)

    mz = np.zeros((len(frame_ids), len(TOF_GRID)), dtype=np.float64)
    im = np.zeros((len(frame_ids), len(scan_grid)), dtype=np.float64)
    for i, fid in enumerate(frame_ids):
        mz[i] = np.asarray(ds.tof_to_mz(int(fid), TOF_GRID), dtype=np.float64)
        im[i] = np.asarray(ds.scan_to_inverse_mobility(int(fid), scan_grid), dtype=np.float64)

    meta = {
        "dataset": str(args.dataset),
        "n_frames": n_frames,
        "build": _build_fingerprint(),
        "calibration": _calibration_model_types(args.dataset),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out, frame_ids=frame_ids, tof_grid=TOF_GRID, scan_grid=scan_grid,
             mz=mz, inv_mobility=im, meta=json.dumps(meta))
    print(json.dumps(meta, indent=2))
    print(f"-> {args.out}  ({len(frame_ids)} frames x {len(TOF_GRID)} tof, "
          f"{len(scan_grid)} scans)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
