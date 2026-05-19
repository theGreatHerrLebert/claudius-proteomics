#!/usr/bin/env python3
"""QC gate for the San José runner.

Two modes:

  Single dataset:   qc_gate.py <accession> [--data-dir DIR]
      Exit 0 if the dataset's merged manifest passes QC, 1 if it fails,
      2 if there is no manifest yet.

  Batch halt check: qc_gate.py --batch-check [--data-dir DIR]
      Exit 1 (halt the batch) if the most recent HALT_AFTER completed
      datasets ALL fail QC — the signature of a systemic extraction bug.
      Exit 0 otherwise. Fewer than HALT_AFTER completed datasets never halts.

A dataset fails QC if its manifest's isotope_cosim_median is below COSIM_MIN.
Healthy datasets run ~0.78-0.98; a broken m/z calibration drives it to ~0.
"""
import argparse
import glob
import json
import os
import sys

COSIM_MIN = 0.5   # isotope_cosim_median below this == catastrophic
HALT_AFTER = 3    # consecutive recent failures that halt a batch


def _cosim(manifest):
    return manifest.get("quality_summary", {}).get("isotope_cosim_median")


def _passes(manifest):
    c = _cosim(manifest)
    return c is not None and c >= COSIM_MIN


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("accession", nargs="?", help="dataset accession to check")
    ap.add_argument("--batch-check", action="store_true",
                    help="halt-check the most recent completed datasets")
    ap.add_argument("--data-dir", default="data",
                    help="data root containing merged/ (default: data)")
    args = ap.parse_args()

    if args.batch_check:
        manifests = sorted(
            glob.glob(os.path.join(args.data_dir, "merged", "*", "manifest.json")),
            key=os.path.getmtime,
        )
        recent = manifests[-HALT_AFTER:]
        if len(recent) < HALT_AFTER:
            print(f"qc_gate: only {len(recent)} completed dataset(s) — no halt")
            return 0
        verdicts = []
        for mp in recent:
            with open(mp) as f:
                m = json.load(f)
            acc = m.get("accession", os.path.basename(os.path.dirname(mp)))
            verdicts.append((acc, _cosim(m), _passes(m)))
        for acc, c, ok in verdicts:
            print(f"  {acc}: isotope_cosim_median={c}  {'PASS' if ok else 'FAIL'}")
        if all(not ok for _, _, ok in verdicts):
            print(f"qc_gate: HALT — last {HALT_AFTER} datasets all failed QC "
                  f"(isotope_cosim_median < {COSIM_MIN}); likely a systemic bug")
            return 1
        print("qc_gate: batch OK")
        return 0

    if not args.accession:
        ap.error("provide an accession, or use --batch-check")
    path = os.path.join(args.data_dir, "merged", args.accession, "manifest.json")
    if not os.path.exists(path):
        print(f"qc_gate: no manifest at {path} — cannot evaluate")
        return 2
    with open(path) as f:
        m = json.load(f)
    ok = _passes(m)
    print(f"qc_gate: {args.accession} isotope_cosim_median={_cosim(m)}  "
          f"{'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
