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

COSIM_MIN = 0.5   # isotope_cosim_median below this == catastrophic m/z failure
HQ_MIN = 5.0      # pct_high_quality below this == degraded extraction.
                  # Needed because a boundary-m/z dataset can still show a
                  # misleadingly-moderate isotope_cosim (~0.8) computed from
                  # tiny, wrong-shaped envelopes — pct_high_quality (0-2% when
                  # affected, 12-28% when healthy) is what actually catches it.
HALT_AFTER = 3    # consecutive recent failures that halt a batch


def _cosim(manifest):
    return manifest.get("quality_summary", {}).get("isotope_cosim_median")


def _hq(manifest):
    return manifest.get("quality_summary", {}).get("pct_high_quality")


def _verdict(manifest):
    """Return (passes: bool, summary: str)."""
    c, h = _cosim(manifest), _hq(manifest)
    cok = c is not None and c >= COSIM_MIN
    hok = h is not None and h >= HQ_MIN
    fails = []
    if not cok:
        fails.append(f"isotope_cosim_median={c} < {COSIM_MIN}")
    if not hok:
        fails.append(f"pct_high_quality={h} < {HQ_MIN}")
    summary = (f"isotope_cosim_median={c} pct_high_quality={h}"
               + ("" if (cok and hok) else "  [" + "; ".join(fails) + "]"))
    return (cok and hok), summary


def _passes(manifest):
    return _verdict(manifest)[0]


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
            ok, summary = _verdict(m)
            verdicts.append(ok)
            print(f"  {acc}: {summary}  {'PASS' if ok else 'FAIL'}")
        if all(not ok for ok in verdicts):
            print(f"qc_gate: HALT — last {HALT_AFTER} datasets all failed QC "
                  f"(isotope_cosim_median < {COSIM_MIN} or pct_high_quality < "
                  f"{HQ_MIN}); likely a systemic bug")
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
    ok, summary = _verdict(m)
    print(f"qc_gate: {args.accession} {summary}  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
