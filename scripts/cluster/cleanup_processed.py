#!/usr/bin/env python3
"""Safely reclaim space by removing processed/{accession}/ intermediates.

processed/{acc}/ holds the search-engine outputs + precursor_index. They are
only needed until a dataset is merged and QC-approved; after that they can go.
But if a dataset failed QC it may need a step-4 re-extraction, which needs the
search results — so this tool removes processed/{acc}/ ONLY when the dataset is
fully merged AND passes QC. Anything else is kept.

Dry-run by default; pass --apply to actually delete.

    cleanup_processed.py --data-dir /lustre/.../data            # dry-run
    cleanup_processed.py --data-dir /lustre/.../data --apply
"""
import argparse
import json
import os
import shutil
import sys

COSIM_MIN = 0.5   # keep in sync with qc_gate.py
HQ_MIN = 5.0      # keep in sync with qc_gate.py


def _qc_ok(manifest_path):
    """Return (merged, ok, reason)."""
    if not os.path.exists(manifest_path):
        return False, False, "not merged"
    with open(manifest_path) as f:
        m = json.load(f)
    q = m.get("quality_summary", {})
    c = q.get("isotope_cosim_median")
    h = q.get("pct_high_quality")
    if c is None or h is None:
        return True, False, "merged but manifest missing QC fields"
    if c < COSIM_MIN or h < HQ_MIN:
        return True, False, (f"QC fail (isotope_cosim_median={c}, "
                             f"pct_high_quality={h}) — may need re-extraction")
    return True, True, (f"merged + QC pass (isotope_cosim_median={c}, "
                        f"pct_high_quality={h})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default="data",
                    help="data root containing processed/ and merged/ (default: data)")
    ap.add_argument("--apply", action="store_true",
                    help="actually delete (default: dry-run, list only)")
    args = ap.parse_args()

    processed = os.path.join(args.data_dir, "processed")
    if not os.path.isdir(processed):
        print(f"no processed/ directory at {processed}")
        return 0

    removable, kept = [], []
    for acc in sorted(os.listdir(processed)):
        pdir = os.path.join(processed, acc)
        if not os.path.isdir(pdir):
            continue
        manifest = os.path.join(args.data_dir, "merged", acc, "manifest.json")
        merged, ok, reason = _qc_ok(manifest)
        (removable if ok else kept).append((acc, reason))

    for acc, reason in kept:
        print(f"  KEEP          {acc:14s} {reason}")
    for acc, reason in removable:
        verb = "REMOVED      " if args.apply else "WOULD REMOVE "
        print(f"  {verb} {acc:14s} {reason}")
        if args.apply:
            shutil.rmtree(os.path.join(processed, acc))

    tail = "" if args.apply else "   (dry-run — pass --apply to delete)"
    print(f"\n{len(removable)} removable, {len(kept)} kept{tail}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
