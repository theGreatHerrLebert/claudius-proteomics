#!/usr/bin/env python3
"""
Derive the actually-processed instrument from Bruker analysis.tdf
GlobalMetadata, and diff against the currently-recorded value in
pride_metadata.yaml. Read-only — produces a report, no edits.

Follow-up to audit_pride_instruments.py. The audit flags deposits where
PRIDE's instruments[0] is unlikely to be what we processed; this script
answers definitively by reading the truth from each .d folder's
analysis.tdf SQLite database, exactly as scripts/raw_metadata.py already
does for full-metadata extraction (key=InstrumentName in GlobalMetadata).

Why not just re-run extract_raw_metadata.py? That script does the same
read and overwrites pride_metadata.yaml with the merged value (raw wins
over PRIDE for instrument — see scripts/extract_raw_metadata.py:60-71).
Running it would already fix the bug. This script gives the diff first
so you can review the proposed changes before any file is touched —
especially for the MALDI-fleX cluster and the mega-deposits where the
right call may need a human.

Stdlib only and 3.6-compatible — runs on the cluster login node without
a `module load`. See [[mogon-nhr-access]] for the 3.6.8 quirk.

Usage:
  # Diff for a specific set of accessions:
  python scripts/derive_instrument_from_tdf.py \\
      --raw-root /lustre/project/ki-proanagi/dateschn/data/raw \\
      --meta-root /lustre/project/ki-proanagi/dateschn/data/metadata \\
      PXD073076 PXD050530 PXD058376

  # Pipe accessions from the audit (TSV first column):
  python scripts/audit_pride_instruments.py \\
      /lustre/project/ki-proanagi/dateschn/data/metadata \\
      --tsv --flagged-only \\
    | tail -n +2 | cut -f1 \\
    | python scripts/derive_instrument_from_tdf.py \\
        --raw-root /lustre/project/ki-proanagi/dateschn/data/raw \\
        --meta-root /lustre/project/ki-proanagi/dateschn/data/metadata
"""
import argparse
import re
import sqlite3
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


_YAML_FIELD_RE = re.compile(
    r"^\s{2}([a-z_]+):\s*\n(?:\s{4}.*\n)*?\s{4}value:\s*(.*)$",
    re.MULTILINE,
)


def extract_yaml_field(yaml_text: str, name: str) -> Optional[str]:
    """Lift `fields.<name>.value` out of pride_metadata.yaml without pyyaml."""
    for fname, value in _YAML_FIELD_RE.findall(yaml_text):
        if fname == name:
            v = value.strip()
            if v in ("null", "~", ""):
                return None
            return v.strip('"').strip("'")
    return None


def read_global_metadata(tdf_path: Path) -> Dict[str, str]:
    """Return a dict of Key→Value from analysis.tdf's GlobalMetadata table.

    Read-only sqlite open (URI mode) — safe to call against files that
    other jobs may be reading. Mirrors raw_metadata._extract_global_metadata.
    """
    con = sqlite3.connect("file:%s?mode=ro" % tdf_path, uri=True)
    try:
        rows = con.execute("SELECT Key, Value FROM GlobalMetadata").fetchall()
    finally:
        con.close()
    return {k: v for k, v in rows}


def derive_for_accession(acc: str, raw_root: Path) -> List[Tuple[Path, Optional[str], Optional[str]]]:
    """For acc, return [(d_path, instrument_name, error)] across all its .d folders."""
    acc_dir = raw_root / acc
    out = []  # type: List[Tuple[Path, Optional[str], Optional[str]]]
    if not acc_dir.exists():
        return out
    for tdf in sorted(acc_dir.rglob("*.d/analysis.tdf")):
        # Quarantined .d folders are intentionally excluded (see validate_d.py).
        if "_quarantine" in tdf.parts:
            continue
        try:
            gm = read_global_metadata(tdf)
            out.append((tdf.parent, gm.get("InstrumentName"), None))
        except sqlite3.DatabaseError as e:
            out.append((tdf.parent, None, "sqlite: %s" % e))
        except OSError as e:
            out.append((tdf.parent, None, "io: %s" % e))
    return out


def read_accessions_stdin() -> List[str]:
    """Accept either raw accessions (one per line) or the audit TSV first column."""
    seen = []
    for line in sys.stdin:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        first = line.split("\t")[0].strip()
        if first.startswith("PXD") or first.startswith("SIM"):
            if first not in seen:
                seen.append(first)
    return seen


def classify(recorded: Optional[str], distinct: List[str]) -> str:
    if not distinct:
        return "NO_D_FILES"
    if any(d is None for d in distinct):
        return "READ_ERROR"
    if len(set(distinct)) > 1:
        return "MIXED_WITHIN_DEPOSIT"
    derived = distinct[0]
    if not recorded:
        return "NO_RECORDED"
    if recorded == derived:
        return "MATCH"
    rec_low, der_low = recorded.lower(), derived.lower()
    if "tims" in rec_low and "tims" in der_low:
        return "TIMSTOF_VARIANT_DIFF"
    return "MISMATCH"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("accessions", nargs="*",
                    help="Accessions to check. If none given, read from stdin.")
    ap.add_argument("--raw-root", type=Path, required=True,
                    help="Directory containing per-accession .d folders.")
    ap.add_argument("--meta-root", type=Path, required=True,
                    help="Directory containing per-accession metadata caches.")
    ap.add_argument("--tsv", action="store_true",
                    help="Emit TSV (one row per accession) instead of human report.")
    args = ap.parse_args()

    accs = list(args.accessions) or read_accessions_stdin()
    if not accs:
        print("error: no accessions provided (positional or stdin)", file=sys.stderr)
        return 2

    rows = []
    for acc in accs:
        pmeta = args.meta_root / acc / "pride_metadata.yaml"
        recorded = extract_yaml_field(pmeta.read_text(), "instrument") if pmeta.exists() else None
        results = derive_for_accession(acc, args.raw_root)
        names = [name for _, name, _ in results]
        verdict = classify(recorded, names)
        rows.append({
            "acc": acc,
            "recorded": recorded,
            "results": results,
            "verdict": verdict,
        })

    if args.tsv:
        print("accession\trecorded\tderived\tdistinct\tn_d\tverdict")
        for r in rows:
            names = [n for _, n, _ in r["results"]]
            distinct = sorted({n for n in names if n is not None})
            derived = distinct[0] if len(distinct) == 1 else "|".join(distinct)
            print("\t".join([
                r["acc"], r["recorded"] or "", derived,
                "|".join(distinct), str(len(r["results"])), r["verdict"],
            ]))
    else:
        for r in rows:
            names = [n for _, n, _ in r["results"]]
            distinct = sorted({n for n in names if n is not None})
            print()
            print("%s  [%s]" % (r["acc"], r["verdict"]))
            print("  recorded:  %r" % r["recorded"])
            if distinct:
                print("  derived:   %r" % (distinct[0] if len(distinct) == 1 else distinct))
            print("  n .d files: %d" % len(r["results"]))
            errors = [(d, e) for d, _, e in r["results"] if e]
            for d, e in errors[:3]:
                print("  ERR  %s: %s" % (d.name, e))
            if r["verdict"] == "MIXED_WITHIN_DEPOSIT":
                # Show counts so we can see whether one variant dominates.
                from collections import Counter
                counts = Counter(names)
                for name, n in counts.most_common():
                    print("    %d× %r" % (n, name))

    bad = {"MISMATCH", "MIXED_WITHIN_DEPOSIT", "TIMSTOF_VARIANT_DIFF", "READ_ERROR"}
    return 1 if any(r["verdict"] in bad for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
