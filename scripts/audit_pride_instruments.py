#!/usr/bin/env python3
"""
Audit PRIDE metadata for the instrument-mismatch pattern found in PXD073076.

The bug: pride_metadata.py records `instrument` as `instruments[0].name` from
the PRIDE API. When a deposit lists multiple instruments and the first one is
not the one we actually processed, the recorded label is wrong even though the
.d-folder search ran correctly. PXD073076 is the canonical example —
instruments=[Orbitrap Exploris 480, timsTOF HT], recorded=Orbitrap, processed
data is timsTOF (ion_mobility populated, .d filenames prefixed HT_).

This script flags two patterns across all cached PRIDE metadata:

  MULTI_NON_TIMS    len(instruments) > 1 AND recorded instrument is not timsTOF
                    (= the PXD073076 pattern; almost certainly wrong if .d was
                    processed)

  PASEF_NON_TIMS    acquisition_mode == PASEF AND recorded instrument is not
                    timsTOF (catches single-instrument deposits where someone
                    mislabeled the API entry — PASEF is Bruker-only)

It does not touch .d files or the raw API; it reads only the metadata cache
(raw_api_response.json + pride_metadata.yaml). Safe to run on a login node.

Stdlib only and 3.6-compatible so the cluster's login-node /usr/bin/python3
(Python 3.6.8) can run it without `module load`. For scripts that need a
modern Python, `module load lang/Python/3.12.3-GCCcore-13.3.0` works on the
login node too.

Usage:
  # On the cluster, scanning the full campaign metadata cache:
  python scripts/audit_pride_instruments.py /lustre/project/ki-proanagi/dateschn/data/metadata

  # Locally against the synced collection:
  python scripts/audit_pride_instruments.py /scratch/claudius-proteomics --layout collection

  # TSV output for diffing / piping:
  python scripts/audit_pride_instruments.py data/metadata --tsv > audit.tsv
"""
import argparse
import json
import re
import sys
from pathlib import Path
from typing import List, Optional


def find_metadata_dirs(root: Path, layout: str) -> List[Path]:
    """
    Return per-dataset directories that contain a metadata cache.

    layout=campaign:   root is `data/metadata/`, children are accession dirs.
    layout=collection: root is `/scratch/claudius-proteomics/` style, children
                       hold metadata in `<acc>/metadata/`.
    layout=auto:       try campaign first, fall back to collection.
    """
    if layout in ("campaign", "auto"):
        campaign = [p for p in sorted(root.iterdir())
                    if p.is_dir() and (p / "raw_api_response.json").exists()]
        if campaign or layout == "campaign":
            return campaign
    return [p / "metadata" for p in sorted(root.iterdir())
            if p.is_dir() and (p / "metadata" / "raw_api_response.json").exists()]


_YAML_FIELD_RE = re.compile(
    r"^\s{2}([a-z_]+):\s*\n(?:\s{4}.*\n)*?\s{4}value:\s*(.*)$",
    re.MULTILINE,
)


def extract_yaml_field(yaml_text: str, name: str) -> Optional[str]:
    """
    Pull `fields.<name>.value` out of a pride_metadata.yaml without a YAML lib.
    The schema is fixed (pride_metadata.py writes it), so a regex is enough and
    keeps this script stdlib-only.
    """
    for fname, value in _YAML_FIELD_RE.findall(yaml_text):
        if fname == name:
            v = value.strip()
            if v in ("null", "~", ""):
                return None
            return v.strip('"').strip("'")
    return None


def audit_one(meta_dir: Path) -> Optional[dict]:
    raw = meta_dir / "raw_api_response.json"
    if not raw.exists():
        return None
    try:
        api = json.load(open(raw))
    except json.JSONDecodeError as e:
        return {"accession": meta_dir.parent.name if meta_dir.name == "metadata" else meta_dir.name,
                "error": f"raw_api_response.json: {e}"}

    accession = api.get("accession") or meta_dir.name
    instruments = [i.get("name") for i in (api.get("instruments") or [])]

    recorded_instr = None
    recorded_mode = None
    pmeta = meta_dir / "pride_metadata.yaml"
    if pmeta.exists():
        text = pmeta.read_text()
        recorded_instr = extract_yaml_field(text, "instrument")
        recorded_mode = extract_yaml_field(text, "acquisition_mode")

    flags = []
    rec_low = (recorded_instr or "").lower()
    has_tims = "tims" in rec_low
    if len(instruments) > 1 and not has_tims:
        flags.append("MULTI_NON_TIMS")
    if recorded_mode and recorded_mode.upper() == "PASEF" and not has_tims and recorded_instr:
        flags.append("PASEF_NON_TIMS")

    return {
        "accession": accession,
        "n_instruments": len(instruments),
        "instruments_all": instruments,
        "recorded_instrument": recorded_instr,
        "acquisition_mode": recorded_mode,
        "flags": flags,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path,
                    help="Metadata root (data/metadata) or collection root")
    ap.add_argument("--layout", choices=["auto", "campaign", "collection"],
                    default="auto",
                    help="campaign: <root>/<acc>/raw_api_response.json (cluster). "
                         "collection: <root>/<acc>/metadata/raw_api_response.json "
                         "(/scratch sync). auto: try campaign, then collection.")
    ap.add_argument("--tsv", action="store_true",
                    help="Emit a TSV table instead of the human report (one row "
                         "per dataset, all instruments joined with '|').")
    ap.add_argument("--flagged-only", action="store_true",
                    help="In TSV mode, omit datasets with no flag.")
    args = ap.parse_args()

    if not args.root.exists():
        print(f"error: {args.root} does not exist", file=sys.stderr)
        return 2

    dirs = find_metadata_dirs(args.root, args.layout)
    if not dirs:
        print(f"error: no metadata caches found under {args.root} "
              f"(layout={args.layout})", file=sys.stderr)
        return 2

    rows = [r for r in (audit_one(d) for d in dirs) if r is not None]
    flagged = [r for r in rows if r.get("flags")]

    if args.tsv:
        print("accession\tn_instr\trecorded_instrument\tacquisition_mode\tflags\tinstruments_all")
        for r in rows:
            if args.flagged_only and not r.get("flags"):
                continue
            if "error" in r:
                print(f"{r['accession']}\t?\t?\t?\tERROR\t{r['error']}")
                continue
            print("\t".join([
                r["accession"],
                str(r["n_instruments"]),
                r["recorded_instrument"] or "",
                r["acquisition_mode"] or "",
                ",".join(r["flags"]) or "ok",
                "|".join(i or "" for i in r["instruments_all"]),
            ]))
        return 0 if not flagged else 1

    print(f"Scanned {len(rows)} datasets under {args.root} (layout={args.layout}).")
    multi = [r for r in rows if r.get("n_instruments", 0) > 1]
    print(f"  multi-instrument deposits: {len(multi)}")
    print(f"  flagged:                   {len(flagged)}")
    print()

    if flagged:
        print("== FLAGGED ==")
        for r in flagged:
            print(f"\n{r['accession']}  [{', '.join(r['flags'])}]")
            print(f"  recorded:        {r['recorded_instrument']!r}")
            print(f"  acquisition:     {r['acquisition_mode']!r}")
            print(f"  PRIDE instruments[]: {r['instruments_all']}")

    multi_not_flagged = [r for r in multi if not r.get("flags")]
    if multi_not_flagged:
        print("\n== Multi-instrument deposits where instruments[0] is timsTOF ==")
        print("(probably fine, but record which timsTOF variant we actually used)")
        for r in multi_not_flagged:
            print(f"  {r['accession']}: recorded={r['recorded_instrument']!r} "
                  f"all={r['instruments_all']}")

    return 0 if not flagged else 1


if __name__ == "__main__":
    sys.exit(main())
