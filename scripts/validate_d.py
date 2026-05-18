#!/usr/bin/env python3
"""
Pre-flight validation of Bruker timsTOF .d folders.

A .d folder is valid when analysis.tdf is a readable SQLite database with a
non-empty Frames table and analysis.tdf_bin exists and is not truncated (its
size exceeds the largest per-frame offset). This catches incomplete / corrupt
downloads — e.g. a truncated .tdf_bin — which otherwise abort a whole dataset's
search: MSFragger hard-fails on the first unreadable file.

quarantine_bad_d() moves bad .d folders into raw_dir/_quarantine/ (with a
quarantine.json record) so the search engines only ever see valid files.

Usage:
  python scripts/validate_d.py data/raw/PXD046777               # report only
  python scripts/validate_d.py --quarantine data/raw/PXD046777  # move bad ones
"""
import argparse
import json
import shutil
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path


def validate_d_folder(d_path: Path) -> tuple:
    """Return (ok: bool, reason: str) for a Bruker .d folder."""
    d_path = Path(d_path)
    tdf = d_path / "analysis.tdf"
    tdf_bin = d_path / "analysis.tdf_bin"
    if not tdf.exists():
        return False, "missing analysis.tdf"
    if not tdf_bin.exists():
        return False, "missing analysis.tdf_bin"

    # Core: analysis.tdf must be a readable SQLite db with frames.
    try:
        con = sqlite3.connect(f"file:{tdf}?mode=ro", uri=True)
        try:
            n_frames = con.execute("SELECT COUNT(*) FROM Frames").fetchone()[0]
        finally:
            con.close()
    except Exception as e:
        return False, f"analysis.tdf unreadable: {str(e)[:120]}"
    if not n_frames:
        return False, "Frames table is empty"

    bin_size = tdf_bin.stat().st_size
    if bin_size == 0:
        return False, "analysis.tdf_bin is empty"

    # Best-effort truncation check: the .tdf_bin must extend past the last
    # frame's offset. Skipped silently if the schema lacks TimsId.
    try:
        con = sqlite3.connect(f"file:{tdf}?mode=ro", uri=True)
        try:
            max_off = con.execute("SELECT MAX(TimsId) FROM Frames").fetchone()[0]
        finally:
            con.close()
        if max_off is not None and bin_size <= max_off:
            return False, (f"analysis.tdf_bin truncated "
                           f"({bin_size} bytes <= last frame offset {max_off})")
    except Exception:
        pass

    return True, "ok"


def quarantine_bad_d(raw_dir: Path) -> dict:
    """Validate every *.d in raw_dir; move bad ones into _quarantine/.

    Returns {checked, valid, quarantined: [{name, reason}]}.
    """
    raw_dir = Path(raw_dir)
    result = {"checked": 0, "valid": 0, "quarantined": []}
    if not raw_dir.is_dir():
        return result

    for d in sorted(raw_dir.glob("*.d")):
        if not d.is_dir():
            continue
        result["checked"] += 1
        ok, reason = validate_d_folder(d)
        if ok:
            result["valid"] += 1
            continue
        qdir = raw_dir / "_quarantine"
        qdir.mkdir(exist_ok=True)
        dest = qdir / d.name
        if dest.exists():
            shutil.rmtree(dest, ignore_errors=True)
        shutil.move(str(d), str(dest))
        result["quarantined"].append({"name": d.name, "reason": reason})

    if result["quarantined"]:
        rec = {
            "raw_dir": str(raw_dir),
            "quarantined_at_utc": datetime.now(timezone.utc)
            .strftime("%Y-%m-%dT%H:%M:%SZ"),
            "checked": result["checked"],
            "valid": result["valid"],
            "quarantined": result["quarantined"],
        }
        (raw_dir / "_quarantine" / "quarantine.json").write_text(
            json.dumps(rec, indent=2))
    return result


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("raw_dir", type=Path, help="data/raw/{accession} directory")
    ap.add_argument("--quarantine", action="store_true",
                    help="move bad .d folders into _quarantine/ (default: report only)")
    args = ap.parse_args()

    if args.quarantine:
        r = quarantine_bad_d(args.raw_dir)
        for q in r["quarantined"]:
            print(f"  QUARANTINED {q['name']}: {q['reason']}")
        print(f"{r['valid']}/{r['checked']} valid, {len(r['quarantined'])} quarantined")
    else:
        checked = valid = 0
        for d in sorted(Path(args.raw_dir).glob("*.d")):
            if not d.is_dir():
                continue
            checked += 1
            ok, reason = validate_d_folder(d)
            valid += int(ok)
            if not ok:
                print(f"  BAD  {d.name}: {reason}")
        print(f"{valid}/{checked} valid")
    return 0


if __name__ == "__main__":
    sys.exit(main())
