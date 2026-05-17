#!/usr/bin/env python3
"""Extract individual .d folders from a remote PRIDE zip via HTTP range requests.

PRIDE distributes some timsTOF datasets (e.g. PXD019086) as huge per-organism
zip bundles (60-160 GB) with no individual-file granularity, and its FTP is
often blocked on HPC. This pulls only the bytes for the selected .d folders
over HTTPS, so a 1-file test costs ~1-3 GB instead of the whole bundle.

Usage:
    python fetch_pride_member.py <zip-url> --out <dir> --max-d 1
"""
import argparse
import sys
from pathlib import Path

from remotezip import RemoteZip


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("url", help="HTTPS URL of the PRIDE .zip bundle")
    ap.add_argument("--out", required=True, type=Path, help="Output directory")
    ap.add_argument("--max-d", type=int, default=1, help="Number of .d folders to extract")
    ap.add_argument("--contains", default="", help="Only .d folders whose name contains this substring")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    with RemoteZip(args.url) as rz:
        names = rz.namelist()
        d_roots = sorted({
            n[: n.index(".d/") + 2]
            for n in names
            if ".d/" in n
        })
        if args.contains:
            d_roots = [d for d in d_roots if args.contains in d]
        print(f"Found {len(d_roots)} matching .d folder(s) in the bundle")
        if not d_roots:
            print("ERROR: no .d folders found", file=sys.stderr)
            return 1

        selected = d_roots[: args.max_d]
        for d in selected:
            members = [n for n in names if n.startswith(d + "/")]
            total = sum(rz.getinfo(m).file_size for m in members)
            print(f"Extracting {d}  ({len(members)} members, {total/1e9:.2f} GB)...")
            for m in members:
                rz.extract(m, args.out)
        print(f"\nExtracted {len(selected)} .d folder(s) to {args.out}")
        for d in selected:
            print(f"  {args.out / d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
