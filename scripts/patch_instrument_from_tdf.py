#!/usr/bin/env python3
"""
Surgically patch the `instrument` field in pride_metadata.yaml using the
truth from Bruker analysis.tdf GlobalMetadata.

This is the third step in the PXD073076 cleanup chain:
  1. audit_pride_instruments.py     — flag deposits where instruments[0] is wrong
  2. derive_instrument_from_tdf.py  — read truth from .d/analysis.tdf, diff
  3. patch_instrument_from_tdf.py   — apply the fix (this script)

Why not just re-run scripts/extract_raw_metadata.py? That script does the
full PRIDE+raw merge and rewrites gradient_length, acquisition_mode, and
lc_system in addition to instrument. Useful in general but a bigger blast
than the audit authorized; this script touches only `instrument`. The
systemic fix (wire extract_raw_metadata.py into the campaign pipeline so
new deposits never have this bug) is a separate piece of work.

Behavior:
  * Default is DRY-RUN — prints the proposed before/after for each file.
  * `--apply` writes the change. The pre-edit file is preserved as
    pride_metadata.yaml.bak.<UTC timestamp> next to the original.
  * `--apply` refuses to overwrite if the recorded value already matches
    the .d truth, or if no .d files are available for the accession.
  * The patched block records source = raw_data.GlobalMetadata.InstrumentName,
    so a re-run of audit_pride_instruments.py will no longer flag the
    deposit (audit only flags `pride_api.instruments[0].name`).

Stdlib only, 3.6-compatible (see [[mogon-nhr-access]]).

Usage:
  # Dry-run for the 8 definitive cases:
  python scripts/patch_instrument_from_tdf.py \\
      --raw-root  /lustre/project/ki-proanagi/dateschn/data/raw \\
      --meta-root /lustre/project/ki-proanagi/dateschn/data/metadata \\
      PXD050530 PXD053881 PXD056978 PXD058376 \\
      PXD070632 PXD073076 PXD046755 PXD069593

  # Apply the changes:
  python scripts/patch_instrument_from_tdf.py --apply \\
      --raw-root  /lustre/project/ki-proanagi/dateschn/data/raw \\
      --meta-root /lustre/project/ki-proanagi/dateschn/data/metadata \\
      PXD050530 PXD053881 PXD056978 PXD058376 \\
      PXD070632 PXD073076 PXD046755 PXD069593
"""
import argparse
import datetime as _dt
import re
import sqlite3
import sys
from pathlib import Path
from typing import List, Optional, Tuple


_INSTRUMENT_BLOCK_RE = re.compile(
    r"^  instrument:\n(?:    .*\n)*",
    re.MULTILINE,
)


def read_truth(acc: str, raw_root: Path) -> Tuple[Optional[str], List[str], int]:
    """Return (single_truth, distinct_values, n_d_files_read).

    single_truth is None if no .d files, or if the deposit is genuinely
    multi-instrument (distinct values > 1) — both are refused as auto-patchable.
    """
    acc_dir = raw_root / acc
    if not acc_dir.exists():
        return None, [], 0
    names = []
    for tdf in sorted(acc_dir.rglob("*.d/analysis.tdf")):
        if "_quarantine" in tdf.parts:
            continue
        try:
            con = sqlite3.connect("file:%s?mode=ro" % tdf, uri=True)
            try:
                row = con.execute(
                    "SELECT Value FROM GlobalMetadata WHERE Key = 'InstrumentName'"
                ).fetchone()
            finally:
                con.close()
        except sqlite3.DatabaseError:
            continue
        if row and row[0]:
            names.append(row[0].strip())
    distinct = sorted(set(names))
    if len(distinct) == 1:
        return distinct[0], distinct, len(names)
    return None, distinct, len(names)


def patch_yaml(text: str, new_value: str) -> str:
    """Rewrite the `instrument:` block in pride_metadata.yaml."""
    if not _INSTRUMENT_BLOCK_RE.search(text):
        raise ValueError("No `instrument:` block found in metadata yaml")
    block = (
        "  instrument:\n"
        "    value: %s\n"
        "    status: auto\n"
        "    source: raw_data.GlobalMetadata.InstrumentName\n"
    ) % new_value
    return _INSTRUMENT_BLOCK_RE.sub(block, text, count=1)


def current_value(text: str) -> Optional[str]:
    """Pull the existing `instrument.value` out of the YAML."""
    m = _INSTRUMENT_BLOCK_RE.search(text)
    if not m:
        return None
    block = m.group(0)
    vm = re.search(r"^    value:\s*(.*)$", block, re.MULTILINE)
    if not vm:
        return None
    v = vm.group(1).strip()
    if v in ("null", "~", ""):
        return None
    return v.strip('"').strip("'")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("accessions", nargs="+", help="Accessions to patch.")
    ap.add_argument("--raw-root", type=Path, required=True)
    ap.add_argument("--meta-root", type=Path, required=True)
    ap.add_argument("--apply", action="store_true",
                    help="Actually write the changes. Without this flag, "
                         "the script only prints proposed before/after.")
    args = ap.parse_args()

    stamp = _dt.datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")

    n_changed = 0
    n_skipped = 0
    n_no_op = 0
    for acc in args.accessions:
        yaml_path = args.meta_root / acc / "pride_metadata.yaml"
        if not yaml_path.exists():
            print("\n%s  SKIP: no pride_metadata.yaml at %s" % (acc, yaml_path))
            n_skipped += 1
            continue
        text = yaml_path.read_text()
        current = current_value(text)

        truth, distinct, n_files = read_truth(acc, args.raw_root)
        if truth is None:
            reason = "no .d files" if n_files == 0 else "ambiguous: %s" % distinct
            print("\n%s  SKIP: %s (recorded=%r, %d .d read)" %
                  (acc, reason, current, n_files))
            n_skipped += 1
            continue

        if current == truth:
            print("\n%s  NOOP: already %r (%d .d agree)" % (acc, truth, n_files))
            n_no_op += 1
            continue

        new_text = patch_yaml(text, truth)
        print()
        print("%s  PATCH: %r -> %r  (%d .d agree)" % (acc, current, truth, n_files))
        # Show the precise diff of the instrument block for review.
        old_block = _INSTRUMENT_BLOCK_RE.search(text).group(0)
        new_block = _INSTRUMENT_BLOCK_RE.search(new_text).group(0)
        print("  --- before ---")
        for line in old_block.rstrip().split("\n"):
            print("  - " + line)
        print("  --- after ---")
        for line in new_block.rstrip().split("\n"):
            print("  + " + line)

        if args.apply:
            backup = yaml_path.with_suffix(".yaml.bak." + stamp)
            backup.write_text(text)
            yaml_path.write_text(new_text)
            print("  -> wrote %s (backup at %s)" % (yaml_path, backup.name))
            n_changed += 1

    print()
    if args.apply:
        print("Done: %d patched, %d already correct, %d skipped." %
              (n_changed, n_no_op, n_skipped))
    else:
        print("Dry run: would patch %d, %d already correct, %d skipped. "
              "Re-run with --apply to write." %
              (sum(1 for _ in args.accessions) - n_no_op - n_skipped,
               n_no_op, n_skipped))
    return 0


if __name__ == "__main__":
    sys.exit(main())
