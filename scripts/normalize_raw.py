#!/usr/bin/env python3
"""
Normalise downloaded raw-data layout.

A retroactive pass over data/raw/{accession}/ for datasets fetched before the
download extraction was fixed (commit 9cb13ed). Per dataset:
  - delete archives that are redundant (their .d content is already extracted),
  - extract any archive whose .d content is still missing (.zip/.7z/.tar*/.rar),
  - flatten .d folders nested inside a wrapper directory,
  - drop broken / empty .d folders,
so every data/raw/{acc}/ holds clean, flat, processable .d folders and no
redundant archives. Idempotent and safe to re-run.

Usage:
  python scripts/normalize_raw.py --raw-root /lustre/.../data/raw --dry-run
  python scripts/normalize_raw.py --raw-root /lustre/.../data/raw
  python scripts/normalize_raw.py --raw-root ... PXD073076 PXD069027
"""
import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from download_pride import extract_archives, _flatten_d_wrappers  # noqa: E402

ARCHIVE_GLOBS = ("*.zip", "*.7z", "*.rar", "*.tar", "*.tar.*")
_D_EXTS = (".d.zip", ".d.7z", ".d.rar", ".d.tar.gz", ".d.tar.bz2",
           ".d.tar.xz", ".d.tar")


def _is_real_d(p: Path) -> bool:
    return p.is_dir() and ((p / "analysis.tdf").exists()
                           or (p / "analysis.tdf_bin").exists())


def _has_inner_d(p: Path) -> bool:
    return p.is_dir() and any(
        c.is_dir() and c.name.endswith(".d") for c in p.iterdir())


def _d_name_for(archive: Path) -> str:
    """The .d folder name an archive should produce."""
    n = archive.name
    for ext in _D_EXTS:
        if n.lower().endswith(ext):
            return n[: -(len(ext) - 2)]  # strip archive suffix, keep ".d"
    return archive.stem


def normalize_dataset(d_root: Path, dry_run: bool) -> dict:
    s = {"freed": 0, "deleted": 0, "to_extract": 0, "flatten": 0, "broken": 0}
    archives = [a for g in ARCHIVE_GLOBS for a in d_root.glob(g)]

    # 1. redundant archives — their .d content is already on disk
    remaining = []
    for arc in archives:
        d_dir = d_root / _d_name_for(arc)
        if _is_real_d(d_dir) or _has_inner_d(d_dir):
            s["freed"] += arc.stat().st_size
            s["deleted"] += 1
            if not dry_run:
                arc.unlink()
        else:
            remaining.append(arc)
    s["to_extract"] = len(remaining)

    # 2. nested .d wrappers to flatten (count before; fix below)
    s["flatten"] = sum(
        1 for sub in d_root.iterdir()
        if sub.is_dir() and not _is_real_d(sub) and _has_inner_d(sub))

    if not dry_run:
        if remaining:
            extract_archives(d_root, cleanup=True)  # extracts + flattens + cleans
        else:
            _flatten_d_wrappers(d_root)

    # 3. broken / empty top-level .d folders
    for sub in list(d_root.glob("*.d")):
        if sub.is_dir() and not _is_real_d(sub) and not _has_inner_d(sub):
            s["broken"] += 1
            if not dry_run:
                shutil.rmtree(sub, ignore_errors=True)
    return s


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("accessions", nargs="*")
    ap.add_argument("--raw-root", type=Path, required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    dirs = ([args.raw_root / a for a in args.accessions] if args.accessions
            else sorted(d for d in args.raw_root.glob("PXD*") if d.is_dir()))

    tot = {"freed": 0, "deleted": 0, "to_extract": 0, "flatten": 0, "broken": 0}
    for d in dirs:
        if not d.is_dir():
            print(f"  {d.name}: missing")
            continue
        s = normalize_dataset(d, args.dry_run)
        if s["deleted"] or s["to_extract"] or s["flatten"] or s["broken"]:
            print(f"  {d.name}: -{s['deleted']} archives ({s['freed']/1e9:.1f} GB), "
                  f"extract {s['to_extract']}, flatten {s['flatten']}, "
                  f"broken .d {s['broken']}")
        for k in tot:
            tot[k] += s[k]

    tag = "[dry-run] would: " if args.dry_run else "done: "
    print(f"\n{tag}{len(dirs)} datasets | delete {tot['deleted']} archives "
          f"({tot['freed']/1e12:.2f} TB freed) | extract {tot['to_extract']} | "
          f"flatten {tot['flatten']} | remove {tot['broken']} broken .d")
    return 0


if __name__ == "__main__":
    sys.exit(main())
