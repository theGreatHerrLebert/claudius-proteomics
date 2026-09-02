#!/usr/bin/env python3
"""preflight_raw.py — normalise and gate a dataset's raw dir before searching.

Array job 569993 (2026-08-24) burned three node allocations on datasets that
could never have searched: the failures only surfaced after the download and
the engine launch, ~35 min in. Every one of them was detectable in seconds.

Runs after download, before step 2. Two jobs:

  NORMALISE (repairs, in place)
    - flatten a doubly-nested `X.d/X.d/` (some PRIDE zips wrap the .d again)
    - drop `__MACOSX` sidecar dirs
    - rename an extension-less dir that contains analysis.tdf to `<name>.d`
    - restore anything in `_quarantine/` that is valid after the above

  GATE (fails fast, non-zero)
    - at least one valid `.d` (analysis.tdf readable, .tdf_bin not truncated)
      must remain at depth 1
    - the resolved organism must have a usable FASTA source in the config

Usage:
  python scripts/cluster/preflight_raw.py <ACCESSION> [--config config/config.mogon.yaml]
                                          [--data-dir DIR] [--normalise-only]

Exit codes: 0 ok · 2 no usable raw · 3 no FASTA source · 4 bad invocation.
"""
import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from validate_d import validate_d_folder  # noqa: E402

OK, NO_RAW, NO_FASTA, BAD_ARGS = 0, 2, 3, 4


def _log(msg):
    print(f"  preflight: {msg}", flush=True)


def normalise(raw_dir: Path) -> dict:
    """Repair the common packaging variants. Returns a summary dict."""
    acts = {"flattened": [], "macosx": [], "renamed": [], "unquarantined": []}
    if not raw_dir.is_dir():
        return acts

    for junk in list(raw_dir.rglob("__MACOSX")):
        if junk.is_dir():
            shutil.rmtree(junk, ignore_errors=True)
            acts["macosx"].append(junk.name)

    # `X.d/X.d/` -> `X.d/`   (zip wrapped the .d folder a second time)
    for d in sorted(raw_dir.glob("*.d")):
        if not d.is_dir() or (d / "analysis.tdf").exists():
            continue
        inner = d / d.name
        if (inner / "analysis.tdf").exists():
            tmp = raw_dir / f"_flatten_{d.name}"
            shutil.move(str(inner), str(tmp))
            shutil.rmtree(d, ignore_errors=True)
            shutil.move(str(tmp), str(d))
            acts["flattened"].append(d.name)

    # extension-less dir holding a .d payload -> add the suffix
    for d in sorted(p for p in raw_dir.iterdir() if p.is_dir()):
        if d.name.startswith("_") or d.suffix == ".d":
            continue
        if (d / "analysis.tdf").exists():
            dest = d.with_name(d.name + ".d")
            if not dest.exists():
                shutil.move(str(d), str(dest))
                acts["renamed"].append(d.name)

    # anything quarantined that is valid after the repairs comes back
    qdir = raw_dir / "_quarantine"
    if qdir.is_dir():
        for d in sorted(qdir.glob("*.d")):
            inner = d / d.name
            src = inner if (inner / "analysis.tdf").exists() else d
            if not (src / "analysis.tdf").exists():
                continue
            dest = raw_dir / d.name
            if dest.exists():
                continue
            shutil.move(str(src), str(dest))
            shutil.rmtree(d, ignore_errors=True)
            acts["unquarantined"].append(d.name)
    return acts


def fasta_source(config: dict, organism: str, project_root: Path) -> str | None:
    """Describe where this organism's FASTA would come from, or None."""
    entry = (config.get("organisms") or {}).get(organism) or {}
    if entry.get("local_fasta") and Path(entry["local_fasta"]).exists():
        return f"local_fasta {entry['local_fasta']}"
    if entry.get("proteome_id"):
        return f"UniProt {entry['proteome_id']}"
    prebuilt = project_root / "resources" / "fasta" / "search_db" / f"{organism}_decoys.fasta"
    if prebuilt.exists():
        return f"prebuilt {prebuilt.name}"
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("accession")
    ap.add_argument("--config", type=Path, default=Path("config/config.mogon.yaml"))
    ap.add_argument("--data-dir", type=Path, default=Path("data"))
    ap.add_argument("--normalise-only", action="store_true")
    args = ap.parse_args()

    raw_dir = args.data_dir / "raw" / args.accession
    acts = normalise(raw_dir)
    for key, label in (("flattened", "flattened nested .d"),
                       ("macosx", "removed __MACOSX"),
                       ("renamed", "added .d suffix"),
                       ("unquarantined", "restored from _quarantine")):
        if acts[key]:
            _log(f"{label}: {len(acts[key])} ({', '.join(acts[key][:3])}"
                 f"{'...' if len(acts[key]) > 3 else ''})")

    valid, bad = [], []
    for d in sorted(raw_dir.glob("*.d")) if raw_dir.is_dir() else []:
        (valid if validate_d_folder(d)[0] else bad).append(d.name)
    _log(f"{len(valid)} valid .d, {len(bad)} unusable")

    if args.normalise_only:
        return OK
    if not valid:
        _log(f"FAIL: no usable .d in {raw_dir} — nothing to search. "
             f"The deposit may not ship timsTOF .d (check the RAW file names on "
             f"PRIDE) or the archives did not extract.")
        return NO_RAW

    import yaml
    config = yaml.safe_load(args.config.read_text())
    organism = ((config.get("dataset_metadata") or {})
                .get(args.accession, {}).get("organism"))
    if organism:
        src = fasta_source(config, organism, args.config.resolve().parents[1])
        if src is None:
            _log(f"FAIL: organism '{organism}' has no FASTA source. Set "
                 f"organisms.{organism}.proteome_id or .local_fasta, or drop the "
                 f"dataset. Note: not every species has a UniProt reference "
                 f"proteome.")
            return NO_FASTA
        _log(f"FASTA source for '{organism}': {src}")
    else:
        _log(f"no organism in dataset_metadata for {args.accession}; "
             f"groups will be inferred and the FASTA checked at search time")

    _log("OK")
    return OK


if __name__ == "__main__":
    sys.exit(main())
