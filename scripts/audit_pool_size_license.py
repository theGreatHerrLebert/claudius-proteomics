#!/usr/bin/env python3
"""audit_pool_size_license.py — size + CC0 audit for a pool of PRIDE accessions.

The discovery sweep's stage-3 file check is the first thing to get API-throttled,
so most of the catalog carries no `bruker_size_gb` and no verified license. Both
are hard gates before a processing wave: CC0 because the corpus guarantees it,
size because the July 2026 quota wall showed what an unmeasured wave costs.

For each accession: pull the project record (license, PTMs, instrument, organism)
and the file listing (Bruker .d/.zip count + bytes). Rows are appended as they
land, so an interrupted run resumes from the TSV instead of re-fetching.

    python scripts/audit_pool_size_license.py <accessions.txt> -o <out.tsv>

Deliberately gentle: one dataset at a time, sleep between calls, exponential
backoff on 429/503. This walks a public API on someone else's budget.
"""
import argparse
import csv
import fnmatch
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

API = "https://www.ebi.ac.uk/pride/ws/archive/v2"
API3 = "https://www.ebi.ac.uk/pride/ws/archive/v3"   # v2 files/byProject returns an empty body
BRUKER_PATTERNS = ["*.d", "*.d.zip", "*.d.rar", "*.tdf", "*.tdf_bin", "*brukertimstof*"]
FIELDS = ["accession", "cc0", "license", "n_bruker_files", "bruker_size_gb",
          "n_raw_files", "raw_size_gb",
          "instrument", "organism", "ptms", "publication_date", "title"]


class FetchError(Exception):
    """A request did not produce an answer. NOT the same as an empty answer.

    This distinction is the whole point: PRIDE's v2 files/byProject returns
    HTTP 200 with an empty body, which an earlier version of this script read
    as "this dataset has no raw files" and silently sized 191/544 datasets at
    0 GB. A failure must never be representable as data.
    """


def get(url, sleep, max_retries=4):
    """Return the decoded JSON, or raise FetchError. Never returns None to mean
    'failed' — callers cannot distinguish that from a legitimately empty list."""
    last = None
    for attempt in range(max_retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = json.load(r)
            time.sleep(sleep)
            return data
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (429, 503):
                wait = 2 ** (attempt + 2)
                print(f"    rate limited ({e.code}), backing off {wait}s", flush=True)
                time.sleep(wait)
                continue
            if e.code == 404:
                raise FetchError(f"404 {url}") from e
            time.sleep(2 ** (attempt + 1))
        except Exception as e:
            last = e
            if attempt == max_retries - 1:
                break
            time.sleep(2 ** (attempt + 1))
    raise FetchError(f"{type(last).__name__}: {last} ({url})")


def _failed_row(acc, why):
    """A row that cannot be mistaken for a measurement: sizes stay empty, and
    the reason is carried so a failed fetch is never read as 'no files'."""
    row = {k: "" for k in FIELDS}
    row["accession"], row["cc0"], row["license"] = acc, "FETCH_FAIL", str(why)[:110]
    return row


def audit(acc, sleep):
    try:
        proj = get(f"{API}/projects/{acc}", sleep)
    except FetchError as e:
        return _failed_row(acc, e)
    lic = proj.get("license") or ""
    files, page = [], 0
    while True:   # v3 paginates; a few deposits run to thousands of files
        try:
            chunk = get(f"{API3}/projects/{acc}/files?pageSize=1000&page={page}", sleep)
        except FetchError as e:
            # Do NOT treat a failed page as end-of-listing: that is exactly the
            # bug this script exists to catch. Fail the row instead.
            return _failed_row(acc, f"file listing page {page}: {e}")
        if not isinstance(chunk, list):
            return _failed_row(acc, f"file listing page {page}: unexpected {type(chunk).__name__}")
        files.extend(chunk)
        if len(chunk) < 1000:
            break
        page += 1
    bruker = [f for f in files
              if any(fnmatch.fnmatch((f.get("fileName") or "").lower(), p) for p in BRUKER_PATTERNS)]
    size_gb = sum(f.get("fileSizeBytes", 0) for f in bruker) / (1024 ** 3)
    # Filename globs miss archive-wrapped raw (raw.zip, *.tar.gz, Topology.rar), so
    # size on PRIDE's own RAW category too - the catalog is already timsTOF-only,
    # so RAW here means Bruker .d in some wrapper.
    raw = [f for f in files if ((f.get("fileCategory") or {}).get("value") or "").upper() == "RAW"]
    raw_gb = sum(f.get("fileSizeBytes", 0) for f in raw) / (1024 ** 3)
    return {
        "accession": acc,
        "cc0": "Y" if ("CC0" in lic.upper() or "public domain" in lic.lower()) else "N",
        "license": lic,
        "n_bruker_files": len(bruker),
        "bruker_size_gb": f"{size_gb:.2f}",
        "n_raw_files": len(raw),
        "raw_size_gb": f"{raw_gb:.2f}",
        "instrument": ",".join(i.get("name", "") for i in proj.get("instruments", []))[:60],
        "organism": ",".join(o.get("name", "") for o in proj.get("organisms", []))[:60],
        "ptms": ";".join(p.get("name", "") for p in (proj.get("identifiedPTMStrings") or []))[:120],
        "publication_date": (proj.get("publicationDate") or "")[:10],
        "title": (proj.get("title") or "").replace("\t", " ")[:110],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accessions", type=Path)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--sleep", type=float, default=0.6, help="seconds between API calls")
    args = ap.parse_args()

    accs = [l.strip() for l in args.accessions.read_text().splitlines() if l.strip().startswith("PXD")]
    done = set()
    if args.out.exists():
        with args.out.open() as f:
            done = {r["accession"] for r in csv.DictReader(f, delimiter="\t")}
        print(f"resuming: {len(done)} already audited")
    todo = [a for a in accs if a not in done]
    print(f"{len(accs)} accessions, {len(todo)} to fetch (~{len(todo) * args.sleep * 2 / 60:.0f} min)", flush=True)

    new_file = not args.out.exists()
    with args.out.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        if new_file:
            w.writeheader()
        for i, acc in enumerate(todo, 1):
            row = audit(acc, args.sleep)
            w.writerow(row)
            f.flush()
            if i % 25 == 0 or i == len(todo):
                print(f"  [{i}/{len(todo)}] {acc} cc0={row['cc0']} "
                      f"{row['bruker_size_gb']}GB", flush=True)

    with args.out.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    cc0 = [r for r in rows if r["cc0"] == "Y"]
    sizes = [float(r["bruker_size_gb"]) for r in cc0 if r["bruker_size_gb"]]
    print(f"\naudited {len(rows)}: CC0 {len(cc0)}, non-CC0 "
          f"{sum(1 for r in rows if r['cc0'] == 'N')}, failed "
          f"{sum(1 for r in rows if r['cc0'] == 'FETCH_FAIL')}")
    if sizes:
        sizes.sort()
        print(f"CC0 raw volume: total {sum(sizes) / 1024:.1f} TB, median "
              f"{sizes[len(sizes) // 2]:.0f} GB, no-files {sum(1 for s in sizes if s == 0)}")
    print(f"-> {args.out}")


if __name__ == "__main__":
    sys.exit(main())
