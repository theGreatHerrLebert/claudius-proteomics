#!/usr/bin/env python3
"""
Batch paper-metadata extraction — runs OFF-CLUSTER.

Mogon NHR has no Anthropic API access, so paper-based metadata extraction runs
here (locally), and the resulting paper_extraction.yaml files are shipped to
Mogon, where the runner merges them with no API needed. See docs/RUN_LEDGER.md.

For each accession this:
  1. fetches PRIDE metadata (for the DOI) into pride_metadata.yaml, if absent;
  2. downloads the open-access PDF via Unpaywall;
  3. extracts LC-MS fields with the Claude API;
  4. writes data/metadata/{accession}/paper_extraction.yaml.

Required in the environment (never committed to the repo):
  ANTHROPIC_API_KEY   Claude API key
  UNPAYWALL_EMAIL     contact email for the Unpaywall API

Usage:
  ANTHROPIC_API_KEY=... UNPAYWALL_EMAIL=you@example.org \\
      python scripts/extract_papers_batch.py --list scripts/cluster/download_accessions.txt
  python scripts/extract_papers_batch.py PXD019086 PXD051790

Then ship the results to Mogon, e.g.:
  rsync -av data/metadata/ mogon-nhr:/lustre/project/ki-proanagi/dateschn/data/metadata/
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from pride_metadata import DatasetMetadata          # noqa: E402
from paper_metadata import run_paper_extraction      # noqa: E402

MODEL = "claude-sonnet-4-6"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("accessions", nargs="*", help="PRIDE accessions to process")
    ap.add_argument("--list", type=Path,
                    help="File with one accession per line (e.g. download_accessions.txt)")
    ap.add_argument("--metadata-root", type=Path, default=Path("data/metadata"),
                    help="Per-dataset metadata root (default: data/metadata)")
    ap.add_argument("--limit", type=int, default=0, help="Process at most N accessions")
    args = ap.parse_args()

    accessions = list(args.accessions)
    if args.list and args.list.exists():
        accessions += [
            ln.split()[0] for ln in args.list.read_text().splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")
        ]
    seen: set = set()
    accessions = [a for a in accessions if not (a in seen or seen.add(a))]
    if args.limit > 0:
        accessions = accessions[:args.limit]
    if not accessions:
        ap.error("no accessions given (positional args or --list)")

    if not os.environ.get("ANTHROPIC_API_KEY"):
        print("WARNING: ANTHROPIC_API_KEY not set — LLM extraction will be skipped.")
    if not os.environ.get("UNPAYWALL_EMAIL"):
        print("WARNING: UNPAYWALL_EMAIL not set — PDF download will be skipped.")

    config = {"paper_extraction": {
        "enabled": True,
        "skip_if_exists": True,
        "max_pdf_chars": 80_000,
        "model": MODEL,
    }}

    counts: dict = {}
    for i, acc in enumerate(accessions, 1):
        print(f"\n[{i}/{len(accessions)}] {acc}")
        meta_dir = args.metadata_root / acc
        meta_dir.mkdir(parents=True, exist_ok=True)

        # Pick up any manually-downloaded PDFs from the committed paper inbox
        # (paper_inbox/{acc}/) into the dataset's working papers/ folder.
        inbox = Path("paper_inbox") / acc
        if inbox.is_dir():
            papers_dir = meta_dir / "papers"
            papers_dir.mkdir(parents=True, exist_ok=True)
            for pdf in list(inbox.glob("*.pdf")) + list(inbox.glob("*.PDF")):
                dest = papers_dir / pdf.name
                if not dest.exists():
                    shutil.copy2(pdf, dest)
                    print(f"  picked up manual PDF: {pdf.name}")

        # Re-process datasets whose prior extraction was unsuccessful
        # (no_pdf / no_text / no_api_key / llm_failed): only a 'success' result
        # is final. Lets a manual-PDF backfill or a later LLM run re-run cleanly.
        prev = meta_dir / "paper_extraction.yaml"
        if prev.exists():
            try:
                prev_status = (yaml.safe_load(prev.read_text()) or {}).get("status")
            except Exception:
                prev_status = None
            if prev_status and prev_status != "success":
                prev.unlink()

        # run_paper_extraction needs the DOI from pride_metadata.yaml — fetch if absent.
        pride_yaml = meta_dir / "pride_metadata.yaml"
        if not pride_yaml.exists():
            try:
                DatasetMetadata.from_pride_api(acc).to_yaml(pride_yaml)
                print("  fetched PRIDE metadata")
            except Exception as e:
                print(f"  PRIDE metadata fetch failed: {e}")
                counts["error"] = counts.get("error", 0) + 1
                continue

        try:
            result = run_paper_extraction(acc, meta_dir, config)
            status = result.get("status", "unknown") if result else "skipped"
        except Exception as e:
            print(f"  extraction failed: {e}")
            status = "error"
        counts[status] = counts.get(status, 0) + 1

    print("\n=== summary ===")
    for status, n in sorted(counts.items()):
        print(f"  {status}: {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
