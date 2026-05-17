#!/usr/bin/env python3
"""
Build the manual paper inbox.

For every dataset whose publication PDF was NOT auto-downloaded (Unpaywall found
no open-access copy, or the publisher blocked the download), create
paper_inbox/{accession}/ with a committed SOURCE.md recording the DOI / link.

The folder structure and SOURCE.md files are committed; PDFs dropped into them
are gitignored (.gitignore: paper_inbox/**/*.pdf). Workflow: clone main, read
each SOURCE.md, download the paper via institutional access, drop the .pdf into
the matching folder, copy the tree to the extraction machine — extract_papers_
batch.py then picks those PDFs up automatically.

Usage:
  python scripts/build_paper_inbox.py [--metadata-root data/metadata]
"""
import argparse
import sys
from pathlib import Path

import yaml

INBOX = Path("paper_inbox")


def _load(p: Path) -> dict:
    try:
        return yaml.safe_load(p.read_text()) or {}
    except Exception:
        return {}


def _field(d: dict, key: str):
    """Pull a value from a pride_metadata.yaml field (dict-with-value or scalar)."""
    v = d.get(key)
    return v.get("value") if isinstance(v, dict) else v


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--metadata-root", type=Path, default=Path("data/metadata"))
    args = ap.parse_args()

    INBOX.mkdir(exist_ok=True)
    needed = []

    for meta_dir in sorted(args.metadata_root.glob("PXD*")):
        acc = meta_dir.name
        extraction = _load(meta_dir / "paper_extraction.yaml")
        if not extraction:
            continue  # not processed by the download sweep yet
        if extraction.get("pdf_sources"):
            continue  # a PDF was obtained — no manual entry needed

        pride = _load(meta_dir / "pride_metadata.yaml")
        title = _field(pride, "title")
        doi = _field(pride, "publication_doi")
        status = extraction.get("status", "unknown")

        d = INBOX / acc
        d.mkdir(exist_ok=True)
        lines = [f"# {acc} — publication", "", f"- **Title:** {title or '(unknown)'}"]
        if doi:
            lines += [f"- **DOI:** {doi}", f"- **Link:** https://doi.org/{doi}"]
        else:
            lines.append("- **DOI:** not recorded in PRIDE metadata — find the paper "
                         "by title search (the publication usually still exists)")
        lines += [
            f"- **Auto-download:** failed (`{status}`)",
            "",
            "Download the PDF (institutional access if paywalled) and place it in",
            "this folder as a `.pdf`. PDFs here are gitignored — never pushed.",
            "",
        ]
        (d / "SOURCE.md").write_text("\n".join(lines))
        needed.append((acc, doi))

    readme = [
        "# Paper inbox", "",
        f"{len(needed)} dataset(s) need a manually-downloaded PDF — the auto-download "
        "(Unpaywall) found no open-access copy, or the publisher blocked it.", "",
        "Each subfolder has a `SOURCE.md` with the DOI/link. Download the paper, drop "
        "the `.pdf` into the folder (it is gitignored), and copy the tree to the "
        "extraction machine.", "",
        "Regenerate with `python scripts/build_paper_inbox.py`.", "",
    ]
    readme += [
        f"- [{a}]({a}/SOURCE.md)" + (f" — {doi}" if doi else " — (no DOI)")
        for a, doi in needed
    ]
    (INBOX / "README.md").write_text("\n".join(readme) + "\n")

    print(f"paper_inbox: {len(needed)} dataset(s) need a manual PDF")
    return 0


if __name__ == "__main__":
    sys.exit(main())
