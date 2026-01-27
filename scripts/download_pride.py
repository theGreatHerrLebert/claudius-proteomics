#!/usr/bin/env python3
"""
Download raw data from PRIDE Archive.

Downloads timsTOF .d folders for a given accession number.
"""

import argparse
import subprocess
import sys
from pathlib import Path


def download_pride(accession: str, output_dir: Path) -> None:
    """
    Download raw files from PRIDE Archive.

    Args:
        accession: PRIDE accession number (e.g., PXD019086)
        output_dir: Directory to save downloaded files
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Downloading {accession} to {output_dir}")

    # Option 1: Use pride-archive-downloader (if available)
    # Option 2: Use aspera (faster for large files)
    # Option 3: Use FTP/HTTP fallback

    # For now, use pride_inspector or wget
    # TODO: Implement proper PRIDE API download

    pride_ftp_base = f"ftp://ftp.pride.ebi.ac.uk/pride/data/archive"

    # Get file list from PRIDE API
    # For POC, we'll document manual download steps

    print(f"TODO: Implement automated download for {accession}")
    print(f"Manual download instructions:")
    print(f"  1. Go to https://www.ebi.ac.uk/pride/archive/projects/{accession}")
    print(f"  2. Download raw files (.d folders for timsTOF)")
    print(f"  3. Place them in {output_dir}/")

    # Create placeholder to indicate download location
    placeholder = output_dir / ".download_instructions.txt"
    placeholder.write_text(
        f"Download raw files from:\n"
        f"https://www.ebi.ac.uk/pride/archive/projects/{accession}\n"
        f"\n"
        f"For timsTOF data, download the .d folders.\n"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Download raw data from PRIDE Archive"
    )
    parser.add_argument(
        "--accession", "-a",
        required=True,
        help="PRIDE accession number (e.g., PXD019086)"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        type=Path,
        help="Output directory"
    )

    args = parser.parse_args()

    download_pride(args.accession, args.output)


if __name__ == "__main__":
    main()
