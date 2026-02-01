#!/usr/bin/env python3
"""
Download raw data from PRIDE Archive with metadata extraction.

This script:
1. Fetches metadata from PRIDE REST API
2. Parses protocol text for gradient/column/mode hints
3. Downloads raw files using pridepy
4. Creates metadata YAML with validation status

Usage:
    # Metadata only (no download)
    python download_pride.py --accession PXD019086 --metadata-only

    # Full download
    python download_pride.py --accession PXD019086 --output data/raw/PXD019086
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional

import requests

from pride_metadata import (
    DatasetMetadata,
    MetadataField,
    fetch_pride_metadata,
    save_raw_response,
)

# PRIDE REST API v2 base URL
PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v2"


def get_file_list(accession: str) -> List[dict]:
    """
    Get list of files for a PRIDE accession.

    Returns list of file info dicts with fileName, fileSize, ftpLink, etc.
    """
    url = f"{PRIDE_API_BASE}/files/byProject?accession={accession}"

    response = requests.get(url, timeout=60)
    response.raise_for_status()

    return response.json()


def filter_raw_files(
    files: List[dict],
    patterns: Optional[List[str]] = None,
    max_files: int = 0,
) -> List[dict]:
    """
    Filter files to download based on patterns.

    Args:
        files: List of file info dicts from PRIDE API
        patterns: Glob patterns to match (e.g., ["*.d.zip", "*.d.tar"])
        max_files: Maximum number of files (0 = unlimited)

    Returns:
        Filtered list of files
    """
    if patterns is None:
        # Default: timsTOF raw data patterns
        patterns = ["*.d.zip", "*.d.tar", "*.d", "*.raw"]

    import fnmatch

    filtered = []
    for f in files:
        filename = f.get("fileName", "")
        for pattern in patterns:
            if fnmatch.fnmatch(filename.lower(), pattern.lower()):
                filtered.append(f)
                break

    # Sort by filename for consistent ordering
    filtered.sort(key=lambda x: x.get("fileName", ""))

    # Apply max files limit
    if max_files > 0:
        filtered = filtered[:max_files]

    return filtered


def download_with_pridepy(
    accession: str,
    output_dir: Path,
    file_list: Optional[List[dict]] = None,
    protocol: str = "ftp",
    retry_count: int = 3,
) -> bool:
    """
    Download files using pridepy CLI.

    Args:
        accession: PRIDE accession
        output_dir: Directory to save files
        file_list: Optional list of specific files to download
        protocol: Download protocol ("ftp" or "aspera")
        retry_count: Number of retries for failed downloads

    Returns:
        True if download succeeded
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if pridepy is available
    try:
        subprocess.run(
            ["pridepy", "--help"],
            capture_output=True,
            check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: pridepy is not installed. Install with: pip install pridepy")
        return False

    # Build pridepy command
    cmd = [
        "pridepy",
        "download-all-raw-files",
        "-a", accession,
        "-o", str(output_dir),
        "-p", protocol,
    ]

    print(f"Running: {' '.join(cmd)}")

    for attempt in range(retry_count):
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600 * 4,  # 4 hour timeout for large datasets
            )

            if result.returncode == 0:
                print("Download completed successfully")
                return True

            print(f"Download attempt {attempt + 1} failed: {result.stderr}")

        except subprocess.TimeoutExpired:
            print(f"Download attempt {attempt + 1} timed out")
        except Exception as e:
            print(f"Download attempt {attempt + 1} error: {e}")

        if attempt < retry_count - 1:
            wait_time = 30 * (attempt + 1)
            print(f"Retrying in {wait_time} seconds...")
            time.sleep(wait_time)

    return False


def download_with_wget(
    files: List[dict],
    output_dir: Path,
    retry_count: int = 3,
) -> bool:
    """
    Fallback download using wget for individual files.

    Args:
        files: List of file info dicts with ftpLink
        output_dir: Directory to save files
        retry_count: Number of retries per file

    Returns:
        True if all downloads succeeded
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    success_count = 0
    fail_count = 0

    for f in files:
        url = f.get("ftpLink") or f.get("publicFileLocations", [{}])[0].get("value")
        filename = f.get("fileName", "")

        if not url:
            print(f"Warning: No download URL for {filename}")
            fail_count += 1
            continue

        output_path = output_dir / filename

        # Skip if already downloaded
        if output_path.exists():
            expected_size = f.get("fileSize", 0)
            actual_size = output_path.stat().st_size
            if actual_size == expected_size:
                print(f"Skipping {filename} (already downloaded)")
                success_count += 1
                continue

        print(f"Downloading {filename}...")

        cmd = [
            "wget",
            "-c",  # Continue partial downloads
            "-q",
            "--show-progress",
            "-O", str(output_path),
            url,
        ]

        for attempt in range(retry_count):
            try:
                result = subprocess.run(cmd, timeout=3600)
                if result.returncode == 0:
                    success_count += 1
                    break
            except subprocess.TimeoutExpired:
                print(f"Download timed out for {filename}")
            except Exception as e:
                print(f"Download error for {filename}: {e}")

            if attempt == retry_count - 1:
                fail_count += 1

    print(f"Downloaded {success_count}/{len(files)} files ({fail_count} failed)")
    return fail_count == 0


def extract_archives(output_dir: Path) -> None:
    """
    Extract downloaded .zip and .tar archives.

    Args:
        output_dir: Directory containing downloaded archives
    """
    output_dir = Path(output_dir)

    # Extract .zip files
    for zip_file in output_dir.glob("*.zip"):
        print(f"Extracting {zip_file.name}...")
        try:
            import zipfile
            with zipfile.ZipFile(zip_file, "r") as zf:
                zf.extractall(output_dir)
            # Optionally remove archive after extraction
            # zip_file.unlink()
        except Exception as e:
            print(f"Failed to extract {zip_file.name}: {e}")

    # Extract .tar files
    for tar_file in output_dir.glob("*.tar"):
        print(f"Extracting {tar_file.name}...")
        try:
            import tarfile
            with tarfile.open(tar_file, "r") as tf:
                tf.extractall(output_dir)
        except Exception as e:
            print(f"Failed to extract {tar_file.name}: {e}")

    # Extract .tar.gz files
    for tar_gz in output_dir.glob("*.tar.gz"):
        print(f"Extracting {tar_gz.name}...")
        try:
            import tarfile
            with tarfile.open(tar_gz, "r:gz") as tf:
                tf.extractall(output_dir)
        except Exception as e:
            print(f"Failed to extract {tar_gz.name}: {e}")


def download_pride(
    accession: str,
    output_dir: Path,
    metadata_dir: Optional[Path] = None,
    protocol: str = "ftp",
    retry_count: int = 3,
    max_files: int = 0,
    file_patterns: Optional[List[str]] = None,
    extract: bool = True,
) -> bool:
    """
    Download raw files from PRIDE Archive.

    Args:
        accession: PRIDE accession number (e.g., PXD019086)
        output_dir: Directory to save downloaded files
        metadata_dir: Directory for metadata files (default: data/metadata/{accession})
        protocol: Download protocol ("ftp" or "aspera")
        retry_count: Number of retries for failed downloads
        max_files: Maximum number of files (0 = unlimited, for testing)
        file_patterns: Glob patterns for file selection
        extract: Whether to extract archives after download

    Returns:
        True if download succeeded
    """
    output_dir = Path(output_dir)
    if metadata_dir is None:
        metadata_dir = Path(f"data/metadata/{accession}")

    print(f"=" * 60)
    print(f"PRIDE Download: {accession}")
    print(f"=" * 60)

    # Step 1: Fetch metadata
    print(f"\n1. Fetching PRIDE metadata...")
    metadata = fetch_pride_metadata(accession)

    # Save PRIDE metadata
    metadata_dir.mkdir(parents=True, exist_ok=True)
    pride_meta_file = metadata_dir / "pride_metadata.yaml"
    metadata.to_yaml(pride_meta_file)
    print(f"   Saved PRIDE metadata to {pride_meta_file}")

    # Save raw API response
    try:
        raw_file = save_raw_response(accession, metadata_dir)
        print(f"   Saved raw API response to {raw_file}")
    except Exception as e:
        print(f"   Warning: Could not save raw API response: {e}")

    # Step 2: Get file list
    print(f"\n2. Getting file list...")
    try:
        all_files = get_file_list(accession)
        print(f"   Found {len(all_files)} total files")
    except Exception as e:
        print(f"   Error getting file list: {e}")
        return False

    # Filter files
    files = filter_raw_files(all_files, file_patterns, max_files)
    print(f"   Selected {len(files)} files for download")

    if not files:
        print("   No matching files found!")
        return False

    # Show files to download
    total_size = sum(f.get("fileSize", 0) for f in files)
    print(f"   Total size: {total_size / (1024**3):.2f} GB")
    for f in files[:5]:
        print(f"     - {f.get('fileName')}")
    if len(files) > 5:
        print(f"     ... and {len(files) - 5} more")

    # Step 3: Download files
    print(f"\n3. Downloading files (protocol: {protocol})...")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Try pridepy first
    success = download_with_pridepy(
        accession, output_dir, files, protocol, retry_count
    )

    if not success:
        print("   pridepy failed, trying wget fallback...")
        success = download_with_wget(files, output_dir, retry_count)

    if not success:
        print("   Download failed!")
        return False

    # Step 4: Extract archives
    if extract:
        print(f"\n4. Extracting archives...")
        extract_archives(output_dir)

    # Step 5: Create download complete flag
    flag_file = output_dir / ".download_complete"
    flag_file.write_text(f"Downloaded at: {time.strftime('%Y-%m-%dT%H:%M:%S')}\n")
    print(f"\n5. Download complete!")
    print(f"   Files saved to: {output_dir}")
    print(f"   Metadata saved to: {metadata_dir}")

    return True


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
        type=Path,
        help="Output directory for raw files"
    )
    parser.add_argument(
        "--metadata-dir",
        type=Path,
        help="Directory for metadata files (default: data/metadata/{accession})"
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Only fetch metadata, don't download files"
    )
    parser.add_argument(
        "--protocol",
        choices=["ftp", "aspera"],
        default="ftp",
        help="Download protocol (default: ftp)"
    )
    parser.add_argument(
        "--retry",
        type=int,
        default=3,
        help="Number of download retries (default: 3)"
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Maximum number of files to download (0 = all)"
    )
    parser.add_argument(
        "--no-extract",
        action="store_true",
        help="Don't extract archives after download"
    )

    args = parser.parse_args()

    # Set default output directory
    if args.output is None:
        args.output = Path(f"data/raw/{args.accession}")

    if args.metadata_dir is None:
        args.metadata_dir = Path(f"data/metadata/{args.accession}")

    if args.metadata_only:
        # Just fetch and save metadata
        print(f"Fetching metadata for {args.accession}...")
        metadata = fetch_pride_metadata(args.accession)

        args.metadata_dir.mkdir(parents=True, exist_ok=True)
        pride_meta_file = args.metadata_dir / "pride_metadata.yaml"
        metadata.to_yaml(pride_meta_file)
        print(f"Saved PRIDE metadata to {pride_meta_file}")

        # Also save raw API response
        try:
            raw_file = save_raw_response(args.accession, args.metadata_dir)
            print(f"Saved raw API response to {raw_file}")
        except Exception as e:
            print(f"Warning: Could not save raw API response: {e}")

        # Print summary
        summary = metadata.validation_summary()
        print(f"\nValidation summary:")
        print(f"  Auto fields: {summary['auto_fields']}")
        print(f"  Inferred fields: {summary['inferred_fields']}")
        print(f"  Missing fields: {summary['missing_fields']}")
        if summary["missing_field_names"]:
            print(f"  Missing: {', '.join(summary['missing_field_names'])}")

        return 0

    # Full download
    success = download_pride(
        accession=args.accession,
        output_dir=args.output,
        metadata_dir=args.metadata_dir,
        protocol=args.protocol,
        retry_count=args.retry,
        max_files=args.max_files,
        extract=not args.no_extract,
    )

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
