#!/usr/bin/env python3
"""
Submit job results back to central storage.

Supports:
- S3 bucket upload
- REST API upload
- Local copy (for testing)
"""

import argparse
import os
import shutil
import sys
from pathlib import Path
from urllib.parse import urlparse

import yaml


def upload_to_s3(results_dir: Path, url: str, job: dict) -> None:
    """Upload results to S3 bucket."""
    import boto3

    parsed = urlparse(url)
    bucket = parsed.netloc
    prefix = parsed.path.lstrip("/") or "results/"

    accession = job["accession"]
    job_id = job["job_id"]

    s3 = boto3.client("s3")

    # Upload all files in results directory
    for file_path in results_dir.rglob("*"):
        if file_path.is_file():
            rel_path = file_path.relative_to(results_dir)
            s3_key = f"{prefix}{accession}/{rel_path}"

            print(f"Uploading: {rel_path} -> s3://{bucket}/{s3_key}")
            s3.upload_file(str(file_path), bucket, s3_key)

    # Mark job as completed
    completed_key = f"jobs/completed/{job_id}.yaml"
    s3.put_object(
        Bucket=bucket,
        Key=completed_key,
        Body=yaml.dump(job).encode(),
    )

    print(f"Results uploaded to s3://{bucket}/{prefix}{accession}/")


def upload_to_api(results_dir: Path, url: str, job: dict) -> None:
    """Upload results via REST API."""
    import requests

    accession = job["accession"]
    job_id = job["job_id"]

    # Create tarball
    archive_path = results_dir.parent / f"{accession}.tar.gz"
    shutil.make_archive(
        str(archive_path).replace(".tar.gz", ""),
        "gztar",
        results_dir.parent,
        results_dir.name,
    )

    # Upload
    with open(archive_path, "rb") as f:
        response = requests.post(
            f"{url}/upload",
            files={"results": f},
            data={"job_id": job_id, "accession": accession},
            timeout=300,
        )
        response.raise_for_status()

    # Clean up
    archive_path.unlink()

    print(f"Results uploaded to {url}")


def copy_to_directory(results_dir: Path, dest_path: str, job: dict) -> None:
    """Copy results to local directory."""
    accession = job["accession"]
    dest = Path(dest_path) / accession

    if dest.exists():
        shutil.rmtree(dest)

    shutil.copytree(results_dir, dest)
    print(f"Results copied to {dest}")


def main():
    parser = argparse.ArgumentParser(description="Submit job results")
    parser.add_argument("--job", required=True, help="Job manifest YAML")
    parser.add_argument("--results", required=True, help="Results directory")
    parser.add_argument("--url", required=True, help="Upload destination (s3://, http://, or file://)")
    args = parser.parse_args()

    with open(args.job) as f:
        job = yaml.safe_load(f)

    results_dir = Path(args.results)
    accession = job["accession"]

    # Find the actual results directory
    if (results_dir / accession).exists():
        results_dir = results_dir / accession
    elif not (results_dir / "COMPLETED").exists():
        print(f"ERROR: Results not found or incomplete: {results_dir}", file=sys.stderr)
        sys.exit(1)

    parsed = urlparse(args.url)

    try:
        if parsed.scheme == "s3":
            upload_to_s3(results_dir, args.url, job)
        elif parsed.scheme in ("http", "https"):
            upload_to_api(results_dir, args.url, job)
        elif parsed.scheme == "file" or not parsed.scheme:
            path = parsed.path or args.url
            copy_to_directory(results_dir, path, job)
        else:
            print(f"Unknown destination type: {parsed.scheme}", file=sys.stderr)
            sys.exit(1)

        print("Results submitted successfully")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
