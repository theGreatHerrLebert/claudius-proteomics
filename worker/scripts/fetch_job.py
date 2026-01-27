#!/usr/bin/env python3
"""
Fetch next available job from queue.

Supports multiple queue backends:
- S3 bucket (jobs/ prefix)
- Simple REST API
- Local directory (for testing)
"""

import argparse
import json
import os
import sys
from pathlib import Path
from urllib.parse import urlparse

import yaml


def fetch_from_s3(url: str, worker_id: str, output: Path) -> Path | None:
    """Fetch job from S3 bucket."""
    import boto3

    parsed = urlparse(url)
    bucket = parsed.netloc
    prefix = parsed.path.lstrip("/") or "jobs/pending/"

    s3 = boto3.client("s3")

    # List pending jobs
    response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix, MaxKeys=1)

    if "Contents" not in response:
        return None

    # Get first job
    job_key = response["Contents"][0]["Key"]

    # Download job manifest
    s3.download_file(bucket, job_key, str(output))

    # Move to processing (claim the job)
    new_key = job_key.replace("/pending/", f"/processing/{worker_id}/")
    s3.copy_object(
        Bucket=bucket,
        CopySource={"Bucket": bucket, "Key": job_key},
        Key=new_key,
    )
    s3.delete_object(Bucket=bucket, Key=job_key)

    print(f"Claimed job: {job_key}")
    return output


def fetch_from_api(url: str, worker_id: str, output: Path) -> Path | None:
    """Fetch job from REST API."""
    import requests

    response = requests.post(
        f"{url}/claim",
        json={"worker_id": worker_id},
        timeout=30,
    )

    if response.status_code == 204:
        # No jobs available
        return None

    response.raise_for_status()
    job = response.json()

    with open(output, "w") as f:
        yaml.dump(job, f)

    print(f"Claimed job: {job.get('job_id', 'unknown')}")
    return output


def fetch_from_directory(path: str, worker_id: str, output: Path) -> Path | None:
    """Fetch job from local directory (for testing)."""
    pending_dir = Path(path) / "pending"
    processing_dir = Path(path) / "processing" / worker_id

    if not pending_dir.exists():
        return None

    # Get first pending job
    jobs = list(pending_dir.glob("*.yaml"))
    if not jobs:
        return None

    job_file = jobs[0]

    # Move to processing
    processing_dir.mkdir(parents=True, exist_ok=True)
    dest = processing_dir / job_file.name
    job_file.rename(dest)

    # Copy to output
    with open(dest) as f:
        job = yaml.safe_load(f)
    with open(output, "w") as f:
        yaml.dump(job, f)

    print(f"Claimed job: {job_file.name}")
    return output


def main():
    parser = argparse.ArgumentParser(description="Fetch job from queue")
    parser.add_argument("--url", required=True, help="Queue URL (s3://, http://, or file://)")
    parser.add_argument("--worker", required=True, help="Worker ID")
    parser.add_argument("--output", required=True, help="Output file for job manifest")
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    parsed = urlparse(args.url)

    try:
        if parsed.scheme == "s3":
            result = fetch_from_s3(args.url, args.worker, output)
        elif parsed.scheme in ("http", "https"):
            result = fetch_from_api(args.url, args.worker, output)
        elif parsed.scheme == "file" or not parsed.scheme:
            path = parsed.path or args.url
            result = fetch_from_directory(path, args.worker, output)
        else:
            print(f"Unknown queue type: {parsed.scheme}", file=sys.stderr)
            sys.exit(1)

        if result:
            print(str(result))
        else:
            sys.exit(0)  # No jobs, but not an error

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
