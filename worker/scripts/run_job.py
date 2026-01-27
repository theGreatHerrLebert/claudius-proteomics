#!/usr/bin/env python3
"""
Execute a San José job.

Orchestrates the full pipeline:
1. Download FASTA (if needed)
2. Download raw files from PRIDE
3. Run FragPipe search
4. Run DIA-NN search
5. Extract raw features
6. Package results
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import requests
import yaml
from tqdm import tqdm


def load_job(job_path: str) -> dict:
    """Load job manifest from YAML."""
    with open(job_path) as f:
        return yaml.safe_load(f)


def download_file(url: str, dest: Path, desc: str = None) -> Path:
    """Download file with progress bar."""
    response = requests.get(url, stream=True, timeout=30)
    response.raise_for_status()

    total_size = int(response.headers.get("content-length", 0))
    desc = desc or dest.name

    with open(dest, "wb") as f:
        with tqdm(total=total_size, unit="B", unit_scale=True, desc=desc) as pbar:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                pbar.update(len(chunk))

    return dest


def prepare_fasta(job: dict, fasta_dir: Path) -> Path:
    """Download and prepare FASTA database."""
    fasta_config = job["fasta"]
    organism = fasta_config["organism"]

    # Check if already exists
    final_fasta = fasta_dir / f"{organism}_crap_decoy.fasta"
    if final_fasta.exists():
        print(f"FASTA exists: {final_fasta}")
        return final_fasta

    print(f"Preparing FASTA for {organism}...")

    # Download organism FASTA
    fasta_url = fasta_config["url"]
    gz_path = fasta_dir / f"{organism}.fasta.gz"

    if not gz_path.exists():
        download_file(fasta_url, gz_path, f"{organism}.fasta.gz")

    # Decompress
    fasta_path = fasta_dir / f"{organism}.fasta"
    if not fasta_path.exists():
        print("Decompressing FASTA...")
        with gzip.open(gz_path, "rb") as f_in:
            with open(fasta_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    # Add contaminants
    if fasta_config.get("include_contaminants", True):
        crap_url = "https://ftp.thegpm.org/fasta/cRAP/crap.fasta"
        crap_path = fasta_dir / "crap.fasta"
        if not crap_path.exists():
            download_file(crap_url, crap_path, "crap.fasta")

        combined_path = fasta_dir / f"{organism}_crap.fasta"
        if not combined_path.exists():
            print("Adding contaminants...")
            with open(combined_path, "w") as out:
                with open(fasta_path) as f:
                    out.write(f.read())
                with open(crap_path) as f:
                    out.write(f.read())
            fasta_path = combined_path

    # Generate decoys
    if fasta_config.get("generate_decoys", True):
        print("Generating decoy sequences...")
        decoy_fasta = generate_decoys(fasta_path, final_fasta)
        return decoy_fasta

    return fasta_path


def generate_decoys(input_fasta: Path, output_fasta: Path) -> Path:
    """Generate reversed decoy sequences."""
    with open(input_fasta) as f_in, open(output_fasta, "w") as f_out:
        # Write original sequences
        sequences = []
        current_header = None
        current_seq = []

        for line in f_in:
            if line.startswith(">"):
                if current_header:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line.strip()
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_header:
            sequences.append((current_header, "".join(current_seq)))

        # Write targets
        for header, seq in sequences:
            f_out.write(f"{header}\n")
            # Write sequence in 60-char lines
            for i in range(0, len(seq), 60):
                f_out.write(f"{seq[i:i+60]}\n")

        # Write decoys (reversed)
        for header, seq in sequences:
            decoy_header = header.replace(">", ">DECOY_")
            f_out.write(f"{decoy_header}\n")
            reversed_seq = seq[::-1]
            for i in range(0, len(reversed_seq), 60):
                f_out.write(f"{reversed_seq[i:i+60]}\n")

    return output_fasta


def download_raw_files(job: dict, data_dir: Path) -> list[Path]:
    """Download raw files from PRIDE."""
    raw_files = job["raw_files"]
    downloaded = []

    for file_info in raw_files:
        name = file_info["name"]
        url = file_info["url"]
        dest = data_dir / name

        if dest.exists():
            print(f"Already downloaded: {name}")
        else:
            print(f"Downloading: {name}")
            download_file(url, dest, name)

            # Unzip if needed
            if name.endswith(".d.zip"):
                print(f"Extracting: {name}")
                shutil.unpack_archive(dest, data_dir)
                dest = data_dir / name.replace(".zip", "")

        downloaded.append(dest)

    return downloaded


def run_fragpipe(
    job: dict,
    raw_files: list[Path],
    fasta: Path,
    fragpipe_path: str,
    output_dir: Path,
) -> Path:
    """Run FragPipe search."""
    print("\n" + "=" * 50)
    print("Running FragPipe")
    print("=" * 50)

    fp_config = job["search_config"]["fragpipe"]
    fp_output = output_dir / "fragpipe"
    fp_output.mkdir(parents=True, exist_ok=True)

    # Create manifest file
    manifest_path = fp_output / "manifest.fp-manifest"
    with open(manifest_path, "w") as f:
        for raw_file in raw_files:
            # Format: path \t experiment \t biorep \t techrep \t data_type
            f.write(f"{raw_file}\texperiment\t1\t1\tDDA\n")

    # Create workflow file path
    workflow = fp_config.get("workflow", "LFQ-MBR")
    workflow_path = Path(fragpipe_path) / "workflows" / f"{workflow}.workflow"

    # Run FragPipe headless
    cmd = [
        "java", "-jar",
        str(Path(fragpipe_path) / "lib" / "fragpipe-*.jar"),
        "--headless",
        "--workflow", str(workflow_path),
        "--manifest", str(manifest_path),
        "--workdir", str(fp_output),
        "--db", str(fasta),
        "--threads", str(fp_config.get("threads", 16)),
        "--ram", str(fp_config.get("memory_gb", 64)),
    ]

    # Use shell expansion for jar wildcard
    cmd_str = " ".join(cmd)
    print(f"Command: {cmd_str}")

    result = subprocess.run(
        cmd_str,
        shell=True,
        cwd=fragpipe_path,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(f"FragPipe STDERR: {result.stderr}")
        raise RuntimeError(f"FragPipe failed with code {result.returncode}")

    print("FragPipe completed")
    return fp_output


def run_diann(
    job: dict,
    raw_files: list[Path],
    fasta: Path,
    diann_path: str,
    output_dir: Path,
) -> Path:
    """Run DIA-NN search."""
    print("\n" + "=" * 50)
    print("Running DIA-NN")
    print("=" * 50)

    dn_config = job["search_config"]["diann"]
    dn_output = output_dir / "diann"
    dn_output.mkdir(parents=True, exist_ok=True)

    report_path = dn_output / "report.parquet"

    # Build command
    cmd = [
        diann_path,
        "--fasta", str(fasta),
        "--out", str(report_path),
        "--threads", str(dn_config.get("threads", 16)),
        "--qvalue", str(dn_config.get("qvalue", 0.01)),
        "--min-pep-len", str(dn_config.get("min_peptide_length", 7)),
        "--max-pep-len", str(dn_config.get("max_peptide_length", 30)),
        "--predictor",  # Use deep learning
        "--smart-profiling",
        "--no-ifs-removal",
    ]

    # Add raw files
    for raw_file in raw_files:
        cmd.extend(["--f", str(raw_file.absolute())])

    print(f"Command: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"DIA-NN STDERR: {result.stderr}")
        raise RuntimeError(f"DIA-NN failed with code {result.returncode}")

    print("DIA-NN completed")
    return dn_output


def package_results(job: dict, output_dir: Path) -> Path:
    """Package results into final output structure."""
    print("\n" + "=" * 50)
    print("Packaging results")
    print("=" * 50)

    accession = job["accession"]
    results_dir = output_dir / accession
    results_dir.mkdir(parents=True, exist_ok=True)

    # Copy FragPipe PSMs
    fp_psm = output_dir / "fragpipe" / "psm.tsv"
    if fp_psm.exists():
        shutil.copy(fp_psm, results_dir / "fragpipe_psm.tsv")

    # Copy DIA-NN report
    dn_report = output_dir / "diann" / "report.parquet"
    if dn_report.exists():
        shutil.copy(dn_report, results_dir / "diann_report.parquet")

    # Write job manifest
    with open(results_dir / "job.yaml", "w") as f:
        yaml.dump(job, f)

    # Write completion marker
    (results_dir / "COMPLETED").touch()

    print(f"Results packaged: {results_dir}")
    return results_dir


def main():
    parser = argparse.ArgumentParser(description="Execute San José job")
    parser.add_argument("--job", required=True, help="Job manifest YAML")
    parser.add_argument("--fragpipe", required=True, help="FragPipe installation path")
    parser.add_argument("--diann", required=True, help="DIA-NN binary path")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--fasta-dir", default="/fastas", help="FASTA cache directory")
    parser.add_argument("--data-dir", default="/data", help="Raw data directory")
    parser.add_argument("--skip-fragpipe", action="store_true", help="Skip FragPipe")
    parser.add_argument("--skip-diann", action="store_true", help="Skip DIA-NN")
    args = parser.parse_args()

    job = load_job(args.job)
    print(f"Job: {job['job_id']}")
    print(f"Accession: {job['accession']}")
    print(f"Organism: {job['organism']}")
    print(f"Raw files: {job['n_raw_files']}")

    output_dir = Path(args.output)
    fasta_dir = Path(args.fasta_dir)
    data_dir = Path(args.data_dir)

    for d in [output_dir, fasta_dir, data_dir]:
        d.mkdir(parents=True, exist_ok=True)

    try:
        # Step 1: Prepare FASTA
        fasta = prepare_fasta(job, fasta_dir)

        # Step 2: Download raw files
        raw_files = download_raw_files(job, data_dir)

        # Step 3: Run searches
        if not args.skip_fragpipe:
            run_fragpipe(job, raw_files, fasta, args.fragpipe, output_dir)

        if not args.skip_diann:
            run_diann(job, raw_files, fasta, args.diann, output_dir)

        # Step 4: Package results
        results = package_results(job, output_dir)

        print("\n" + "=" * 50)
        print("JOB COMPLETED SUCCESSFULLY")
        print(f"Results: {results}")
        print("=" * 50)

    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
