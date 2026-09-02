#!/usr/bin/env python3
"""
Generate a job manifest for a PRIDE dataset.

Queries PRIDE API to get dataset metadata and file list,
determines organism, and creates a complete job manifest.
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import requests
import yaml

# Organism mapping based on common PRIDE tags
ORGANISM_MAP = {
    "homo sapiens": "human",
    "human": "human",
    "mus musculus": "mouse",
    "mouse": "mouse",
    "saccharomyces cerevisiae": "yeast",
    "yeast": "yeast",
    "escherichia coli": "ecoli",
    "e. coli": "ecoli",
}

# FASTA URLs per organism (UniProt reference proteomes)
FASTA_SOURCES = {
    "human": {
        "proteome_id": "UP000005640",
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz",
    },
    "mouse": {
        "proteome_id": "UP000000589",
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000589/UP000000589_10090.fasta.gz",
    },
    "yeast": {
        "proteome_id": "UP000002311",
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002311/UP000002311_559292.fasta.gz",
    },
    "ecoli": {
        "proteome_id": "UP000000625",
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000625/UP000000625_83333.fasta.gz",
    },
}

PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v3"


def get_pride_metadata(accession: str) -> dict:
    """Fetch dataset metadata from PRIDE API."""
    url = f"{PRIDE_API_BASE}/projects/{accession}"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    return response.json()


def get_pride_files(accession: str) -> list[dict]:
    """Fetch the full file list from the PRIDE API, following pagination.

    v3 returns a bare JSON list per page, not the HAL ``_embedded.files``
    envelope the retired v2 files/byProject endpoint used.
    """
    files: list[dict] = []
    page = 0
    page_size = 1000
    while True:
        url = f"{PRIDE_API_BASE}/projects/{accession}/files"
        params = {"pageSize": page_size, "page": page}
        response = requests.get(url, params=params, timeout=60)
        response.raise_for_status()
        batch = response.json()
        files.extend(batch)
        if len(batch) < page_size:
            break
        page += 1
    return files


def detect_organism(metadata: dict) -> str:
    """Detect organism from PRIDE metadata."""
    organisms = metadata.get("organisms", [])
    if organisms:
        org_name = organisms[0].get("name", "").lower()
        for key, value in ORGANISM_MAP.items():
            if key in org_name:
                return value

    # Fallback: check title and description
    text = (metadata.get("title", "") + " " + metadata.get("projectDescription", "")).lower()
    for key, value in ORGANISM_MAP.items():
        if key in text:
            return value

    return "unknown"


def filter_raw_files(files: list[dict]) -> list[dict]:
    """Filter for timsTOF .d files or other raw formats."""
    raw_files = []
    for f in files:
        name = f.get("fileName", "")
        # timsTOF .d folders are often zipped
        if name.endswith(".d.zip") or name.endswith(".d.tar"):
            raw_files.append(f)
        # Thermo .raw
        elif name.endswith(".raw"):
            raw_files.append(f)
        # Bruker .d (sometimes listed as folder)
        elif name.endswith(".d"):
            raw_files.append(f)
    return raw_files


def generate_job_manifest(accession: str) -> dict:
    """Generate complete job manifest for a dataset."""
    print(f"Fetching metadata for {accession}...")
    metadata = get_pride_metadata(accession)

    print(f"Fetching file list...")
    files = get_pride_files(accession)
    raw_files = filter_raw_files(files)

    if not raw_files:
        raise ValueError(f"No raw files found for {accession}")

    organism = detect_organism(metadata)
    print(f"Detected organism: {organism}")

    if organism == "unknown":
        print("WARNING: Could not detect organism, defaulting to human")
        organism = "human"

    # Build file URLs
    file_entries = []
    for f in raw_files:
        # PRIDE FTP structure
        ftp_url = f.get("publicFileLocations", [{}])[0].get("value", "")
        if ftp_url:
            file_entries.append({
                "name": f.get("fileName"),
                "url": ftp_url,
                "size_mb": f.get("fileSizeBytes", 0) / 1024 / 1024,
            })

    # Get FASTA info
    fasta_info = FASTA_SOURCES.get(organism, FASTA_SOURCES["human"])

    job = {
        "job_id": f"{accession}_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        "accession": accession,
        "title": metadata.get("title", ""),
        "organism": organism,
        "instrument": extract_instrument(metadata),
        "created_at": datetime.now().isoformat(),

        "fasta": {
            "organism": organism,
            "proteome_id": fasta_info["proteome_id"],
            "url": fasta_info["url"],
            "include_contaminants": True,
            "generate_decoys": True,
        },

        "raw_files": file_entries,
        "n_raw_files": len(file_entries),

        "search_config": {
            "fragpipe": {
                "workflow": "LFQ-MBR",
                "threads": 16,
                "memory_gb": 64,
            },
            "diann": {
                "threads": 16,
                "qvalue": 0.01,
                "min_peptide_length": 7,
                "max_peptide_length": 30,
            },
        },

        "extraction_config": {
            "mz_tolerance_ppm": 20,
            "rt_tolerance_sec": 30,
            "mobility_tolerance": 0.05,
            "extract_xic": True,
            "extract_mobilogram": True,
            "extract_isotopes": True,
        },

        "output": {
            "format": "parquet",
            "compress": True,
        },
    }

    return job


def extract_instrument(metadata: dict) -> str:
    """Extract instrument type from metadata."""
    instruments = metadata.get("instruments", [])
    if instruments:
        return instruments[0].get("name", "unknown")
    return "unknown"


def main():
    parser = argparse.ArgumentParser(description="Generate job manifest for PRIDE dataset")
    parser.add_argument("--accession", required=True, help="PRIDE accession (e.g., PXD019086)")
    parser.add_argument("--output", required=True, help="Output YAML file path")
    parser.add_argument("--max-files", type=int, default=0, help="Limit number of raw files (0=all)")
    args = parser.parse_args()

    try:
        job = generate_job_manifest(args.accession)

        # Optionally limit files for testing
        if args.max_files > 0:
            job["raw_files"] = job["raw_files"][:args.max_files]
            job["n_raw_files"] = len(job["raw_files"])

        # Write manifest
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as f:
            yaml.dump(job, f, default_flow_style=False, sort_keys=False)

        print(f"Job manifest written to: {output_path}")
        print(f"  Accession: {job['accession']}")
        print(f"  Organism: {job['organism']}")
        print(f"  Raw files: {job['n_raw_files']}")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
