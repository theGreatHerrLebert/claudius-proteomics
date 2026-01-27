"""
Rules for creating database snapshots for model training.

Snapshots are versioned, immutable exports of the database
that can be used for reproducible model training.
"""

from datetime import datetime
import json


rule create_training_snapshot:
    """
    Create a versioned snapshot from the database for training.

    Combines all datasets, applies filters, and creates train/val/test splits.
    """
    input:
        metadata="database/metadata.json",
        peptides=lambda wildcards: expand(
            "database/peptides/{accession}/peptides.parquet",
            accession=config.get("datasets", ["PXD019086"])
        )
    output:
        training_set="database/snapshots/{version}/training_set.parquet",
        manifest="database/snapshots/{version}/manifest.json"
    params:
        test_split=config["training"]["test_split"],
        min_ccs=200,  # Filter range
        max_ccs=1500,
    singularity:
        config["containers"]["imspy"]
    resources:
        mem_mb=32000,
        time="1:00:00"
    log:
        "logs/snapshot/{version}.log"
    script:
        "../scripts/create_snapshot.py"


rule list_snapshot_contents:
    """
    Show contents of a specific snapshot.
    """
    input:
        manifest="database/snapshots/{version}/manifest.json"
    run:
        with open(input.manifest) as f:
            manifest = json.load(f)

        print(f"\n=== Snapshot {wildcards.version} ===")
        print(f"Created: {manifest.get('created', 'unknown')}")
        print(f"Source datasets: {manifest.get('source_datasets', [])}")
        print(f"Total PSMs: {manifest.get('n_psms', 0):,}")
        print(f"Unique peptides: {manifest.get('n_unique_peptides', 0):,}")
        print(f"Train/Val/Test split: {manifest.get('split_ratio', 'unknown')}")
        print(f"CCS range: {manifest.get('ccs_min', '?')} - {manifest.get('ccs_max', '?')} Å²")
        print("==============================\n")
