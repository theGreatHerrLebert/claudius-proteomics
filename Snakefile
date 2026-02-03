"""
CLAUDIUS-PROTEOMICS: Peptide Property Prediction Pipeline

Two-phase architecture:
  Phase 1: DATABASE - Build extensible peptide property database
  Phase 2: MODELING - Train models on database snapshots

Usage:
  # Add a dataset to the database
  snakemake --profile profiles/mogon2 add_dataset --config accession=PXD019086

  # Create a database snapshot for training
  snakemake --profile profiles/mogon2 create_snapshot --config version=v1.0

  # Train a model on a snapshot
  snakemake --profile profiles/mogon2 train_model --config snapshot=v1.0 model=ccs_v1
"""

from pathlib import Path
import json

configfile: "config/config.yaml"

# =============================================================================
# PHASE 1: DATABASE BUILDING
# =============================================================================

include: "rules/fasta.smk"
include: "rules/download.smk"
include: "rules/fragpipe.smk"
include: "rules/diann.smk"
include: "rules/sage.smk"
include: "rules/extract.smk"
include: "rules/database.smk"

# =============================================================================
# PHASE 2: MODEL TRAINING
# =============================================================================

include: "rules/snapshot.smk"
include: "rules/train.smk"


# =============================================================================
# HIGH-LEVEL TARGETS
# =============================================================================

# Default: show help
rule help:
    run:
        print("""
CLAUDIUS-PROTEOMICS Pipeline

Database Building (Phase 1):
  snakemake add_dataset --config accession=PXD019086
      Add a single dataset to the database

  snakemake add_all_datasets
      Add all configured datasets to the database

  snakemake create_snapshot --config version=v1.0
      Create a versioned snapshot for training

Processing:
  snakemake process_full
      Run extraction + all search engines (RECOMMENDED)
      Extraction first (IM calibration), then triple search

  snakemake extract_only
      Run raw feature extraction with IM calibration (imspy)

Search Methods (orthogonal validation):
  snakemake process_only
      Run FragPipe/MSFragger (spectrum-centric)

  snakemake process_diann
      Run DIA-NN (peptide-centric)

  snakemake process_sage
      Run Sage (fast Rust-based search)

  snakemake process_both
      Run FragPipe + DIA-NN

  snakemake process_all
      Run all three engines (triple orthogonal validation)

Model Training (Phase 2):
  snakemake train_model --config snapshot=v1.0 model=ccs_v1
      Train a model on a database snapshot

  snakemake evaluate_model --config snapshot=v1.0 model=ccs_v1
      Evaluate a trained model

Utilities:
  snakemake db_stats
      Show database statistics

  snakemake list_snapshots
      List available snapshots

  snakemake prepare_fasta --config accession=PXD019086
      Prepare FASTA database for a dataset (organism + contaminants)

  snakemake list_fasta_databases
      List available FASTA databases

Use --profile profiles/mogon2 for HPC execution.
        """)


# -----------------------------------------------------------------------------
# Database targets
# -----------------------------------------------------------------------------

rule add_dataset:
    """Add a single dataset to the database."""
    input:
        lambda wildcards: f"database/peptides/{config.get('accession', 'PXD019086')}/manifest.json"


rule add_all_datasets:
    """Add all configured datasets to the database."""
    input:
        expand(
            "database/peptides/{accession}/manifest.json",
            accession=config.get("datasets", ["PXD019086"])
        )


rule create_snapshot:
    """Create a versioned snapshot for training."""
    input:
        lambda wildcards: f"database/snapshots/{config.get('version', 'v1.0')}/training_set.parquet"


rule db_stats:
    """Show database statistics."""
    input:
        "database/metadata.json"
    run:
        with open(input[0]) as f:
            metadata = json.load(f)
        print("\n=== Database Statistics ===")
        print(f"Version: {metadata.get('version', 'unknown')}")
        print(f"Datasets: {metadata.get('n_datasets', 0)}")
        print(f"Total peptides: {metadata.get('n_peptides', 0)}")
        print(f"Total PSMs: {metadata.get('n_psms', 0)}")
        print(f"Last updated: {metadata.get('last_updated', 'unknown')}")
        print("===========================\n")


rule list_snapshots:
    """List available database snapshots."""
    run:
        snapshot_dir = Path("database/snapshots")
        if snapshot_dir.exists():
            snapshots = [d.name for d in snapshot_dir.iterdir() if d.is_dir()]
            print("\n=== Available Snapshots ===")
            for s in sorted(snapshots):
                manifest = snapshot_dir / s / "manifest.json"
                if manifest.exists():
                    with open(manifest) as f:
                        info = json.load(f)
                    print(f"  {s}: {info.get('n_peptides', '?')} peptides, {info.get('created', '?')}")
                else:
                    print(f"  {s}: (no manifest)")
            print("===========================\n")
        else:
            print("No snapshots found.")


# -----------------------------------------------------------------------------
# Model targets
# -----------------------------------------------------------------------------

rule train_model:
    """Train a model on a database snapshot."""
    input:
        lambda wildcards: f"models/{config.get('model', 'ccs_v1')}/metrics.json"


rule evaluate_model:
    """Evaluate a trained model."""
    input:
        lambda wildcards: f"models/{config.get('model', 'ccs_v1')}/evaluation_report.html"


# -----------------------------------------------------------------------------
# Utility targets
# -----------------------------------------------------------------------------

rule clean_temp:
    """Remove temporary files (keeps database and models)."""
    shell:
        """
        rm -rf data/raw/*/
        rm -rf data/processed/*/
        rm -rf data/extracted/*/
        rm -rf .snakemake/
        echo "Temporary files cleaned."
        """


rule download_only:
    """Download raw data without processing."""
    input:
        expand(
            "data/raw/{accession}",
            accession=config.get("datasets", ["PXD019086"])
        )


rule process_only:
    """Run FragPipe on all datasets."""
    input:
        expand(
            "data/processed/{accession}/combined_peptide.tsv",
            accession=config.get("datasets", ["PXD019086"])
        )


rule process_diann:
    """Run DIA-NN on all datasets (peptide-centric search)."""
    input:
        expand(
            "data/processed/{accession}/diann/report.parquet",
            accession=config.get("datasets", ["PXD019086"])
        )


rule process_both:
    """Run both FragPipe and DIA-NN for orthogonal validation."""
    input:
        # Spectrum-centric (FragPipe/MSFragger)
        expand(
            "data/processed/{accession}/combined_peptide.tsv",
            accession=config.get("datasets", ["PXD019086"])
        ),
        # Peptide-centric (DIA-NN)
        expand(
            "data/processed/{accession}/diann/report.parquet",
            accession=config.get("datasets", ["PXD019086"])
        )


rule process_sage:
    """Run Sage on all datasets (fast Rust-based search)."""
    input:
        expand(
            "data/processed/{accession}/sage/results.sage.parquet",
            accession=config.get("datasets", ["PXD019086"])
        )


rule process_all:
    """Run all three search engines for triple orthogonal validation."""
    input:
        # Spectrum-centric (FragPipe/MSFragger)
        expand(
            "data/processed/{accession}/combined_peptide.tsv",
            accession=config.get("datasets", ["PXD019086"])
        ),
        # Peptide-centric (DIA-NN)
        expand(
            "data/processed/{accession}/diann/report.parquet",
            accession=config.get("datasets", ["PXD019086"])
        ),
        # Fast Rust-based (Sage)
        expand(
            "data/processed/{accession}/sage/results.sage.parquet",
            accession=config.get("datasets", ["PXD019086"])
        )


rule extract_only:
    """Run raw feature extraction with IM calibration (imspy)."""
    input:
        expand(
            "data/extracted/{accession}/raw_features.parquet",
            accession=config.get("datasets", ["PXD019086"])
        )


rule process_full:
    """
    Run extraction + all search engines (recommended workflow).

    Extraction runs first (uses Bruker SDK for accurate IM calibration),
    then all three search engines run for triple orthogonal validation.
    """
    input:
        # Raw extraction with IM calibration (runs first)
        expand(
            "data/extracted/{accession}/raw_features.parquet",
            accession=config.get("datasets", ["PXD019086"])
        ),
        # Spectrum-centric (FragPipe/MSFragger)
        expand(
            "data/processed/{accession}/combined_peptide.tsv",
            accession=config.get("datasets", ["PXD019086"])
        ),
        # Peptide-centric (DIA-NN)
        expand(
            "data/processed/{accession}/diann/report.parquet",
            accession=config.get("datasets", ["PXD019086"])
        ),
        # Fast Rust-based (Sage)
        expand(
            "data/processed/{accession}/sage/results.sage.parquet",
            accession=config.get("datasets", ["PXD019086"])
        )


rule prepare_fasta:
    """Prepare FASTA database for a dataset."""
    input:
        lambda wildcards: f"resources/fasta/search_db/{config.get('accession', 'PXD019086')}.fasta"
