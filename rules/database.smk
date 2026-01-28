"""
Rules for managing the peptide property database.

The database is the core "resource" - a collection of processed datasets
that can be extended over time and used for model training.
"""

from datetime import datetime
import json


rule merge_to_database:
    """
    Merge extracted data and insert into database.

    This combines FragPipe PSMs with raw extracted features
    and stores them in the database structure.
    """
    input:
        peptides="data/processed/{accession}/combined_peptide.tsv",
        raw_features="data/extracted/{accession}/raw_features.parquet"
    output:
        peptides="database/peptides/{accession}/peptides.parquet",
        manifest="database/peptides/{accession}/manifest.json"
    singularity:
        config["containers"]["imspy"]
    resources:
        mem_mb=16000,
        time="1:00:00"
    log:
        "logs/database/{accession}.log"
    script:
        "../scripts/insert_to_database.py"


rule update_database_metadata:
    """
    Update database-level metadata after adding a dataset.

    Aggregates statistics from all datasets.
    """
    input:
        manifests=lambda wildcards: expand(
            "database/peptides/{accession}/manifest.json",
            accession=config.get("datasets", ["PXD019086"])
        )
    output:
        metadata="database/metadata.json"
    run:
        import json
        from pathlib import Path
        from datetime import datetime

        # Aggregate stats from all manifests
        total_peptides = 0
        total_psms = 0
        datasets = []

        for manifest_path in input.manifests:
            with open(manifest_path) as f:
                manifest = json.load(f)
            total_peptides += manifest.get("n_unique_peptides", 0)
            total_psms += manifest.get("n_psms", 0)
            datasets.append(manifest.get("accession", Path(manifest_path).parent.name))

        metadata = {
            "version": "1.0",
            "last_updated": datetime.now().isoformat(),
            "n_datasets": len(datasets),
            "datasets": datasets,
            "n_peptides": total_peptides,
            "n_psms": total_psms,
            "schema_version": "1.0",
        }

        Path(output.metadata).parent.mkdir(parents=True, exist_ok=True)
        with open(output.metadata, "w") as f:
            json.dump(metadata, f, indent=2)


rule init_database:
    """
    Initialize empty database structure.
    """
    output:
        schema="database/schema.json"
    run:
        import json
        from pathlib import Path

        schema = {
            "version": "1.0",
            "columns": {
                "sequence": {"type": "string", "description": "Peptide sequence (unmodified)"},
                "modified_sequence": {"type": "string", "description": "Peptide sequence with modifications"},
                "charge": {"type": "int32", "description": "Charge state"},
                "mz": {"type": "float64", "description": "Precursor m/z"},
                "mass": {"type": "float64", "description": "Peptide mass"},
                "retention_time": {"type": "float64", "description": "Retention time (seconds)"},
                "mobility": {"type": "float64", "description": "Ion mobility (1/K0)"},
                "ccs": {"type": "float64", "description": "Collisional cross section (Å²)"},
                "raw_file": {"type": "string", "description": "Source raw file"},
                "accession": {"type": "string", "description": "PRIDE accession"},
                # Raw features
                "xic_rt": {"type": "list[float64]", "description": "XIC retention times"},
                "xic_intensity": {"type": "list[float64]", "description": "XIC intensities"},
                "mobilogram_mobility": {"type": "list[float64]", "description": "Mobilogram mobility values"},
                "mobilogram_intensity": {"type": "list[float64]", "description": "Mobilogram intensities"},
                "isotope_mz": {"type": "list[float64]", "description": "Isotope m/z values"},
                "isotope_intensity": {"type": "list[float64]", "description": "Isotope intensities"},
            }
        }

        Path(output.schema).parent.mkdir(parents=True, exist_ok=True)
        with open(output.schema, "w") as f:
            json.dump(schema, f, indent=2)
