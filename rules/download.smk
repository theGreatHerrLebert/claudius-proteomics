"""
Rules for downloading data from PRIDE Archive with metadata extraction.

Multi-stage workflow:
1. fetch_pride_metadata - Get metadata from PRIDE API
2. download_raw_files - Download raw files using pridepy
3. extract_raw_metadata - Extract metadata from .d folders (fills gaps)
4. validate_metadata - Human checkpoint for incomplete metadata
"""


rule fetch_pride_metadata:
    """
    Fetch initial metadata from PRIDE REST API.

    This is the first step - can proceed independently of download.
    Extracts: organism, instrument, lab_id, title, DOI, and parses protocol
    text for gradient/column/mode hints.
    """
    output:
        metadata="data/metadata/{accession}/pride_metadata.yaml",
        raw_response="data/metadata/{accession}/raw_api_response.json"
    params:
        accession="{accession}"
    log:
        "logs/download/{accession}_metadata.log"
    shell:
        """
        python scripts/download_pride.py \
            --accession {params.accession} \
            --metadata-only \
            --metadata-dir data/metadata/{params.accession} \
            2>&1 | tee {log}
        """


rule download_raw_files:
    """
    Download raw files from PRIDE (can proceed with partial metadata).

    Uses pridepy for FTP/Aspera download with retry logic.
    Downloads to data/raw/{accession}/.
    """
    input:
        pride_meta="data/metadata/{accession}/pride_metadata.yaml"
    output:
        directory("data/raw/{accession}"),
        flag="data/raw/{accession}/.download_complete"
    params:
        accession="{accession}",
        local_path=lambda wildcards: config.get("local_data", {}).get(wildcards.accession, None),
        protocol=config.get("download", {}).get("protocol", "ftp"),
        retry_count=config.get("download", {}).get("retry_count", 3),
        max_files=config.get("test_mode", {}).get("max_files", 0) if config.get("test_mode", {}).get("enabled", False) else 0
    log:
        "logs/download/{accession}_download.log"
    run:
        import os
        from pathlib import Path

        output_path = Path(output[0])
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if params.local_path and Path(params.local_path).exists():
            # Symlink to local data
            local = Path(params.local_path)
            print(f"Linking local data: {local} -> {output_path}")

            if output_path.exists() or output_path.is_symlink():
                os.remove(output_path)

            os.symlink(local, output_path)
            print(f"Created symlink to {local}")

            # Create download complete flag for symlinked data
            flag_path = Path(output[1])
            flag_path.write_text(f"Symlinked from: {local}\n")

        else:
            # Download from PRIDE
            print(f"Downloading {params.accession} from PRIDE...")
            max_files_arg = f"--max-files {params.max_files}" if params.max_files > 0 else ""
            shell(f"""
                python scripts/download_pride.py \
                    --accession {params.accession} \
                    --output {output_path} \
                    --metadata-dir data/metadata/{params.accession} \
                    --protocol {params.protocol} \
                    --retry {params.retry_count} \
                    {max_files_arg} \
                    2>&1 | tee {log}
            """)


rule extract_raw_metadata:
    """
    Extract metadata from downloaded .d folders (fills gaps from PRIDE API).

    Reads:
    - analysis.tdf/GlobalMetadata: instrument, operator, method, datetime
    - SampleInfo.xml: LC/MS method names (gradient length, mode)
    - Frames table: actual run duration
    - chromatography-data.sqlite: LC system type

    Merges with PRIDE metadata, keeping highest-priority values.
    """
    input:
        raw_dir="data/raw/{accession}",
        pride_meta="data/metadata/{accession}/pride_metadata.yaml"
    output:
        metadata="data/metadata/{accession}/metadata.yaml"
    log:
        "logs/download/{accession}_raw_metadata.log"
    shell:
        """
        python scripts/extract_raw_metadata.py \
            --raw-dir {input.raw_dir} \
            --pride-meta {input.pride_meta} \
            --output {output.metadata} \
            2>&1 | tee {log}
        """


rule validate_metadata:
    """
    Human checkpoint - warns if required fields still missing.

    Checks metadata.yaml and reports:
    - Complete: all required fields populated
    - Incomplete: lists missing fields with hints

    Critical missing fields will raise an error to block pipeline.
    Optional missing fields (like column_type) generate warnings only.
    """
    input:
        metadata="data/metadata/{accession}/metadata.yaml"
    output:
        touch("data/metadata/{accession}/.validated")
    run:
        import yaml
        from pathlib import Path

        metadata_path = Path(input.metadata)
        with open(metadata_path) as f:
            metadata = yaml.safe_load(f)

        validation = metadata.get("validation", {})
        complete = validation.get("complete", False)
        missing = validation.get("missing_field_names", [])

        print(f"\n{'=' * 60}")
        print(f"Metadata Validation: {wildcards.accession}")
        print(f"{'=' * 60}")

        if complete:
            print("Status: COMPLETE - All required fields populated")
        else:
            print(f"Status: INCOMPLETE - {len(missing)} field(s) missing")

        # Report field status
        fields = metadata.get("fields", {})
        print(f"\nField Status:")
        for field_name, field_data in fields.items():
            status = field_data.get("status", "unknown")
            value = field_data.get("value")
            icon = {"auto": "✓", "inferred": "~", "manual": "✎", "missing": "✗"}.get(status, "?")
            print(f"  {icon} {field_name}: {value if value else f'[{status}]'}")

        # Critical fields that block pipeline
        critical_fields = ["organism", "instrument"]
        critical_missing = [f for f in critical_fields if f in missing]

        if critical_missing:
            print(f"\n⚠️  CRITICAL FIELDS MISSING: {', '.join(critical_missing)}")
            print(f"   Edit {input.metadata} to add these fields before continuing.")
            raise ValueError(f"Critical metadata fields missing: {critical_missing}")

        # Optional fields that generate warnings
        if missing:
            print(f"\n⚠️  Optional fields missing: {', '.join(missing)}")
            print(f"   Consider editing {input.metadata} to improve metadata quality.")

            # Show hints for missing fields
            for field_name in missing:
                field_data = fields.get(field_name, {})
                hint = field_data.get("hint")
                if hint:
                    print(f"   - {field_name}: {hint}")

        print(f"\nValidation checkpoint passed.")


rule download_or_link:
    """
    DEPRECATED: Use download_raw_files instead.

    This rule is kept for backwards compatibility but redirects to the new workflow.
    """
    output:
        directory("data/raw_legacy/{accession}")
    params:
        accession="{accession}",
        local_path=lambda wildcards: config.get("local_data", {}).get(wildcards.accession, None)
    log:
        "logs/download/{accession}_legacy.log"
    run:
        import os
        from pathlib import Path

        output_path = Path(output[0])
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if params.local_path and Path(params.local_path).exists():
            local = Path(params.local_path)
            print(f"Linking local data: {local} -> {output_path}")

            if output_path.exists() or output_path.is_symlink():
                os.remove(output_path)

            os.symlink(local, output_path)
        else:
            print(f"DEPRECATED: Use 'snakemake download_raw_files --config accession={params.accession}' instead")
            shell(f"""
                python scripts/download_pride.py \
                    --accession {params.accession} \
                    --output {output_path} \
                    2>&1 | tee {log}
            """)


rule list_raw_files:
    """
    List all .d files available for an accession.
    """
    input:
        "data/raw/{accession}"
    output:
        "data/raw/{accession}/.file_list.txt"
    params:
        max_files=config.get("test_mode", {}).get("max_files", 0) if config.get("test_mode", {}).get("enabled", False) else 0
    run:
        from pathlib import Path

        raw_dir = Path(input[0])
        d_files = sorted(raw_dir.glob("*.d"))

        # Apply test mode limit if set
        if params.max_files > 0:
            d_files = d_files[:params.max_files]
            print(f"Test mode: limiting to {params.max_files} files")

        print(f"Found {len(d_files)} .d files")

        with open(output[0], 'w') as f:
            for d_file in d_files:
                f.write(f"{d_file}\n")


# Convenience targets
rule prepare_dataset:
    """
    Prepare a dataset for processing: metadata + download + extract + validate.

    This rule prepares raw data and metadata but doesn't run search engines.
    Use add_dataset in main Snakefile for full pipeline.

    Usage:
        snakemake prepare_dataset --config accession=PXD019086
    """
    input:
        raw="data/raw/{accession}",
        metadata="data/metadata/{accession}/metadata.yaml",
        validated="data/metadata/{accession}/.validated"


rule fetch_all_metadata:
    """
    Fetch PRIDE metadata for all configured datasets.
    """
    input:
        expand("data/metadata/{accession}/pride_metadata.yaml",
               accession=config.get("datasets", []))


rule download_all_datasets:
    """
    Download raw files for all configured datasets.
    """
    input:
        expand("data/raw/{accession}/.download_complete",
               accession=config.get("datasets", []))


rule extract_all_metadata:
    """
    Extract and merge metadata for all datasets with downloaded raw files.
    """
    input:
        expand("data/metadata/{accession}/metadata.yaml",
               accession=config.get("datasets", []))


rule validate_all_metadata:
    """
    Validate metadata for all datasets.
    """
    input:
        expand("data/metadata/{accession}/.validated",
               accession=config.get("datasets", []))
