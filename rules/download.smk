"""
Rules for downloading data from PRIDE Archive or linking local data.
"""


rule download_or_link:
    """
    Download raw timsTOF data from PRIDE Archive, or symlink local data.

    If local_data path is configured for the accession, creates a symlink.
    Otherwise, downloads from PRIDE.
    """
    output:
        directory("data/raw/{accession}")
    params:
        accession="{accession}",
        local_path=lambda wildcards: config.get("local_data", {}).get(wildcards.accession, None)
    log:
        "logs/download/{accession}.log"
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

        else:
            # Download from PRIDE
            print(f"Downloading {params.accession} from PRIDE...")
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
