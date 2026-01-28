"""
Rules for Sage search engine processing of timsTOF DDA data.

Sage is a fast, open-source proteomics search engine written in Rust.
Uses timsrust internally for reading Bruker .d files.
"""

import json
from pathlib import Path


rule sage_search:
    """
    Run Sage on timsTOF DDA data.

    Sage reads .d files directly via timsrust.
    Generates decoys internally, supports LFQ quantification.
    """
    input:
        raw_dir="data/raw/{accession}",
        fasta="resources/fasta/search_db/{accession}.fasta"
    output:
        results="data/processed/{accession}/sage/results.sage.parquet",
        done="data/processed/{accession}/sage/.sage_done"
    params:
        sage_path=config["sage"]["path"],
        threads=config["sage"].get("threads", 16),
        config_template=config["sage"].get("config", "config/sage_config.json"),
        max_files=config.get("test_mode", {}).get("max_files", 0) if config.get("test_mode", {}).get("enabled", False) else 0,
        outdir="data/processed/{accession}/sage"
    resources:
        mem_mb=32000,
        time="2:00:00",
        cpus_per_task=config["sage"].get("threads", 16)
    log:
        "logs/sage/{accession}.log"
    run:
        import subprocess

        # Create output directory
        outdir = Path(params.outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        # Find raw files (.d directories)
        raw_dir = Path(input.raw_dir)
        raw_files = sorted(raw_dir.glob("*.d"))

        if params.max_files > 0:
            raw_files = raw_files[:params.max_files]
            print(f"Test mode: limiting to {params.max_files} files")

        print(f"Found {len(raw_files)} raw files")

        # Load config template and customize
        with open(params.config_template) as f:
            sage_config = json.load(f)

        # Set FASTA path (absolute)
        fasta_abs = Path(input.fasta).resolve()
        sage_config["database"]["fasta"] = str(fasta_abs)

        # Set input files (absolute paths)
        sage_config["mzml_paths"] = [str(f.resolve()) for f in raw_files]

        # Write customized config (use absolute path)
        config_path = (outdir / "sage_config.json").resolve()
        with open(config_path, 'w') as f:
            json.dump(sage_config, f, indent=2)

        # Build Sage command
        outdir_abs = outdir.resolve()
        cmd = [
            params.sage_path,
            str(config_path),  # Absolute path
            "-o", str(outdir_abs),
            "--parquet",  # Output in parquet format
            "--write-pin",  # Write percolator input for rescoring
        ]

        print(f"Running Sage:")
        print(f"  Config: {config_path}")
        print(f"  Output: {outdir_abs}")
        print(f"  Files: {len(raw_files)}")

        # Run Sage
        with open(log[0], 'w') as log_file:
            log_file.write(f"Command: {' '.join(cmd)}\n")
            log_file.write(f"Config:\n{json.dumps(sage_config, indent=2)}\n\n")
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                cwd=str(outdir)
            )
            log_file.write(result.stdout)
            log_file.write(f"\nReturn code: {result.returncode}\n")

        if result.returncode != 0:
            print(f"Sage failed with exit code {result.returncode}")
            print(f"Check log: {log[0]}")
            raise RuntimeError(f"Sage failed with exit code {result.returncode}")

        # Sage outputs results.sage.parquet by default
        expected_output = outdir / "results.sage.parquet"
        if not expected_output.exists():
            # Check for alternative naming
            parquet_files = list(outdir.glob("*.parquet"))
            if parquet_files:
                print(f"Found parquet files: {parquet_files}")
                # Rename first one to expected name
                parquet_files[0].rename(expected_output)
            else:
                raise RuntimeError(f"Sage did not produce expected output: {expected_output}")

        # Mark completion
        Path(output.done).touch()
        print("Sage completed successfully!")
