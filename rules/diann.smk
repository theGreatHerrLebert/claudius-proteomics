"""
Rules for DIA-NN processing of timsTOF DDA data.

DIA-NN 2.3+ supports DDA mode with --dda flag.
Can be used as an alternative or complement to FragPipe.
"""


rule diann_search:
    """
    Run DIA-NN on timsTOF DDA data.

    DIA-NN 2.3+ can process DDA data in library-free mode.
    Generates decoys internally, no need for pre-generated decoys.
    """
    input:
        raw_dir="data/raw/{accession}",
        fasta="resources/fasta/search_db/{accession}.fasta"
    output:
        report="data/processed/{accession}/diann/report.tsv",
        done="data/processed/{accession}/diann/.diann_done"
    params:
        diann_path=config["diann"]["path"],
        threads=config["diann"]["threads"],
        qvalue=config["diann"]["qvalue"],
        max_files=config.get("test_mode", {}).get("max_files", 0) if config.get("test_mode", {}).get("enabled", False) else 0,
        outdir="data/processed/{accession}/diann"
    resources:
        mem_mb=64000,
        time="4:00:00",
        cpus_per_task=config["diann"]["threads"]
    log:
        "logs/diann/{accession}.log"
    run:
        from pathlib import Path
        import subprocess

        # Create output directory
        outdir = Path(params.outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        # Find raw files
        raw_dir = Path(input.raw_dir)
        raw_files = sorted(raw_dir.glob("*.d"))

        if params.max_files > 0:
            raw_files = raw_files[:params.max_files]
            print(f"Test mode: limiting to {params.max_files} files")

        print(f"Found {len(raw_files)} raw files")

        # Build DIA-NN command
        cmd = [
            params.diann_path,
            "--fasta", input.fasta,
            "--out", output.report,
            "--threads", str(params.threads),
            "--qvalue", str(params.qvalue),
            "--dda",  # DDA mode
            "--fasta-search",  # Library-free mode
            "--predictor",  # Use deep learning predictor
            "--met-excision",  # Methionine excision
            "--cut", "K*,R*",  # Trypsin
            "--missed-cleavages", "2",
            "--min-pep-len", "7",
            "--max-pep-len", "30",
            "--var-mods", "2",
        ]

        # Add input files
        for raw_file in raw_files:
            cmd.extend(["--f", str(raw_file)])

        print(f"Running DIA-NN:")
        print(f"  Command: {' '.join(cmd[:20])}...")  # Truncate for display

        # Run DIA-NN
        with open(log[0], 'w') as log_file:
            log_file.write(f"Command: {' '.join(cmd)}\n\n")
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
            print(f"DIA-NN failed with exit code {result.returncode}")
            print(f"Check log: {log[0]}")
            raise RuntimeError(f"DIA-NN failed with exit code {result.returncode}")

        # Mark completion
        Path(output.done).touch()
        print("DIA-NN completed successfully!")


rule diann_search_single:
    """
    Run DIA-NN on a single .d file.
    """
    input:
        raw_file="data/raw/{accession}/{sample}.d",
        fasta="resources/fasta/search_db/{accession}.fasta"
    output:
        report="data/processed/{accession}/{sample}/diann/report.tsv",
        done="data/processed/{accession}/{sample}/diann/.diann_done"
    params:
        diann_path=config["diann"]["path"],
        threads=config["diann"]["threads"],
        qvalue=config["diann"]["qvalue"],
        outdir="data/processed/{accession}/{sample}/diann"
    resources:
        mem_mb=32000,
        time="1:00:00",
        cpus_per_task=config["diann"]["threads"]
    log:
        "logs/diann/{accession}/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}

        {params.diann_path} \
            --f {input.raw_file} \
            --fasta {input.fasta} \
            --out {output.report} \
            --threads {params.threads} \
            --qvalue {params.qvalue} \
            --dda \
            --fasta-search \
            --predictor \
            --met-excision \
            --cut "K*,R*" \
            --missed-cleavages 2 \
            --min-pep-len 7 \
            --max-pep-len 30 \
            --var-mods 2 \
            2>&1 | tee {log}

        touch {output.done}
        """
