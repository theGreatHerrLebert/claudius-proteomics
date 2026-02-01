"""
Rules for FragPipe processing of timsTOF data.

FragPipe runs in headless mode to identify PSMs from DDA timsTOF data.
Uses LFQ-MBR workflow by default for label-free quantification with
match-between-runs enabled.
"""


rule fragpipe_search:
    """
    Run FragPipe on timsTOF .d files to generate PSM identifications.

    FragPipe is mounted from user-provided path due to licensing.
    Uses the per-accession FASTA from the fasta.smk rules (organism + contaminants).
    """
    input:
        raw_dir="data/raw/{accession}",
        fasta="resources/fasta/search_db/{accession}_decoys.fasta"
    output:
        peptide="data/processed/{accession}/combined_peptide.tsv",
        ion="data/processed/{accession}/combined_ion.tsv",
        protein="data/processed/{accession}/combined_protein.tsv",
        done="data/processed/{accession}/.fragpipe_done"
    params:
        fragpipe_path=config["fragpipe"]["path"],
        workflow=config["fragpipe"]["workflow"],
        threads=config["fragpipe"]["threads"],
        ram=config["fragpipe"]["memory_gb"],
        max_files=config.get("test_mode", {}).get("max_files", 0) if config.get("test_mode", {}).get("enabled", False) else 0,
        outdir="data/processed/{accession}",
        temp_dir=config.get("temp_dir", "/scratch/timsim/claudius-proteomics/tmp")
    resources:
        mem_mb=config["fragpipe"]["memory_gb"] * 1000,
        time="8:00:00",
        cpus_per_task=config["fragpipe"]["threads"]
    log:
        "logs/fragpipe/{accession}.log"
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}

        # Run FragPipe via Python wrapper
        # San José mode: disable FDR filtering to report ALL PSMs
        FDR_FLAG=""
        if [ "{config[san_jose][report_all]}" != "True" ]; then
            FDR_FLAG="--enable-fdr-filter"
        fi

        python3 scripts/run_fragpipe.py \
            --fragpipe {params.fragpipe_path} \
            --input {input.raw_dir} \
            --output {params.outdir} \
            --fasta {input.fasta} \
            --workflow {params.workflow} \
            --threads {params.threads} \
            --ram {params.ram} \
            --max-files {params.max_files} \
            --temp-dir {params.temp_dir} \
            $FDR_FLAG \
            2>&1 | tee {log}

        # Mark completion
        touch {output.done}
        """


rule fragpipe_search_single:
    """
    Run FragPipe on a single .d file.

    Useful for testing or processing individual files.
    """
    input:
        raw_file="data/raw/{accession}/{sample}.d",
        fasta="resources/fasta/search_db/{accession}.fasta"
    output:
        psm="data/processed/{accession}/{sample}/psm.tsv",
        done="data/processed/{accession}/{sample}/.fragpipe_done"
    params:
        fragpipe_path=config["fragpipe"]["path"],
        workflow=config["fragpipe"]["workflow"],
        threads=config["fragpipe"]["threads"],
        ram=config["fragpipe"]["memory_gb"],
        outdir="data/processed/{accession}/{sample}"
    resources:
        mem_mb=config["fragpipe"]["memory_gb"] * 1000,
        time="2:00:00",
        cpus_per_task=config["fragpipe"]["threads"]
    log:
        "logs/fragpipe/{accession}/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}

        # Create single-file manifest
        echo -e "{input.raw_file}\t{wildcards.sample}\t1\t1\tDDA" > {params.outdir}/manifest.fp-manifest

        # Get workflow path
        WORKFLOW_PATH="{params.fragpipe_path}/workflows/{params.workflow}.workflow"

        # Run FragPipe
        {params.fragpipe_path}/bin/fragpipe \
            --headless \
            --workflow "$WORKFLOW_PATH" \
            --manifest {params.outdir}/manifest.fp-manifest \
            --workdir {params.outdir} \
            --threads {params.threads} \
            --ram {params.ram} \
            2>&1 | tee {log}

        touch {output.done}
        """
