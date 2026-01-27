"""
Rules for FragPipe processing of timsTOF data
"""


rule fragpipe_search:
    """
    Run FragPipe on timsTOF .d files to generate PSM identifications.

    FragPipe is mounted from user-provided path due to licensing.
    The container only provides Java and dependencies.

    Uses the per-accession FASTA from the fasta.smk rules (organism + contaminants).
    """
    input:
        raw_dir="data/raw/{accession}",
        fasta="resources/fasta/search_db/{accession}.fasta"
    output:
        psm="data/processed/{accession}/psm.tsv",
        ion="data/processed/{accession}/ion.tsv",
        done="data/processed/{accession}/.fragpipe_done"
    params:
        fragpipe_path=config["fragpipe"]["path"],
        workflow=config["fragpipe"]["workflow"],
        threads=config["fragpipe"]["threads"],
        outdir="data/processed/{accession}"
    singularity:
        config["containers"]["fragpipe_base"]
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

        # Run FragPipe via mounted path
        bash scripts/run_fragpipe.sh \
            --fragpipe {params.fragpipe_path} \
            --input {input.raw_dir} \
            --output {params.outdir} \
            --fasta {input.fasta} \
            --workflow {params.workflow} \
            --threads {params.threads} \
            2>&1 | tee {log}

        # Mark completion
        touch {output.done}
        """
