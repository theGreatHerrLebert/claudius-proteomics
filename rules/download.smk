"""
Rules for downloading data from PRIDE Archive
"""


rule download_pride:
    """
    Download raw timsTOF data from PRIDE Archive.

    Uses pride-archive-downloader or similar tool to fetch .d folders.
    """
    output:
        directory("data/raw/{accession}")
    params:
        accession="{accession}"
    singularity:
        config["containers"]["download"]
    resources:
        mem_mb=4000,
        time="8:00:00",  # Downloads can take a while
        disk_mb=500000   # Reserve disk space for large downloads
    log:
        "logs/download/{accession}.log"
    shell:
        """
        python scripts/download_pride.py \
            --accession {params.accession} \
            --output {output} \
            2>&1 | tee {log}
        """
