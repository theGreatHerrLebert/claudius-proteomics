"""
Rules for merging FragPipe PSMs with raw extracted features
"""


rule merge_psm_raw:
    """
    Merge FragPipe PSM identifications with raw extracted features.

    Aligns by precursor ID to create a unified dataset with:
    - Peptide sequence and modifications
    - CCS values
    - Raw signal features (chromatograms, mobilograms, isotopes)
    """
    input:
        psm="data/processed/{accession}/psm.tsv",
        raw_features="data/extracted/{accession}/raw_features.parquet"
    output:
        merged="data/merged/{accession}/peptides_full.parquet"
    singularity:
        config["containers"]["imspy"]
    resources:
        mem_mb=16000,
        time="1:00:00"
    log:
        "logs/merge/{accession}.log"
    script:
        "../scripts/merge_psm_raw.py"


rule combine_datasets:
    """
    Combine multiple merged datasets into a single training set.

    Handles deduplication and alignment across studies.
    """
    input:
        datasets=expand(
            "data/merged/{accession}/peptides_full.parquet",
            accession=config["datasets"]
        )
    output:
        combined="data/datasets/combined_training.parquet"
    singularity:
        config["containers"]["imspy"]
    resources:
        mem_mb=32000,
        time="1:00:00"
    log:
        "logs/combine/combined.log"
    script:
        "../scripts/combine_datasets.py"
