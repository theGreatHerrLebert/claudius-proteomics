"""
Rules for raw feature extraction via imspy/rustdf

This extracts signal-level features directly from timsTOF .d files:
- Ion distributions
- Chromatographic peaks
- Mobilograms
- Isotope patterns
"""


rule extract_raw_features:
    """
    Extract raw signal features from timsTOF .d files using imspy.

    This goes beyond what FragPipe exports, giving access to:
    - Full ion distributions (not just top-k)
    - Raw chromatographic peak shapes
    - Complete mobilograms
    - Isotope envelopes
    """
    input:
        raw_dir="data/raw/{accession}",
        psm="data/processed/{accession}/psm.tsv"  # Need PSMs to know what to extract
    output:
        features="data/extracted/{accession}/raw_features.parquet"
    params:
        mz_tol=config["extraction"]["mz_tolerance_ppm"],
        rt_tol=config["extraction"]["rt_tolerance_sec"],
        mob_tol=config["extraction"]["mobility_tolerance"],
        extract_chrom=config["extraction"]["extract_chromatogram"],
        extract_mob=config["extraction"]["extract_mobilogram"],
        extract_iso=config["extraction"]["extract_isotopes"]
    singularity:
        config["containers"]["imspy"]
    resources:
        mem_mb=32000,
        time="4:00:00",
        cpus_per_task=8
    log:
        "logs/extract/{accession}.log"
    script:
        "../scripts/extract_raw_features.py"
