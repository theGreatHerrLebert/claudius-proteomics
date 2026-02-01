"""
Rules for raw feature extraction via imspy

Extracts precursor+fragment data directly from timsTOF .d files:
- Full 4D fragment spectra (merged if re-fragmented)
- Precursor metadata (m/z, charge, isolation window)
- Ion mobility (1/K0 apex) - with accurate Bruker SDK calibration via lookup tables
- Collision energies

Output: index.parquet + blobs.bin per raw file

Ion Mobility Calibration:
The Bruker SDK provides accurate scan→1/K0 conversion but is not thread-safe.
We extract calibration lookup tables once per raw file, enabling fast parallel
extraction with accurate mobility values.
"""

import os


def get_raw_files(wildcards):
    """Get all .d folders for an accession."""
    raw_dir = f"data/raw/{wildcards.accession}"
    if os.path.exists(raw_dir):
        return [f for f in os.listdir(raw_dir) if f.endswith(".d")]
    return []


rule extract_im_calibration:
    """
    Extract ion mobility calibration lookup table from a .d file.

    Uses the Bruker SDK (single-threaded) to probe the exact scan→1/K0 mapping,
    then saves it as a numpy array for fast parallel extraction later.

    This is a one-time extraction per raw file (~1-2 seconds).
    The resulting lookup table is small (~8KB for 1000 scans).
    """
    input:
        raw_file="data/raw/{accession}/{raw_file}"
    output:
        calibration="data/raw/{accession}/{raw_file}.im_calibration.npy"
    params:
        # Configuration from config.yaml
        use_bruker_sdk=config.get("extraction", {}).get("use_bruker_sdk", True),
    resources:
        mem_mb=4000,
        time="0:05:00",
        cpus_per_task=1  # Single-threaded for Bruker SDK
    log:
        "logs/calibration/{accession}/{raw_file}.log"
    run:
        import sys
        sys.path.insert(0, "scripts")

        from pathlib import Path
        from extract_calibration import extract_and_save_calibration

        if not params.use_bruker_sdk:
            # If not using SDK, create a dummy calibration file to satisfy dependencies
            # The extraction will use linear interpolation instead
            import numpy as np
            print(f"Bruker SDK disabled - creating dummy calibration file")
            np.save(output.calibration, np.array([]))
        else:
            print(f"Extracting IM calibration for: {input.raw_file}")
            extract_and_save_calibration(
                str(input.raw_file),
                str(output.calibration),
                verbose=True,
            )
            print(f"Calibration saved to: {output.calibration}")


rule extract_all_precursors:
    """
    Extract precursors from all raw files for an accession.
    """
    input:
        lambda wildcards: expand(
            "data/extracted/{accession}/{raw_file}/index.parquet",
            accession=wildcards.accession,
            raw_file=get_raw_files(wildcards),
        )
    output:
        touch("data/extracted/{accession}/.done")


rule extract_precursors:
    """
    Extract precursors with fragment spectra from a single .d file using imspy.

    Uses pre-computed IM calibration for accurate mobility values with fast parallel extraction.

    Outputs:
    - index.parquet: Lightweight precursor metadata (m/z, charge, RT, mobility, n_peaks)
    - blobs.bin: Serialized TimsFrame data for each precursor (compressed)
    """
    input:
        raw_file="data/raw/{accession}/{raw_file}",
        calibration="data/raw/{accession}/{raw_file}.im_calibration.npy"
    output:
        index="data/extracted/{accession}/{raw_file}/index.parquet",
        blobs="data/extracted/{accession}/{raw_file}/blobs.bin"
    params:
        output_dir="data/extracted/{accession}/{raw_file}",
        threads=config.get("extraction", {}).get("threads", 4),
        use_calibration=config.get("extraction", {}).get("use_bruker_sdk", True),
    resources:
        mem_mb=32000,
        time="2:00:00",
        cpus_per_task=8
    log:
        "logs/extract/{accession}/{raw_file}.log"
    run:
        import sys
        sys.path.insert(0, "scripts")

        from pathlib import Path
        from extract_precursors import (
            TimsDatasetDDA,
            extract_precursors,
            write_extraction,
            setup_logging,
        )

        logger = setup_logging(Path(log[0]))

        logger.info("=" * 60)
        logger.info("Precursor extraction with imspy")
        logger.info("=" * 60)
        logger.info(f"Input: {input.raw_file}")
        logger.info(f"Output: {params.output_dir}")
        logger.info(f"Using calibration: {params.use_calibration}")

        # Load dataset with calibration if available
        import numpy as np
        logger.info("Loading dataset...")

        if params.use_calibration and Path(input.calibration).exists():
            cal_data = np.load(input.calibration)
            if len(cal_data) > 0:
                # Use calibrated dataset (accurate + thread-safe)
                from imspy_connector.py_dda import PyTimsDatasetDDA as PyTimsDatasetDDARust
                rust_dataset = PyTimsDatasetDDARust.with_calibration(
                    str(input.raw_file), False, cal_data.tolist()
                )
                logger.info("  Using pre-computed IM calibration (accurate + fast)")
                # Also load Python wrapper for metadata access
                dataset = TimsDatasetDDA(str(input.raw_file), in_memory=False, use_bruker_sdk=False)
            else:
                # Dummy calibration file - use linear interpolation
                logger.info("  No calibration data - using linear interpolation")
                dataset = TimsDatasetDDA(str(input.raw_file), in_memory=False, use_bruker_sdk=False)
                rust_dataset = None
        else:
            # Fall back to SDK or linear interpolation
            dataset = TimsDatasetDDA(str(input.raw_file), in_memory=False, use_bruker_sdk=False)
            rust_dataset = None

        logger.info(f"  Frames: {dataset.frame_count}")
        logger.info(f"  Fragmented precursors: {len(dataset.fragmented_precursors)}")

        # Extract
        raw_file_name = Path(input.raw_file).name
        precursors = extract_precursors(
            dataset=dataset,
            raw_file_name=raw_file_name,
            num_threads=params.threads,
            logger=logger,
        )

        # Write
        stats = write_extraction(
            precursors=precursors,
            output_dir=Path(params.output_dir),
            write_blobs=True,
            logger=logger,
        )

        logger.info("=" * 60)
        logger.info(f"Extraction complete: {stats['n_precursors']} precursors")
        logger.info(f"Blob size: {stats['blob_size_bytes'] / 1024 / 1024:.1f} MB")
        logger.info("=" * 60)


rule merge_extracted_precursors:
    """
    Merge all extracted precursor indices into a single parquet per accession.
    """
    input:
        indices=lambda wildcards: expand(
            "data/extracted/{accession}/{raw_file}/index.parquet",
            accession=wildcards.accession,
            raw_file=get_raw_files(wildcards),
        )
    output:
        merged="data/extracted/{accession}/precursors.parquet"
    run:
        import pandas as pd

        dfs = []
        for index_file in input.indices:
            df = pd.read_parquet(index_file)
            dfs.append(df)

        if dfs:
            merged = pd.concat(dfs, ignore_index=True)
            merged.to_parquet(output.merged, index=False)
            print(f"Merged {len(merged)} precursors from {len(dfs)} files")
        else:
            # Empty output
            pd.DataFrame().to_parquet(output.merged, index=False)


rule build_training_spectra:
    """
    Build training-ready spectra with normalization and fragment annotation.

    Takes a precursor store (with fragment_mz/fragment_intensity columns) and adds:
    - Normalized intensities (configurable method)
    - Fragment ion annotations (b/y ions matched to theoretical)
    - Annotation quality metrics (intensity explained, coverage)

    Uses build_training_spectra.py which combines:
    - spectrum_normalization.py for intensity normalization
    - annotate_spectrum.py for b/y ion annotation
    - fragment_matching.py (imspy) for theoretical fragment generation
    """
    input:
        store="data/processed/{accession}/precursors.parquet"
    output:
        training="data/processed/{accession}/training_spectra.parquet"
    params:
        normalization=config.get("training", {}).get("normalization", "base_peak"),
        mz_tolerance=config.get("training", {}).get("mz_tolerance_ppm", 20.0),
        min_intensity=config.get("training", {}).get("min_intensity", 0.0),
        batch_size=config.get("training", {}).get("batch_size", 1000),
    resources:
        mem_mb=16000,
        time="1:00:00",
        cpus_per_task=4
    log:
        "logs/training/{accession}/build_training_spectra.log"
    run:
        import sys
        sys.path.insert(0, "scripts")

        from build_training_spectra import build_training_spectra

        print(f"Building training spectra for {wildcards.accession}")
        print(f"  Input:  {input.store}")
        print(f"  Output: {output.training}")
        print(f"  Normalization: {params.normalization}")
        print(f"  m/z tolerance: {params.mz_tolerance} ppm")

        stats = build_training_spectra(
            input_path=input.store,
            output_path=output.training,
            normalization=params.normalization,
            mz_tolerance_ppm=params.mz_tolerance,
            min_intensity=params.min_intensity,
            batch_size=params.batch_size,
        )

        print(f"\nTraining spectra built:")
        print(f"  Precursors: {stats.get('n_precursors', 0):,}")
        print(f"  Annotated:  {stats.get('n_annotated', 0):,}")
        print(f"  Mean intensity explained: {stats.get('mean_intensity_explained', 0):.1%}")
