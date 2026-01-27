#!/usr/bin/env python3
"""
Extract raw signal features from timsTOF .d files using imspy/rustdf.

This extracts features beyond what FragPipe provides:
- Full ion distributions
- Chromatographic peak shapes (XICs)
- Mobilograms
- Isotope patterns
"""

import sys
from pathlib import Path
from typing import List, Dict, Any, Optional
import logging

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# imspy imports for timsTOF data handling
from imspy.timstof import TimsDatasetDDA
from imspy.timstof.frame import TimsFrame
from imspy.core import MzSpectrum

# Snakemake provides these variables
raw_dir = Path(snakemake.input.raw_dir)
psm_file = Path(snakemake.input.psm)
output_file = Path(snakemake.output.features)
log_file = Path(snakemake.log[0])

# Parameters from config
mz_tolerance_ppm = snakemake.params.mz_tol
rt_tolerance_sec = snakemake.params.rt_tol
mobility_tolerance = snakemake.params.mob_tol
extract_chromatogram = snakemake.params.extract_chrom
extract_mobilogram = snakemake.params.extract_mob
extract_isotopes = snakemake.params.extract_iso


def setup_logging(log_path: Path) -> logging.Logger:
    """Set up logging to file and stdout."""
    log_path.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # File handler
    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.INFO)

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger


def load_psms(psm_path: Path) -> pd.DataFrame:
    """
    Load PSM identifications from FragPipe output.

    Normalizes column names for consistent access.
    """
    logger.info(f"Loading PSMs from {psm_path}")
    psms = pd.read_csv(psm_path, sep='\t')

    # Normalize column names (FragPipe output varies by version)
    column_mapping = {
        'Peptide': 'sequence',
        'peptide': 'sequence',
        'Modified Sequence': 'modified_sequence',
        'Charge': 'charge',
        'charge': 'charge',
        'Retention': 'retention_time',
        'RT': 'retention_time',
        'Calibrated Retention Time': 'retention_time',
        'Ion Mobility': 'mobility',
        'IonMobility': 'mobility',
        '1/K0': 'inv_mobility',
        'CCS': 'ccs',
        'Calculated M/Z': 'mz',
        'Observed M/Z': 'mz_observed',
        'm/z': 'mz',
        'Calculated Peptide Mass': 'mass',
        'Spectrum': 'spectrum_id',
        'Spectrum File': 'raw_file',
    }

    psms = psms.rename(columns={k: v for k, v in column_mapping.items() if k in psms.columns})

    logger.info(f"Loaded {len(psms)} PSMs")
    logger.info(f"Columns: {list(psms.columns)}")

    return psms


def find_raw_files(raw_dir: Path) -> Dict[str, TimsDatasetDDA]:
    """
    Find and load all .d folders in the raw directory.

    Returns a dict mapping raw file stem to TimsDatasetDDA object.
    """
    d_folders = list(raw_dir.glob("*.d"))
    logger.info(f"Found {len(d_folders)} .d folders in {raw_dir}")

    datasets = {}
    for d_folder in d_folders:
        try:
            logger.info(f"Loading {d_folder.name}...")
            dataset = TimsDatasetDDA(str(d_folder))
            datasets[d_folder.stem] = dataset
            logger.info(f"  Loaded: {dataset.num_frames} frames")
        except Exception as e:
            logger.error(f"  Failed to load {d_folder.name}: {e}")

    return datasets


def extract_xic(
    dataset: TimsDatasetDDA,
    mz: float,
    rt_center: float,
    charge: int,
    mz_tol_ppm: float = 20.0,
    rt_window: float = 30.0,
) -> Dict[str, List[float]]:
    """
    Extract chromatographic ion current (XIC) for a precursor.

    Args:
        dataset: TimsDatasetDDA object
        mz: Target m/z value
        rt_center: Retention time center (seconds)
        charge: Charge state
        mz_tol_ppm: m/z tolerance in ppm
        rt_window: RT window in seconds (total width)

    Returns:
        Dict with 'rt' and 'intensity' arrays
    """
    mz_tol = mz * mz_tol_ppm / 1e6
    rt_start = rt_center - rt_window / 2
    rt_end = rt_center + rt_window / 2

    rts = []
    intensities = []

    # Get frames within RT window
    for frame in dataset.get_frames_by_rt_range(rt_start, rt_end):
        # Sum intensity within m/z tolerance
        mz_mask = (frame.mz >= mz - mz_tol) & (frame.mz <= mz + mz_tol)
        intensity = frame.intensity[mz_mask].sum() if mz_mask.any() else 0.0

        rts.append(frame.rt)
        intensities.append(float(intensity))

    return {
        'xic_rt': rts,
        'xic_intensity': intensities,
    }


def extract_mobilogram(
    dataset: TimsDatasetDDA,
    mz: float,
    rt_center: float,
    mobility_center: float,
    mz_tol_ppm: float = 20.0,
    rt_window: float = 10.0,
    mobility_window: float = 0.1,
) -> Dict[str, List[float]]:
    """
    Extract ion mobility profile (mobilogram) for a precursor.

    Args:
        dataset: TimsDatasetDDA object
        mz: Target m/z value
        rt_center: Retention time center (seconds)
        mobility_center: Ion mobility center (1/K0)
        mz_tol_ppm: m/z tolerance in ppm
        rt_window: RT window in seconds
        mobility_window: Mobility window (1/K0 units)

    Returns:
        Dict with 'mobility' and 'intensity' arrays
    """
    mz_tol = mz * mz_tol_ppm / 1e6
    rt_start = rt_center - rt_window / 2
    rt_end = rt_center + rt_window / 2

    # Aggregate across frames in RT window
    mobility_dict = {}

    for frame in dataset.get_frames_by_rt_range(rt_start, rt_end):
        mz_mask = (frame.mz >= mz - mz_tol) & (frame.mz <= mz + mz_tol)

        if mz_mask.any():
            for mob, inten in zip(frame.mobility[mz_mask], frame.intensity[mz_mask]):
                # Bin mobility values
                mob_bin = round(mob, 3)
                mobility_dict[mob_bin] = mobility_dict.get(mob_bin, 0) + inten

    # Sort by mobility
    sorted_mobs = sorted(mobility_dict.keys())

    return {
        'mobilogram_mobility': sorted_mobs,
        'mobilogram_intensity': [mobility_dict[m] for m in sorted_mobs],
    }


def extract_isotope_pattern(
    dataset: TimsDatasetDDA,
    mz: float,
    rt_center: float,
    mobility_center: float,
    charge: int,
    mz_tol_ppm: float = 20.0,
    rt_window: float = 5.0,
    n_isotopes: int = 6,
) -> Dict[str, List[float]]:
    """
    Extract isotope distribution for a precursor.

    Args:
        dataset: TimsDatasetDDA object
        mz: Monoisotopic m/z
        rt_center: Retention time center
        mobility_center: Ion mobility center
        charge: Charge state
        mz_tol_ppm: m/z tolerance
        rt_window: RT window
        n_isotopes: Number of isotope peaks to extract

    Returns:
        Dict with 'isotope_mz' and 'isotope_intensity' arrays
    """
    mz_tol = mz * mz_tol_ppm / 1e6
    rt_start = rt_center - rt_window / 2
    rt_end = rt_center + rt_window / 2

    # Isotope spacing based on charge
    isotope_spacing = 1.003355 / charge  # Neutron mass / charge

    isotope_mzs = [mz + i * isotope_spacing for i in range(n_isotopes)]
    isotope_intensities = [0.0] * n_isotopes

    # Aggregate across frames
    for frame in dataset.get_frames_by_rt_range(rt_start, rt_end):
        for i, iso_mz in enumerate(isotope_mzs):
            mz_mask = (frame.mz >= iso_mz - mz_tol) & (frame.mz <= iso_mz + mz_tol)
            if mz_mask.any():
                isotope_intensities[i] += frame.intensity[mz_mask].sum()

    return {
        'isotope_mz': isotope_mzs,
        'isotope_intensity': isotope_intensities,
    }


def extract_features_for_psm(
    psm: pd.Series,
    dataset: TimsDatasetDDA,
    mz_tol_ppm: float,
    rt_window: float,
    mobility_window: float,
    do_xic: bool,
    do_mobilogram: bool,
    do_isotopes: bool,
) -> Dict[str, Any]:
    """
    Extract all raw features for a single PSM.
    """
    features = {
        'psm_id': psm.name,
        'sequence': psm.get('sequence', ''),
        'modified_sequence': psm.get('modified_sequence', psm.get('sequence', '')),
        'charge': int(psm.get('charge', 0)),
        'mz': float(psm.get('mz', 0)),
        'retention_time': float(psm.get('retention_time', 0)),
        'ccs': float(psm.get('ccs', 0)) if 'ccs' in psm else None,
    }

    # Get mobility (handle different column names)
    mobility = psm.get('mobility', psm.get('inv_mobility', 0))
    features['mobility'] = float(mobility) if mobility else 0.0

    mz = features['mz']
    rt = features['retention_time']
    mob = features['mobility']
    charge = features['charge']

    # Extract chromatogram (XIC)
    if do_xic and mz > 0 and rt > 0:
        try:
            xic = extract_xic(dataset, mz, rt, charge, mz_tol_ppm, rt_window)
            features.update(xic)
        except Exception as e:
            logger.debug(f"XIC extraction failed for PSM {psm.name}: {e}")
            features['xic_rt'] = []
            features['xic_intensity'] = []
    else:
        features['xic_rt'] = []
        features['xic_intensity'] = []

    # Extract mobilogram
    if do_mobilogram and mz > 0 and rt > 0 and mob > 0:
        try:
            mobilogram = extract_mobilogram(
                dataset, mz, rt, mob, mz_tol_ppm, rt_window / 3, mobility_window
            )
            features.update(mobilogram)
        except Exception as e:
            logger.debug(f"Mobilogram extraction failed for PSM {psm.name}: {e}")
            features['mobilogram_mobility'] = []
            features['mobilogram_intensity'] = []
    else:
        features['mobilogram_mobility'] = []
        features['mobilogram_intensity'] = []

    # Extract isotope pattern
    if do_isotopes and mz > 0 and rt > 0 and charge > 0:
        try:
            isotopes = extract_isotope_pattern(
                dataset, mz, rt, mob, charge, mz_tol_ppm, rt_window / 3
            )
            features.update(isotopes)
        except Exception as e:
            logger.debug(f"Isotope extraction failed for PSM {psm.name}: {e}")
            features['isotope_mz'] = []
            features['isotope_intensity'] = []
    else:
        features['isotope_mz'] = []
        features['isotope_intensity'] = []

    return features


def match_psm_to_raw_file(psm: pd.Series, datasets: Dict[str, TimsDatasetDDA]) -> Optional[str]:
    """
    Match a PSM to its source raw file.

    Uses the 'raw_file' or 'spectrum_id' column to identify the source.
    """
    # Try direct raw_file column
    if 'raw_file' in psm.index:
        raw_file = str(psm['raw_file'])
        # Extract stem (remove path and extension)
        stem = Path(raw_file).stem.replace('.d', '')
        if stem in datasets:
            return stem

    # Try spectrum_id (format: rawfile.scan.scan.charge)
    if 'spectrum_id' in psm.index:
        spectrum_id = str(psm['spectrum_id'])
        parts = spectrum_id.split('.')
        if len(parts) >= 1:
            stem = parts[0].replace('.d', '')
            if stem in datasets:
                return stem

    # Fallback: return first dataset if only one
    if len(datasets) == 1:
        return list(datasets.keys())[0]

    return None


def main():
    global logger
    logger = setup_logging(log_file)

    logger.info("=" * 60)
    logger.info("Starting raw feature extraction with imspy")
    logger.info("=" * 60)
    logger.info(f"Raw directory: {raw_dir}")
    logger.info(f"PSM file: {psm_file}")
    logger.info(f"Output: {output_file}")
    logger.info(f"Parameters:")
    logger.info(f"  mz_tolerance_ppm: {mz_tolerance_ppm}")
    logger.info(f"  rt_tolerance_sec: {rt_tolerance_sec}")
    logger.info(f"  mobility_tolerance: {mobility_tolerance}")
    logger.info(f"  extract_chromatogram: {extract_chromatogram}")
    logger.info(f"  extract_mobilogram: {extract_mobilogram}")
    logger.info(f"  extract_isotopes: {extract_isotopes}")

    # Load PSMs
    psms = load_psms(psm_file)

    # Load raw datasets
    datasets = find_raw_files(raw_dir)

    if not datasets:
        logger.error("No .d folders could be loaded!")
        sys.exit(1)

    # Extract features for each PSM
    all_features = []
    n_success = 0
    n_failed = 0

    logger.info(f"Extracting features for {len(psms)} PSMs...")

    for idx, (psm_idx, psm) in enumerate(psms.iterrows()):
        if idx % 1000 == 0:
            logger.info(f"Processing PSM {idx}/{len(psms)} ({n_success} success, {n_failed} failed)")

        # Find matching raw file
        raw_file_key = match_psm_to_raw_file(psm, datasets)

        if raw_file_key is None:
            logger.debug(f"Could not match PSM {psm_idx} to a raw file")
            n_failed += 1
            continue

        dataset = datasets[raw_file_key]

        try:
            features = extract_features_for_psm(
                psm,
                dataset,
                mz_tolerance_ppm,
                rt_tolerance_sec,
                mobility_tolerance,
                extract_chromatogram,
                extract_mobilogram,
                extract_isotopes,
            )
            features['raw_file'] = raw_file_key
            all_features.append(features)
            n_success += 1

        except Exception as e:
            logger.debug(f"Failed to extract features for PSM {psm_idx}: {e}")
            n_failed += 1

    logger.info(f"Extraction complete: {n_success} success, {n_failed} failed")

    # Convert to DataFrame
    logger.info("Converting to DataFrame...")
    df = pd.DataFrame(all_features)

    # Save to parquet
    logger.info(f"Saving to {output_file}...")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_file, index=False)

    logger.info(f"Saved {len(df)} feature records")
    logger.info("Raw feature extraction completed successfully")


if __name__ == "__main__":
    main()
