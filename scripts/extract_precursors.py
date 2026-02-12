#!/usr/bin/env python3
"""
Extract precursors with fragment spectra AND MS1 signal from timsTOF DDA data.

This script:
1. Loads a timsTOF .d dataset
2. Extracts all PASEF precursor-fragment pairs (via Rust backend)
3. Extracts MS1 precursor signal (XIC, mobilogram, isotopes) in parallel (via Rust)
4. Calculates moment statistics for all dimensions
5. Outputs combined metadata parquet + serialized blobs

Uses imspy's optimized Rust backend for both extractions.
"""

import sys
import io
import logging
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Optional, Iterator, Dict, Any
import json

import numpy as np
from numpy.typing import NDArray
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# imspy imports
from imspy_core.timstof import TimsDatasetDDA
from imspy_core.timstof.frame import TimsFrame
from imspy_connector import py_dda, py_chemistry

# scipy for Gaussian fitting
from scipy.optimize import curve_fit

# Optional zstd compression (fallback to gzip if not available)
try:
    import zstandard as zstd
    HAS_ZSTD = True
except ImportError:
    import gzip
    HAS_ZSTD = False


def calculate_moments(coords: np.ndarray, intensities: np.ndarray) -> Dict[str, float]:
    """Calculate statistical moments from coordinate and intensity arrays."""
    if len(coords) == 0 or np.sum(intensities) == 0:
        return {
            "mean": 0.0, "variance": 0.0, "skewness": 0.0,
            "apex": 0.0, "fwhm": 0.0, "total_intensity": 0.0
        }

    total = np.sum(intensities)
    weights = intensities / total

    mean = np.sum(coords * weights)
    variance = np.sum(weights * (coords - mean) ** 2)
    std = np.sqrt(variance) if variance > 0 else 1e-10
    skewness = np.sum(weights * ((coords - mean) / std) ** 3)

    apex_idx = np.argmax(intensities)
    apex = coords[apex_idx]

    half_max = intensities[apex_idx] / 2
    above_half = coords[intensities >= half_max]
    if len(above_half) >= 2:
        fwhm = above_half[-1] - above_half[0]
    else:
        fwhm = 2.355 * std

    return {
        "mean": float(mean),
        "variance": float(variance),
        "skewness": float(skewness),
        "apex": float(apex),
        "fwhm": float(fwhm),
        "total_intensity": float(total),
    }


def gaussian(x: np.ndarray, amplitude: float, mu: float, sigma: float) -> np.ndarray:
    """Gaussian function for curve fitting."""
    return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def fit_gaussian(coords: np.ndarray, intensities: np.ndarray) -> Dict[str, float]:
    """Fit Gaussian to 1D signal, return apex, sigma, r2."""
    if len(coords) < 3 or np.sum(intensities) == 0:
        return {"apex": 0.0, "sigma": 0.0, "r2": 0.0}

    total = np.sum(intensities)
    weights = intensities / total
    mu_init = np.sum(coords * weights)
    var_init = np.sum(weights * (coords - mu_init) ** 2)
    sigma_init = np.sqrt(var_init) if var_init > 0 else 0.1
    amp_init = np.max(intensities)

    try:
        popt, _ = curve_fit(
            gaussian, coords, intensities,
            p0=[amp_init, mu_init, sigma_init],
            bounds=([0, coords.min(), 1e-6],
                    [np.inf, coords.max(), coords.max() - coords.min() + 1e-6]),
            maxfev=1000
        )
        amplitude, mu, sigma = popt

        y_pred = gaussian(coords, amplitude, mu, sigma)
        ss_res = np.sum((intensities - y_pred) ** 2)
        ss_tot = np.sum((intensities - np.mean(intensities)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        r2 = max(0.0, min(1.0, r2))

        return {"apex": float(mu), "sigma": float(sigma), "r2": float(r2)}
    except (RuntimeError, ValueError):
        return {"apex": float(mu_init), "sigma": float(sigma_init), "r2": 0.0}


def compute_isotope_cosine_similarity(
    mass: float, charge: int, observed: List[float], n_isotopes: int = 5
) -> float:
    """Compute cosine similarity between observed and theoretical isotope envelope."""
    if charge <= 0 or mass <= 0 or not observed:
        return 0.0

    theo = py_chemistry.generate_precursor_spectrum(
        mass=mass, charge=charge, min_intensity=1,
        k=n_isotopes, resolution=3, centroid=True
    )
    theo_int = np.array(theo.intensity[:n_isotopes])
    obs = np.array(observed[:n_isotopes])

    # Pad to same length
    max_len = max(len(theo_int), len(obs))
    theo_int = np.pad(theo_int, (0, max_len - len(theo_int)))
    obs = np.pad(obs, (0, max_len - len(obs)))

    # Normalize and compute cosine
    theo_norm = np.linalg.norm(theo_int)
    obs_norm = np.linalg.norm(obs)
    if theo_norm == 0 or obs_norm == 0:
        return 0.0

    return float(np.dot(theo_int, obs) / (theo_norm * obs_norm))


def compute_sage_fragment_cosine(
    obs_mz: np.ndarray, obs_int: np.ndarray,
    sage_frags: List[dict], ppm_tol: float = 20.0,
) -> float:
    """Cosine similarity between observed fragment spectrum and Sage matched ions.

    For each Sage fragment, finds the closest observed peak within ppm_tol
    using binary search and builds matched intensity vectors for cosine.

    Returns NaN if fewer than 2 fragments match.
    """
    if len(obs_mz) == 0 or not sage_frags:
        return float('nan')

    # Sort observed spectrum by m/z for binary search
    sort_idx = np.argsort(obs_mz)
    obs_mz_sorted = obs_mz[sort_idx]
    obs_int_sorted = obs_int[sort_idx]

    matched_obs = []
    matched_sage = []

    for frag in sage_frags:
        sage_mz = frag['mz_experimental']
        sage_int = frag['intensity']
        tol = sage_mz * ppm_tol / 1e6

        # Binary search for closest peak
        idx = np.searchsorted(obs_mz_sorted, sage_mz)
        best_idx = None
        best_delta = tol

        for candidate in [idx - 1, idx]:
            if 0 <= candidate < len(obs_mz_sorted):
                delta = abs(obs_mz_sorted[candidate] - sage_mz)
                if delta < best_delta:
                    best_delta = delta
                    best_idx = candidate

        if best_idx is not None:
            matched_obs.append(obs_int_sorted[best_idx])
            matched_sage.append(sage_int)

    if len(matched_obs) < 2:
        return float('nan')

    a = np.array(matched_obs)
    b = np.array(matched_sage)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)

    if norm_a == 0 or norm_b == 0:
        return float('nan')

    return float(np.dot(a, b) / (norm_a * norm_b))


@dataclass
class ExtractedPrecursor:
    """Complete extracted precursor with fragment spectrum and MS1 signal."""

    # Identifiers
    precursor_id: int
    raw_file: str

    # Precursor properties (from Bruker metadata)
    charge: Optional[int]
    mono_mz: Optional[float]
    largest_peak_mz: float
    average_mz: float
    isolation_mz: float
    isolation_width: float
    precursor_intensity: float

    # Temporal/mobility from PASEF
    rt_seconds: float
    mobility: float  # 1/K0 from fragment signal

    # Fragment spectrum
    fragment_frame: TimsFrame
    n_fragments_merged: int
    collision_energies: List[float]

    # Fragment signal moments
    fragment_mz_mean: float
    fragment_mz_var: float
    fragment_im_mean: float
    fragment_im_var: float
    fragment_im_apex: float
    fragment_im_fwhm: float

    # MS1 precursor signal (from Rust extraction)
    ms1_rt_apex: float
    ms1_rt_fwhm: float
    ms1_rt_skew: float
    ms1_im_apex: float
    ms1_im_fwhm: float
    ms1_im_skew: float
    ms1_total_intensity: float
    ms1_isotope_intensities: List[float]

    # Raw MS1 signal arrays (for visualization)
    ms1_rt_coords: np.ndarray
    ms1_rt_intensities: np.ndarray
    ms1_im_coords: np.ndarray
    ms1_im_intensities: np.ndarray

    # Raw 4D MS1 data (full point cloud for RT vs IM heatmap)
    raw_rt: np.ndarray
    raw_mz: np.ndarray
    raw_mobility: np.ndarray
    raw_intensity: np.ndarray

    # Gaussian fit quality metrics
    ms1_rt_sigma: float
    ms1_rt_r2: float
    ms1_im_sigma: float
    ms1_im_r2: float

    # Isotope envelope quality
    isotope_cosim: float

    # Sage fragment cosine similarity
    sage_cosine: float  # Cosine sim: observed fragment spectrum vs Sage matched ions

    @property
    def mz(self) -> float:
        if self.mono_mz is not None and not np.isnan(self.mono_mz):
            return self.mono_mz
        return self.largest_peak_mz

    @property
    def collision_energy(self) -> float:
        """Average collision energy across all PASEF events for this precursor."""
        if not self.collision_energies:
            return 0.0
        return float(np.mean(self.collision_energies))

    @property
    def n_peaks(self) -> int:
        return len(self.fragment_frame.mz)

    @property
    def fragment_total_intensity(self) -> float:
        return float(np.sum(self.fragment_frame.intensity))

    def to_index_dict(self) -> dict:
        """Return dict for index parquet."""
        return {
            "precursor_id": self.precursor_id,
            "raw_file": self.raw_file,
            "charge": self.charge if self.charge is not None else 0,
            "mono_mz": self.mono_mz,
            "mz": self.mz,
            "isolation_mz": self.isolation_mz,
            "isolation_width": self.isolation_width,
            "precursor_intensity": self.precursor_intensity,
            "rt_seconds": self.rt_seconds,
            "mobility": self.mobility,
            "collision_energy": self.collision_energy,  # Average CE across PASEF events
            "n_fragments_merged": self.n_fragments_merged,
            "n_peaks": self.n_peaks,
            "fragment_total_intensity": self.fragment_total_intensity,
            # Fragment moments
            "fragment_mz_mean": self.fragment_mz_mean,
            "fragment_mz_var": self.fragment_mz_var,
            "fragment_im_mean": self.fragment_im_mean,
            "fragment_im_var": self.fragment_im_var,
            "fragment_im_apex": self.fragment_im_apex,
            "fragment_im_fwhm": self.fragment_im_fwhm,
            # MS1 moments
            "ms1_rt_apex": self.ms1_rt_apex,
            "ms1_rt_fwhm": self.ms1_rt_fwhm,
            "ms1_rt_skew": self.ms1_rt_skew,
            "ms1_im_apex": self.ms1_im_apex,
            "ms1_im_fwhm": self.ms1_im_fwhm,
            "ms1_im_skew": self.ms1_im_skew,
            "ms1_total_intensity": self.ms1_total_intensity,
            # Isotope pattern (first 5)
            "ms1_iso_0": self.ms1_isotope_intensities[0] if len(self.ms1_isotope_intensities) > 0 else 0.0,
            "ms1_iso_1": self.ms1_isotope_intensities[1] if len(self.ms1_isotope_intensities) > 1 else 0.0,
            "ms1_iso_2": self.ms1_isotope_intensities[2] if len(self.ms1_isotope_intensities) > 2 else 0.0,
            "ms1_iso_3": self.ms1_isotope_intensities[3] if len(self.ms1_isotope_intensities) > 3 else 0.0,
            "ms1_iso_4": self.ms1_isotope_intensities[4] if len(self.ms1_isotope_intensities) > 4 else 0.0,
            # Gaussian fit quality metrics
            "ms1_rt_sigma": self.ms1_rt_sigma,
            "ms1_rt_r2": self.ms1_rt_r2,
            "ms1_im_sigma": self.ms1_im_sigma,
            "ms1_im_r2": self.ms1_im_r2,
            # Isotope envelope quality
            "isotope_cosim": self.isotope_cosim,
            # Sage fragment cosine similarity
            "sage_cosine": self.sage_cosine,
        }

    def serialize_blob(self) -> bytes:
        """Serialize fragment frame + MS1 signal to compressed bytes."""
        frame = self.fragment_frame

        # Fragment arrays - use float32 for signal data (sufficient precision, 50% size reduction)
        fragment_arrays = {
            "frag_frame_id": np.array([frame.frame_id], dtype=np.int32),
            "frag_ms_type": np.array([frame.ms_type], dtype=np.int32),
            "frag_rt": np.array([frame.retention_time], dtype=np.float32),
            "frag_scan": frame.scan.astype(np.int32),
            "frag_mobility": frame.mobility.astype(np.float32),
            "frag_tof": frame.tof.astype(np.int32),
            "frag_mz": frame.mz.astype(np.float32),
            "frag_intensity": frame.intensity.astype(np.float32),
            # MS1 signal arrays - float32 sufficient for visualization/training
            "ms1_rt_coords": self.ms1_rt_coords.astype(np.float32),
            "ms1_rt_intensities": self.ms1_rt_intensities.astype(np.float32),
            "ms1_im_coords": self.ms1_im_coords.astype(np.float32),
            "ms1_im_intensities": self.ms1_im_intensities.astype(np.float32),
            # Raw 4D MS1 point cloud (for RT vs IM heatmap)
            "raw_rt": self.raw_rt.astype(np.float32),
            "raw_mz": self.raw_mz.astype(np.float32),
            "raw_mobility": self.raw_mobility.astype(np.float32),
            "raw_intensity": self.raw_intensity.astype(np.float32),
        }

        # Metadata
        metadata = {
            "precursor_id": self.precursor_id,
            "collision_energies": self.collision_energies,
            "n_fragments_merged": self.n_fragments_merged,
            "ms1_isotope_intensities": self.ms1_isotope_intensities,
        }
        metadata_bytes = json.dumps(metadata).encode("utf-8")

        # Pack arrays
        npz_buffer = io.BytesIO()
        np.savez_compressed(npz_buffer, **fragment_arrays)
        npz_bytes = npz_buffer.getvalue()

        # Combine: [metadata_len (4 bytes)][metadata][npz]
        metadata_len = len(metadata_bytes).to_bytes(4, "little")
        combined = metadata_len + metadata_bytes + npz_bytes

        # Compress
        if HAS_ZSTD:
            cctx = zstd.ZstdCompressor(level=3)
            return cctx.compress(combined)
        else:
            return gzip.compress(combined)

    @classmethod
    def deserialize_blob(cls, data: bytes) -> dict:
        """Deserialize blob to dict with arrays."""
        if HAS_ZSTD:
            dctx = zstd.ZstdDecompressor()
            combined = dctx.decompress(data)
        else:
            combined = gzip.decompress(data)

        metadata_len = int.from_bytes(combined[:4], "little")
        metadata_bytes = combined[4:4+metadata_len]
        metadata = json.loads(metadata_bytes.decode("utf-8"))

        npz_bytes = combined[4+metadata_len:]
        npz_buffer = io.BytesIO(npz_bytes)
        arrays = dict(np.load(npz_buffer))

        return {"metadata": metadata, "arrays": arrays}


def extract_precursors(
    dataset: TimsDatasetDDA,
    raw_file_name: str,
    num_threads: int = 16,
    rt_window_sec: float = 15.0,
    mz_tol_ppm: float = 20.0,
    im_window: float = 0.2,  # Doubled from 0.1 for better coverage
    n_isotopes: int = 5,
    calibration: Optional[NDArray[np.float64]] = None,
    sage_fragments_by_id: Optional[Dict[int, List[dict]]] = None,
    logger: Optional[logging.Logger] = None,
) -> List[ExtractedPrecursor]:
    """
    Extract all precursors with fragment spectra and MS1 signal.

    Args:
        dataset: Loaded TimsDatasetDDA object
        raw_file_name: Name of the .d file
        num_threads: Threads for parallel processing
        rt_window_sec: RT window for MS1 extraction (seconds)
        mz_tol_ppm: m/z tolerance for MS1 extraction (ppm)
        im_window: IM window for MS1 extraction (1/K0), default 0.2 (doubled for coverage)
        n_isotopes: Number of isotope peaks to extract
        calibration: Optional IM calibration array (scan → 1/K0 lookup).
                     If provided, uses LookupIndexConverter for fast + accurate extraction.
                     If None, uses the dataset's default converter.
        sage_fragments_by_id: Optional dict mapping precursor_id to list of Sage fragment dicts
        logger: Optional logger

    Returns:
        List of ExtractedPrecursor objects
    """
    if logger:
        logger.info(f"Extracting from {raw_file_name}")

    # Step 1: Get PASEF fragments
    if logger:
        logger.info("  Getting PASEF fragments...")
    fragments_df = dataset.get_pasef_fragments(num_threads=num_threads)
    if logger:
        logger.info(f"    Raw fragments: {len(fragments_df)}")

    # Group by precursor_id
    grouped = fragments_df.groupby('precursor_id').agg({
        'frame_id': 'first',
        'time': 'first',
        'raw_data': 'sum',
        'scan_begin': 'first',
        'scan_end': 'first',
        'isolation_mz': 'first',
        'isolation_width': 'first',
        'collision_energy': list,
        'largest_peak_mz': 'first',
        'average_mz': 'first',
        'monoisotopic_mz': 'first',
        'charge': 'first',
        'average_scan': 'first',
        'intensity': 'first',
        'parent_id': 'first',
    })
    if logger:
        logger.info(f"    Unique precursors: {len(grouped)}")

    # Step 2: Build coordinates for MS1 extraction
    if logger:
        logger.info("  Preparing MS1 extraction coordinates...")

    coords = []
    precursor_data = []  # Store fragment data for later

    for precursor_id, row in grouped.iterrows():
        fragment_frame = row['raw_data']
        mobility = fragment_frame.get_inverse_mobility_along_scan_marginal()

        charge = int(row['charge']) if pd.notna(row['charge']) else None
        mono_mz = float(row['monoisotopic_mz']) if pd.notna(row['monoisotopic_mz']) else None
        largest_peak_mz = float(row['largest_peak_mz'])
        rt_sec = float(row['time']) * 60.0

        # Get IM bounds from fragment selection (scan_begin/scan_end)
        # These are the IM values used for PASEF selection, useful for plotting
        scan_begin = int(row['scan_begin'])
        scan_end = int(row['scan_end'])
        # Convert scans to 1/K0 using the dataset (scan_end has lower mobility number)
        im_values = dataset.scan_to_inverse_mobility(int(row['frame_id']), [scan_begin, scan_end])
        im_start = float(im_values[0])  # Higher 1/K0 (lower scan)
        im_end = float(im_values[1])    # Lower 1/K0 (higher scan)

        # mono_mz: pass as 0.0 if not available (Rust will use largest_peak_mz as fallback)
        mono_mz_value = mono_mz if mono_mz is not None and not np.isnan(mono_mz) else 0.0

        # Build coordinate for Rust extraction
        coord = py_dda.PyPrecursorCoord(
            precursor_id=int(precursor_id),
            mz=float(largest_peak_mz),       # Fallback m/z if mono_mz is 0
            mono_mz=float(mono_mz_value),    # 0 if unknown, triggers single-peak extraction
            rt_seconds=float(rt_sec),
            mobility=float(mobility),
            im_start=float(im_start),        # Fragment selection IM bounds (for plotting)
            im_end=float(im_end),
            charge=int(charge) if charge is not None else 2,
        )
        coords.append(coord)

        # Store fragment data
        precursor_data.append({
            'precursor_id': int(precursor_id),
            'charge': charge,
            'mono_mz': mono_mz,
            'largest_peak_mz': largest_peak_mz,
            'average_mz': float(row['average_mz']),
            'isolation_mz': float(row['isolation_mz']),
            'isolation_width': float(row['isolation_width']),
            'precursor_intensity': float(row['intensity']),
            'rt_seconds': rt_sec,
            'mobility': float(mobility),
            'im_start': im_start,
            'im_end': im_end,
            'fragment_frame': fragment_frame,
            'n_fragments_merged': len(row['collision_energy']),
            'collision_energies': [float(ce) for ce in row['collision_energy']],
        })

    # Step 3: Extract MS1 signals in parallel (Rust)
    if logger:
        logger.info(f"  Extracting MS1 signals ({len(coords)} precursors)...")

    # Use calibrated dataset if calibration provided (fast + accurate via LookupIndexConverter)
    if calibration is not None and len(calibration) > 0:
        if logger:
            logger.info(f"    Using IM calibration ({len(calibration)} scans) for accurate extraction")
        rust_dataset = py_dda.PyTimsDatasetDDA.with_calibration(
            dataset.data_path, False, calibration.tolist()
        )
    else:
        rust_dataset = dataset.get_py_ptr()

    ms1_signals = rust_dataset.extract_precursor_ms1_signals(
        coords,
        rt_window_sec=rt_window_sec,
        mz_tol_ppm=mz_tol_ppm,
        im_window=im_window,
        n_isotopes=n_isotopes,
        num_threads=num_threads,
    )

    if logger:
        logger.info(f"    MS1 extraction complete")

    # Step 4: Combine fragment and MS1 data
    if logger:
        logger.info("  Combining fragment and MS1 data...")

    # Build lookup for MS1 signals
    ms1_lookup = {s.precursor_id: s for s in ms1_signals}

    results = []
    for pdata in precursor_data:
        pid = pdata['precursor_id']
        ms1 = ms1_lookup.get(pid)

        # Calculate fragment moments
        frag_frame = pdata['fragment_frame']
        frag_mz_moments = calculate_moments(frag_frame.mz, frag_frame.intensity)
        frag_im_moments = calculate_moments(frag_frame.mobility, frag_frame.intensity)

        # Get MS1 data
        if ms1 is not None:
            ms1_rt_apex = ms1.rt_moments.apex
            ms1_rt_fwhm = ms1.rt_moments.fwhm
            ms1_rt_skew = ms1.rt_moments.skewness
            ms1_im_apex = ms1.im_moments.apex
            ms1_im_fwhm = ms1.im_moments.fwhm
            ms1_im_skew = ms1.im_moments.skewness
            ms1_total_intensity = ms1.rt_moments.total_intensity
            ms1_isotope_intensities = list(ms1.isotope_intensity)
            ms1_rt_coords = np.array(ms1.rt_coords)
            ms1_rt_intensities = np.array(ms1.rt_intensities)
            ms1_im_coords = np.array(ms1.im_coords)
            ms1_im_intensities = np.array(ms1.im_intensities)

            # Raw 4D MS1 point cloud
            raw_rt = np.array(ms1.raw_rt)
            raw_mz = np.array(ms1.raw_mz)
            raw_mobility = np.array(ms1.raw_mobility)
            raw_intensity = np.array(ms1.raw_intensity)

            # Gaussian fit quality metrics
            rt_fit = fit_gaussian(ms1_rt_coords, ms1_rt_intensities)
            im_fit = fit_gaussian(ms1_im_coords, ms1_im_intensities)
            ms1_rt_sigma = rt_fit["sigma"]
            ms1_rt_r2 = rt_fit["r2"]
            ms1_im_sigma = im_fit["sigma"]
            ms1_im_r2 = im_fit["r2"]

            # Isotope cosine similarity
            precursor_mz = pdata['mono_mz'] if pdata['mono_mz'] is not None else pdata['largest_peak_mz']
            charge = pdata['charge'] if pdata['charge'] is not None else 2
            mass = (precursor_mz - 1.007276) * charge
            isotope_cosim = compute_isotope_cosine_similarity(
                mass, charge, ms1_isotope_intensities, n_isotopes=5
            )
        else:
            ms1_rt_apex = 0.0
            ms1_rt_fwhm = 0.0
            ms1_rt_skew = 0.0
            ms1_im_apex = 0.0
            ms1_im_fwhm = 0.0
            ms1_im_skew = 0.0
            ms1_total_intensity = 0.0
            ms1_isotope_intensities = [0.0] * n_isotopes
            ms1_rt_coords = np.array([])
            ms1_rt_intensities = np.array([])
            ms1_im_coords = np.array([])
            ms1_im_intensities = np.array([])
            raw_rt = np.array([])
            raw_mz = np.array([])
            raw_mobility = np.array([])
            raw_intensity = np.array([])
            ms1_rt_sigma = 0.0
            ms1_rt_r2 = 0.0
            ms1_im_sigma = 0.0
            ms1_im_r2 = 0.0
            isotope_cosim = 0.0

        # Sage fragment cosine similarity
        sage_cosine = float('nan')
        if sage_fragments_by_id is not None:
            sage_frags = sage_fragments_by_id.get(pid)
            if sage_frags and len(frag_frame.mz) > 0:
                sage_cosine = compute_sage_fragment_cosine(
                    frag_frame.mz, frag_frame.intensity, sage_frags
                )

        precursor = ExtractedPrecursor(
            precursor_id=pid,
            raw_file=raw_file_name,
            charge=pdata['charge'],
            mono_mz=pdata['mono_mz'],
            largest_peak_mz=pdata['largest_peak_mz'],
            average_mz=pdata['average_mz'],
            isolation_mz=pdata['isolation_mz'],
            isolation_width=pdata['isolation_width'],
            precursor_intensity=pdata['precursor_intensity'],
            rt_seconds=pdata['rt_seconds'],
            mobility=pdata['mobility'],
            fragment_frame=pdata['fragment_frame'],
            n_fragments_merged=pdata['n_fragments_merged'],
            collision_energies=pdata['collision_energies'],
            fragment_mz_mean=frag_mz_moments['mean'],
            fragment_mz_var=frag_mz_moments['variance'],
            fragment_im_mean=frag_im_moments['mean'],
            fragment_im_var=frag_im_moments['variance'],
            fragment_im_apex=frag_im_moments['apex'],
            fragment_im_fwhm=frag_im_moments['fwhm'],
            ms1_rt_apex=ms1_rt_apex,
            ms1_rt_fwhm=ms1_rt_fwhm,
            ms1_rt_skew=ms1_rt_skew,
            ms1_im_apex=ms1_im_apex,
            ms1_im_fwhm=ms1_im_fwhm,
            ms1_im_skew=ms1_im_skew,
            ms1_total_intensity=ms1_total_intensity,
            ms1_isotope_intensities=ms1_isotope_intensities,
            ms1_rt_coords=ms1_rt_coords,
            ms1_rt_intensities=ms1_rt_intensities,
            ms1_im_coords=ms1_im_coords,
            ms1_im_intensities=ms1_im_intensities,
            raw_rt=raw_rt,
            raw_mz=raw_mz,
            raw_mobility=raw_mobility,
            raw_intensity=raw_intensity,
            ms1_rt_sigma=ms1_rt_sigma,
            ms1_rt_r2=ms1_rt_r2,
            ms1_im_sigma=ms1_im_sigma,
            ms1_im_r2=ms1_im_r2,
            isotope_cosim=isotope_cosim,
            sage_cosine=sage_cosine,
        )
        results.append(precursor)

    if logger:
        logger.info(f"  Extracted {len(results)} precursors")

        # Report collision energy statistics
        ce_values = [p.collision_energy for p in results if p.collision_energy > 0]
        if ce_values:
            ce_arr = np.array(ce_values)
            logger.info(f"  Collision energy: mean={np.mean(ce_arr):.1f} eV, "
                       f"range=[{np.min(ce_arr):.1f}, {np.max(ce_arr):.1f}] eV")

            # Check for CE variance in re-fragmented precursors
            n_refragmented = sum(1 for p in results if p.n_fragments_merged > 1)
            if n_refragmented > 0:
                ce_ranges = []
                for p in results:
                    if p.n_fragments_merged > 1 and len(p.collision_energies) > 1:
                        ce_range = max(p.collision_energies) - min(p.collision_energies)
                        ce_ranges.append(ce_range)
                if ce_ranges:
                    max_range = max(ce_ranges)
                    avg_range = np.mean(ce_ranges)
                    logger.info(f"  Re-fragmented precursors: {n_refragmented}, "
                               f"CE variance: avg={avg_range:.2f} eV, max={max_range:.2f} eV")
                    if max_range > 5.0:
                        logger.warning(f"    Note: Some precursors have CE spread > 5 eV - using average")

    return results


def write_extraction(
    precursors: List[ExtractedPrecursor],
    output_dir: Path,
    write_blobs: bool = True,
    logger: Optional[logging.Logger] = None,
) -> dict:
    """Write extracted precursors to index parquet + blob file."""
    output_dir.mkdir(parents=True, exist_ok=True)

    index_path = output_dir / "index.parquet"
    blobs_path = output_dir / "blobs.bin"

    index_records = []
    current_offset = 0

    blob_file = open(blobs_path, "wb") if write_blobs else None

    try:
        for precursor in precursors:
            record = precursor.to_index_dict()

            if write_blobs:
                blob_bytes = precursor.serialize_blob()
                blob_size = len(blob_bytes)
                record["blob_offset"] = current_offset
                record["blob_size"] = blob_size
                blob_file.write(blob_bytes)
                current_offset += blob_size

            index_records.append(record)
    finally:
        if blob_file:
            blob_file.close()

    # Write index
    if index_records:
        index_df = pd.DataFrame(index_records)
        index_df.to_parquet(index_path, index=False)

        if logger:
            logger.info(f"Wrote index: {index_path} ({len(index_df)} precursors)")
            if write_blobs:
                logger.info(f"Wrote blobs: {blobs_path} ({current_offset / 1024 / 1024:.1f} MB)")

    return {
        "n_precursors": len(index_records),
        "blob_size_bytes": current_offset if write_blobs else 0,
    }


def extract_precursors_batched(
    dataset: TimsDatasetDDA,
    raw_file_name: str,
    output_dir: Path,
    batch_size: int = 10000,
    num_threads: int = 16,
    rt_window_sec: float = 15.0,
    mz_tol_ppm: float = 20.0,
    im_window: float = 0.2,  # Doubled from 0.1 for better coverage
    n_isotopes: int = 5,
    calibration: Optional[NDArray[np.float64]] = None,
    sage_fragments_by_id: Optional[Dict[int, List[dict]]] = None,
    logger: Optional[logging.Logger] = None,
) -> dict:
    """
    Extract precursors in batches to reduce peak memory usage.

    Instead of loading all 200k+ precursors into memory, processes in chunks
    and writes directly to disk. Reduces memory from ~65GB to ~10GB per file.

    Args:
        dataset: Loaded TimsDatasetDDA object
        raw_file_name: Name of the .d file
        output_dir: Output directory for index.parquet and blobs.bin
        batch_size: Number of precursors per batch (default 10k)
        num_threads: Threads for parallel MS1 extraction
        rt_window_sec: RT window for MS1 extraction (seconds)
        mz_tol_ppm: m/z tolerance for MS1 extraction (ppm)
        im_window: IM window for MS1 extraction (1/K0)
        n_isotopes: Number of isotope peaks to extract
        calibration: Optional IM calibration array (scan → 1/K0 lookup)
        sage_fragments_by_id: Optional dict mapping precursor_id to list of Sage fragment dicts
        logger: Optional logger

    Returns:
        Dict with extraction statistics
    """
    import gc
    from math import ceil

    output_dir.mkdir(parents=True, exist_ok=True)
    index_path = output_dir / "index.parquet"
    blobs_path = output_dir / "blobs.bin"

    if logger:
        logger.info(f"Batched extraction from {raw_file_name}")
        logger.info(f"  Batch size: {batch_size}")

    # Step 1: Get PASEF fragments (this is relatively lightweight)
    if logger:
        logger.info("  Getting PASEF fragments...")
    fragments_df = dataset.get_pasef_fragments(num_threads=num_threads)
    if logger:
        logger.info(f"    Raw fragments: {len(fragments_df)}")

    # Group by precursor_id - keep only metadata, not the heavy raw_data yet
    grouped = fragments_df.groupby('precursor_id').agg({
        'frame_id': 'first',
        'time': 'first',
        'raw_data': 'sum',  # Merge fragment frames
        'scan_begin': 'first',
        'scan_end': 'first',
        'isolation_mz': 'first',
        'isolation_width': 'first',
        'collision_energy': list,
        'largest_peak_mz': 'first',
        'average_mz': 'first',
        'monoisotopic_mz': 'first',
        'charge': 'first',
        'average_scan': 'first',
        'intensity': 'first',
        'parent_id': 'first',
    })

    precursor_ids = list(grouped.index)
    n_precursors = len(precursor_ids)
    n_batches = ceil(n_precursors / batch_size)

    if logger:
        logger.info(f"    Unique precursors: {n_precursors}")
        logger.info(f"    Batches: {n_batches} x {batch_size}")

    # Create calibrated Rust dataset once (reused across batches)
    if calibration is not None and len(calibration) > 0:
        if logger:
            logger.info(f"  Using IM calibration ({len(calibration)} scans)")
        rust_dataset = py_dda.PyTimsDatasetDDA.with_calibration(
            dataset.data_path, False, calibration.tolist()
        )
    else:
        rust_dataset = dataset.get_py_ptr()

    # Process batches
    all_records = []
    current_offset = 0

    with open(blobs_path, "wb") as blob_file:
        for batch_idx in range(n_batches):
            batch_start = batch_idx * batch_size
            batch_end = min(batch_start + batch_size, n_precursors)
            batch_ids = precursor_ids[batch_start:batch_end]

            if logger:
                logger.info(f"  Batch {batch_idx + 1}/{n_batches}: precursors {batch_start}-{batch_end}")

            # Get grouped data for this batch
            batch_grouped = grouped.loc[batch_ids]

            # Build coords and precursor data for batch
            coords = []
            batch_data = []

            for precursor_id, row in batch_grouped.iterrows():
                fragment_frame = row['raw_data']
                mobility = fragment_frame.get_inverse_mobility_along_scan_marginal()

                charge = int(row['charge']) if pd.notna(row['charge']) else None
                mono_mz = float(row['monoisotopic_mz']) if pd.notna(row['monoisotopic_mz']) else None
                largest_peak_mz = float(row['largest_peak_mz'])
                rt_sec = float(row['time']) * 60.0

                # Get IM bounds from fragment selection (scan_begin/scan_end)
                scan_begin = int(row['scan_begin'])
                scan_end = int(row['scan_end'])
                im_values = dataset.scan_to_inverse_mobility(int(row['frame_id']), [scan_begin, scan_end])
                im_start = float(im_values[0])
                im_end = float(im_values[1])

                # mono_mz: pass as 0.0 if not available (Rust will use largest_peak_mz as fallback)
                mono_mz_value = mono_mz if mono_mz is not None and not np.isnan(mono_mz) else 0.0

                coord = py_dda.PyPrecursorCoord(
                    precursor_id=int(precursor_id),
                    mz=float(largest_peak_mz),       # Fallback m/z if mono_mz is 0
                    mono_mz=float(mono_mz_value),    # 0 if unknown, triggers single-peak extraction
                    rt_seconds=float(rt_sec),
                    mobility=float(mobility),
                    im_start=float(im_start),        # Fragment selection IM bounds
                    im_end=float(im_end),
                    charge=int(charge) if charge is not None else 2,
                )
                coords.append(coord)

                batch_data.append({
                    'precursor_id': int(precursor_id),
                    'charge': charge,
                    'mono_mz': mono_mz,
                    'largest_peak_mz': largest_peak_mz,
                    'average_mz': float(row['average_mz']),
                    'isolation_mz': float(row['isolation_mz']),
                    'isolation_width': float(row['isolation_width']),
                    'precursor_intensity': float(row['intensity']),
                    'rt_seconds': rt_sec,
                    'mobility': float(mobility),
                    'im_start': im_start,
                    'im_end': im_end,
                    'fragment_frame': fragment_frame,
                    'n_fragments_merged': len(row['collision_energy']),
                    'collision_energies': [float(ce) for ce in row['collision_energy']],
                })

            # Extract MS1 signals for batch
            ms1_signals = rust_dataset.extract_precursor_ms1_signals(
                coords,
                rt_window_sec=rt_window_sec,
                mz_tol_ppm=mz_tol_ppm,
                im_window=im_window,
                n_isotopes=n_isotopes,
                num_threads=num_threads,
            )

            # Build lookup for MS1 signals
            ms1_lookup = {s.precursor_id: s for s in ms1_signals}

            # Combine and write batch
            for pdata in batch_data:
                pid = pdata['precursor_id']
                ms1 = ms1_lookup.get(pid)

                # Calculate fragment moments
                frag_frame = pdata['fragment_frame']
                frag_mz_moments = calculate_moments(frag_frame.mz, frag_frame.intensity)
                frag_im_moments = calculate_moments(frag_frame.mobility, frag_frame.intensity)

                # Get MS1 data
                if ms1 is not None:
                    ms1_rt_apex = ms1.rt_moments.apex
                    ms1_rt_fwhm = ms1.rt_moments.fwhm
                    ms1_rt_skew = ms1.rt_moments.skewness
                    ms1_im_apex = ms1.im_moments.apex
                    ms1_im_fwhm = ms1.im_moments.fwhm
                    ms1_im_skew = ms1.im_moments.skewness
                    ms1_total_intensity = ms1.rt_moments.total_intensity
                    ms1_isotope_intensities = list(ms1.isotope_intensity)
                    ms1_rt_coords = np.array(ms1.rt_coords)
                    ms1_rt_intensities = np.array(ms1.rt_intensities)
                    ms1_im_coords = np.array(ms1.im_coords)
                    ms1_im_intensities = np.array(ms1.im_intensities)

                    # Raw 4D MS1 point cloud
                    raw_rt = np.array(ms1.raw_rt)
                    raw_mz = np.array(ms1.raw_mz)
                    raw_mobility = np.array(ms1.raw_mobility)
                    raw_intensity = np.array(ms1.raw_intensity)

                    # Gaussian fit quality metrics
                    rt_fit = fit_gaussian(ms1_rt_coords, ms1_rt_intensities)
                    im_fit = fit_gaussian(ms1_im_coords, ms1_im_intensities)
                    ms1_rt_sigma = rt_fit["sigma"]
                    ms1_rt_r2 = rt_fit["r2"]
                    ms1_im_sigma = im_fit["sigma"]
                    ms1_im_r2 = im_fit["r2"]

                    # Isotope cosine similarity
                    precursor_mz = pdata['mono_mz'] if pdata['mono_mz'] is not None else pdata['largest_peak_mz']
                    charge = pdata['charge'] if pdata['charge'] is not None else 2
                    mass = (precursor_mz - 1.007276) * charge
                    isotope_cosim = compute_isotope_cosine_similarity(
                        mass, charge, ms1_isotope_intensities, n_isotopes=5
                    )
                else:
                    ms1_rt_apex = 0.0
                    ms1_rt_fwhm = 0.0
                    ms1_rt_skew = 0.0
                    ms1_im_apex = 0.0
                    ms1_im_fwhm = 0.0
                    ms1_im_skew = 0.0
                    ms1_total_intensity = 0.0
                    ms1_isotope_intensities = [0.0] * n_isotopes
                    ms1_rt_coords = np.array([])
                    ms1_rt_intensities = np.array([])
                    ms1_im_coords = np.array([])
                    ms1_im_intensities = np.array([])
                    raw_rt = np.array([])
                    raw_mz = np.array([])
                    raw_mobility = np.array([])
                    raw_intensity = np.array([])
                    ms1_rt_sigma = 0.0
                    ms1_rt_r2 = 0.0
                    ms1_im_sigma = 0.0
                    ms1_im_r2 = 0.0
                    isotope_cosim = 0.0

                # Sage fragment cosine similarity
                sage_cosine = float('nan')
                if sage_fragments_by_id is not None:
                    sage_frags = sage_fragments_by_id.get(pid)
                    if sage_frags and len(frag_frame.mz) > 0:
                        sage_cosine = compute_sage_fragment_cosine(
                            frag_frame.mz, frag_frame.intensity, sage_frags
                        )

                precursor = ExtractedPrecursor(
                    precursor_id=pid,
                    raw_file=raw_file_name,
                    charge=pdata['charge'],
                    mono_mz=pdata['mono_mz'],
                    largest_peak_mz=pdata['largest_peak_mz'],
                    average_mz=pdata['average_mz'],
                    isolation_mz=pdata['isolation_mz'],
                    isolation_width=pdata['isolation_width'],
                    precursor_intensity=pdata['precursor_intensity'],
                    rt_seconds=pdata['rt_seconds'],
                    mobility=pdata['mobility'],
                    fragment_frame=pdata['fragment_frame'],
                    n_fragments_merged=pdata['n_fragments_merged'],
                    collision_energies=pdata['collision_energies'],
                    fragment_mz_mean=frag_mz_moments['mean'],
                    fragment_mz_var=frag_mz_moments['variance'],
                    fragment_im_mean=frag_im_moments['mean'],
                    fragment_im_var=frag_im_moments['variance'],
                    fragment_im_apex=frag_im_moments['apex'],
                    fragment_im_fwhm=frag_im_moments['fwhm'],
                    ms1_rt_apex=ms1_rt_apex,
                    ms1_rt_fwhm=ms1_rt_fwhm,
                    ms1_rt_skew=ms1_rt_skew,
                    ms1_im_apex=ms1_im_apex,
                    ms1_im_fwhm=ms1_im_fwhm,
                    ms1_im_skew=ms1_im_skew,
                    ms1_total_intensity=ms1_total_intensity,
                    ms1_isotope_intensities=ms1_isotope_intensities,
                    ms1_rt_coords=ms1_rt_coords,
                    ms1_rt_intensities=ms1_rt_intensities,
                    ms1_im_coords=ms1_im_coords,
                    ms1_im_intensities=ms1_im_intensities,
                    raw_rt=raw_rt,
                    raw_mz=raw_mz,
                    raw_mobility=raw_mobility,
                    raw_intensity=raw_intensity,
                    ms1_rt_sigma=ms1_rt_sigma,
                    ms1_rt_r2=ms1_rt_r2,
                    ms1_im_sigma=ms1_im_sigma,
                    ms1_im_r2=ms1_im_r2,
                    isotope_cosim=isotope_cosim,
                    sage_cosine=sage_cosine,
                )

                # Write blob and record offset
                record = precursor.to_index_dict()
                blob_bytes = precursor.serialize_blob()
                record["blob_offset"] = current_offset
                record["blob_size"] = len(blob_bytes)
                blob_file.write(blob_bytes)
                current_offset += len(blob_bytes)
                all_records.append(record)

            # Clear batch memory
            del coords, batch_data, ms1_signals, ms1_lookup
            gc.collect()

            if logger:
                logger.info(f"    Batch complete, blob size: {current_offset / 1024 / 1024:.1f} MB")

    # Write final index
    if all_records:
        index_df = pd.DataFrame(all_records)
        index_df.to_parquet(index_path, index=False)

        if logger:
            logger.info(f"Wrote index: {index_path} ({len(index_df)} precursors)")
            logger.info(f"Wrote blobs: {blobs_path} ({current_offset / 1024 / 1024:.1f} MB)")

            # Report collision energy statistics
            if "collision_energy" in index_df.columns:
                ce_values = index_df["collision_energy"].dropna()
                ce_values = ce_values[ce_values > 0]
                if len(ce_values) > 0:
                    logger.info(f"Collision energy: mean={ce_values.mean():.1f} eV, "
                               f"range=[{ce_values.min():.1f}, {ce_values.max():.1f}] eV")

    return {
        "n_precursors": len(all_records),
        "blob_size_bytes": current_offset,
        "n_batches": n_batches,
    }


def setup_logging(log_path: Optional[Path] = None) -> logging.Logger:
    """Set up logging."""
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_path)
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Extract precursors from timsTOF DDA data")
    parser.add_argument("input", type=Path, help="Path to .d folder")
    parser.add_argument("output", type=Path, help="Output directory")
    parser.add_argument("--threads", type=int, default=16, help="Number of threads")
    parser.add_argument("--no-blobs", action="store_true", help="Skip blob writing")
    parser.add_argument("--no-bruker-sdk", action="store_true", help="Don't use Bruker SDK")
    parser.add_argument("--rt-window", type=float, default=15.0, help="RT window for MS1 (seconds, ±7.5s around apex)")
    parser.add_argument("--mz-tol", type=float, default=20.0, help="m/z tolerance for MS1 (ppm)")
    parser.add_argument("--im-window", type=float, default=0.1, help="IM window for MS1 (1/K0)")
    parser.add_argument("--log", type=Path, help="Log file path")

    args = parser.parse_args()

    logger = setup_logging(args.log)

    logger.info("=" * 60)
    logger.info("Precursor + MS1 extraction with imspy")
    logger.info("=" * 60)
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    logger.info(f"Threads: {args.threads}")
    logger.info(f"RT window: {args.rt_window} sec")
    logger.info(f"m/z tol: {args.mz_tol} ppm")
    logger.info(f"IM window: {args.im_window}")

    # Load dataset
    logger.info("Loading dataset...")
    dataset = TimsDatasetDDA(
        str(args.input),
        in_memory=False,
        use_bruker_sdk=not args.no_bruker_sdk,
    )
    logger.info(f"  Frames: {dataset.frame_count}")
    logger.info(f"  Fragmented precursors: {len(dataset.fragmented_precursors)}")

    # Extract
    precursors = extract_precursors(
        dataset=dataset,
        raw_file_name=args.input.name,
        num_threads=args.threads,
        rt_window_sec=args.rt_window,
        mz_tol_ppm=args.mz_tol,
        im_window=args.im_window,
        logger=logger,
    )

    # Write
    stats = write_extraction(
        precursors=precursors,
        output_dir=args.output,
        write_blobs=not args.no_blobs,
        logger=logger,
    )

    logger.info("=" * 60)
    logger.info(f"Extraction complete: {stats['n_precursors']} precursors")
    if stats['blob_size_bytes'] > 0:
        logger.info(f"Blob size: {stats['blob_size_bytes'] / 1024 / 1024:.1f} MB")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
