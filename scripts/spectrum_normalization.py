#!/usr/bin/env python3
"""
Spectrum Normalization Utilities

Provides multiple normalization methods for MS2 spectra, preserving raw intensities
while computing normalized values for training.

Normalization methods:
- base_peak: I / max(I) * 100 (standard, preserves relative intensities)
- tic: I / sum(I) * 100 (total ion current, accounts for variable spectrum complexity)
- sqrt: sqrt(I / max(I)) * 100 (compressed dynamic range, reduces dominance of high peaks)
- log: log10(I + 1) / log10(max(I) + 1) * 100 (log scale, useful for wide dynamic range)

Usage:
    from spectrum_normalization import normalize_spectrum, NormalizationMethod

    mz = np.array([100.0, 200.0, 300.0])
    intensity = np.array([1000.0, 5000.0, 500.0])

    result = normalize_spectrum(mz, intensity, NormalizationMethod.BASE_PEAK)
    print(result.normalized_intensity)  # [20.0, 100.0, 10.0]
    print(result.normalization_factor)  # 5000.0
"""

from dataclasses import dataclass
from enum import Enum
from typing import Optional
import numpy as np


class NormalizationMethod(Enum):
    """Supported spectrum normalization methods."""
    BASE_PEAK = "base_peak"  # I / max(I) * 100
    TIC = "tic"              # I / sum(I) * 100
    SQRT = "sqrt"            # sqrt(I / max(I)) * 100
    LOG = "log"              # log10(I + 1) / log10(max(I) + 1) * 100

    @classmethod
    def from_string(cls, s: str) -> 'NormalizationMethod':
        """Create from string value."""
        s = s.lower().strip()
        for method in cls:
            if method.value == s:
                return method
        raise ValueError(f"Unknown normalization method: {s}. "
                        f"Valid options: {[m.value for m in cls]}")


@dataclass
class NormalizedSpectrum:
    """Result of spectrum normalization."""
    mz: np.ndarray                  # Original m/z values
    raw_intensity: np.ndarray       # Original intensities (preserved)
    normalized_intensity: np.ndarray  # Normalized intensities (0-100 scale)
    normalization_factor: float     # Factor used for normalization (for reversibility)
    method: NormalizationMethod     # Method used

    @property
    def n_peaks(self) -> int:
        """Number of peaks in spectrum."""
        return len(self.mz)

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            'mz': self.mz.tolist(),
            'raw_intensity': self.raw_intensity.tolist(),
            'normalized_intensity': self.normalized_intensity.tolist(),
            'normalization_factor': self.normalization_factor,
            'method': self.method.value,
        }


def normalize_spectrum(
    mz: np.ndarray,
    intensity: np.ndarray,
    method: NormalizationMethod = NormalizationMethod.BASE_PEAK,
    min_intensity: float = 0.0,
) -> NormalizedSpectrum:
    """
    Normalize a mass spectrum.

    Args:
        mz: Array of m/z values
        intensity: Array of intensity values
        method: Normalization method to use
        min_intensity: Minimum intensity threshold (peaks below are removed)

    Returns:
        NormalizedSpectrum with raw and normalized intensities
    """
    # Ensure numpy arrays
    mz = np.asarray(mz, dtype=np.float64)
    intensity = np.asarray(intensity, dtype=np.float64)

    # Handle empty spectra
    if len(mz) == 0 or len(intensity) == 0:
        return NormalizedSpectrum(
            mz=np.array([], dtype=np.float64),
            raw_intensity=np.array([], dtype=np.float64),
            normalized_intensity=np.array([], dtype=np.float64),
            normalization_factor=0.0,
            method=method,
        )

    # Filter by minimum intensity
    if min_intensity > 0:
        mask = intensity >= min_intensity
        mz = mz[mask]
        intensity = intensity[mask]

        if len(mz) == 0:
            return NormalizedSpectrum(
                mz=np.array([], dtype=np.float64),
                raw_intensity=np.array([], dtype=np.float64),
                normalized_intensity=np.array([], dtype=np.float64),
                normalization_factor=0.0,
                method=method,
            )

    # Compute normalized intensities based on method
    max_intensity = np.max(intensity)
    sum_intensity = np.sum(intensity)

    if method == NormalizationMethod.BASE_PEAK:
        # Standard base peak normalization
        if max_intensity > 0:
            normalized = (intensity / max_intensity) * 100.0
            factor = max_intensity
        else:
            normalized = np.zeros_like(intensity)
            factor = 1.0

    elif method == NormalizationMethod.TIC:
        # Total ion current normalization
        if sum_intensity > 0:
            normalized = (intensity / sum_intensity) * 100.0
            factor = sum_intensity
        else:
            normalized = np.zeros_like(intensity)
            factor = 1.0

    elif method == NormalizationMethod.SQRT:
        # Square root normalization (compressed dynamic range)
        if max_intensity > 0:
            normalized = np.sqrt(intensity / max_intensity) * 100.0
            factor = max_intensity
        else:
            normalized = np.zeros_like(intensity)
            factor = 1.0

    elif method == NormalizationMethod.LOG:
        # Log normalization
        if max_intensity > 0:
            log_max = np.log10(max_intensity + 1)
            if log_max > 0:
                normalized = (np.log10(intensity + 1) / log_max) * 100.0
            else:
                normalized = np.zeros_like(intensity)
            factor = max_intensity
        else:
            normalized = np.zeros_like(intensity)
            factor = 1.0

    else:
        raise ValueError(f"Unknown normalization method: {method}")

    return NormalizedSpectrum(
        mz=mz.copy(),
        raw_intensity=intensity.copy(),
        normalized_intensity=normalized,
        normalization_factor=factor,
        method=method,
    )


def denormalize_spectrum(
    normalized_intensity: np.ndarray,
    normalization_factor: float,
    method: NormalizationMethod = NormalizationMethod.BASE_PEAK,
) -> np.ndarray:
    """
    Reverse normalization to recover approximate raw intensities.

    Note: For SQRT and LOG methods, this is not perfectly reversible
    due to the non-linear transformations.

    Args:
        normalized_intensity: Normalized intensities (0-100 scale)
        normalization_factor: Factor from original normalization
        method: Method used for normalization

    Returns:
        Approximate raw intensities
    """
    normalized = np.asarray(normalized_intensity, dtype=np.float64)

    if method == NormalizationMethod.BASE_PEAK:
        return (normalized / 100.0) * normalization_factor

    elif method == NormalizationMethod.TIC:
        return (normalized / 100.0) * normalization_factor

    elif method == NormalizationMethod.SQRT:
        # Reverse: I_norm = sqrt(I/max)*100 => I = (I_norm/100)^2 * max
        return np.square(normalized / 100.0) * normalization_factor

    elif method == NormalizationMethod.LOG:
        # Reverse: I_norm = log10(I+1)/log10(max+1)*100 => I = 10^(I_norm/100*log10(max+1)) - 1
        log_max = np.log10(normalization_factor + 1)
        return np.power(10, (normalized / 100.0) * log_max) - 1

    else:
        raise ValueError(f"Unknown normalization method: {method}")


def normalize_batch(
    spectra: list,
    method: NormalizationMethod = NormalizationMethod.BASE_PEAK,
    min_intensity: float = 0.0,
) -> list:
    """
    Normalize a batch of spectra.

    Args:
        spectra: List of (mz_array, intensity_array) tuples
        method: Normalization method
        min_intensity: Minimum intensity threshold

    Returns:
        List of NormalizedSpectrum objects
    """
    results = []
    for mz, intensity in spectra:
        result = normalize_spectrum(mz, intensity, method, min_intensity)
        results.append(result)
    return results


# ============================================================================
# Demo / Testing
# ============================================================================

def demo():
    """Demonstrate normalization methods."""
    print("Spectrum Normalization Demo")
    print("=" * 60)

    # Example spectrum
    mz = np.array([100.0, 200.0, 300.0, 400.0, 500.0])
    intensity = np.array([1000.0, 5000.0, 500.0, 2000.0, 100.0])

    print(f"\nOriginal spectrum:")
    print(f"  m/z:       {mz}")
    print(f"  intensity: {intensity}")
    print(f"  max: {np.max(intensity)}, sum: {np.sum(intensity)}")

    for method in NormalizationMethod:
        result = normalize_spectrum(mz, intensity, method)
        print(f"\n{method.value.upper()} normalization:")
        print(f"  normalized: {np.round(result.normalized_intensity, 2)}")
        print(f"  factor:     {result.normalization_factor}")

        # Test reversibility
        recovered = denormalize_spectrum(
            result.normalized_intensity,
            result.normalization_factor,
            method,
        )
        error = np.max(np.abs(recovered - intensity))
        print(f"  max recovery error: {error:.6f}")

    # Test empty spectrum
    print("\n" + "=" * 60)
    print("Empty spectrum handling:")
    empty = normalize_spectrum(np.array([]), np.array([]), NormalizationMethod.BASE_PEAK)
    print(f"  n_peaks: {empty.n_peaks}")
    print(f"  factor: {empty.normalization_factor}")


if __name__ == "__main__":
    demo()
