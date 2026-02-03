#!/usr/bin/env python3
"""
Quality Metrics for Raw Data Extraction

Provides statistical measures and quality scores for MS1 signal:
- Moment statistics (mean, variance, skewness, FWHM)
- Gaussian fitting with R² quality metric
- Isotope envelope cosine similarity

These metrics are used in Step 4 (Raw Data Extraction) and
displayed in the Dashboard (Step 6).
"""

from typing import Dict, List, Optional
import numpy as np
from numpy.typing import NDArray


def calculate_moments(coords: np.ndarray, intensities: np.ndarray) -> Dict[str, float]:
    """Calculate statistical moments from coordinate and intensity arrays.

    Args:
        coords: 1D array of coordinates (RT, m/z, or mobility)
        intensities: 1D array of intensities

    Returns:
        Dict with: mean, variance, skewness, apex, fwhm, total_intensity
    """
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
        fwhm = 2.355 * std  # Gaussian FWHM approximation

    return {
        "mean": float(mean),
        "variance": float(variance),
        "skewness": float(skewness),
        "apex": float(apex),
        "fwhm": float(fwhm),
        "total_intensity": float(total),
    }


def gaussian(x: np.ndarray, amplitude: float, mu: float, sigma: float) -> np.ndarray:
    """Gaussian function for curve fitting.

    Args:
        x: Input coordinates
        amplitude: Peak height
        mu: Center position
        sigma: Standard deviation

    Returns:
        Gaussian values at input coordinates
    """
    return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def fit_gaussian(coords: np.ndarray, intensities: np.ndarray) -> Dict[str, float]:
    """Fit Gaussian to 1D signal, return apex, sigma, r2.

    Uses scipy.optimize.curve_fit with bounded parameters.
    Falls back to moment-based estimates if fitting fails.

    Args:
        coords: 1D array of coordinates
        intensities: 1D array of intensities

    Returns:
        Dict with: apex (mu), sigma, r2 (goodness of fit)
    """
    if len(coords) < 3 or np.sum(intensities) == 0:
        return {"apex": 0.0, "sigma": 0.0, "r2": 0.0}

    # Initial estimates from moments
    total = np.sum(intensities)
    weights = intensities / total
    mu_init = np.sum(coords * weights)
    var_init = np.sum(weights * (coords - mu_init) ** 2)
    sigma_init = np.sqrt(var_init) if var_init > 0 else 0.1
    amp_init = np.max(intensities)

    try:
        from scipy.optimize import curve_fit

        popt, _ = curve_fit(
            gaussian, coords, intensities,
            p0=[amp_init, mu_init, sigma_init],
            bounds=([0, coords.min(), 1e-6],
                    [np.inf, coords.max(), coords.max() - coords.min() + 1e-6]),
            maxfev=1000
        )
        amplitude, mu, sigma = popt

        # Compute R² (coefficient of determination)
        y_pred = gaussian(coords, amplitude, mu, sigma)
        ss_res = np.sum((intensities - y_pred) ** 2)
        ss_tot = np.sum((intensities - np.mean(intensities)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        r2 = max(0.0, min(1.0, r2))  # Clamp to [0, 1]

        return {"apex": float(mu), "sigma": float(sigma), "r2": float(r2)}

    except (RuntimeError, ValueError, ImportError):
        # Fitting failed - return moment-based estimates with r2=0
        return {"apex": float(mu_init), "sigma": float(sigma_init), "r2": 0.0}


def compute_isotope_cosine_similarity(
    mass: float,
    charge: int,
    observed: List[float],
    n_isotopes: int = 5
) -> float:
    """Compute cosine similarity between observed and theoretical isotope envelope.

    Uses imspy's theoretical isotope pattern generator based on averagine model.

    Args:
        mass: Precursor neutral mass (Da)
        charge: Precursor charge state
        observed: Observed isotope intensities (normalized or raw)
        n_isotopes: Number of isotope peaks to consider

    Returns:
        Cosine similarity [0, 1], where 1 = perfect match
    """
    if charge <= 0 or mass <= 0 or not observed:
        return 0.0

    try:
        from imspy_connector import py_chemistry

        # Generate theoretical spectrum
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

    except ImportError:
        # Fall back to simple pattern comparison if imspy not available
        return _simple_isotope_score(observed, n_isotopes)


def _simple_isotope_score(observed: List[float], n_isotopes: int = 5) -> float:
    """Simple isotope pattern quality score without theoretical comparison.

    Checks for reasonable decreasing pattern typical of peptide isotopes.

    Args:
        observed: Observed isotope intensities
        n_isotopes: Number of peaks to check

    Returns:
        Score [0, 1] based on pattern quality
    """
    if len(observed) < 2:
        return 0.0

    obs = np.array(observed[:n_isotopes])
    if np.sum(obs) == 0:
        return 0.0

    # Normalize
    obs = obs / np.max(obs)

    # Score based on expected pattern:
    # 1. First peak should be among highest
    # 2. Generally decreasing after first 2-3 peaks
    # 3. No large gaps

    score = 0.0

    # First peak reasonable (>0.5 of max)
    if obs[0] > 0.5:
        score += 0.3

    # Second peak present and reasonable
    if len(obs) > 1 and obs[1] > 0.2:
        score += 0.3

    # Pattern roughly decreasing after peak
    if len(obs) > 2:
        peak_idx = np.argmax(obs)
        after_peak = obs[peak_idx + 1:] if peak_idx + 1 < len(obs) else []
        if len(after_peak) > 0:
            # Check if generally decreasing
            decreasing = sum(1 for i in range(len(after_peak) - 1)
                           if after_peak[i] >= after_peak[i + 1])
            if decreasing >= len(after_peak) - 1:
                score += 0.4

    return score


def compute_quality_summary(
    rt_r2: float,
    im_r2: float,
    isotope_cosim: float,
    thresholds: Optional[Dict[str, float]] = None,
) -> Dict[str, any]:
    """Compute quality summary from individual metrics.

    Args:
        rt_r2: RT Gaussian fit R²
        im_r2: IM Gaussian fit R²
        isotope_cosim: Isotope cosine similarity
        thresholds: Custom thresholds (default: r2>=0.8, cosim>=0.9)

    Returns:
        Dict with is_high_quality, individual flags, combined score
    """
    if thresholds is None:
        thresholds = {
            "rt_r2_min": 0.8,
            "im_r2_min": 0.8,
            "isotope_cosim_min": 0.9,
        }

    rt_good = rt_r2 >= thresholds["rt_r2_min"]
    im_good = im_r2 >= thresholds["im_r2_min"]
    iso_good = isotope_cosim >= thresholds["isotope_cosim_min"]

    # Combined score (weighted average)
    combined = (rt_r2 * 0.3 + im_r2 * 0.3 + isotope_cosim * 0.4)

    return {
        "is_high_quality": rt_good and im_good and iso_good,
        "rt_quality": "good" if rt_good else "poor",
        "im_quality": "good" if im_good else "poor",
        "isotope_quality": "good" if iso_good else "poor",
        "combined_score": round(combined, 3),
    }


if __name__ == "__main__":
    # Test examples
    print("Testing quality metrics:\n")

    # Test moment calculation
    coords = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    intensities = np.array([10, 50, 100, 50, 10])
    moments = calculate_moments(coords, intensities)
    print(f"Moments test:")
    print(f"  Mean: {moments['mean']:.2f} (expected ~3.0)")
    print(f"  Apex: {moments['apex']:.2f} (expected 3.0)")
    print(f"  FWHM: {moments['fwhm']:.2f}")

    # Test Gaussian fit
    x = np.linspace(0, 10, 50)
    y = gaussian(x, 100, 5.0, 1.0) + np.random.normal(0, 2, len(x))
    fit_result = fit_gaussian(x, y)
    print(f"\nGaussian fit test:")
    print(f"  Fitted apex: {fit_result['apex']:.2f} (expected ~5.0)")
    print(f"  Fitted sigma: {fit_result['sigma']:.2f} (expected ~1.0)")
    print(f"  R²: {fit_result['r2']:.3f}")

    # Test isotope score (simple fallback)
    observed_isotopes = [0.8, 1.0, 0.6, 0.3, 0.1]
    score = _simple_isotope_score(observed_isotopes)
    print(f"\nIsotope pattern test:")
    print(f"  Observed: {observed_isotopes}")
    print(f"  Simple score: {score:.2f}")

    # Test quality summary
    summary = compute_quality_summary(rt_r2=0.95, im_r2=0.88, isotope_cosim=0.92)
    print(f"\nQuality summary:")
    print(f"  High quality: {summary['is_high_quality']}")
    print(f"  Combined score: {summary['combined_score']}")
