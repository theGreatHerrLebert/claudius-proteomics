"""Sum-integration fragment intensity readout for raw PASEF point-cloud blobs.

The blobs.bin "merged fragment spectrum" is a raw, un-centroided point cloud:
each ion is spread over several m/z points and the intact precursor survives
into the MS2 region. `extract_fragment_peaks._find_peak` returns the single
m/z-closest point, which severely under-reads every ion (~60x vs Sage).

`integrate_peaks_for_psm` instead sums every point intensity within tolerance
of each theoretical ion m/z -- the correct integration of a point cloud --
which recovers the true monoisotopic ion intensity.
"""
from __future__ import annotations

import numpy as np


def integrate_peak(mz: np.ndarray, inten: np.ndarray, target_mz: float,
                   tol_ppm: float) -> tuple[float, float, int]:
    """Sum all point intensities within +/- tol_ppm of target_mz.

    mz must be sorted ascending. Returns (summed_intensity,
    intensity-weighted ppm error, n_points)."""
    if mz.size == 0:
        return 0.0, np.nan, 0
    tol_da = target_mz * tol_ppm * 1e-6
    lo = np.searchsorted(mz, target_mz - tol_da)
    hi = np.searchsorted(mz, target_mz + tol_da, side="right")
    if lo >= hi:
        return 0.0, np.nan, 0
    w_mz, w_int = mz[lo:hi], inten[lo:hi]
    tot = float(w_int.sum())
    if tot > 0:
        err = float(np.sum((w_mz - target_mz) / target_mz * 1e6 * w_int) / tot)
    else:
        err = np.nan
    return tot, err, int(hi - lo)


def integrate_peaks_for_psm(frag_mz, frag_intensity, theoretical: dict,
                            tol_ppm: float = 20.0, precursor_mz: float | None = None,
                            precursor_excl_da: float = 2.0) -> dict:
    """Integrate every theoretical b/y ion against a raw point-cloud spectrum.

    theoretical: dict from FragmentMatcher.generate_theoretical_fragments,
                 {ion_type: [(ion_number, charge, mz, seq), ...]}.
    precursor_mz: if given, points within precursor_excl_da are dropped from
                  the total-intensity denominator (the intact precursor would
                  otherwise dominate intensity_explained)."""
    order = np.argsort(frag_mz)
    mz = np.asarray(frag_mz, dtype=float)[order]
    inten = np.asarray(frag_intensity, dtype=float)[order]

    if precursor_mz:
        keep = np.abs(mz - precursor_mz) > precursor_excl_da
        total = float(inten[keep].sum())
    else:
        total = float(inten.sum())

    ion_type, ion_number, ion_charge = [], [], []
    intensity, error_ppm, n_points = [], [], []
    matched_sum = 0.0
    nb = ny = nbm = nym = 0
    for it, ions in theoretical.items():
        for num, chg, theo_mz, _seq in ions:
            tot, err, n = integrate_peak(mz, inten, theo_mz, tol_ppm)
            ion_type.append(it); ion_number.append(num); ion_charge.append(chg)
            intensity.append(tot); error_ppm.append(err); n_points.append(n)
            if it == "b":
                nb += 1
            else:
                ny += 1
            if tot > 0:
                matched_sum += tot
                if it == "b":
                    nbm += 1
                else:
                    nym += 1
    return {
        "ion_type": ion_type, "ion_number": ion_number, "ion_charge": ion_charge,
        "intensity": intensity, "error_ppm": error_ppm, "n_points": n_points,
        "n_theoretical": nb + ny,
        "coverage_b": nbm / nb if nb else 0.0,
        "coverage_y": nym / ny if ny else 0.0,
        "intensity_explained": matched_sum / total if total > 0 else 0.0,
    }
