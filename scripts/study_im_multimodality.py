#!/usr/bin/env python3
"""Decisive multi-modality test on raw mobilograms (frac01 blobs), NOISE-ROBUST
and intensity-stratified — IM traces are noisy, so naive peak-count / Gaussian-r2
inflate 'multimodality' at low intensity.

Robust criterion per mobilogram:
  - estimate noise floor sd from the high-freq residual (MAD of y - smooth(y));
  - a peak is REAL only if its prominence >= max(0.20*global_max, 5*noise_sd);
  - 'bimodal' = >=2 real peaks AND the valley between the top two drops at least
    5*noise_sd below the smaller of the two (a significant separation, not a shoulder).
Stratify the bimodal fraction by charge x intensity tercile. The readout that
matters: the HIGH-intensity (high-SNR) tercile — does bimodality persist there
(real conformers) and is it 3+-enriched, or does it collapse (it was noise)?
"""
import sys
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.signal import find_peaks

DATA = "/scratch/claudius-proteomics"
BASE = f"{DATA}/_extracted_poc/PXD019086"
RAW = "20190504_TIMS1_FlMe_SA_HeLa_frac01_A10_1_93"
sys.path.insert(0, "/home/administrator/Documents/promotion/claudius-proteomics")
from sanjose.blob import BlobReader  # noqa: E402

N_PER_Z = 5000
MIN_PTS = 11


def is_bimodal(inten):
    y = np.asarray(inten, float)
    if len(y) < MIN_PTS or y.max() <= 0:
        return None
    ys = np.convolve(y, np.ones(3) / 3, mode="same")
    noise = 1.4826 * np.median(np.abs(y - ys)) + 1e-9      # robust hi-freq noise sd
    thr = max(0.20 * ys.max(), 5 * noise)
    peaks, props = find_peaks(ys, prominence=thr, distance=2)
    if len(peaks) < 2:
        return False
    h = ys[peaks]; order = np.argsort(h)[::-1]
    p1, p2 = peaks[order[0]], peaks[order[1]]
    lo, hi = min(p1, p2), max(p1, p2)
    valley = ys[lo:hi + 1].min()
    return bool(valley <= h[order[1]] - 5 * noise)          # significant dip between top-2


def main():
    rf = pq.read_table(f"{BASE}/raw_features.parquet",
                       columns=["precursor_id", "blob_offset", "blob_size",
                                "charge", "precursor_intensity"]).to_pandas()
    reader = BlobReader(BASE)
    print(f"{'z':>3} {'int_tercile':>11}{'n':>8}{'bimodal%':>10}")
    table = {}
    for z in [1, 2, 3, 4]:
        sub = rf[(rf.charge == z) & (rf.precursor_intensity > 0)].copy()
        if len(sub) < 300:
            continue
        sub["it"] = pd.qcut(sub.precursor_intensity, 3, labels=["low", "mid", "high"])
        for tier in ["low", "mid", "high"]:
            s = sub[sub.it == tier]
            if len(s) > N_PER_Z:
                s = s.iloc[np.linspace(0, len(s) - 1, N_PER_Z).astype(int)]
            sigs = reader.read_batch(RAW, s.blob_offset.tolist(), s.blob_size.tolist())
            bm = [is_bimodal(x.ms1_signal.mobilogram_intensity) for x in sigs if x is not None]
            bm = [b for b in bm if b is not None]
            if bm:
                frac = np.mean(bm)
                table[(z, tier)] = frac
                print(f"{z:>3} {tier:>11}{len(bm):>8}{frac*100:>9.0f}%")
    print("\n=== high-SNR readout (the real test) ===")
    for z in [1, 2, 3, 4]:
        if (z, "high") in table:
            print(f"  charge {z}: bimodal(high-intensity) = {table[(z,'high')]*100:.0f}%")
    if (3, "high") in table and (2, "high") in table:
        print(f"  => at high SNR, 3+ vs 2+: {table[(3,'high')]*100:.0f}% vs {table[(2,'high')]*100:.0f}% "
              f"(3+ enriched: {table[(3,'high')] > table[(2,'high')]})")
        print(f"  => noise check: does bimodal% fall low->high intensity? "
              f"z2 {table.get((2,'low'),0)*100:.0f}->{table.get((2,'high'),0)*100:.0f}%, "
              f"z3 {table.get((3,'low'),0)*100:.0f}->{table.get((3,'high'),0)*100:.0f}%")


if __name__ == "__main__":
    main()
