#!/usr/bin/env python3
"""Probe how reliable the EMG sigma estimate is under intensity / SNR shifts.

Two views:
  A) SYNTHETIC ground-truth: inject EMG XICs at a realistic MS1 frame sampling,
     sweep peak amplitude (= intensity) with a Poisson+baseline noise model, refit,
     and measure sigma recovery error + sigma_ok rate vs SNR. Tells us where the
     width estimate becomes unreliable, with a known true sigma.
  B) EMPIRICAL: bin the real frac01 precursors by measured precursor intensity and
     report sigma_ok rate / EMG-R2 / n_pts / sigma — does the synthetic story hold?

Sample-type axis is deferred: only frac01 (one HeLa sample) has local blobs;
needs a selective multi-dataset blob pull once raw_features identifies sample types.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import exponnorm

REPO = Path("/home/administrator/Documents/promotion/claudius-proteomics")
DATA = Path("/scratch/claudius-proteomics")
BASE = DATA / "_extracted_poc" / "PXD019086"
# emg_refit_poc.py lives on the feat/sim-shape-predictor branch (not the currently
# checked-out fix branch's working tree); materialised to /tmp/ssprobe by the runner.
sys.path.insert(0, "/tmp/ssprobe")
sys.path.insert(0, str(REPO))  # for sanjose import inside emg_refit_poc
from emg_refit_poc import fit_emg  # reuse the hardened fitter  # noqa: E402

RNG = np.random.default_rng(0)


def synthetic_intensity_sweep():
    """Known sigma; sweep amplitude (intensity) under Poisson+baseline noise."""
    print("=== A) SYNTHETIC: sigma recovery vs intensity (true_sigma=3.5s) ===")
    true_sig, true_lam = 3.5, 0.15
    K = 1.0 / (true_lam * true_sig)
    mu = 30.0
    dt = 1.1                      # ~ real frac01 MS1 cycle (rt_span/n_pts ~ 1s)
    x = np.arange(mu - 6, mu + 12, dt)   # ~16 pts spanning the peak (real: 7-13)
    shape = exponnorm.pdf(x, K, loc=mu, scale=true_sig)
    shape = shape / shape.max()
    baseline = 50.0              # constant chemical-noise floor (counts)
    N = 400
    print(f"{'apex_intensity':>14}{'~SNR':>7}{'sig_relerr_med':>15}{'within20%':>10}{'sigma_ok%':>10}")
    for apex in [200, 500, 1000, 3000, 10000, 50000]:
        rel, oks, win = [], [], []
        for _ in range(N):
            clean = apex * shape
            # Poisson shot noise on signal + Gaussian baseline noise
            noisy = RNG.poisson(np.clip(clean, 0, None)) + RNG.normal(0, baseline, size=x.shape)
            f = fit_emg(x, noisy)
            if f is None:
                oks.append(False); continue
            oks.append(f["sigma_ok"])
            if f["sigma_ok"]:
                e = abs(f["sigma"] - true_sig) / true_sig
                rel.append(e); win.append(e <= 0.2)
        snr = apex / baseline
        rm = np.median(rel) if rel else np.nan
        wn = 100 * np.mean(win) if win else np.nan
        print(f"{apex:>14}{snr:>7.0f}{rm:>15.3f}{wn:>9.0f}%{100*np.mean(oks):>9.0f}%")
    print("  (sigma_relerr over sigma_ok fits only; within20% = |sig-true|/true<=0.2)")


def empirical_by_intensity():
    refit = BASE / "emg_refit_frac01.parquet"
    if not refit.exists():
        print("\n(skip B: emg_refit_frac01.parquet not found)"); return
    print("\n=== B) EMPIRICAL frac01: reliability vs measured precursor intensity ===")
    em = pq.read_table(refit).to_pandas()
    rf = pq.read_table(BASE / "raw_features.parquet",
                       columns=["precursor_id", "precursor_intensity", "ms1_total_intensity"]).to_pandas()
    df = em.merge(rf, on="precursor_id", how="left")
    df = df[df.emg_r2.notna()].copy()
    df["int_bin"] = pd.qcut(df.precursor_intensity, 6, labels=False, duplicates="drop")
    print(f"{'int_decile':>10}{'med_intensity':>14}{'n':>8}{'sigma_ok%':>10}{'emg_r2_med':>11}{'n_pts_med':>10}{'sigma_med':>10}")
    for b, g in df.groupby("int_bin"):
        print(f"{int(b):>10}{g.precursor_intensity.median():>14.0f}{len(g):>8}"
              f"{100*g.sigma_ok.mean():>9.0f}%{g.emg_r2.median():>11.3f}"
              f"{g.n_pts.median():>10.0f}{g[g.sigma_ok].emg_sigma.median():>10.2f}")
    # correlation summary
    r = np.corrcoef(np.log10(df.precursor_intensity.clip(lower=1)), df.emg_r2)[0, 1]
    print(f"  corr(log10 intensity, EMG-R2) = {r:.3f}   "
          f"sigma_ok rate low vs high decile: "
          f"{100*df[df.int_bin==df.int_bin.min()].sigma_ok.mean():.0f}% -> "
          f"{100*df[df.int_bin==df.int_bin.max()].sigma_ok.mean():.0f}%")


if __name__ == "__main__":
    synthetic_intensity_sweep()
    empirical_by_intensity()
