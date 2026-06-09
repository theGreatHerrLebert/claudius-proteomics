#!/usr/bin/env python3
"""Step 1-2 of plan v2: fit the empirical FWHM(G) target curve on the 102-dataset
corpus, cluster-safe, with dataset-clustered LODO + gradient-proxy sensitivity.

Target is the MODEL-FREE ms1_rt_fwhm (half-max width), not Gaussian/EMG sigma.
SINGLE pass over the 37GB corpus (all gradient proxies computed together).
Output = the real FWHM(G) law timsim's prior must hit (step 4 tunes timsim's sigma
by simulation to match this; that needs timsim runs and is separate).
"""
import glob
import re
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = "/scratch/claudius-proteomics/_mogon_raw_features"
R2_MIN = 0.8
MIN_GOOD = 200
QUANTS = [95, 97.5, 99, 100]


def per_run_table():
    """ONE pass: per run, G at several rt-quantiles + median model-free FWHM."""
    rows = []
    files = sorted(glob.glob(f"{ROOT}/**/raw_features.parquet", recursive=True))
    for i, f in enumerate(files):
        m = re.search(r"(PXD\d+)", f); ds = m.group(1) if m else f
        try:
            t = pq.read_table(f, columns=["raw_file", "rt_seconds", "ms1_rt_fwhm", "ms1_rt_r2"])
        except Exception:
            continue
        def col(n): return t.column(n).combine_chunks().to_numpy(zero_copy_only=False)
        rf = col("raw_file"); rt = col("rt_seconds").astype(float)
        fw = col("ms1_rt_fwhm").astype(float); r2 = col("ms1_rt_r2").astype(float)
        for run in np.unique(rf):
            mm = rf == run
            good = mm & (r2 >= R2_MIN) & np.isfinite(fw) & (fw > 0)
            if good.sum() < MIN_GOOD:
                continue
            qs = np.nanpercentile(rt[mm], QUANTS)
            if not np.isfinite(qs[2]) or qs[2] <= 0:
                continue
            rows.append((ds, qs[0], qs[1], qs[2], qs[3], float(np.median(fw[good]))))
        if i % 20 == 0:
            print(f"  read {i}/{len(files)}", flush=True)
    return pd.DataFrame(rows, columns=["dataset", "G95", "G975", "G99", "Gmax", "fwhm"])


def fit_and_lodo(ds_tbl, label):
    G = ds_tbl.G.values; F = ds_tbl.fwhm.values
    b, a = np.polyfit(np.log(G), np.log(F), 1)
    r = np.corrcoef(np.log(G), np.log(F))[0, 1]
    e_pow, e_lin, e_const = [], [], []
    for i in range(len(ds_tbl)):
        o = np.arange(len(ds_tbl)) != i
        bp, ap = np.polyfit(np.log(G[o]), np.log(F[o]), 1)
        e_pow.append(abs(np.exp(ap) * G[i] ** bp - F[i]) / F[i])
        ml, cl = np.polyfit(G[o], F[o], 1)
        e_lin.append(abs(ml * G[i] + cl - F[i]) / F[i])
        e_const.append(abs(np.median(F[o]) - F[i]) / F[i])
    print(f"[{label}] n_ds={len(ds_tbl)}  FWHM = {np.exp(a):.3f}·G^{b:.3f}  (log-log r={r:.3f})")
    print(f"   LODO median rel-err  power={np.median(e_pow):.3f}  linear={np.median(e_lin):.3f}  const={np.median(e_const):.3f}")
    return b, np.exp(a), np.median(e_pow), np.median(e_const)


def main():
    df = per_run_table()
    print(f"\nruns={len(df)} datasets={df.dataset.nunique()}  "
          f"FWHM {df.fwhm.min():.1f}-{df.fwhm.max():.1f}s  G99 {df.G99.min():.0f}-{df.G99.max():.0f}s")

    print("\n=== primary fit (G = p99 rt, dataset-median, cluster-safe) ===")
    ds = df.groupby("dataset").agg(G=("G99", "median"), fwhm=("fwhm", "median")).reset_index()
    b, a, e_pow, e_const = fit_and_lodo(ds, "p99")
    print(f"   => power-law beats constant OOS: {e_pow < e_const} ({e_const/e_pow:.2f}x)  "
          f"exponent b={b:.2f} (0.5≈√G, 1.0=linear)")

    print("\n=== gradient-proxy sensitivity (exponent stability) ===")
    for col in ["G95", "G975", "G99", "Gmax"]:
        d2 = df.groupby("dataset").agg(G=(col, "median"), fwhm=("fwhm", "median")).reset_index()
        bb, aa = np.polyfit(np.log(d2.G), np.log(d2.fwhm), 1)
        print(f"   {col}: b={bb:.3f}  a={np.exp(aa):.3f}")

    print("\n=== target FWHM(G) at representative gradients (p99 fit) ===")
    for mins in (30, 60, 90, 120, 197):
        Gs = mins * 60
        print(f"   {mins:>3} min (G={Gs}s): FWHM ≈ {a*Gs**b:.2f} s")
    ds.to_csv(f"{ROOT}/_fwhm_gradient.csv", index=False)
    print("\nwrote _fwhm_gradient.csv  (FWHM->timsim-σ needs the EMG λ; step 4 sim-matches it)")


if __name__ == "__main__":
    main()
