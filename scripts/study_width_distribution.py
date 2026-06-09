#!/usr/bin/env python3
"""Pilot distributional study (within-run, local frac01) — learn what RT peak
width sigma actually depends on, before fixing a model target.

Uses the existing frac01 EMG fits (no blobs needed) joined with raw_features
(rt_seconds, intensity, charge) + precursor_index (peptide). Decomposes sigma
variance across factors via eta^2 (= SS_between / SS_total per factor):
  - RT-POSITION in the gradient  (Codex: local slope/position may beat total duration)
  - intensity
  - charge
  - peptide identity
Between-dataset (gradient, sample type) needs the multi-dataset blob subset (MOGON).
"""
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

DATA = "/scratch/claudius-proteomics"
BASE = f"{DATA}/_extracted_poc/PXD019086"
RAW = "20190504_TIMS1_FlMe_SA_HeLa_frac01_A10_1_93"


def eta_sq(y, groups):
    """Fraction of variance of y explained by a categorical grouping."""
    y = np.asarray(y, float)
    g = pd.Series(groups).values
    grand = y.mean()
    ss_tot = ((y - grand) ** 2).sum()
    ss_between = 0.0
    for _, idx in pd.Series(range(len(y))).groupby(g):
        yi = y[idx.values]
        ss_between += len(yi) * (yi.mean() - grand) ** 2
    return ss_between / ss_tot if ss_tot > 0 else np.nan


def main():
    em = pq.read_table(f"{BASE}/emg_refit_frac01.parquet").to_pandas()
    rf = pq.read_table(f"{BASE}/raw_features.parquet",
                       columns=["precursor_id", "rt_seconds", "precursor_intensity", "charge"]).to_pandas()
    pi = pq.read_table(f"{DATA}/PXD019086/precursor_index.parquet").to_pandas()
    pi = pi[pi.raw_file == RAW][["precursor_id", "sage_modified", "sage_qvalue"]]

    df = em.merge(rf, on="precursor_id", how="left").merge(pi, on="precursor_id", how="left")
    sig = df[df.sigma_ok].copy()
    print(f"sigma_ok precursors: {len(sig)}  (RT-position span {sig.rt_seconds.min():.0f}-{sig.rt_seconds.max():.0f}s)")

    # factor bins
    sig["rt_pos"] = pd.qcut(sig.rt_seconds, 10, labels=False, duplicates="drop")
    sig["int_bin"] = pd.qcut(sig.precursor_intensity, 10, labels=False, duplicates="drop")
    conf = sig[sig.sage_qvalue <= 0.01].dropna(subset=["sage_modified"]).copy()
    conf["pep"] = conf.sage_modified.astype(str) + "/" + conf.charge.astype(int).astype(str)

    print("\n=== eta^2: fraction of sigma variance explained by each factor ===")
    print(f"  RT-position (decile): {eta_sq(sig.emg_sigma, sig.rt_pos):.3f}")
    print(f"  intensity   (decile): {eta_sq(sig.emg_sigma, sig.int_bin):.3f}")
    print(f"  charge             : {eta_sq(sig.emg_sigma, sig.charge):.3f}")
    # peptide eta^2 on confident-ID subset with repeats
    rep = conf.groupby("pep").filter(lambda g: len(g) >= 2)
    print(f"  peptide (>=2 obs, n={rep.pep.nunique()} peps / {len(rep)} obs): "
          f"{eta_sq(rep.emg_sigma, rep.pep):.3f}")

    print("\n=== sigma vs RT-position in the gradient (within one run) ===")
    g = sig.groupby("rt_pos")
    tab = pd.DataFrame({"med_rt_s": g.rt_seconds.median(), "n": g.size(),
                        "sigma_med": g.emg_sigma.median(), "sigma_iqr": g.emg_sigma.quantile(.75) - g.emg_sigma.quantile(.25)})
    print(tab.to_string(formatters={"med_rt_s": "{:.0f}".format, "sigma_med": "{:.2f}".format, "sigma_iqr": "{:.2f}".format}))
    early, late = tab.sigma_med.iloc[0], tab.sigma_med.iloc[-1]
    print(f"  early->late gradient sigma: {early:.2f} -> {late:.2f} s  ({late/early:.2f}x)")

    print("\n=== lambda (tail) vs RT-position [lambda_ok subset] ===")
    lam = df[df.lambda_ok].copy()
    lam["rt_pos"] = pd.qcut(lam.rt_seconds, 10, labels=False, duplicates="drop")
    lg = lam.groupby("rt_pos").emg_lambda.median()
    print(f"  lambda eta^2 (RT-position): {eta_sq(lam.emg_lambda, lam.rt_pos):.3f}   "
          f"median lambda early->late: {lg.iloc[0]:.3f} -> {lg.iloc[-1]:.3f}")


if __name__ == "__main__":
    main()
