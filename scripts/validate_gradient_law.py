#!/usr/bin/env python3
"""Multi-gradient validation: does RT peak-width sigma scale with gradient length?

Uses the synced corpus-wide raw_features.parquet (102 datasets). For each LC run
(raw_file): gradient proxy G = p99(rt_seconds); width = median ms1_rt_sigma over
well-fit peaks (ms1_rt_r2>=0.8). Aggregated per dataset.

Tests (closes the n=4 caveat from the Codex review):
  - corr(G, sigma) and the log-log slope (sigma ∝ G => slope ~ 1).
  - dispersion of absolute sigma vs gradient-normalised sigma/G (normalise wins?).
  - LEAVE-ONE-DATASET-OUT: fit k = median(sigma/G) on the other datasets, predict
    the held-out dataset's sigma = k*G_held; compare to a no-gradient baseline
    (predict the global median sigma). If G-based beats baseline out-of-sample,
    the normalise-and-rescale assumption the RT head rests on is justified.

Caveat: G = p99(rt_seconds) is a proxy (tracks ID density, not the programmed
gradient); confounded with instrument/column/sample. Real LC metadata would be
better — flagged, not yet available.
"""
import glob
from pathlib import Path
import numpy as np
import pyarrow.parquet as pq

ROOT = Path("/scratch/claudius-proteomics/_mogon_raw_features")
R2_MIN = 0.8
MIN_GOOD = 200   # min well-fit peaks per run to trust its median sigma


def dataset_of(path: str) -> str:
    import re
    m = re.search(r"(PXD\d+)", path)
    return m.group(1) if m else Path(path).parent.name


def main():
    files = sorted(glob.glob(str(ROOT / "**" / "raw_features.parquet"), recursive=True))
    print(f"raw_features files: {len(files)}")
    rows = []  # (dataset, raw_file, G, sigma_med, n_good)
    for i, f in enumerate(files):
        ds = dataset_of(f)
        try:
            tbl = pq.read_table(f, columns=["raw_file", "rt_seconds", "ms1_rt_sigma", "ms1_rt_r2"])
        except Exception as e:
            print(f"  skip {ds}: {e}", flush=True); continue
        # to_numpy (near zero-copy for numerics) — NOT to_pydict, which is O(rows) slow
        def col(name):
            return tbl.column(name).combine_chunks().to_numpy(zero_copy_only=False)
        rf = col("raw_file")
        rt = col("rt_seconds").astype(float)
        sig = col("ms1_rt_sigma").astype(float)
        r2 = col("ms1_rt_r2").astype(float)
        if i % 10 == 0:
            print(f"  read {i}/{len(files)} ({ds})", flush=True)
        for run in np.unique(rf):
            m = rf == run
            rts = rt[m]; sigs = sig[m]; r2s = r2[m]
            G = np.nanpercentile(rts, 99)
            good = (r2s >= R2_MIN) & np.isfinite(sigs) & (sigs > 0)
            if good.sum() < MIN_GOOD or not np.isfinite(G) or G <= 0:
                continue
            rows.append((ds, str(run), float(G), float(np.median(sigs[good])), int(good.sum())))

    if not rows:
        print("no usable runs"); return
    import pandas as pd
    df = pd.DataFrame(rows, columns=["dataset", "raw_file", "G", "sigma", "n_good"])
    # aggregate to per-dataset (median across runs)
    ds = df.groupby("dataset").agg(G=("G", "median"), sigma=("sigma", "median"),
                                   n_runs=("raw_file", "size")).reset_index()
    print(f"runs used: {len(df)}   datasets: {len(ds)}")
    print(f"gradient range (p99 rt_s): {ds.G.min():.0f}–{ds.G.max():.0f} s "
          f"({ds.G.min()/60:.0f}–{ds.G.max()/60:.0f} min)")

    G = ds.G.values; S = ds.sigma.values
    print("\n=== scaling ===")
    print(f"corr(G, sigma)        = {np.corrcoef(G, S)[0,1]:.3f}")
    print(f"corr(log G, log sigma)= {np.corrcoef(np.log(G), np.log(S))[0,1]:.3f}")
    slope, intercept = np.polyfit(np.log(G), np.log(S), 1)
    print(f"log-log slope         = {slope:.3f}  (1.0 = perfectly proportional)")
    cv_abs = np.std(S) / np.mean(S)
    cv_norm = np.std(S / G) / np.mean(S / G)
    print(f"CV(sigma)={cv_abs:.3f}   CV(sigma/G)={cv_norm:.3f}   "
          f"=> normalisation {'REDUCES' if cv_norm < cv_abs else 'does NOT reduce'} dispersion "
          f"({cv_abs/cv_norm:.2f}x)")

    # leave-one-dataset-out
    print("\n=== leave-one-dataset-out (out-of-sample prediction of per-dataset sigma) ===")
    err_grad, err_base = [], []
    for i in range(len(ds)):
        others = np.arange(len(ds)) != i
        k = np.median(S[others] / G[others])          # proportional law from others
        pred_grad = k * G[i]
        pred_base = np.median(S[others])               # no-gradient baseline
        err_grad.append(abs(pred_grad - S[i]) / S[i])
        err_base.append(abs(pred_base - S[i]) / S[i])
    eg, eb = np.median(err_grad), np.median(err_base)
    print(f"  LODO median rel-error  gradient-law: {eg:.3f}   no-gradient baseline: {eb:.3f}")
    print(f"  => gradient law {'BEATS' if eg < eb else 'does NOT beat'} baseline out-of-sample "
          f"({eb/eg:.2f}x better)" if eg>0 else "")

    # show the spread
    ds["sigma_norm"] = ds.sigma / ds.G
    ds = ds.sort_values("G")
    print("\nper-dataset (sorted by gradient):")
    print(ds.to_string(index=False, formatters={"G": "{:.0f}".format, "sigma": "{:.2f}".format,
                                                 "sigma_norm": "{:.5f}".format}))
    out = "/scratch/claudius-proteomics/_mogon_raw_features/_gradient_law.csv"
    ds.to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
