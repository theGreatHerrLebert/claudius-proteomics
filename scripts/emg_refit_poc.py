#!/usr/bin/env python3
"""EMG refit of RT XIC profiles for the SIM shape predictor (POC).

Refits the raw RT XIC (stored in blobs.bin) with an exponentially-modified
Gaussian (EMG) — the shape timsim consumes — instead of the symmetric Gaussian
the extractor stored. Parameterized in scipy's `exponnorm`, which maps directly
to timsim's (sigma, lambda, k):

    exponnorm.pdf(x, K, loc=mu, scale=sigma)
    sigma = scale ;  k = K (== timsim rt_k) ;  lambda = 1 / (K * sigma)

Outputs per precursor: EMG (mu, sigma, lambda, k, r2) alongside the stored
Gaussian (sigma, r2) so we can confirm the EMG actually improves the fit on
skewed peaks (and is a valid label source, vs the invalid moment-conversion).

Usage:
    python scripts/emg_refit_poc.py --n 3000        # sample
    python scripts/emg_refit_poc.py --all           # full frac01 (~214k)
"""
import argparse, sys, time
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from scipy.optimize import curve_fit
from scipy.stats import exponnorm

REPO = Path("/home/administrator/Documents/promotion/claudius-proteomics")
DATA = Path("/scratch/claudius-proteomics")
EXTRACTED = DATA / "_extracted_poc" / "PXD019086"
# index co-located with this blobs.bin (matches its offsets); NOT the top-level
# merged precursor_store.parquet, whose offsets index a different/newer blob file.
STORE = EXTRACTED / "raw_features.parquet"
RAW_FILE = "20190504_TIMS1_FlMe_SA_HeLa_frac01_A10_1_93"
OUT = DATA / "_extracted_poc" / "PXD019086" / "emg_refit_frac01.parquet"

sys.path.insert(0, str(REPO))
from sanjose.blob import BlobReader  # noqa: E402


MIN_PTS = 8        # below this an EMG (4 params) overfits sparse XICs
K_LO, K_HI = 1e-2, 50.0   # exponnorm shape; rails flag tail non-identifiability
R2_OK = 0.5        # minimum EMG R² for a usable shape label


def emg_model(x, amp, K, mu, sigma):
    return amp * exponnorm.pdf(x, K, loc=mu, scale=sigma)


def fit_emg(coords, inten):
    """Fit EMG to a 1D XIC with physical bounds + validity flags.

    Returns dict with mu/sigma/lambda/k/r2/n and booleans sigma_ok, lambda_ok.
    sigma_ok: enough points, converged, decent R², sigma not at a rail.
    lambda_ok: sigma_ok AND the exponential tail is identifiable (K interior) —
    sparse XICs frequently fail this (lambda → degenerate), so it gates whether
    lambda is a usable label vs left to timsim's sampler.
    """
    if coords is None or len(coords) < MIN_PTS:
        return None
    x = np.asarray(coords, float)
    y = np.asarray(inten, float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    if y.max() <= 0:
        return None
    span = float(x.max() - x.min())
    w = y / y.sum()
    mean = float((w * x).sum())
    var = float((w * (x - mean) ** 2).sum())
    sd = np.sqrt(var) if var > 0 else (span / 6 + 1e-6)
    sigma_max = max(0.5 * span, 1e-2)            # a peak can't fill the window
    sigma_init = min(max(sd, 1e-2), sigma_max)
    amp0 = y.max() * sigma_init * 2.5
    p0 = [amp0, 1.0, x[np.argmax(y)], sigma_init]
    bounds = ([0, K_LO, x.min() - 3 * sd, 1e-2],
              [np.inf, K_HI, x.max() + 3 * sd, sigma_max])
    try:
        popt, _ = curve_fit(emg_model, x, y, p0=p0, bounds=bounds, maxfev=4000)
    except Exception:
        return None
    amp, K, mu, sigma = popt
    yhat = emg_model(x, *popt)
    ss_res = float(((y - yhat) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    r2 = float(max(0.0, min(1.0, r2)))
    lam = 1.0 / (K * sigma)
    sigma_at_rail = sigma <= 1.05e-2 or sigma >= 0.97 * sigma_max
    K_at_rail = K <= 1.2 * K_LO or K >= 0.8 * K_HI
    sigma_ok = (r2 >= R2_OK) and (not sigma_at_rail)
    lambda_ok = sigma_ok and (not K_at_rail)
    return dict(mu=float(mu), sigma=float(sigma), lam=float(lam), k=float(K),
                r2=r2, n=len(x), sigma_ok=bool(sigma_ok), lambda_ok=bool(lambda_ok))


def _fit_chunk(payload):
    """Worker: open own BlobReader, read+fit a chunk of precursors."""
    extracted, raw_file, recs = payload
    reader = BlobReader(extracted)
    offs = [r["blob_offset"] for r in recs]
    szs = [r["blob_size"] for r in recs]
    sigs = reader.read_batch(raw_file, offs, szs)
    out = []
    for r, sig in zip(recs, sigs):
        rec = dict(precursor_id=r["precursor_id"],
                   gauss_sigma=r["ms1_rt_sigma"], gauss_r2=r["ms1_rt_r2"],
                   gauss_skew=r["ms1_rt_skew"], im_sigma=r["ms1_im_sigma"],
                   im_r2=r["ms1_im_r2"], emg_sigma=np.nan, emg_lambda=np.nan,
                   emg_k=np.nan, emg_mu=np.nan, emg_r2=np.nan, n_pts=0,
                   sigma_ok=False, lambda_ok=False)
        if sig is not None:
            f = fit_emg(sig.ms1_signal.xic_rt, sig.ms1_signal.xic_intensity)
            if f:
                rec.update(emg_sigma=f["sigma"], emg_lambda=f["lam"], emg_k=f["k"],
                           emg_mu=f["mu"], emg_r2=f["r2"], n_pts=f["n"],
                           sigma_ok=f["sigma_ok"], lambda_ok=f["lambda_ok"])
        out.append(rec)
    return out


def synthetic_recovery(rng_seed_terms=range(6)):
    """Inject XICs with known (sigma, lambda) + noise, refit, report recovery.
    Deterministic noise (no RNG: Date/random unavailable); vary by index."""
    print("\n=== synthetic EMG recovery (known sigma/lambda + structured noise) ===")
    print(f"{'true_sig':>9}{'true_lam':>9}{'rec_sig':>9}{'rec_lam':>9}{'r2':>7}")
    dt = 0.7  # ~frame spacing seconds
    for i in rng_seed_terms:
        true_sig = 2.0 + i * 1.5
        true_lam = 0.05 + 0.04 * i
        K = 1.0 / (true_lam * true_sig)
        x = np.arange(0, 60 + 8 * true_sig, dt)
        mu = 20.0
        y = 1000 * exponnorm.pdf(x, K, loc=mu, scale=true_sig)
        # structured multiplicative + additive perturbation (deterministic)
        y = y * (1 + 0.05 * np.sin(0.5 * x + i)) + 0.5 * (i + 1)
        f = fit_emg(x, y)
        if f:
            print(f"{true_sig:9.2f}{true_lam:9.3f}{f['sigma']:9.2f}{f['lam']:9.3f}{f['r2']:7.3f}")
        else:
            print(f"{true_sig:9.2f}{true_lam:9.3f}    FIT-FAILED")


def main():
    import pandas as pd
    from concurrent.futures import ProcessPoolExecutor

    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=3000)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--no-synth", action="store_true")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    if not args.no_synth:
        synthetic_recovery()

    print(f"\n=== loading raw_features for {RAW_FILE} ===")
    df = pq.read_table(STORE, columns=["precursor_id", "blob_offset", "blob_size",
                                       "ms1_rt_sigma", "ms1_rt_r2", "ms1_rt_skew",
                                       "ms1_im_sigma", "ms1_im_r2"]).to_pandas().reset_index(drop=True)
    print(f"precursors with local blobs: {len(df)}")
    if not args.all:
        idx = np.unique(np.linspace(0, len(df) - 1, min(args.n, len(df))).astype(int))
        df = df.iloc[idx].reset_index(drop=True)
    print(f"refitting {len(df)} precursors with {args.workers} workers ...")

    recs = df.to_dict("records")
    for r in recs:  # native types for pickling
        r["precursor_id"] = int(r["precursor_id"]); r["blob_offset"] = int(r["blob_offset"])
        r["blob_size"] = int(r["blob_size"])
    CHUNK = 1500
    chunks = [(str(EXTRACTED), RAW_FILE, recs[i:i + CHUNK]) for i in range(0, len(recs), CHUNK)]
    t0 = time.time(); rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for j, res in enumerate(ex.map(_fit_chunk, chunks)):
            rows.extend(res)
            if j % 10 == 0:
                print(f"  {len(rows)}/{len(df)}  ({time.time()-t0:.0f}s)")

    out = pd.DataFrame(rows)
    fit = out[out.emg_r2.notna()]

    def med(x):
        x = np.asarray(x, float); x = x[np.isfinite(x)]
        return f"{np.median(x):.3f}/{np.percentile(x,90):.3f}" if len(x) else "NA"

    n = len(out)
    print(f"\n=== results: {len(fit)}/{n} fit, in {time.time()-t0:.0f}s ===")
    print(f"sigma_ok:  {out.sigma_ok.sum()} ({100*out.sigma_ok.mean():.0f}%)   "
          f"lambda_ok: {out.lambda_ok.sum()} ({100*out.lambda_ok.mean():.0f}%)")
    sok = out[out.sigma_ok]
    print(f"EMG   r2 med/p90: {med(fit.emg_r2)}    Gauss r2 med/p90: {med(fit.gauss_r2)}")
    delta = (fit.emg_r2 - fit.gauss_r2).values
    print(f"EMG-Gauss r2 median: {np.median(delta):.3f}   frac EMG>=Gauss: {100*np.mean(delta>=-1e-6):.0f}%")
    print(f"[sigma_ok] EMG sigma med/p90: {med(sok.emg_sigma)} s   (Gauss sigma: {med(fit.gauss_sigma)})")
    print(f"[lambda_ok] EMG lambda med/p90: {med(out[out.lambda_ok].emg_lambda)}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT)
    print(f"\nwrote {OUT}  ({len(out)} rows; sigma label = emg_sigma where sigma_ok)")


if __name__ == "__main__":
    main()
