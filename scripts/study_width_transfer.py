#!/usr/bin/env python3
"""Cross-LC transfer test: does a peptide's RT peak width transfer across
different gradients/datasets, once the per-dataset SCALE is removed?

Pilot on the 6-dataset blob subset (11->197 min gradients). For each dataset:
EMG-fit a sample of confident-ID precursors from its pulled run, get sigma per
(peptide,charge). Then:
  - per-dataset SCALE = median sigma over its peptides.
  - sigma_rel = sigma / scale.
  - For peptides shared across >=2 datasets: is sigma_rel consistent across
    datasets (peptide eta^2 on sigma_rel, BETWEEN datasets)? High => relative
    width is peptide-intrinsic and transfers => supports sigma ~ f(peptide)*scale(LC).
  - Does scale track gradient length G (p99 rt_seconds)?

Decides the modelling target: peptide relative-width + learned setup scale, vs
absolute-from-sequence (no transfer) vs sigma/G (wrong scale law).
"""
import glob
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from concurrent.futures import ProcessPoolExecutor

DATA = "/scratch/claudius-proteomics"
BLOBS = Path(f"{DATA}/_mogon_blobs")
RFROOT = Path(f"{DATA}/_mogon_raw_features")
sys.path.insert(0, "/tmp/ssprobe")           # materialised emg_refit_poc (feat branch)
sys.path.insert(0, "/home/administrator/Documents/promotion/claudius-proteomics")
from emg_refit_poc import _fit_chunk          # noqa: E402

N_PER_DS = 6000   # confident precursors to fit per dataset (speed vs signal)
RFCOLS = ["precursor_id", "raw_file", "blob_offset", "blob_size", "rt_seconds", "charge",
          "ms1_rt_sigma", "ms1_rt_r2", "ms1_rt_skew", "ms1_im_sigma", "ms1_im_r2"]


def find_rf(pxd):
    hits = glob.glob(str(RFROOT / "**" / pxd / "**" / "raw_features.parquet"), recursive=True)
    hits += glob.glob(str(RFROOT / "**" / pxd / "raw_features.parquet"), recursive=True)
    return hits[0] if hits else None


def eta_sq(y, groups):
    y = np.asarray(y, float); g = pd.Series(groups).values
    grand = y.mean(); sst = ((y - grand) ** 2).sum(); ssb = 0.0
    for _, idx in pd.Series(range(len(y))).groupby(g):
        yi = y[idx.values]; ssb += len(yi) * (yi.mean() - grand) ** 2
    return ssb / sst if sst > 0 else np.nan


def fit_dataset(blob_path):
    p = Path(blob_path)
    d_dir = p.parent                       # <raw_file>.d
    raw_file_d = d_dir.name                 # with .d (matches raw_features)
    raw_file = raw_file_d[:-2] if raw_file_d.endswith(".d") else raw_file_d
    extracted_dir = str(d_dir.parent)
    pxd = next(s for s in p.parts if s.startswith("PXD"))

    rf_path = find_rf(pxd)
    if not rf_path:
        print(f"  {pxd}: no raw_features", flush=True); return None
    rf = pq.read_table(rf_path, columns=RFCOLS).to_pandas()
    rf = rf[rf.raw_file == raw_file_d]
    if len(rf) == 0:  # try without .d
        rf = pq.read_table(rf_path, columns=RFCOLS).to_pandas()
        rf = rf[rf.raw_file.str.replace(".d", "", regex=False) == raw_file]
    G = float(np.nanpercentile(rf.rt_seconds, 99))

    # confident peptide IDs for this run
    pidx_path = None
    for cand in (f"{DATA}/_mogon_precursor_index/{pxd}/precursor_index.parquet",
                 f"{DATA}/{pxd}/precursor_index.parquet"):
        if Path(cand).exists():
            pidx_path = cand; break
    if pidx_path is None:
        print(f"  {pxd}: no precursor_index", flush=True); return None
    pidx = pq.read_table(pidx_path).to_pandas()
    pidx = pidx[(pidx.raw_file == raw_file) | (pidx.raw_file == raw_file_d)]
    # schema-agnostic peptide ID: sequence_normalized (present in all variants,
    # and a clean cross-dataset key); confidence via n_engines (>=1 = identified).
    keep = pidx[pidx["sequence_normalized"].notna()].copy()
    if "n_engines" in keep.columns:
        keep = keep[keep.n_engines >= 1]
    pidx = keep[["precursor_id", "sequence_normalized"]].rename(columns={"sequence_normalized": "peptide"})

    m = rf.merge(pidx, on="precursor_id", how="inner")
    if len(m) == 0:
        print(f"  {pxd}: 0 confident precursors matched", flush=True); return None
    if len(m) > N_PER_DS:
        m = m.iloc[np.linspace(0, len(m) - 1, N_PER_DS).astype(int)]
    recs = m.to_dict("records")
    for r in recs:
        r["precursor_id"] = int(r["precursor_id"]); r["blob_offset"] = int(r["blob_offset"]); r["blob_size"] = int(r["blob_size"])
    # fit in chunks (parallel)
    CH = 1500
    chunks = [(extracted_dir, raw_file_d, recs[i:i+CH]) for i in range(0, len(recs), CH)]
    with ProcessPoolExecutor(max_workers=8) as ex:
        fitted = [x for chunk in ex.map(_fit_chunk, chunks) for x in chunk]
    fdf = pd.DataFrame(fitted)
    fdf = fdf[fdf.sigma_ok][["precursor_id", "emg_sigma"]]
    out = m.merge(fdf, on="precursor_id", how="inner")
    out["dataset"] = pxd; out["G"] = G
    out["pep"] = out.peptide.astype(str) + "/" + out.charge.astype(int).astype(str)
    print(f"  {pxd}: G={G:.0f}s, {len(out)} sigma_ok fits, {out.pep.nunique()} peptides", flush=True)
    return out[["dataset", "G", "pep", "emg_sigma"]]


def main():
    blobs = sorted(glob.glob(str(BLOBS / "**" / "blobs.bin"), recursive=True))
    print(f"pilot blob files: {len(blobs)}")
    parts = [fit_dataset(b) for b in blobs]
    parts = [p for p in parts if p is not None and len(p)]
    if len(parts) < 2:
        print("need >=2 datasets fitted"); return
    df = pd.concat(parts, ignore_index=True)

    # per-dataset scale + relative sigma
    scale = df.groupby("dataset").emg_sigma.median().rename("scale")
    df = df.merge(scale, on="dataset")
    df["sigma_rel"] = df.emg_sigma / df.scale
    # per (dataset, peptide) median
    pe = df.groupby(["dataset", "pep"]).agg(sigma=("emg_sigma", "median"),
                                            sigma_rel=("sigma_rel", "median")).reset_index()

    print("\n=== per-dataset scale vs gradient ===")
    sg = df.groupby("dataset").agg(G=("G", "first"), scale=("scale", "first"), n=("pep", "size")).sort_values("G")
    print(sg.to_string(formatters={"G": "{:.0f}".format, "scale": "{:.2f}".format}))
    print(f"  corr(G, scale) = {np.corrcoef(sg.G, sg.scale)[0,1]:.3f}")

    # transfer: peptides shared across >=2 datasets
    shared = pe.groupby("pep").filter(lambda g: g.dataset.nunique() >= 2)
    n_shared = shared.pep.nunique()
    print(f"\n=== TRANSFER TEST: {n_shared} peptides shared across >=2 datasets ===")
    if n_shared >= 20:
        # peptide eta^2 BETWEEN datasets: absolute vs relative sigma
        eta_abs = eta_sq(shared.sigma, shared.pep)
        eta_rel = eta_sq(shared.sigma_rel, shared.pep)
        print(f"  peptide eta^2 on ABSOLUTE sigma : {eta_abs:.3f}")
        print(f"  peptide eta^2 on RELATIVE sigma : {eta_rel:.3f}  "
              f"(higher => peptide relative-width TRANSFERS across LC)")
        # within-peptide CV of relative sigma (lower = better transfer)
        cv = shared.groupby("pep").sigma_rel.apply(lambda x: x.std() / x.mean() if x.mean() else np.nan)
        print(f"  within-peptide CV of sigma_rel (median): {np.nanmedian(cv):.3f}")
    else:
        print("  too few shared peptides for a stable test (diverse samples).")
    df.to_parquet(f"{DATA}/_mogon_blobs/_transfer_fits.parquet")
    print(f"\nwrote {DATA}/_mogon_blobs/_transfer_fits.parquet")


if __name__ == "__main__":
    main()
