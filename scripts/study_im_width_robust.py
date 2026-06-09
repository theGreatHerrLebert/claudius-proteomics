#!/usr/bin/env python3
"""Settle the original hypothesis: is IM peak WIDTH sequence-determined?
Robust version — uses the model-free ms1_im_fwhm (NOT the poor Gaussian sigma),
gated to HIGH-INTENSITY observations (where mobilograms are reliable), with the
precursor_index merge de-duplicated. Cross-run peptide replication on PXD061039.
Reports the gated result vs ungated, and vs RT's ~20% reference.
"""
import glob
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

DATA = "/scratch/claudius-proteomics"
PXD = "PXD061039"


def eta_sq(y, groups):
    y = np.asarray(y, float); g = pd.Series(groups).values
    grand = y.mean(); sst = ((y - grand) ** 2).sum(); ssb = 0.0
    for _, idx in pd.Series(range(len(y))).groupby(g):
        yi = y[idx.values]; ssb += len(yi) * (yi.mean() - grand) ** 2
    return ssb / sst if sst > 0 else np.nan


def replication(g, label):
    g = g.copy()
    g["run"] = g.raw_file
    g["pep"] = g.sequence_normalized.astype(str) + "/" + g.charge.astype(int).astype(str)
    g["rel"] = g.ms1_im_fwhm / g.groupby("run").ms1_im_fwhm.transform("median")
    pr = g.groupby(["pep", "run"]).rel.median().reset_index()
    seen = pr.groupby("pep").run.nunique()
    multi = pr[pr.pep.isin(seen[seen >= 2].index)]
    if multi.pep.nunique() < 30:
        print(f"  [{label}] too few cross-run peptides"); return
    eta = eta_sq(multi.rel, multi.pep)
    rng = np.random.default_rng(0)
    base = np.mean([eta_sq(multi.rel.values, rng.permutation(multi.pep.values)) for _ in range(5)])
    genuine = (eta - base) / (1 - base)
    print(f"  [{label}] {multi.pep.nunique()} peps in>=2 runs | eta^2={eta:.3f} base={base:.3f} "
          f"=> genuine peptide-explained ~{genuine:.2f}")
    return genuine


def main():
    rf_path = (glob.glob(f"{DATA}/_mogon_raw_features/**/{PXD}/**/raw_features.parquet", recursive=True) +
               glob.glob(f"{DATA}/_mogon_raw_features/**/{PXD}/raw_features.parquet", recursive=True))[0]
    rf = pq.read_table(rf_path, columns=["precursor_id", "raw_file", "charge",
                                         "ms1_im_fwhm", "precursor_intensity"]).to_pandas()
    pi = pq.read_table(f"{DATA}/_mogon_precursor_index/{PXD}/precursor_index.parquet").to_pandas()
    pi = pi[(pi.sequence_normalized.notna()) & (pi.sage_qvalue <= 0.01)]
    pi = pi[["precursor_id", "sequence_normalized"]].drop_duplicates("precursor_id")   # DEDUP
    d = rf.merge(pi, on="precursor_id", how="inner")
    d = d[(d.ms1_im_fwhm > 0)]
    thr = d.precursor_intensity.quantile(2 / 3)     # high-intensity = top tercile
    print(f"{PXD}: {len(d)} confident im-fwhm obs; high-intensity thr={thr:.0f}")

    print("\n=== IM-FWHM cross-run sequence-determinism (model-free width) ===")
    replication(d, "ALL intensity")
    replication(d[d.precursor_intensity >= thr], "HIGH intensity only")
    print("  (RT reference ~0.20; IM-sigma-Gaussian earlier ~0.18 [unreliable])")

    # by charge, high intensity
    print("\n=== high-intensity, by charge ===")
    hi = d[d.precursor_intensity >= thr]
    for z in [2, 3]:
        replication(hi[hi.charge == z], f"charge {z}")


if __name__ == "__main__":
    main()
