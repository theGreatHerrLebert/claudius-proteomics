#!/usr/bin/env python3
"""IM peak-width study (analog of the RT one). Hypotheses:
  (1) IM width is more SEQUENCE-determined than RT (mobility ~ molecular structure)
      -> cross-run peptide replication eta^2 should beat RT's ~20%.
  (2) multi-modality (gas-phase conformers) is charge-dependent, esp. 3+
      -> single-Gaussian fit quality im_r2 worse / skew higher for charge 3+.

On PXD061039 (12 runs, ~981k confident IDs, full schema) — single file, local, fast.
Uses raw_features (ms1_im_*) + precursor_index. Raw-mobilogram bimodality (decisive
multimodality) is a separate blob-based step.
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


def main():
    rf_path = (glob.glob(f"{DATA}/_mogon_raw_features/**/{PXD}/**/raw_features.parquet", recursive=True) +
               glob.glob(f"{DATA}/_mogon_raw_features/**/{PXD}/raw_features.parquet", recursive=True))[0]
    rf = pq.read_table(rf_path, columns=["precursor_id", "raw_file", "charge",
                                         "ms1_im_sigma", "ms1_im_r2", "ms1_im_skew",
                                         "ms1_rt_r2"]).to_pandas()
    pi = pq.read_table(f"{DATA}/_mogon_precursor_index/{PXD}/precursor_index.parquet").to_pandas()
    pi = pi[(pi.sequence_normalized.notna()) & (pi.sage_qvalue <= 0.01)][["precursor_id", "sequence_normalized"]]
    d = rf.merge(pi, on="precursor_id", how="inner")
    g = d[(d.ms1_im_r2 >= 0.8) & (d.ms1_im_sigma > 0)].copy()      # quality-gated for the width analysis
    print(f"{PXD}: {len(rf)} precursors, {len(d)} confident, {len(g)} im-quality-gated; {d.raw_file.nunique()} runs")

    # (1) cross-run replication of im_sigma
    g["run"] = g.raw_file
    g["pep"] = g.sequence_normalized.astype(str) + "/" + g.charge.astype(int).astype(str)
    g["sig_rel"] = g.ms1_im_sigma / g.groupby("run").ms1_im_sigma.transform("median")
    pr = g.groupby(["pep", "run"]).sig_rel.median().reset_index()
    seen = pr.groupby("pep").run.nunique()
    multi = pr[pr.pep.isin(seen[seen >= 2].index)]
    eta = eta_sq(multi.sig_rel, multi.pep)
    rng = np.random.default_rng(0)
    base = np.mean([eta_sq(multi.sig_rel.values, rng.permutation(multi.pep.values)) for _ in range(5)])
    genuine = (eta - base) / (1 - base)
    print(f"\n(1) IM-sigma cross-run replication: {multi.pep.nunique()} peptides in >=2 runs")
    print(f"    eta^2={eta:.3f}  chance-baseline={base:.3f}  => genuine peptide-explained ~{genuine:.2f}")
    print(f"    (RT was ~0.20; IM {'HIGHER -> more sequence-determined' if genuine > 0.25 else 'similar/lower'})")

    # (2) multi-modality proxy by charge: single-Gaussian fit quality
    print("\n(2) single-Gaussian IM fit quality by charge (low r2 / high |skew| => multi-modal):")
    print(f"    {'z':>3}{'n':>9}{'im_r2_med':>11}{'frac_r2<0.8':>13}{'|im_skew|_med':>15}{'im_sigma_med':>14}")
    d["zc"] = d.charge.clip(1, 5)
    for z, gg in d.groupby("zc"):
        gg2 = gg[gg.ms1_im_sigma > 0]
        print(f"    {int(z):>3}{len(gg):>9}{gg.ms1_im_r2.median():>11.3f}"
              f"{(gg.ms1_im_r2 < 0.8).mean()*100:>12.0f}%{gg.ms1_im_skew.abs().median():>15.3f}"
              f"{gg2.ms1_im_sigma.median():>14.4f}")
    print("    (compare: same table for RT fit quality)")
    print(f"    {'z':>3}{'rt_r2_med':>11}{'frac rt_r2<0.8':>16}")
    for z, gg in d.groupby("zc"):
        print(f"    {int(z):>3}{gg.ms1_rt_r2.median():>11.3f}{(gg.ms1_rt_r2 < 0.8).mean()*100:>15.0f}%")


if __name__ == "__main__":
    main()
