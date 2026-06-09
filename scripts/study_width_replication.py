#!/usr/bin/env python3
"""DECISIVE test: does a peptide's RT peak width REPLICATE across runs?

If sigma is genuinely peptide-determined, the same (peptide,charge) measured in
different runs of one dataset should get a consistent width (after removing the
per-run scale). If not, sigma is mostly noise/LC and per-peptide prediction has a
low ceiling. Removes the within-run shared-elution confound that inflated the
earlier 62%.

Uses raw_features (all runs, Gaussian ms1_rt_sigma — biased absolute but fine for
a consistency test) + precursor_index (peptide IDs). No blobs / no MOGON needed.
Calibrates eta^2 against a label-shuffled baseline.
"""
import glob
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

DATA = "/scratch/claudius-proteomics"
RFROOT = f"{DATA}/_mogon_raw_features"


def eta_sq(y, groups):
    y = np.asarray(y, float); g = pd.Series(groups).values
    grand = y.mean(); sst = ((y - grand) ** 2).sum(); ssb = 0.0
    for _, idx in pd.Series(range(len(y))).groupby(g):
        yi = y[idx.values]; ssb += len(yi) * (yi.mean() - grand) ** 2
    return ssb / sst if sst > 0 else np.nan


def run_dataset(pxd):
    rf_hits = glob.glob(f"{RFROOT}/**/{pxd}/**/raw_features.parquet", recursive=True) + \
              glob.glob(f"{RFROOT}/**/{pxd}/raw_features.parquet", recursive=True)
    pi_path = next((p for p in (f"{DATA}/_mogon_precursor_index/{pxd}/precursor_index.parquet",
                                 f"{DATA}/{pxd}/precursor_index.parquet") if glob.glob(p)), None)
    if not rf_hits or pi_path is None:
        print(f"{pxd}: missing inputs"); return
    rf = pq.read_table(rf_hits[0], columns=["precursor_id", "raw_file", "ms1_rt_sigma",
                                            "ms1_rt_r2", "charge"]).to_pandas()
    pi = pq.read_table(pi_path).to_pandas()
    if "sage_qvalue" in pi.columns:
        pi = pi[(pi.sequence_normalized.notna()) & (pi.sage_qvalue <= 0.01)]
    else:
        pi = pi[pi.sequence_normalized.notna()]
    pi = pi[["precursor_id", "sequence_normalized"]]
    d = rf.merge(pi, on="precursor_id", how="inner")
    d = d[(d.ms1_rt_r2 >= 0.8) & (d.ms1_rt_sigma > 0)].copy()
    d["run"] = d.raw_file
    d["pep"] = d.sequence_normalized.astype(str) + "/" + d.charge.astype(int).astype(str)
    # remove per-run scale
    d["sig_rel"] = d.ms1_rt_sigma / d.groupby("run").ms1_rt_sigma.transform("median")

    nruns = d.run.nunique()
    # collapse to one value per (pep, run) then keep peptides seen in >=2 RUNS
    pr = d.groupby(["pep", "run"]).sig_rel.median().reset_index()
    seen = pr.groupby("pep").run.nunique()
    multi = pr[pr.pep.isin(seen[seen >= 2].index)]
    npep = multi.pep.nunique()
    print(f"\n=== {pxd}: {nruns} runs, {len(d)} good fits, "
          f"{npep} peptides in >=2 runs ({len(multi)} obs) ===")
    if npep < 30:
        print("  too few cross-run peptides (likely fractionated, not replicate runs)."); return

    eta = eta_sq(multi.sig_rel, multi.pep)
    # shuffled baseline: permute peptide labels, recompute eta^2 (chance level)
    rng = np.random.default_rng(0)
    shuf = [eta_sq(multi.sig_rel.values, rng.permutation(multi.pep.values)) for _ in range(5)]
    # within-peptide CV across runs
    cv = multi.groupby("pep").sig_rel.apply(lambda x: x.std() / x.mean() if x.mean() else np.nan)
    print(f"  peptide eta^2 on sig_rel (cross-run) : {eta:.3f}")
    print(f"  shuffled-baseline eta^2 (chance)     : {np.mean(shuf):.3f} +/- {np.std(shuf):.3f}")
    print(f"  within-peptide CV of sig_rel (median): {np.nanmedian(cv):.3f}")
    verdict = ("REPLICATES (peptide signal real)" if eta > np.mean(shuf) + 5 * (np.std(shuf) + 1e-6) and eta > 0.15
               else "DOES NOT replicate beyond chance -> sigma ~ noise/LC, low per-peptide ceiling")
    print(f"  => {verdict}")


def main():
    for pxd in ["PXD061039", "PXD040716", "PXD019086"]:
        run_dataset(pxd)


if __name__ == "__main__":
    main()
