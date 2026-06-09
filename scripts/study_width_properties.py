#!/usr/bin/env python3
"""Property-mediated width study: does RT peak width sigma track peptide LENGTH,
consistently across datasets, and does it survive controlling for RT-position?

Combines 3 datasets at different gradients (no shared peptides needed):
  - frac01 PXD019086 (110 min): full per-precursor data incl. rt_seconds.
  - PXD040716 (19 min) + PXD061039 (60 min): from the pilot _transfer_fits.

Each dataset is self-normalised (sigma_rel = sigma / per-dataset median). Tests:
  1. corr(length, sigma_rel) per dataset — and is it CONSISTENT across datasets?
     (consistency => relative-width-vs-length transfers via the property.)
  2. control: does length still explain sigma after RT-position is removed
     (frac01 only, has rt) — i.e. is it length or just elution position?
"""
import re
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

DATA = "/scratch/claudius-proteomics"
BASE = f"{DATA}/_extracted_poc/PXD019086"
RAW = "20190504_TIMS1_FlMe_SA_HeLa_frac01_A10_1_93"


def plen(s):
    s = re.sub(r"\[[^\]]*\]", "", str(s))   # strip [UNIMOD:x] / [+57]
    s = re.sub(r"\([^)]*\)", "", s)
    return sum(c.isalpha() and c.isupper() for c in s)


def frac01_df():
    em = pq.read_table(f"{BASE}/emg_refit_frac01.parquet").to_pandas()
    rf = pq.read_table(f"{BASE}/raw_features.parquet",
                       columns=["precursor_id", "rt_seconds"]).to_pandas()
    pi = pq.read_table(f"{DATA}/PXD019086/precursor_index.parquet").to_pandas()
    pi = pi[pi.raw_file == RAW][["precursor_id", "sage_modified", "sage_qvalue"]]
    d = em.merge(rf, on="precursor_id").merge(pi, on="precursor_id")
    d = d[d.sigma_ok & (d.sage_qvalue <= 0.01)].dropna(subset=["sage_modified"]).copy()
    d["length"] = d.sage_modified.map(plen)
    d = d.rename(columns={"emg_sigma": "sigma"})
    d["dataset"] = "PXD019086(110min)"
    return d[["dataset", "sigma", "length", "rt_seconds"]]


def transfer_df():
    t = pq.read_table(f"{DATA}/_mogon_blobs/_transfer_fits.parquet").to_pandas()
    t["length"] = t.pep.map(lambda p: plen(p.rsplit("/", 1)[0]))
    g = {"PXD040716": "PXD040716(19min)", "PXD061039": "PXD061039(60min)"}
    t["dataset"] = t.dataset.map(lambda x: g.get(x, x))
    t = t.rename(columns={"emg_sigma": "sigma"})
    t["rt_seconds"] = np.nan
    return t[["dataset", "sigma", "length", "rt_seconds"]]


def main():
    df = pd.concat([frac01_df(), transfer_df()], ignore_index=True)
    df = df[(df.length >= 5) & (df.length <= 40) & np.isfinite(df.sigma)]
    # self-normalise each dataset
    df["sigma_rel"] = df.sigma / df.groupby("dataset").sigma.transform("median")

    print("=== sigma_rel vs peptide length, per dataset ===")
    print(f"{'dataset':<20}{'n':>7}{'corr(len,sig_rel)':>18}{'slope/aa':>10}{'sig_rel @len10':>14}{'@len25':>9}")
    slopes = {}
    for ds, g in df.groupby("dataset"):
        r = np.corrcoef(g.length, g.sigma_rel)[0, 1]
        sl, intc = np.polyfit(g.length, g.sigma_rel, 1)
        slopes[ds] = sl
        p10 = sl * 10 + intc; p25 = sl * 25 + intc
        print(f"{ds:<20}{len(g):>7}{r:>18.3f}{sl:>10.4f}{p10:>14.2f}{p25:>9.2f}")
    sv = np.array(list(slopes.values()))
    print(f"\nslope consistency across datasets: mean={sv.mean():.4f} sd={sv.std():.4f} "
          f"(CV={sv.std()/abs(sv.mean()):.2f}; all same sign={bool(np.all(sv>0) or np.all(sv<0))})")

    # control for RT-position (frac01 only)
    f = df[df.dataset.str.startswith("PXD019086")].copy()
    f["rt_bin"] = pd.qcut(f.rt_seconds, 10, labels=False, duplicates="drop")
    raw_r = np.corrcoef(f.length, f.sigma_rel)[0, 1]
    # within-rt-bin (partial): residualise sigma_rel and length on rt_bin means, then corr
    f["sig_res"] = f.sigma_rel - f.groupby("rt_bin").sigma_rel.transform("mean")
    f["len_res"] = f.length - f.groupby("rt_bin").length.transform("mean")
    part_r = np.corrcoef(f.len_res, f.sig_res)[0, 1]
    print("\n=== control for RT-position (frac01 PXD019086) ===")
    print(f"  raw corr(length, sigma_rel)            = {raw_r:.3f}")
    print(f"  partial corr | RT-position removed     = {part_r:.3f}  "
          f"({'mostly elution-position' if abs(part_r) < 0.5*abs(raw_r) else 'length effect survives'})")
    # length->sigma curve (frac01)
    print("\n  sigma_rel by length decile (frac01):")
    f["lbin"] = pd.qcut(f.length, 8, labels=False, duplicates="drop")
    for b, gg in f.groupby("lbin"):
        print(f"    len~{gg.length.median():.0f}: sigma_rel={gg.sigma_rel.median():.2f}  (n={len(gg)})")


if __name__ == "__main__":
    main()
