#!/usr/bin/env python3
"""Generate pre-computed ion cloud data for the San José landing page visualization.

Loads the ionmob CCS dataset from HuggingFace, stratified-samples 50K peptide ions
(charges 2-4), and writes a compact columnar JSON for the Canvas scatter renderer.

Usage:
    python3 scripts/generate_ion_cloud_data.py

Output:
    dashboard/frontend/public/data/ion_cloud.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from datasets import load_dataset


N_SAMPLE = 50_000
CHARGES = [2, 3, 4]
RANDOM_SEED = 42

OUT_PATH = (
    Path(__file__).resolve().parent.parent
    / "dashboard"
    / "frontend"
    / "public"
    / "data"
    / "ion_cloud.json"
)


def main():
    print("Loading ionmob dataset from HuggingFace...")
    ds = load_dataset("theGreatHerrLebert/ionmob", split="train")
    df = ds.to_pandas()
    print(f"  Loaded {len(df):,} rows, columns: {list(df.columns)}")

    # Filter to charges 2-4
    df = df[df["charge"].isin(CHARGES)].copy()
    print(f"  After charge filter (2-4): {len(df):,} rows")

    # Stratified sample preserving charge distribution
    sample_idx = df.groupby("charge", group_keys=False).apply(
        lambda g: g.sample(
            n=min(len(g), int(N_SAMPLE * len(g) / len(df))),
            random_state=RANDOM_SEED,
        ),
    ).index
    sampled = df.loc[sample_idx]

    # Top up if rounding caused fewer than N_SAMPLE
    if len(sampled) < N_SAMPLE:
        remaining = df.drop(sampled.index)
        extra = remaining.sample(
            min(N_SAMPLE - len(sampled), len(remaining)),
            random_state=RANDOM_SEED,
        )
        sampled = pd.concat([sampled, extra])

    sampled = sampled.head(N_SAMPLE)
    print(f"  Sampled {len(sampled):,} points")

    # Normalize RT per dataset_origin to [0, 1] using quantile rank
    # (robust to outliers and vastly different gradient lengths)
    sampled["rt_norm"] = sampled.groupby("dataset_origin")["rt"].transform(
        lambda x: x.rank(pct=True) if len(x) > 1 else 0.5
    )
    n_origins = sampled["dataset_origin"].nunique()
    print(f"  Normalized RT (quantile rank) across {n_origins} dataset origins")

    # Round for compactness
    mz_arr = np.round(sampled["mz"].values, 1).tolist()
    inv_mob_arr = np.round(sampled["inv_ion_mob"].values, 4).tolist()
    rt_arr = np.round(sampled["rt_norm"].values, 4).tolist()
    charge_arr = sampled["charge"].astype(int).values.tolist()

    # Write columnar JSON
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "mz": mz_arr,
        "inv_mob": inv_mob_arr,
        "rt": rt_arr,
        "charge": charge_arr,
    }
    with open(OUT_PATH, "w") as f:
        json.dump(payload, f, separators=(",", ":"))

    size_kb = OUT_PATH.stat().st_size / 1024
    print(f"\nWrote {len(mz_arr):,} points to {OUT_PATH}")
    print(f"  File size: {size_kb:.0f} KB")
    for c in CHARGES:
        n = charge_arr.count(c)
        print(f"  Charge {c}+: {n:,} ({100 * n / len(charge_arr):.1f}%)")

    # Print data ranges for reference
    print(f"\n  m/z range: {min(mz_arr):.1f} – {max(mz_arr):.1f}")
    print(f"  1/K0 range: {min(inv_mob_arr):.4f} – {max(inv_mob_arr):.4f}")
    print(f"  RT range: {min(rt_arr):.2f} – {max(rt_arr):.2f}")


if __name__ == "__main__":
    main()
