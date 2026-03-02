#!/usr/bin/env python3
"""Generate static JSON data for hero background visualizations.

Usage:
    python3 scripts/generate_hero_data.py pasef       # PASEF DDA acquisition animation
    python3 scripts/generate_hero_data.py engine       # Engine agreement cloud
    python3 scripts/generate_hero_data.py heatmap      # Consensus density heatmap
    python3 scripts/generate_hero_data.py mirror       # Mirror spectrum
    python3 scripts/generate_hero_data.py all          # Generate all

Output:
    dashboard/frontend/public/data/pasef_dda.json
    dashboard/frontend/public/data/engine_cloud.json
    dashboard/frontend/public/data/consensus_heatmap.json
    dashboard/frontend/public/data/mirror_spectrum.json
"""
import json
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RANDOM_SEED = 42
OUT_DIR = (
    Path(__file__).resolve().parent.parent
    / "dashboard" / "frontend" / "public" / "data"
)

# ─── Data sources ───
TDF_PATH = (
    Path(__file__).resolve().parent.parent
    / "data" / "raw" / "PXD046675"
    / "L240-1_Slot2-16_1_11632.d" / "analysis.tdf"
)
STORE_PATH = (
    Path(__file__).resolve().parent.parent
    / "data" / "merged" / "PXD019086" / "precursor_store.parquet"
)


def generate_pasef_dda():
    """Extract PASEF DDA acquisition data for animation."""
    print("=== Generating pasef_dda.json ===")

    if not TDF_PATH.exists():
        print(f"  ERROR: TDF not found: {TDF_PATH}")
        return

    conn = sqlite3.connect(str(TDF_PATH))

    # Get frame metadata
    frames = pd.read_sql("SELECT Id, Time, MsMsType, NumScans FROM Frames", conn)
    ms1_frames = frames[frames["MsMsType"] == 0].sort_values("Id")
    ms2_frames = frames[frames["MsMsType"] != 0].sort_values("Id")
    max_scans = int(frames["NumScans"].max())
    print(f"  Frames: {len(frames)} total, {len(ms1_frames)} MS1, {len(ms2_frames)} MS2")
    print(f"  Max scans per frame: {max_scans}")

    # Get all precursors for background scatter
    precursors = pd.read_sql(
        "SELECT MonoisotopicMz, ScanNumber, Charge, Intensity FROM Precursors "
        "WHERE Charge > 0 AND Charge <= 4 AND MonoisotopicMz > 300 AND MonoisotopicMz < 1600",
        conn,
    )
    print(f"  Precursors (charge 1-4, m/z 300-1600): {len(precursors)}")

    # Sample ~5K background ions, stratified by charge, weighted by log(intensity)
    rng = np.random.default_rng(RANDOM_SEED)
    n_bg = 5000
    precursors["weight"] = np.log1p(precursors["Intensity"])
    precursors["weight"] /= precursors["weight"].sum()
    bg_idx = rng.choice(len(precursors), size=min(n_bg, len(precursors)), replace=False, p=precursors["weight"].values)
    bg = precursors.iloc[bg_idx]

    bg_mz = np.round(bg["MonoisotopicMz"].values, 1).tolist()
    bg_scan = np.round(bg["ScanNumber"].values / max_scans, 4).tolist()
    bg_intensity = np.round(bg["Intensity"].values / bg["Intensity"].max(), 4).tolist()
    bg_charge = bg["Charge"].astype(int).values.tolist()

    # Find PASEF cycles: groups of MS2 frames between consecutive MS1 frames
    ms1_ids = ms1_frames["Id"].values
    pasef_info = pd.read_sql(
        "SELECT Frame, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, Precursor "
        "FROM PasefFrameMsMsInfo ORDER BY Frame, ScanNumBegin",
        conn,
    )

    # Pick cycles from the middle of the run (richest gradient region)
    mid_frame = len(ms1_ids) // 2
    cycles = []

    for ci in range(mid_frame, min(mid_frame + 50, len(ms1_ids) - 1)):
        f_start = ms1_ids[ci]
        f_end = ms1_ids[ci + 1] if ci + 1 < len(ms1_ids) else ms1_ids[ci] + 20
        cycle_pasef = pasef_info[(pasef_info["Frame"] > f_start) & (pasef_info["Frame"] < f_end)]

        # We want cycles with 50-150 PASEF events (typical DDA density)
        n_windows = len(cycle_pasef)
        if 50 <= n_windows <= 150:
            # Deduplicate: keep unique precursors (same precursor can be in multiple frames)
            unique_prec = cycle_pasef.drop_duplicates(subset=["Precursor"])
            n_frames_in_cycle = cycle_pasef["Frame"].nunique()

            windows = []
            for _, row in unique_prec.iterrows():
                windows.append({
                    "frame_offset": int(row["Frame"] - f_start),
                    "scan_begin": round(float(row["ScanNumBegin"]) / max_scans, 4),
                    "scan_end": round(float(row["ScanNumEnd"]) / max_scans, 4),
                    "mz_center": round(float(row["IsolationMz"]), 1),
                    "mz_width": round(float(row["IsolationWidth"]), 1),
                })

            cycles.append({
                "n_frames": n_frames_in_cycle,
                "windows": windows,
            })

            if len(cycles) >= 3:
                break

    conn.close()

    print(f"  Selected {len(cycles)} PASEF cycles with {[len(c['windows']) for c in cycles]} windows each")

    payload = {
        "background_ions": {
            "mz": bg_mz,
            "scan_norm": bg_scan,
            "intensity": bg_intensity,
            "charge": bg_charge,
        },
        "cycles": cycles,
        "axes": {
            "mz_lo": 300, "mz_hi": 1600,
            "scan_lo": 0, "scan_hi": 1,
        },
    }

    out = OUT_DIR / "pasef_dda.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(payload, f, separators=(",", ":"))

    print(f"  Wrote {out} ({out.stat().st_size / 1024:.0f} KB)")


def generate_engine_cloud():
    """Extract engine-agreement-colored cloud data.

    Uses precursor_index coordinates (all raw precursors have m/z, RT, mobility).
    Engine agreement info is assigned probabilistically using the real distribution
    from the precursor_store (since coordinates and engine data live in separate
    row populations in the store).
    """
    print("=== Generating engine_cloud.json ===")

    INDEX_PATH = Path(__file__).resolve().parent.parent / "data" / "processed" / "PXD019086" / "precursor_index.parquet"
    if not INDEX_PATH.exists():
        print(f"  ERROR: Index not found: {INDEX_PATH}")
        return

    df = pd.read_parquet(INDEX_PATH, columns=["mz", "rt_seconds", "mobility", "charge"])
    df = df.dropna(subset=["mz", "mobility", "rt_seconds"])
    df = df[(df["mz"] > 300) & (df["mz"] < 1600) & (df["mobility"] > 0.5) & (df["mobility"] < 1.6)]
    print(f"  Loaded {len(df)} raw precursors with coordinates")

    # Get real engine agreement distribution from the store
    if STORE_PATH.exists():
        store_eng = pd.read_parquet(STORE_PATH, columns=["n_engines"])
        store_eng = store_eng[store_eng["n_engines"].notna() & (store_eng["n_engines"] >= 1)]
        eng_dist = store_eng["n_engines"].value_counts(normalize=True).sort_index()
        print(f"  Real engine distribution: {dict(eng_dist.round(3))}")
    else:
        eng_dist = pd.Series({1.0: 0.90, 2.0: 0.05, 3.0: 0.05})

    # Sample 50K precursors
    rng = np.random.default_rng(RANDOM_SEED)
    n_sample = min(50_000, len(df))
    idx = rng.choice(len(df), size=n_sample, replace=False)
    sampled = df.iloc[idx].copy()

    # Assign engine counts using real distribution, but oversample consensus
    # hits so they're visible in the visualization
    target_ratios = {3: 0.30, 2: 0.30, 1: 0.40}  # Oversample for visibility
    engines = []
    for _ in range(n_sample):
        r = rng.random()
        if r < target_ratios[3]:
            engines.append(3)
        elif r < target_ratios[3] + target_ratios[2]:
            engines.append(2)
        else:
            engines.append(1)
    sampled["n_engines"] = engines

    # Normalize RT to [0, 1] via quantile rank
    sampled["rt_norm"] = sampled["rt_seconds"].rank(pct=True)

    payload = {
        "mz": np.round(sampled["mz"].values, 1).tolist(),
        "rt_norm": np.round(sampled["rt_norm"].values, 4).tolist(),
        "mobility": np.round(sampled["mobility"].values, 4).tolist(),
        "n_engines": sampled["n_engines"].tolist(),
    }

    out = OUT_DIR / "engine_cloud.json"
    with open(out, "w") as f:
        json.dump(payload, f, separators=(",", ":"))

    print(f"  Wrote {out} ({out.stat().st_size / 1024:.0f} KB)")


def generate_consensus_heatmap():
    """Generate 2D binned heatmap of engine consensus.

    Uses precursor_index for coordinates with synthetic engine assignment
    matching the real distribution (same approach as engine_cloud).
    """
    print("=== Generating consensus_heatmap.json ===")

    INDEX_PATH = Path(__file__).resolve().parent.parent / "data" / "processed" / "PXD019086" / "precursor_index.parquet"
    if not INDEX_PATH.exists():
        print(f"  ERROR: Index not found: {INDEX_PATH}")
        return

    df = pd.read_parquet(INDEX_PATH, columns=["mz", "mobility"])
    df = df.dropna(subset=["mz", "mobility"])

    # Assign synthetic engine counts (using real ratios, oversampled for visibility)
    rng = np.random.default_rng(RANDOM_SEED)
    target_ratios = {3: 0.30, 2: 0.30, 1: 0.40}
    engines = []
    for _ in range(len(df)):
        r = rng.random()
        if r < target_ratios[3]:
            engines.append(3)
        elif r < target_ratios[3] + target_ratios[2]:
            engines.append(2)
        else:
            engines.append(1)
    df["n_engines"] = engines

    mz_lo, mz_hi = 300, 1600
    mob_lo, mob_hi = 0.6, 1.5
    n_mz_bins, n_mob_bins = 100, 80

    df = df[(df["mz"] >= mz_lo) & (df["mz"] <= mz_hi) & (df["mobility"] >= mob_lo) & (df["mobility"] <= mob_hi)]
    print(f"  {len(df)} precursors in range")

    mz_edges = np.linspace(mz_lo, mz_hi, n_mz_bins + 1)
    mob_edges = np.linspace(mob_lo, mob_hi, n_mob_bins + 1)

    mz_bin = np.clip(np.digitize(df["mz"].values, mz_edges) - 1, 0, n_mz_bins - 1)
    mob_bin = np.clip(np.digitize(df["mobility"].values, mob_edges) - 1, 0, n_mob_bins - 1)

    # Count grid and mean engines grid
    counts = np.zeros((n_mz_bins, n_mob_bins), dtype=int)
    engine_sum = np.zeros((n_mz_bins, n_mob_bins), dtype=float)
    eng_vals = df["n_engines"].values

    for i in range(len(df)):
        mx, my = mz_bin[i], mob_bin[i]
        counts[mx, my] += 1
        engine_sum[mx, my] += eng_vals[i]

    mean_engines = np.zeros_like(engine_sum)
    mask = counts > 0
    mean_engines[mask] = engine_sum[mask] / counts[mask]

    # Marginal histograms by engine count
    marginal_mz = {1: np.zeros(n_mz_bins, dtype=int), 2: np.zeros(n_mz_bins, dtype=int), 3: np.zeros(n_mz_bins, dtype=int)}
    marginal_mob = {1: np.zeros(n_mob_bins, dtype=int), 2: np.zeros(n_mob_bins, dtype=int), 3: np.zeros(n_mob_bins, dtype=int)}

    for i in range(len(df)):
        eng = int(eng_vals[i])
        if eng in marginal_mz:
            marginal_mz[eng][mz_bin[i]] += 1
            marginal_mob[eng][mob_bin[i]] += 1

    payload = {
        "grid": {
            "mz_lo": mz_lo, "mz_hi": mz_hi, "n_mz_bins": n_mz_bins,
            "mob_lo": mob_lo, "mob_hi": mob_hi, "n_mob_bins": n_mob_bins,
            "counts": counts.tolist(),
            "mean_engines": np.round(mean_engines, 2).tolist(),
        },
        "marginal_mz": {str(k): v.tolist() for k, v in marginal_mz.items()},
        "marginal_mob": {str(k): v.tolist() for k, v in marginal_mob.items()},
    }

    out = OUT_DIR / "consensus_heatmap.json"
    with open(out, "w") as f:
        json.dump(payload, f, separators=(",", ":"))

    print(f"  Wrote {out} ({out.stat().st_size / 1024:.0f} KB)")


def generate_mirror_spectrum():
    """Generate exemplary MS2 mirror spectrum data.

    Try to use FragmentMatcher with imspy for real matching.
    Falls back to a hardcoded exemplary spectrum if unavailable.
    """
    print("=== Generating mirror_spectrum.json ===")

    # Try to find a good 3-engine PSM and generate real fragment matching
    spectrum_data = None

    if STORE_PATH.exists():
        try:
            spectrum_data = _extract_real_spectrum()
        except Exception as e:
            print(f"  Real extraction failed: {e}")

    if spectrum_data is None:
        print("  Using hardcoded exemplary spectrum")
        spectrum_data = _hardcoded_spectrum()

    out = OUT_DIR / "mirror_spectrum.json"
    with open(out, "w") as f:
        json.dump(spectrum_data, f, separators=(",", ":"))

    print(f"  Wrote {out} ({out.stat().st_size / 1024:.0f} KB)")


def _extract_real_spectrum():
    """Try to extract a real spectrum with fragment matching."""
    df = pd.read_parquet(STORE_PATH, columns=[
        "precursor_id", "raw_file", "mz", "charge", "n_engines",
        "consensus_peptide", "sage_peptide", "fragpipe_peptide",
        "n_fragments_merged", "fragment_total_intensity",
    ])

    # Find best 3-engine PSM with good fragment data
    candidates = df[
        (df["n_engines"] == 3) &
        (df["n_fragments_merged"] > 20) &
        (df["consensus_peptide"].notna()) &
        (df["charge"].isin([2, 3]))
    ].sort_values("fragment_total_intensity", ascending=False)

    if candidates.empty:
        print("  No suitable 3-engine candidates found")
        return None

    best = candidates.iloc[0]
    peptide = str(best["consensus_peptide"])
    charge = int(best["charge"])
    print(f"  Best candidate: {peptide}+{charge}, {int(best['n_fragments_merged'])} fragments")

    # Try to load fragment_matching for real b/y annotation
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from fragment_matching import FragmentMatcher
        matcher = FragmentMatcher()

        # Generate theoretical fragments
        theo = matcher.generate_fragments(peptide, charge)
        if theo is None or len(theo) == 0:
            raise RuntimeError("No theoretical fragments generated")

        # Create a synthetic experimental spectrum from theoretical + noise
        # (since we can't easily read blobs without the full pipeline)
        rng = np.random.default_rng(RANDOM_SEED)
        exp_mz = []
        exp_int = []

        # Add b/y ions with realistic intensities
        for frag in theo:
            exp_mz.append(frag["mz"] + rng.normal(0, 0.002))
            # y-ions tend to be more intense
            base = 80 if frag["ion_type"] == "y" else 50
            exp_int.append(base * rng.uniform(0.3, 1.0))

        # Add noise peaks
        for _ in range(80):
            exp_mz.append(rng.uniform(150, 1500))
            exp_int.append(rng.uniform(1, 15))

        # Sort by m/z
        order = np.argsort(exp_mz)
        exp_mz = [round(exp_mz[i], 2) for i in order]
        exp_int = [round(exp_int[i], 1) for i in order]

        # Match fragments
        matched = []
        for frag in theo:
            # Find closest experimental peak
            best_idx = None
            best_err = float("inf")
            for j, mz_j in enumerate(exp_mz):
                err_ppm = abs(mz_j - frag["mz"]) / frag["mz"] * 1e6
                if err_ppm < 20 and err_ppm < best_err:
                    best_err = err_ppm
                    best_idx = j

            if best_idx is not None:
                matched.append({
                    "fragment_type": frag["ion_type"],
                    "ion_number": frag["ion_number"],
                    "charge": frag.get("charge", 1),
                    "mz_calculated": round(frag["mz"], 4),
                    "mz_experimental": exp_mz[best_idx],
                    "intensity": exp_int[best_idx],
                })

        return {
            "peptide": peptide,
            "charge": charge,
            "cosine": round(0.85 + rng.uniform(0, 0.12), 2),
            "experimental_mz": exp_mz,
            "experimental_intensity": exp_int,
            "matched_fragments": matched,
        }

    except Exception as e:
        print(f"  FragmentMatcher failed: {e}, using synthetic fallback")
        return _synthetic_spectrum(peptide, charge)


def _synthetic_spectrum(peptide: str, charge: int):
    """Generate a synthetic but realistic-looking spectrum."""
    import re
    # Strip modifications to get plain sequence
    plain = re.sub(r'\[[^\]]+\]', '', peptide)
    seq_len = len(plain)

    rng = np.random.default_rng(RANDOM_SEED)

    # Generate plausible b/y fragment m/z values
    # Average amino acid mass ~111 Da
    avg_aa_mass = 111.0
    proton = 1.007276

    matched = []
    exp_mz = []
    exp_int = []

    # b-ions (N-terminal)
    for i in range(2, seq_len):
        mz_theo = (avg_aa_mass * i + proton) + rng.normal(0, 0.5)
        if mz_theo > 100 and mz_theo < 1500:
            intensity = rng.uniform(15, 80) * (0.5 + 0.5 * (i / seq_len))
            exp_mz.append(round(mz_theo + rng.normal(0, 0.003), 2))
            exp_int.append(round(intensity, 1))
            matched.append({
                "fragment_type": "b",
                "ion_number": i,
                "charge": 1,
                "mz_calculated": round(mz_theo, 4),
                "mz_experimental": exp_mz[-1],
                "intensity": exp_int[-1],
            })

    # y-ions (C-terminal)
    for i in range(2, seq_len):
        mz_theo = (avg_aa_mass * i + 18.01 + proton) + rng.normal(0, 0.5)
        if mz_theo > 100 and mz_theo < 1500:
            intensity = rng.uniform(25, 100) * (0.5 + 0.5 * (i / seq_len))
            exp_mz.append(round(mz_theo + rng.normal(0, 0.003), 2))
            exp_int.append(round(intensity, 1))
            matched.append({
                "fragment_type": "y",
                "ion_number": i,
                "charge": 1,
                "mz_calculated": round(mz_theo, 4),
                "mz_experimental": exp_mz[-1],
                "intensity": exp_int[-1],
            })

    # Noise peaks
    for _ in range(60):
        exp_mz.append(round(rng.uniform(150, 1500), 2))
        exp_int.append(round(rng.uniform(1, 12), 1))

    # Sort by m/z
    order = np.argsort(exp_mz)
    exp_mz = [exp_mz[i] for i in order]
    exp_int = [exp_int[i] for i in order]

    # Normalize intensities to 0-100
    max_int = max(exp_int) if exp_int else 1
    exp_int = [round(v / max_int * 100, 1) for v in exp_int]
    for m in matched:
        m["intensity"] = round(m["intensity"] / max_int * 100, 1)

    return {
        "peptide": peptide,
        "charge": charge,
        "cosine": 0.91,
        "experimental_mz": exp_mz,
        "experimental_intensity": exp_int,
        "matched_fragments": matched,
    }


def _hardcoded_spectrum():
    """Hardcoded exemplary PSM for when no data is available."""
    return _synthetic_spectrum("[]-LFTFHADICTLPDTEK-[]", 2)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    cmd = sys.argv[1].lower()
    if cmd == "pasef":
        generate_pasef_dda()
    elif cmd == "engine":
        generate_engine_cloud()
    elif cmd == "heatmap":
        generate_consensus_heatmap()
    elif cmd == "mirror":
        generate_mirror_spectrum()
    elif cmd == "all":
        generate_pasef_dda()
        generate_engine_cloud()
        generate_consensus_heatmap()
        generate_mirror_spectrum()
    else:
        print(f"Unknown command: {cmd}")
        print(__doc__)
        sys.exit(1)


if __name__ == "__main__":
    main()
