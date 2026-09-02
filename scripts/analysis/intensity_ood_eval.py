#!/usr/bin/env python3
"""Out-of-distribution evaluation for the production intensity-predictor weights.

Apply every available weights file to deposits that were **deliberately
excluded from training** and measure held-out spectral angle. Tests
generalization beyond the 3 clean training deposits.

Probes:
  PXD019086 (Meier 2021) — early-timsTOF, noisy by today's standards
                            (see [[predictor-training-era]] memory)
  SIM_SMOKE              — synthetic; near-perfect ceiling on the training
                            set, useful as an upper-bound sanity check

For each (probe, weights) pair we score under the paradigm-matched
inference rule:
  - "uniform" rule for v1 and pretrained: predict at a single
    pretrained-calibrated NCE on this probe
  - "per_psm_optimal" rule for v4 variants: predict at the per-PSM
    NCE found by a pretrained sweep on this probe

Each weights file is loaded from disk; pretrained uses
load_deep_intensity_predictor() (hub-managed defaults).

Output:
  notebook/analysis/production_finetune/intensity_ood_eval_<YYYY-MM-DD>.json
"""
from __future__ import annotations

import json
import os
import sys
import time
import warnings
from datetime import date
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc
import torch

from sagepy_rescore.sage_loader import load_psms_from_sage_parquets
from sagepy_rescore.predictors import apply_imspy_patches
from sagepy.utility import psm_collection_to_pandas
from imspy_predictors.intensity.predictors import (
    DeepPeptideIntensityPredictor,
    calibrate_nce,
)
from imspy_predictors.lazy_imports import get_sagepy_fragment_utils

apply_imspy_patches()
associate_fragments, _ = get_sagepy_fragment_utils()


# --- config ----------------------------------------------------------------
DATA_ROOT    = Path(os.environ.get("CLAUDIUS_DATA_ROOT", "/scratch/claudius-proteomics"))
PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR      = PROJECT_ROOT / "notebook" / "analysis" / "production_finetune"
WEIGHTS_DIR  = PROJECT_ROOT / "models" / "production"

PROBES = {
    "PXD019086 (Meier)": {
        "sage":  DATA_ROOT / "PXD019086" / "engines" / "sage",
        "store": DATA_ROOT / "PXD019086" / "precursor_store.parquet",
        "note":  "early-timsTOF (2019), noisy by today's standards",
    },
    "SIM_SMOKE": {
        "sage":  DATA_ROOT / "SIM_SMOKE_human_trypsin" / "engines" / "sage",
        "store": DATA_ROOT / "SIM_SMOKE_human_trypsin" / "precursor_store.parquet",
        "note":  "synthetic; effectively an upper-bound sanity check",
    },
}

# Weights variants. Pretrained uses the default hub load; the others are
# state_dicts from `models/production/`. Variants whose .pt is missing are
# skipped gracefully (helpful while the wider-grid run hasn't landed yet).
WEIGHTS = {
    "pretrained":     None,
    "v1":             WEIGHTS_DIR / "imspy_intensity_finetune_2026-05-20.pt",
    "v4_grid25-65":   WEIGHTS_DIR / "imspy_intensity_finetune_v4_grid25-65_2026-05-20.pt",
    "v4_grid15-75":   WEIGHTS_DIR / "imspy_intensity_finetune_v4_2026-05-20.pt",
}

# Each variant's inference paradigm (matches what it was trained with)
PARADIGM = {
    "pretrained":     "uniform",
    "v1":             "uniform",
    "v4_grid25-65":   "per_psm_optimal",
    "v4_grid15-75":   "per_psm_optimal",
}

PEP_THRESHOLD      = 1e-3
FT_QVALUE          = 0.01
MIN_MATCHED_PEAKS  = 6
EVAL_MAX_PSMS      = 5000           # cap per probe for tractable per-PSM sweep
NCE_GRID           = list(range(15, 76, 2))
RANDOM_SEED        = 0
STAMP              = date.today().isoformat()


# --- helpers ---------------------------------------------------------------
def filter_results(res, *, qval):
    mask = pc.greater_equal(res["matched_peaks"], MIN_MATCHED_PEAKS)
    mask = pc.and_(mask, pc.equal(res["is_decoy"], False))
    mask = pc.and_(mask, pc.less_equal(res["rank"], 1))
    mask = pc.and_(mask, pc.less_equal(res["spectrum_q"], qval))
    return res.filter(mask)


def load_probe_psms(probe_name, meta, n_max=EVAL_MAX_PSMS):
    """Load + filter PSMs for a probe deposit. Caches the filtered subset
    next to the probe's sage outputs so re-runs are fast."""
    sage = meta["sage"]
    cache = OUT_DIR / "ood_cache" / probe_name.replace(" ", "_").replace("(", "").replace(")", "")
    cache.mkdir(parents=True, exist_ok=True)
    rp = cache / "results.parquet"
    fp = cache / "fragments.parquet"

    if not (rp.exists() and fp.exists()):
        print(f"  [{probe_name}] cache miss; filtering...", flush=True)
        res  = pq.read_table(sage / "results.sage.parquet")
        frag = pq.read_table(sage / "matched_fragments.sage.parquet")
        res  = filter_results(res, qval=FT_QVALUE)
        if n_max and res.num_rows > n_max:
            idx = np.sort(np.random.RandomState(RANDOM_SEED).choice(
                res.num_rows, n_max, replace=False))
            res = res.take(pa.array(idx))
        keep = pa.array(list(set(res.column("psm_id").to_pylist())))
        frag = frag.filter(pc.is_in(frag["psm_id"], value_set=keep))
        pq.write_table(res, rp); pq.write_table(frag, fp)

    psms = load_psms_from_sage_parquets(str(rp), str(fp),
                                        default_collision_energy=30.0, max_rank=1)
    print(f"  [{probe_name}] {len(psms):,} PSMs loaded", flush=True)
    return psms


def intensity_mean_sa(model, psms):
    if not psms:
        return float("nan")
    seqs    = [p.sequence_modified for p in psms]
    charges = np.array([p.charge for p in psms])
    ce      = [float(p.collision_energy) for p in psms]
    intens = model.predict_intensities(seqs, charges, ce, batch_size=2048, flatten=True)
    coll = associate_fragments(psms, intens, num_threads=os.cpu_count() or 1)
    for p, c in zip(psms, coll):
        p.prosit_predicted_intensities = c.prosit_predicted_intensities
    df = psm_collection_to_pandas(psms)
    return float(df["spectral_angle_similarity"].mean())


def per_psm_optimal_ce(predictor, psms, nce_grid):
    n = len(psms)
    seqs = [p.sequence_modified for p in psms]
    charges = np.array([p.charge for p in psms])
    sa = np.zeros((n, len(nce_grid)), dtype=np.float32)
    for j, nce in enumerate(nce_grid):
        intens = predictor.predict_intensities(
            seqs, charges, [float(nce)] * n, batch_size=2048, flatten=True)
        coll = associate_fragments(psms, intens, num_threads=os.cpu_count() or 1)
        for p, c in zip(psms, coll):
            p.prosit_predicted_intensities = c.prosit_predicted_intensities
        df = psm_collection_to_pandas(psms)
        sa[:, j] = df["spectral_angle_similarity"].to_numpy(np.float32)
    best_idx = sa.argmax(axis=1)
    return np.array([nce_grid[i] for i in best_idx], dtype=float), sa


def load_predictor(weights_path):
    """Load a fresh predictor and optionally apply a fine-tuned state_dict."""
    p = DeepPeptideIntensityPredictor(verbose=False)
    if weights_path is not None:
        sd = torch.load(str(weights_path), map_location='cpu', weights_only=True)
        p.model.load_state_dict(sd)
        p.model.eval()
    return p


# --- main ------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] OOD intensity eval (stamp={STAMP})", flush=True)

    # Filter out any missing weights files
    variants = [v for v, p in WEIGHTS.items() if p is None or p.exists()]
    missing  = [v for v, p in WEIGHTS.items() if p is not None and not p.exists()]
    if missing:
        print(f"  skipping (file absent): {missing}", flush=True)
    print(f"  variants to evaluate: {variants}", flush=True)

    results = {}
    for probe_name, meta in PROBES.items():
        print(f"\n=== probe: {probe_name} ===", flush=True)
        psms = load_probe_psms(probe_name, meta)

        # One calibration + one per-PSM sweep on the PRETRAINED model;
        # reused across all variants so each variant gets the same CE
        # inputs under its own paradigm.
        print(f"  calibrating CE on pretrained predictor...", flush=True)
        pretrained = load_predictor(None)
        t0 = time.time()
        cal = calibrate_nce(pretrained, psms, nce_grid=NCE_GRID, verbose=False)
        uniform_nce = float(cal["best_nce"])
        print(f"    uniform best NCE = {uniform_nce:.0f}  "
              f"(sweep {time.time()-t0:.1f}s)", flush=True)
        t0 = time.time()
        per_psm_ce, _ = per_psm_optimal_ce(pretrained, psms, NCE_GRID)
        print(f"    per-PSM NCE  median={np.median(per_psm_ce):.0f}  "
              f"std={per_psm_ce.std():.1f}  range=[{per_psm_ce.min():.0f}, "
              f"{per_psm_ce.max():.0f}]  (sweep {time.time()-t0:.1f}s)", flush=True)
        del pretrained

        probe_results = {
            "n_psms":         len(psms),
            "uniform_nce":    uniform_nce,
            "per_psm_nce_stats": {
                "median": float(np.median(per_psm_ce)),
                "mean":   float(per_psm_ce.mean()),
                "std":    float(per_psm_ce.std()),
                "min":    float(per_psm_ce.min()),
                "max":    float(per_psm_ce.max()),
            },
            "variants": {},
        }

        for variant in variants:
            pgrm = PARADIGM[variant]
            if pgrm == "uniform":
                for p in psms:
                    p.collision_energy = uniform_nce
                    p.collision_energy_calibrated = uniform_nce
            else:
                for p, ce in zip(psms, per_psm_ce):
                    p.collision_energy = float(ce)
                    p.collision_energy_calibrated = float(ce)
            model = load_predictor(WEIGHTS[variant])
            sa = intensity_mean_sa(model, psms)
            del model
            print(f"    {variant:<16} ({pgrm:<16})  test SA = {sa:.4f}", flush=True)
            probe_results["variants"][variant] = {
                "paradigm": pgrm,
                "sa_test":  float(sa),
            }

        results[probe_name] = probe_results

    out = {
        "stamp":           STAMP,
        "started_at_utc":  time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(t_start)),
        "elapsed_seconds": round(time.time() - t_start, 1),
        "config": {
            "nce_grid":          NCE_GRID,
            "eval_max_psms":     EVAL_MAX_PSMS,
            "ft_qvalue":         FT_QVALUE,
            "min_matched_peaks": MIN_MATCHED_PEAKS,
            "seed":              RANDOM_SEED,
        },
        "weights": {v: (str(p) if p else "<pretrained>") for v, p in WEIGHTS.items() if v in variants},
        "results": results,
    }
    out_path = OUT_DIR / f"intensity_ood_eval_{STAMP}.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\n[{time.strftime('%H:%M:%S')}] DONE — wrote {out_path}", flush=True)

    # Console summary
    print("\n" + "=" * 72)
    print("OOD held-out test spectral angle, per (probe, variant)")
    print("=" * 72)
    header = f"{'probe':<22}" + " ".join(f"{v:>16}" for v in variants)
    print(header)
    for probe_name, r in results.items():
        row = f"{probe_name:<22}"
        for v in variants:
            sa = r["variants"][v]["sa_test"]
            row += f" {sa:>16.4f}"
        print(row)
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
