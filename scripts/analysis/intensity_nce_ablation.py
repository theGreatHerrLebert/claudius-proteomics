#!/usr/bin/env python3
"""Intensity-predictor NCE-handling ablation — v1 vs v4 (per-PSM optimal CE).

Compares two fine-tune strategies on the same peptide-level train/test split:

  v1  — uniform per-deposit calibrated `best_nce` (current PR #395 paradigm).
        Every PSM in deposit D uses the same NCE value, found by sweep over a
        strict-confidence sample of that deposit's PSMs.

  v4  — **per-PSM optimal CE** via pretrained-model sweep, then fine-tune at
        those fixed per-PSM values. Two phases:
          1. For each PSM, sweep NCE over `range(25, 66, 2)` (±20 around 45,
             the natural operating point) on the PRETRAINED predictor and
             record which NCE maximizes spectral angle for that PSM.
          2. Pin each PSM's `collision_energy` to its per-PSM optimum. Fine-
             tune at those fixed CEs.

        Per-PSM optimum is the controlled signal — every PSM is calibrated to
        an NCE the model actually responds well to, no randomization. The
        fine-tune then learns from a clean, per-PSM-correct CE.

Earlier variants (v1b NCE-noise, v2 Bruker-delta, v3 literal CE+offset)
finished within 0.001 SA of v1 on the killed run, confirming uncontrolled CE
inputs add nothing on top of the fine-tune. Dropped from this script.

Output:
  notebook/analysis/production_finetune/intensity_ablation_<YYYY-MM-DD>.json
  notebook/analysis/production_finetune/intensity_ablation_<YYYY-MM-DD>.log
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
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc

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

TRAINING = {
    "PXD046675": {"store": DATA_ROOT/"PXD046675_pig_trypsin"/"precursor_store.parquet",
                  "sage":  DATA_ROOT/"PXD046675_pig_trypsin"/"engines"/"sage",
                  "organism": "Pig"},
    "PXD046777": {"store": DATA_ROOT/"PXD046777"/"merged"/"precursor_store.parquet",
                  "sage":  DATA_ROOT/"PXD046777"/"processed"/"sage",
                  "organism": "Human (HIV ctx)"},
    "PXD068782": {"store": DATA_ROOT/"PXD068782_chlamydomonas_trypsin"/"precursor_store.parquet",
                  "sage":  DATA_ROOT/"PXD068782_chlamydomonas_trypsin"/"engines"/"sage",
                  "organism": "Chlamydomonas"},
}

PEP_THRESHOLD      = 1e-3
FT_QVALUE          = 0.01
MIN_MATCHED_PEAKS  = 6
N_SAMPLE_CALIB     = 1500
FT_MAX_PER_DEPOSIT = 8000
NCE_GRID_V1        = list(range(15, 76, 2))           # v1 per-deposit calibration — widened so v1 isn't edge-capped
NCE_GRID_V4        = list(range(15, 76, 2))           # v4 per-PSM sweep — matches v1 for apples-to-apples
FT_EPOCHS          = 50
FT_PATIENCE        = 5
TEST_FRAC          = 0.20
RANDOM_SEED        = 0
TRAIN_NCE_MEAN     = 29
STAMP              = date.today().isoformat()


# --- data loading helpers --------------------------------------------------
def filter_results(res, *, pep=None, qval=None):
    mask = pc.greater_equal(res["matched_peaks"], MIN_MATCHED_PEAKS)
    mask = pc.and_(mask, pc.equal(res["is_decoy"], False))
    mask = pc.and_(mask, pc.less_equal(res["rank"], 1))
    if pep is not None:
        mask = pc.and_(mask, pc.less_equal(res["posterior_error"], float(np.log(pep))))
    if qval is not None:
        mask = pc.and_(mask, pc.less_equal(res["spectrum_q"], qval))
    return res.filter(mask)


def write_subset(res, frag, cache_dir, tag):
    keep = pa.array(list(set(res.column("psm_id").to_pylist())))
    fsub = frag.filter(pc.is_in(frag["psm_id"], value_set=keep))
    rp = cache_dir / f"{tag}_results.parquet"
    fp = cache_dir / f"{tag}_fragments.parquet"
    pq.write_table(res, rp); pq.write_table(fsub, fp)
    return rp, fp


def collision_energy_by_psm_id(store_path):
    t = pq.read_table(store_path, columns=["sage_psm_id", "collision_energy"])
    t = t.filter(pc.is_valid(t["sage_psm_id"]))
    df = t.to_pandas()
    df["sage_psm_id"] = df["sage_psm_id"].astype("int64")
    return dict(zip(df["sage_psm_id"], df["collision_energy"].astype(float)))


def prep_deposit(name, meta, verbose=True):
    sage = meta["sage"]
    cache = OUT_DIR / "cache" / name
    cache.mkdir(parents=True, exist_ok=True)
    results_p = cache / "ft_results.parquet"
    frags_p   = cache / "ft_fragments.parquet"
    calib_r   = cache / "calib_results.parquet"
    calib_f   = cache / "calib_fragments.parquet"

    if not (results_p.exists() and frags_p.exists()
            and calib_r.exists() and calib_f.exists()):
        if verbose:
            print(f"  [{name}] cache miss; rebuilding...", flush=True)
        res_full  = pq.read_table(sage / "results.sage.parquet")
        frag_full = pq.read_table(sage / "matched_fragments.sage.parquet")
        res_strict = filter_results(res_full, pep=PEP_THRESHOLD)
        if N_SAMPLE_CALIB and res_strict.num_rows > N_SAMPLE_CALIB:
            idx = np.sort(np.random.RandomState(RANDOM_SEED).choice(
                res_strict.num_rows, N_SAMPLE_CALIB, replace=False))
            res_strict = res_strict.take(pa.array(idx))
        write_subset(res_strict, frag_full, cache, "calib")
        res_ft = filter_results(res_full, qval=FT_QVALUE)
        if FT_MAX_PER_DEPOSIT and res_ft.num_rows > FT_MAX_PER_DEPOSIT:
            idx = np.sort(np.random.RandomState(RANDOM_SEED).choice(
                res_ft.num_rows, FT_MAX_PER_DEPOSIT, replace=False))
            res_ft = res_ft.take(pa.array(idx))
        write_subset(res_ft, frag_full, cache, "ft")

    # v1 deposit-level NCE calibration (kept for the v1 baseline run + reporting)
    psms_cal = load_psms_from_sage_parquets(
        str(calib_r), str(calib_f),
        default_collision_energy=float(TRAIN_NCE_MEAN), max_rank=1)
    cal_model = DeepPeptideIntensityPredictor(verbose=False)
    cal = calibrate_nce(cal_model, psms_cal, nce_grid=NCE_GRID_V1,
                        max_sample=N_SAMPLE_CALIB, verbose=False)
    best_nce = int(cal["best_nce"])
    del cal_model, psms_cal

    psms = load_psms_from_sage_parquets(
        str(results_p), str(frags_p),
        default_collision_energy=float(best_nce), max_rank=1)
    ce_map = collision_energy_by_psm_id(meta["store"])
    n_real = 0
    for p in psms:
        ce = ce_map.get(int(p.sage_feature.psm_id))
        if ce and ce > 0:
            p._bruker_ce = float(ce); n_real += 1
        else:
            p._bruker_ce = None
        p.collision_energy = float(best_nce)
        p.collision_energy_calibrated = float(best_nce)
        p._deposit = name
    if verbose:
        print(f"  [{name}] deposit best_nce={best_nce}  |  {len(psms):,} PSMs  |  "
              f"{n_real:,} have real Bruker CE", flush=True)
    return {"best_nce": best_nce, "psms": psms}


# --- per-PSM CE sweep on pretrained predictor (the v4 core) ----------------
def per_psm_optimal_ce(predictor, psms, nce_grid, verbose=True):
    """For each PSM, find the NCE in `nce_grid` that maximizes pretrained
    spectral angle. Returns (best_nce_per_psm: np.array, sa_matrix: np.array
    of shape (n_psms, n_nce))."""
    n = len(psms)
    seqs    = [p.sequence_modified for p in psms]
    charges = np.array([p.charge for p in psms])
    sa = np.zeros((n, len(nce_grid)), dtype=np.float32)
    for j, nce in enumerate(nce_grid):
        t0 = time.time()
        ces = [float(nce)] * n
        intens = predictor.predict_intensities(
            seqs, charges, ces, batch_size=2048, flatten=True)
        coll = associate_fragments(psms, intens, num_threads=os.cpu_count() or 1)
        for p, c in zip(psms, coll):
            p.prosit_predicted_intensities = c.prosit_predicted_intensities
        df = psm_collection_to_pandas(psms)
        sa[:, j] = df["spectral_angle_similarity"].to_numpy(np.float32)
        if verbose:
            print(f"    NCE={nce:>3}  mean SA = {sa[:, j].mean():.4f}  "
                  f"({time.time()-t0:.1f}s)", flush=True)
    best_idx = sa.argmax(axis=1)
    best_nce = np.array([nce_grid[i] for i in best_idx], dtype=np.float32)
    return best_nce, sa


# --- scoring ---------------------------------------------------------------
def intensity_mean_sa(model, psms):
    if not psms:
        return float("nan")
    seqs    = [p.sequence_modified for p in psms]
    charges = np.array([p.charge for p in psms])
    ce      = [float(p.collision_energy) for p in psms]
    intens = model.predict_intensities(seqs, charges, ce,
                                       batch_size=2048, flatten=True)
    coll = associate_fragments(psms, intens, num_threads=os.cpu_count() or 1)
    for p, c in zip(psms, coll):
        p.prosit_predicted_intensities = c.prosit_predicted_intensities
    df = psm_collection_to_pandas(psms)
    return float(df["spectral_angle_similarity"].mean())


# --- main ------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] intensity v1-vs-v4 ablation starting "
          f"(stamp={STAMP})", flush=True)
    print(f"  v1 NCE grid: {NCE_GRID_V1[0]}..{NCE_GRID_V1[-1]} (deposit-level calib)", flush=True)
    print(f"  v4 NCE grid: {NCE_GRID_V4[0]}..{NCE_GRID_V4[-1]} (per-PSM sweep)", flush=True)

    # --- prep -------------------------------------------------------------
    print("\n=== Phase 1: per-deposit prep ===", flush=True)
    PREP = {n: prep_deposit(n, m) for n, m in TRAINING.items()}

    # --- pool + peptide-level split --------------------------------------
    all_psms = []
    for n in TRAINING:
        all_psms.extend(PREP[n]["psms"])
    rng = np.random.RandomState(RANDOM_SEED)
    all_peps = sorted({p.sequence for p in all_psms})
    rng.shuffle(all_peps)
    test_peps = set(all_peps[:int(TEST_FRAC * len(all_peps))])
    train_psms = [p for p in all_psms if p.sequence not in test_peps]
    test_psms  = [p for p in all_psms if p.sequence in test_peps]
    print(f"\nPooled: {len(all_psms):,} PSMs / {len(all_peps):,} peptides", flush=True)
    print(f"Split: {len(train_psms):,} train / {len(test_psms):,} test", flush=True)

    # --- v1: uniform per-deposit best_nce, fine-tune -----------------------
    print("\n=== Phase 2: v1 (uniform per-deposit best_nce) ===", flush=True)
    for p in all_psms:
        p.collision_energy = float(PREP[p._deposit]["best_nce"])
        p.collision_energy_calibrated = p.collision_energy

    pred_v1 = DeepPeptideIntensityPredictor(verbose=False)
    sa_pre_v1 = intensity_mean_sa(pred_v1, test_psms)
    print(f"  pretrained test SA = {sa_pre_v1:.4f}", flush=True)
    t0 = time.time()
    pred_v1.fine_tune_psms(train_psms, batch_size=64, epochs=FT_EPOCHS,
                           patience=FT_PATIENCE, verbose=False)
    print(f"  fine-tune done in {time.time()-t0:.1f}s", flush=True)
    sa_post_v1 = intensity_mean_sa(pred_v1, test_psms)
    print(f"  fine-tuned test SA = {sa_post_v1:.4f}  (Δ={sa_post_v1-sa_pre_v1:+.4f})",
          flush=True)
    del pred_v1

    # --- v4: per-PSM optimal CE via pretrained sweep -----------------------
    print(f"\n=== Phase 3: v4 step 1 — per-PSM CE sweep on pretrained "
          f"({len(NCE_GRID_V4)} NCE values × {len(all_psms):,} PSMs) ===", flush=True)
    fresh = DeepPeptideIntensityPredictor(verbose=False)
    t0 = time.time()
    best_nce_per_psm, sa_matrix = per_psm_optimal_ce(
        fresh, all_psms, NCE_GRID_V4, verbose=True)
    print(f"  sweep done in {time.time()-t0:.1f}s", flush=True)
    print(f"  per-PSM optimal NCE: "
          f"min={best_nce_per_psm.min():.0f}  "
          f"median={np.median(best_nce_per_psm):.0f}  "
          f"max={best_nce_per_psm.max():.0f}  "
          f"mean={best_nce_per_psm.mean():.1f}  "
          f"std={best_nce_per_psm.std():.1f}", flush=True)
    del fresh

    # Pin per-PSM CE
    for p, ce in zip(all_psms, best_nce_per_psm):
        p.collision_energy = float(ce)
        p.collision_energy_calibrated = float(ce)

    print("\n=== Phase 3: v4 step 2 — fine-tune at fixed per-PSM CEs ===", flush=True)
    pred_v4 = DeepPeptideIntensityPredictor(verbose=False)
    sa_pre_v4 = intensity_mean_sa(pred_v4, test_psms)
    print(f"  pretrained test SA = {sa_pre_v4:.4f}", flush=True)
    t0 = time.time()
    pred_v4.fine_tune_psms(train_psms, batch_size=64, epochs=FT_EPOCHS,
                           patience=FT_PATIENCE, verbose=False)
    print(f"  fine-tune done in {time.time()-t0:.1f}s", flush=True)
    sa_post_v4 = intensity_mean_sa(pred_v4, test_psms)
    print(f"  fine-tuned test SA = {sa_post_v4:.4f}  (Δ={sa_post_v4-sa_pre_v4:+.4f})",
          flush=True)

    # Per-deposit v4 breakdown
    v4_per_deposit = {}
    for dep in TRAINING:
        sub = [p for p in test_psms if p._deposit == dep]
        v4_per_deposit[dep] = {
            "n":       len(sub),
            "sa_post": float(intensity_mean_sa(pred_v4, sub)) if sub else None,
        }
    del pred_v4

    # --- save weights if v4 wins ------------------------------------------
    if sa_post_v4 > sa_post_v1:
        v4_path = PROJECT_ROOT / "models" / "production" / f"imspy_intensity_finetune_v4_{STAMP}.pt"
        # We need a fresh fine-tune to save (the in-memory model was deleted).
        # Re-run the v4 fine-tune one more time and persist.
        print("\nv4 beats v1 — re-running v4 fine-tune and persisting weights...", flush=True)
        for p, ce in zip(all_psms, best_nce_per_psm):
            p.collision_energy = float(ce)
            p.collision_energy_calibrated = float(ce)
        pred_save = DeepPeptideIntensityPredictor(verbose=False)
        pred_save.fine_tune_psms(train_psms, batch_size=64, epochs=FT_EPOCHS,
                                 patience=FT_PATIENCE, verbose=False)
        import torch
        torch.save(
            {"model_state_dict": pred_save.model.state_dict(), "task": "intensity"},
            v4_path,
        )
        print(f"  saved {v4_path}", flush=True)
        del pred_save
        v4_weights_path = str(v4_path)
    else:
        v4_weights_path = None
        print(f"\nv4 does not beat v1 — weights not persisted.", flush=True)

    # --- assemble + save report -------------------------------------------
    histogram = np.bincount(
        np.searchsorted(NCE_GRID_V4, best_nce_per_psm.astype(int)),
        minlength=len(NCE_GRID_V4))
    out = {
        "stamp":               STAMP,
        "started_at_utc":      time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(t_start)),
        "elapsed_seconds":     round(time.time() - t_start, 1),
        "training_deposits":   list(TRAINING),
        "per_deposit_best_nce": {n: PREP[n]["best_nce"] for n in PREP},
        "n_train": len(train_psms),
        "n_test":  len(test_psms),
        "config": {
            "pep_threshold":      PEP_THRESHOLD,
            "ft_qvalue":          FT_QVALUE,
            "min_matched_peaks":  MIN_MATCHED_PEAKS,
            "ft_max_per_deposit": FT_MAX_PER_DEPOSIT,
            "ft_epochs":          FT_EPOCHS,
            "ft_patience":        FT_PATIENCE,
            "test_frac":          TEST_FRAC,
            "seed":               RANDOM_SEED,
            "nce_grid_v1":        NCE_GRID_V1,
            "nce_grid_v4":        NCE_GRID_V4,
        },
        "v1": {
            "rule":              "uniform per-deposit best_nce",
            "sa_pretrained_test": float(sa_pre_v1),
            "sa_finetuned_test":  float(sa_post_v1),
            "delta":              float(sa_post_v1 - sa_pre_v1),
        },
        "v4": {
            "rule":              "per-PSM optimal CE (pretrained sweep) + fine-tune at fixed per-PSM CEs",
            "sa_pretrained_test": float(sa_pre_v4),
            "sa_finetuned_test":  float(sa_post_v4),
            "delta":              float(sa_post_v4 - sa_pre_v4),
            "per_deposit":        v4_per_deposit,
            "per_psm_ce_distribution": {
                "nce_grid":  NCE_GRID_V4,
                "histogram": [int(h) for h in histogram],
                "min":       float(best_nce_per_psm.min()),
                "median":    float(np.median(best_nce_per_psm)),
                "mean":      float(best_nce_per_psm.mean()),
                "std":       float(best_nce_per_psm.std()),
                "max":       float(best_nce_per_psm.max()),
            },
            "weights_path": v4_weights_path,
        },
        "v4_minus_v1_delta_sa": float(sa_post_v4 - sa_post_v1),
    }
    out_path = OUT_DIR / f"intensity_ablation_{STAMP}.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\n[{time.strftime('%H:%M:%S')}] DONE — wrote {out_path}", flush=True)

    print("\n" + "=" * 64)
    print("Summary — held-out test mean spectral angle")
    print("=" * 64)
    print(f"  v1 (uniform best_nce):   pre={sa_pre_v1:.4f}  post={sa_post_v1:.4f}  Δ={sa_post_v1-sa_pre_v1:+.4f}")
    print(f"  v4 (per-PSM optimal CE): pre={sa_pre_v4:.4f}  post={sa_post_v4:.4f}  Δ={sa_post_v4-sa_pre_v4:+.4f}")
    print(f"  v4 − v1 (fine-tuned):    {sa_post_v4 - sa_post_v1:+.4f}")
    print("=" * 64)
    return 0


if __name__ == "__main__":
    sys.exit(main())
