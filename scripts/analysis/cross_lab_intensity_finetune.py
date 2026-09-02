"""Cross-lab fine-tune stability test for the intensity predictor.

Replicates the calibrate-NCE + peptide-level fine-tune of
notebook/intensity_exploration.ipynb (sections 5 and 6), but for every
dataset whose Sage outputs are present, and answers two questions:

  1. Stability  -- does an independent fine-tune on each dataset give a
     comparable held-out spectral-angle gain? (diagonal of the matrix)
  2. Translation -- does a model fine-tuned on dataset A still help on
     dataset B's held-out peptides? (off-diagonal of the matrix)

Sage-only: needs results.sage.parquet + matched_fragments.sage.parquet.
The per-precursor Bruker-CE baseline from the notebook is intentionally
dropped (it needs merged/precursor_store.parquet, not yet produced).
Every dataset is instead evaluated at its own calibrated NCE.

Output: notebook/analysis/intensity_exploration/cross_lab/
"""
from __future__ import annotations

import json
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from sagepy_rescore.sage_loader import load_psms_from_sage_parquets
from sagepy_rescore.predictors import apply_imspy_patches
from sagepy.utility import psm_collection_to_pandas
from imspy_predictors.intensity.predictors import DeepPeptideIntensityPredictor, calibrate_nce
from imspy_predictors.lazy_imports import get_sagepy_fragment_utils

# --- config (mirrors the notebook) ------------------------------------------
DATA_ROOT    = Path(os.environ.get("CLAUDIUS_DATA_ROOT", "/scratch/claudius-proteomics"))
PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR      = PROJECT_ROOT / "notebook" / "analysis" / "intensity_exploration" / "cross_lab"

DATASETS = ["PXD046777", "PXD050342", "PXD051564", "PXD066365", "PXD073076"]
RUN1     = "PXD046777"          # the "first run" the notebook fine-tuned on

PEP_THRESHOLD     = 1e-3        # strict calibration sample
MIN_MATCHED_PEAKS = 6
N_SAMPLE          = 3000        # calibration sweep sample
NCE_GRID          = list(range(12, 51))
FT_QVALUE         = 0.01        # fine-tune set: spectrum_q <= this
FT_MAX            = 20000
RANDOM_SEED       = 0
TEST_FRAC         = 0.20

apply_imspy_patches()
associate_fragments, _ = get_sagepy_fragment_utils()
SIM_COL = "spectral_angle_similarity"


# --- data loading -----------------------------------------------------------
def _sage_paths(ds: str):
    sage = DATA_ROOT / ds / "processed" / "sage"
    return sage / "results.sage.parquet", sage / "matched_fragments.sage.parquet"


def _filter_results(res: pa.Table, *, pep=None, qval=None) -> pa.Table:
    """Target rank-1 PSMs with enough matched peaks, by PEP or q-value."""
    mask = pc.greater_equal(res["matched_peaks"], MIN_MATCHED_PEAKS)
    mask = pc.and_(mask, pc.equal(res["is_decoy"], False))
    mask = pc.and_(mask, pc.less_equal(res["rank"], 1))
    if pep is not None:
        mask = pc.and_(mask, pc.less_equal(res["posterior_error"], float(np.log(pep))))
    if qval is not None:
        mask = pc.and_(mask, pc.less_equal(res["spectrum_q"], qval))
    return res.filter(mask)


def _write_subset(res: pa.Table, frag: pa.Table, cache: Path, tag: str):
    keep = pa.array(list(set(res.column("psm_id").to_pylist())))
    fsub = frag.filter(pc.is_in(frag["psm_id"], value_set=keep))
    rp, fp = cache / f"{tag}_results.parquet", cache / f"{tag}_fragments.parquet"
    pq.write_table(res, rp)
    pq.write_table(fsub, fp)
    return rp, fp


def calibrate_dataset(ds: str, cache: Path) -> int:
    """NCE sweep on the strict PEP<=1e-3 sample -> best aligned NCE."""
    rpath, fpath = _sage_paths(ds)
    res = _filter_results(pq.read_table(rpath), pep=PEP_THRESHOLD)
    n_pass = res.num_rows
    if N_SAMPLE and n_pass > N_SAMPLE:
        idx = np.sort(np.random.RandomState(RANDOM_SEED).choice(n_pass, N_SAMPLE, replace=False))
        res = res.take(pa.array(idx))
    rp, fp = _write_subset(res, pq.read_table(fpath), cache, "calib")
    psms = load_psms_from_sage_parquets(str(rp), str(fp), default_collision_energy=30.0, max_rank=1)
    model = DeepPeptideIntensityPredictor(verbose=False)
    cal = calibrate_nce(model, psms, nce_grid=NCE_GRID, max_sample=N_SAMPLE, verbose=False)
    best = int(cal["best_nce"])
    print(f"  [{ds}] calib sample {len(psms):,} PSMs (of {n_pass:,} PEP<=1e-3)  "
          f"-> best NCE = {best}", flush=True)
    return best


def load_finetune_split(ds: str, best_nce: int, cache: Path):
    """q<=FT_QVALUE set, CE pinned to best_nce, peptide-level 80/20 split."""
    rpath, fpath = _sage_paths(ds)
    res = _filter_results(pq.read_table(rpath), qval=FT_QVALUE)
    if FT_MAX and res.num_rows > FT_MAX:
        idx = np.sort(np.random.RandomState(RANDOM_SEED).choice(res.num_rows, FT_MAX, replace=False))
        res = res.take(pa.array(idx))
    rp, fp = _write_subset(res, pq.read_table(fpath), cache, "ft")
    psms = load_psms_from_sage_parquets(str(rp), str(fp),
                                        default_collision_energy=float(best_nce), max_rank=1)
    for p in psms:
        p.collision_energy = float(best_nce)
        p.collision_energy_calibrated = float(best_nce)
    rng = np.random.RandomState(RANDOM_SEED)
    peps = sorted({p.sequence for p in psms})
    rng.shuffle(peps)
    test_peps = set(peps[:int(TEST_FRAC * len(peps))])
    train = [p for p in psms if p.sequence not in test_peps]
    test  = [p for p in psms if p.sequence in test_peps]
    print(f"  [{ds}] fine-tune set {len(psms):,} PSMs / {len(peps):,} peptides  "
          f"-> {len(train):,} train / {len(test):,} test", flush=True)
    return {"train": train, "test": test, "n_psms": len(psms),
            "n_peptides": len(peps), "train_peptides": set(peps) - test_peps}


# --- scoring -----------------------------------------------------------------
def mean_sa(model: DeepPeptideIntensityPredictor, psms: list, nce: float) -> float:
    """Predict at `nce`, associate to observed fragments, mean spectral angle."""
    if not psms:
        return float("nan")
    seqs    = [p.sequence_modified for p in psms]
    charges = np.array([p.charge for p in psms])
    ce      = [float(nce)] * len(psms)
    intens  = model.predict_intensities(seqs, charges, ce, batch_size=2048, flatten=True)
    coll    = associate_fragments(psms, intens, num_threads=os.cpu_count() or 1)
    for p, c in zip(psms, coll):
        p.prosit_predicted_intensities = c.prosit_predicted_intensities
    df = psm_collection_to_pandas(psms)
    return float(df[SIM_COL].mean())


# --- main --------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    present = [d for d in DATASETS if _sage_paths(d)[0].exists() and _sage_paths(d)[1].exists()]
    missing = [d for d in DATASETS if d not in present]
    if missing:
        print(f"skipping (Sage parquets absent): {missing}", flush=True)
    print(f"datasets: {present}\n", flush=True)

    # Phase 1: calibrate NCE + build fine-tune splits ------------------------
    print("=== phase 1: calibrate + split ===", flush=True)
    info = {}
    for ds in present:
        cache = OUT_DIR / f"{ds}_cache"
        cache.mkdir(parents=True, exist_ok=True)
        best = calibrate_dataset(ds, cache)
        split = load_finetune_split(ds, best, cache)
        info[ds] = {"best_nce": best, **split}

    # Phase 2: pretrained baseline (no fine-tune) ----------------------------
    print("\n=== phase 2: pretrained baseline ===", flush=True)
    pretrained = DeepPeptideIntensityPredictor(verbose=False)
    sa_pretrained = {}
    for ds in present:
        sa = mean_sa(pretrained, info[ds]["test"], info[ds]["best_nce"])
        sa_pretrained[ds] = sa
        print(f"  [{ds}] pretrained held-out SA = {sa:.4f}", flush=True)

    # Phase 3: fine-tune on each dataset, evaluate the full transfer matrix --
    print("\n=== phase 3: fine-tune + transfer matrix ===", flush=True)
    matrix = {}   # matrix[train_ds][test_ds] = held-out mean SA
    history = {}
    for train_ds in present:
        print(f"  fine-tuning on {train_ds} ...", flush=True)
        model = DeepPeptideIntensityPredictor(verbose=False)
        model.fine_tune_psms(info[train_ds]["train"], batch_size=64, epochs=50,
                             patience=5, verbose=False)
        history[train_ds] = getattr(model, "_finetune_history",
                                    {"epochs": [], "train_loss": [], "val_loss": []})
        row = {}
        for test_ds in present:
            sa = mean_sa(model, info[test_ds]["test"], info[test_ds]["best_nce"])
            row[test_ds] = sa
            tag = "self" if test_ds == train_ds else "->"
            print(f"    {train_ds} -> {test_ds:<10} held-out SA = {sa:.4f}  ({tag})", flush=True)
        matrix[train_ds] = row

    # --- assemble report ----------------------------------------------------
    mat_df = pd.DataFrame(matrix).T.loc[present, present]   # rows=train, cols=test
    mat_df.to_csv(OUT_DIR / "transfer_matrix.csv")

    summary = {
        "datasets": present,
        "config": {"pep_threshold": PEP_THRESHOLD, "ft_qvalue": FT_QVALUE,
                   "ft_max": FT_MAX, "min_matched_peaks": MIN_MATCHED_PEAKS,
                   "test_frac": TEST_FRAC, "seed": RANDOM_SEED},
        "best_nce": {d: info[d]["best_nce"] for d in present},
        "n_finetune_psms": {d: info[d]["n_psms"] for d in present},
        "n_finetune_peptides": {d: info[d]["n_peptides"] for d in present},
        "sa_pretrained": sa_pretrained,
        "transfer_matrix": {tr: matrix[tr] for tr in present},
        "stability": {
            d: {"pretrained": sa_pretrained[d],
                "self_finetune": matrix[d][d],
                "gain": matrix[d][d] - sa_pretrained[d]}
            for d in present
        },
        "translation_from_run1": {
            d: {"pretrained": sa_pretrained[d],
                "run1_finetune": matrix[RUN1][d],
                "self_finetune": matrix[d][d],
                "translated_gain": matrix[RUN1][d] - sa_pretrained[d],
                "gain_retained_vs_self":
                    (matrix[RUN1][d] - sa_pretrained[d])
                    / max(matrix[d][d] - sa_pretrained[d], 1e-9)}
            for d in present if d != RUN1
        },
    }
    (OUT_DIR / "cross_lab_summary.json").write_text(json.dumps(summary, indent=2))

    # --- plot ---------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
        im = ax[0].imshow(mat_df.values, cmap="viridis", vmin=0.6, vmax=1.0)
        ax[0].set_xticks(range(len(present)), present, rotation=45, ha="right")
        ax[0].set_yticks(range(len(present)), present)
        ax[0].set_xlabel("evaluated on (held-out test)")
        ax[0].set_ylabel("fine-tuned on")
        ax[0].set_title("Transfer matrix: held-out spectral angle")
        for i in range(len(present)):
            for j in range(len(present)):
                ax[0].text(j, i, f"{mat_df.values[i, j]:.3f}", ha="center", va="center",
                           color="white" if mat_df.values[i, j] < 0.85 else "black", fontsize=9)
        fig.colorbar(im, ax=ax[0], fraction=0.046)

        x = np.arange(len(present))
        ax[1].bar(x - 0.27, [sa_pretrained[d] for d in present], 0.27, label="pretrained")
        ax[1].bar(x,        [matrix[RUN1][d] for d in present], 0.27, label=f"fine-tuned on {RUN1}")
        ax[1].bar(x + 0.27, [matrix[d][d] for d in present], 0.27, label="fine-tuned on self")
        ax[1].set_xticks(x, present, rotation=45, ha="right")
        ax[1].set_ylabel("held-out spectral angle")
        ax[1].set_ylim(0.6, 1.0)
        ax[1].set_title("Translation vs. self fine-tune")
        ax[1].legend()
        fig.tight_layout()
        fig.savefig(OUT_DIR / "cross_lab_finetune.png", dpi=130)
        print(f"\nsaved plot -> {OUT_DIR / 'cross_lab_finetune.png'}", flush=True)
    except Exception as e:
        print(f"plot skipped: {e}", flush=True)

    # --- console summary ----------------------------------------------------
    print("\n" + "=" * 64)
    print("TRANSFER MATRIX (held-out spectral angle; rows=fine-tuned on)")
    print(mat_df.round(4).to_string())
    print("\nSTABILITY (independent fine-tune on each dataset):")
    for d in present:
        s = summary["stability"][d]
        print(f"  {d:<11} {s['pretrained']:.3f} -> {s['self_finetune']:.3f}  "
              f"(gain {s['gain']:+.3f})")
    print(f"\nTRANSLATION (model fine-tuned on {RUN1}, applied zero-shot):")
    for d, t in summary["translation_from_run1"].items():
        print(f"  {d:<11} {t['pretrained']:.3f} -> {t['run1_finetune']:.3f}  "
              f"(gain {t['translated_gain']:+.3f}; "
              f"{t['gain_retained_vs_self']*100:.0f}% of self-finetune gain)")
    print("=" * 64)
    print(f"\nwrote: {OUT_DIR / 'cross_lab_summary.json'}")
    print(f"       {OUT_DIR / 'transfer_matrix.csv'}")


if __name__ == "__main__":
    main()
