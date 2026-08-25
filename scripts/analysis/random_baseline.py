"""Chance-floor baseline for the independent intensity metric.

cosine between two non-negative vectors is substantially positive even when
they are unrelated -- so an independent cosine of 0.64 (pretrained) or 0.80
(fine-tuned) is only interpretable against a null. This builds two nulls on the
PXD046777 held-out test set:

  null A (random spectrum, per the proposal): keep every blob peak's intensity
    but scatter the peaks to random m/z positions on a 0.01-Da grid across the
    spectrum's m/z range, re-integrate the theoretical b/y ions, score the
    prediction against that. Tests the whole extract+score pipeline vs a
    structureless spectrum of identical peak content.

  null B (ion permutation): take the real integrated b/y observed vector and
    shuffle its intensities across ion positions. Same intensity distribution,
    destroyed position<->intensity relationship.

Reported against the real independent cosine, for both the pretrained and the
fine-tuned predictor.

Output: notebook/analysis/intensity_exploration/cross_lab/random_baseline.json + .png
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pyarrow.parquet as pq

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "analysis"))

from sagepy_rescore.sage_loader import load_psms_from_sage_parquets
from sagepy_rescore.predictors import apply_imspy_patches
from imspy_predictors.intensity.predictors import DeepPeptideIntensityPredictor
from fragment_matching import FragmentMatcher, MatchConfig
from extract_fragment_peaks import read_blob, resolve_blob_path
from blob_intensity import integrate_peaks_for_psm

apply_imspy_patches()

DATASET   = "PXD046777"
BEST_NCE  = 47
SEED      = 0
TEST_FRAC = 0.20
TOL_PPM   = 20.0
GRID_DA   = 0.01           # m/z grid for the random-spectrum permutation
N_DRAWS   = 20             # random draws per PSM, per null

DATA_ROOT = Path("/scratch/claudius-proteomics") / DATASET
CACHE     = PROJECT_ROOT / "notebook/analysis/intensity_exploration/cross_lab" / f"{DATASET}_cache"
STORE     = DATA_ROOT / "merged" / "precursor_store.parquet"
BLOB_DIR  = DATA_ROOT / "extracted"
OUT_JSON  = PROJECT_ROOT / "notebook/analysis/intensity_exploration/cross_lab" / "random_baseline.json"
OUT_PNG   = PROJECT_ROOT / "notebook/analysis/intensity_exploration/cross_lab" / "random_baseline.png"

ION_AXIS = {0: 1, 1: 0}    # grid[ordinal-1, ion_axis, charge-1]; y=0, b=1


def cosine(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size < 2 or a.sum() <= 0 or b.sum() <= 0:
        return np.nan
    return float(np.clip(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1))


def grid_value(grid, it, chg, ordn):
    if ordn < 1 or ordn > 29 or chg < 1 or chg > 3:
        return None
    return float(grid[ordn - 1, ION_AXIS[it], chg - 1])


def main():
    print(f"=== random baseline: {DATASET} ===\n", flush=True)

    psms = load_psms_from_sage_parquets(str(CACHE / "ft_results.parquet"),
                                        str(CACHE / "ft_fragments.parquet"),
                                        default_collision_energy=float(BEST_NCE), max_rank=1)
    for p in psms:
        p.collision_energy = p.collision_energy_calibrated = float(BEST_NCE)
    rng = np.random.RandomState(SEED)
    peps = sorted({p.sequence for p in psms})
    rng.shuffle(peps)
    test_peps = set(peps[:int(TEST_FRAC * len(peps))])
    train = [p for p in psms if p.sequence not in test_peps]
    test  = [p for p in psms if p.sequence in test_peps]
    print(f"held-out test: {len(test):,} PSMs", flush=True)

    t = pq.read_table(STORE, columns=["sage_psm_id", "raw_file", "blob_offset",
                                      "blob_size", "mono_mz"])
    df = t.filter(t["sage_psm_id"].is_valid()).to_pandas()
    df["sage_psm_id"] = df["sage_psm_id"].astype("int64")
    bidx = {r.sage_psm_id: (r.raw_file, int(r.blob_offset), int(r.blob_size), float(r.mono_mz))
            for r in df.itertuples()}
    matcher = FragmentMatcher(MatchConfig(mz_tolerance_ppm=TOL_PPM,
                                          ion_types=["b", "y"], max_fragment_charge=2))

    # per-PSM: real integrated observed vector, theoretical-ion keys, blob cloud
    print("integrating real spectra ...", flush=True)
    records = []   # (psm, keys, real_obs[np], blob_mz, blob_int)
    for p in test:
        loc = bidx.get(int(p.sage_feature.psm_id))
        if loc is None:
            continue
        raw_file, off, size, prec_mz = loc
        blob = read_blob(resolve_blob_path(BLOB_DIR, raw_file), off, size)
        if blob is None or len(blob["frag_mz"]) == 0:
            continue
        theo = matcher.generate_theoretical_fragments(p.sequence_modified, p.charge)
        peaks = integrate_peaks_for_psm(blob["frag_mz"], blob["frag_intensity"], theo,
                                        TOL_PPM, precursor_mz=prec_mz)
        keys, real = [], []
        for it, num, chg, inten in zip(peaks["ion_type"], peaks["ion_number"],
                                       peaks["ion_charge"], peaks["intensity"]):
            keys.append((0 if it == "b" else 1, int(chg), int(num)))
            real.append(float(inten))
        records.append((p, theo, keys, np.array(real),
                        np.asarray(blob["frag_mz"], float),
                        np.asarray(blob["frag_intensity"], float)))
    print(f"  {len(records):,} PSMs ready", flush=True)

    def predicted_vectors(model):
        seqs = [r[0].sequence_modified for r in records]
        chs  = [r[0].charge for r in records]
        grids = model.predict_intensities(seqs, chs, [float(BEST_NCE)] * len(records),
                                          batch_size=1024, flatten=False)
        out = []
        for (p, theo, keys, real, bmz, bint), g in zip(records, grids):
            g = np.asarray(g)
            out.append(np.array([grid_value(g, it, chg, ordn) or 0.0
                                 for (it, chg, ordn) in keys]))
        return out

    def evaluate(pred_vecs, tag):
        rs = np.random.RandomState(SEED)
        real_cos, nullA, nullB = [], [], []
        for (p, theo, keys, real, bmz, bint), pred in zip(records, pred_vecs):
            real_cos.append(cosine(pred, real))
            # null A: random-spectrum re-integration
            a = []
            lo, hi = bmz.min(), bmz.max()
            n_grid = max(int((hi - lo) / GRID_DA), 1)
            for _ in range(N_DRAWS):
                rand_mz = lo + GRID_DA * rs.randint(0, n_grid + 1, size=bint.size)
                pk = integrate_peaks_for_psm(rand_mz, bint, theo, TOL_PPM)
                robs = np.array([dict(zip(
                    [(0 if t2 == "b" else 1, int(c2), int(n2)) for t2, n2, c2 in
                     zip(pk["ion_type"], pk["ion_number"], pk["ion_charge"])],
                    pk["intensity"])).get(k, 0.0) for k in keys])
                a.append(cosine(pred, robs))
            nullA.append(np.nanmean(a))
            # null B: ion-permutation of the real observed vector
            b = []
            for _ in range(N_DRAWS):
                b.append(cosine(pred, rs.permutation(real)))
            nullB.append(np.nanmean(b))
        real_cos = np.array(real_cos); nullA = np.array(nullA); nullB = np.array(nullB)
        res = {
            "real_cosine_mean":  float(np.nanmean(real_cos)),
            "null_A_random_spectrum_mean": float(np.nanmean(nullA)),
            "null_B_ion_permutation_mean": float(np.nanmean(nullB)),
            "real_minus_nullA": float(np.nanmean(real_cos) - np.nanmean(nullA)),
            "frac_psms_real_above_nullA": float(np.nanmean(real_cos > nullA)),
        }
        print(f"  [{tag}] real {res['real_cosine_mean']:.4f} | "
              f"null-A(random spectrum) {res['null_A_random_spectrum_mean']:.4f} | "
              f"null-B(ion perm) {res['null_B_ion_permutation_mean']:.4f}", flush=True)
        return res, real_cos, nullA, nullB

    print("\nscoring pretrained ...", flush=True)
    pre = DeepPeptideIntensityPredictor(verbose=False)
    pre_vecs = predicted_vectors(pre)
    pre_res, pre_real, pre_nA, pre_nB = evaluate(pre_vecs, "pretrained")

    print("\nfine-tuning + scoring ...", flush=True)
    ft = DeepPeptideIntensityPredictor(verbose=False)
    ft.fine_tune_psms(train, batch_size=64, epochs=50, patience=5, verbose=False)
    ft_vecs = predicted_vectors(ft)
    ft_res, ft_real, ft_nA, ft_nB = evaluate(ft_vecs, "fine-tuned")

    result = {
        "dataset": DATASET, "n_test_psms": len(records),
        "grid_da": GRID_DA, "n_draws": N_DRAWS,
        "metric": "cosine on raw linear intensities, independent full b/y support",
        "pretrained": pre_res, "finetuned": ft_res,
    }
    OUT_JSON.write_text(json.dumps(result, indent=2))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for a, real, nA, nB, title in [
                (ax[0], pre_real, pre_nA, pre_nB, "pretrained"),
                (ax[1], ft_real, ft_nA, ft_nB, "fine-tuned")]:
            a.hist(nA[~np.isnan(nA)], bins=50, alpha=0.6, color="#9ca3af",
                   label=f"null A random spectrum ({np.nanmean(nA):.3f})")
            a.hist(nB[~np.isnan(nB)], bins=50, alpha=0.6, color="#f59e0b",
                   label=f"null B ion-perm ({np.nanmean(nB):.3f})")
            a.hist(real[~np.isnan(real)], bins=50, alpha=0.6, color="#2563eb",
                   label=f"real ({np.nanmean(real):.3f})")
            a.set_xlabel("independent cosine"); a.set_ylabel("PSMs")
            a.set_title(title); a.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(OUT_PNG, dpi=130)
        print(f"\nsaved plot -> {OUT_PNG}", flush=True)
    except Exception as e:
        print(f"plot skipped: {e}", flush=True)

    print("\n" + "=" * 70)
    print("INDEPENDENT COSINE vs CHANCE FLOOR  (PXD046777 held-out test)")
    for tag, r in [("pretrained", pre_res), ("fine-tuned", ft_res)]:
        print(f"  {tag:11s}: real {r['real_cosine_mean']:.3f}  "
              f"null-A {r['null_A_random_spectrum_mean']:.3f}  "
              f"null-B {r['null_B_ion_permutation_mean']:.3f}  "
              f"(real-nullA {r['real_minus_nullA']:+.3f})")
    print("=" * 70)
    print(f"\nwrote {OUT_JSON}")


if __name__ == "__main__":
    main()
