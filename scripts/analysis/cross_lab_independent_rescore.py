"""Search-engine-independent re-scoring of the intensity fine-tune (v2).

v2 changes vs the first run:
  - blob intensities are SUM-INTEGRATED (blob_intensity.integrate_peaks_for_psm)
    rather than single-closest-point sampled. The blob is a raw point cloud;
    single-point sampling under-read every ion ~60x. This is the fixed readout.
  - metrics reported as raw cosine similarity (primary) AND spectral angle.
    SA is nonlinear and flatters mediocre matches -- cosine is the honest read.

Question unchanged: the fine-tune is trained on Sage matched-ion targets; does
its held-out gain survive on a search-engine-independent support (every
theoretical b/y ion, observed intensity integrated straight from the raw
PASEF blob)? If yes -> real; if no -> Sage matched-ion artifact.

Dataset: PXD046777 (only dataset with raw blobs available locally).
Output: notebook/analysis/intensity_exploration/cross_lab/independent_rescore_integrated.json
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

# --- config ------------------------------------------------------------------
DATASET   = "PXD046777"
BEST_NCE  = 47
SEED      = 0
TEST_FRAC = 0.20
TOL_PPM   = 20.0

DATA_ROOT = Path("/scratch/claudius-proteomics") / DATASET
CACHE     = PROJECT_ROOT / "notebook/analysis/intensity_exploration/cross_lab" / f"{DATASET}_cache"
STORE     = DATA_ROOT / "merged" / "precursor_store.parquet"
BLOB_DIR  = DATA_ROOT / "extracted"
OUT       = PROJECT_ROOT / "notebook/analysis/intensity_exploration/cross_lab" / "independent_rescore_integrated.json"

# (29,2,3) predictor grid: grid[ordinal-1, ion_axis, charge-1]; ion_axis y=0,b=1
ION_AXIS = {0: 1, 1: 0}


# --- metric ------------------------------------------------------------------
def cos_and_sa(obs, pred):
    """Returns (cosine similarity, spectral angle = 1 - 2*arccos(cos)/pi)."""
    obs = np.asarray(obs, float)
    pred = np.asarray(pred, float)
    if obs.size < 2 or obs.sum() <= 0 or pred.sum() <= 0:
        return np.nan, np.nan
    cos = float(np.dot(obs, pred) / (np.linalg.norm(obs) * np.linalg.norm(pred)))
    cos = float(np.clip(cos, -1.0, 1.0))
    return cos, 1.0 - 2.0 * np.arccos(cos) / np.pi


# --- data --------------------------------------------------------------------
def load_split():
    psms = load_psms_from_sage_parquets(
        str(CACHE / "ft_results.parquet"), str(CACHE / "ft_fragments.parquet"),
        default_collision_energy=float(BEST_NCE), max_rank=1)
    for p in psms:
        p.collision_energy = float(BEST_NCE)
        p.collision_energy_calibrated = float(BEST_NCE)
    rng = np.random.RandomState(SEED)
    peps = sorted({p.sequence for p in psms})
    rng.shuffle(peps)
    test_peps = set(peps[:int(TEST_FRAC * len(peps))])
    train = [p for p in psms if p.sequence not in test_peps]
    test  = [p for p in psms if p.sequence in test_peps]
    return psms, train, test


def blob_index():
    """sage_psm_id -> (raw_file, blob_offset, blob_size, precursor_mz)."""
    t = pq.read_table(STORE, columns=["sage_psm_id", "raw_file", "blob_offset",
                                      "blob_size", "mono_mz"])
    df = t.filter(t["sage_psm_id"].is_valid()).to_pandas()
    df["sage_psm_id"] = df["sage_psm_id"].astype("int64")
    return {r.sage_psm_id: (r.raw_file, int(r.blob_offset), int(r.blob_size),
                            float(r.mono_mz))
            for r in df.itertuples()}


def extract_independent(psms, bidx, matcher):
    """psm_id -> {(ion_type,charge,ord): SUM-INTEGRATED observed intensity}."""
    out, stats = {}, {"n": 0, "no_blob": 0, "no_spectrum": 0}
    for p in psms:
        pid = int(p.sage_feature.psm_id)
        loc = bidx.get(pid)
        if loc is None:
            stats["no_blob"] += 1
            continue
        raw_file, off, size, prec_mz = loc
        blob = read_blob(resolve_blob_path(BLOB_DIR, raw_file), off, size)
        if blob is None or len(blob["frag_mz"]) == 0:
            stats["no_spectrum"] += 1
            continue
        theo = matcher.generate_theoretical_fragments(p.sequence_modified, p.charge)
        peaks = integrate_peaks_for_psm(blob["frag_mz"], blob["frag_intensity"], theo,
                                        TOL_PPM, precursor_mz=prec_mz)
        obs = {(0 if it == "b" else 1, int(chg), int(num)): float(inten)
               for it, num, chg, inten in zip(peaks["ion_type"], peaks["ion_number"],
                                              peaks["ion_charge"], peaks["intensity"])}
        out[pid] = obs
        stats["n"] += 1
    return out, stats


# --- scoring -----------------------------------------------------------------
def predict_grids(model, psms):
    seqs = [p.sequence_modified for p in psms]
    charges = [p.charge for p in psms]
    g = model.predict_intensities(seqs, charges, [float(BEST_NCE)] * len(psms),
                                  batch_size=1024, flatten=False)
    return [np.asarray(x) for x in g]


def grid_value(grid, it, chg, ordn):
    if ordn < 1 or ordn > 29 or chg < 1 or chg > 3:
        return None
    return float(grid[ordn - 1, ION_AXIS[it], chg - 1])


def score(model, psms, indep_obs):
    """Mean cosine + spectral angle on the Sage matched-ion support and on the
    search-engine-independent full b/y support."""
    grids = predict_grids(model, psms)
    cos_sage, sa_sage, cos_ind, sa_ind = [], [], [], []
    pred_on_absent = []
    for p, grid in zip(psms, grids):
        pid = int(p.sage_feature.psm_id)

        om = p.observed_fragments_map()
        o, pr = [], []
        for (it, chg, ordn), iv in om.items():
            gv = grid_value(grid, it, chg, ordn)
            if gv is not None:
                o.append(iv); pr.append(gv)
        c, s = cos_and_sa(o, pr)
        cos_sage.append(c); sa_sage.append(s)

        obs = indep_obs.get(pid)
        if obs is None:
            cos_ind.append(np.nan); sa_ind.append(np.nan)
            continue
        o, pr, absent = [], [], 0.0
        for (it, chg, ordn), iv in obs.items():
            gv = grid_value(grid, it, chg, ordn)
            if gv is None:
                continue
            o.append(iv); pr.append(gv)
            if iv <= 0:
                absent += gv
        c, s = cos_and_sa(o, pr)
        cos_ind.append(c); sa_ind.append(s)
        tot = float(np.sum(pr)) if pr else 0.0
        if tot > 0:
            pred_on_absent.append(absent / tot)

    return {
        "cos_sage": float(np.nanmean(cos_sage)), "sa_sage": float(np.nanmean(sa_sage)),
        "cos_indep": float(np.nanmean(cos_ind)), "sa_indep": float(np.nanmean(sa_ind)),
        "n_sage": int(np.sum(~np.isnan(cos_sage))),
        "n_indep": int(np.sum(~np.isnan(cos_ind))),
        "frac_pred_on_absent": float(np.mean(pred_on_absent)) if pred_on_absent else float("nan"),
    }


# --- main --------------------------------------------------------------------
def main():
    print(f"=== independent re-score (sum-integrated): {DATASET} ===\n", flush=True)
    psms, train, test = load_split()
    print(f"fine-tune set {len(psms):,} -> {len(train):,} train / {len(test):,} test", flush=True)

    bidx = blob_index()
    matcher = FragmentMatcher(MatchConfig(mz_tolerance_ppm=TOL_PPM,
                                          ion_types=["b", "y"], max_fragment_charge=2))

    print("integrating search-engine-independent peaks from raw blobs ...", flush=True)
    indep_test, st_te = extract_independent(test, bidx, matcher)
    indep_train, st_tr = extract_independent(train, bidx, matcher)
    print(f"  test {st_te['n']:,} | train {st_tr['n']:,} extracted", flush=True)

    print("\nscoring pretrained ...", flush=True)
    model = DeepPeptideIntensityPredictor(verbose=False)
    pre_test  = score(model, test,  indep_test)
    pre_train = score(model, train, indep_train)
    print(f"  test  Sage cos {pre_test['cos_sage']:.4f} | independent cos {pre_test['cos_indep']:.4f}",
          flush=True)

    print("\nfine-tuning on train (Sage matched-ion targets) ...", flush=True)
    model.fine_tune_psms(train, batch_size=64, epochs=50, patience=5, verbose=False)

    print("scoring fine-tuned ...", flush=True)
    post_test  = score(model, test,  indep_test)
    post_train = score(model, train, indep_train)
    print(f"  test  Sage cos {post_test['cos_sage']:.4f} | independent cos {post_test['cos_indep']:.4f}",
          flush=True)

    def block(pre, post):
        return {
            "sage_support": {
                "cosine": {"pretrained": pre["cos_sage"], "finetuned": post["cos_sage"],
                           "gain": post["cos_sage"] - pre["cos_sage"]},
                "spectral_angle": {"pretrained": pre["sa_sage"], "finetuned": post["sa_sage"],
                                   "gain": post["sa_sage"] - pre["sa_sage"]}},
            "independent": {
                "cosine": {"pretrained": pre["cos_indep"], "finetuned": post["cos_indep"],
                           "gain": post["cos_indep"] - pre["cos_indep"]},
                "spectral_angle": {"pretrained": pre["sa_indep"], "finetuned": post["sa_indep"],
                                   "gain": post["sa_indep"] - pre["sa_indep"]}},
        }

    result = {
        "dataset": DATASET, "best_nce": BEST_NCE, "seed": SEED,
        "n_train": len(train), "n_test": len(test),
        "blob_intensity_readout": "sum-integrated within 20 ppm (blob_intensity.py)",
        "metric_note": "cosine on raw linear intensities; SA = 1-2*arccos(cos)/pi (nonlinear)",
        "held_out_test": block(pre_test, post_test),
        "train": block(pre_train, post_train),
        "pred_intensity_on_absent_ions_test": {
            "pretrained": pre_test["frac_pred_on_absent"],
            "finetuned":  post_test["frac_pred_on_absent"]},
    }
    OUT.write_text(json.dumps(result, indent=2))

    ho = result["held_out_test"]
    print("\n" + "=" * 70)
    print("HELD-OUT TEST -- fine-tune effect, two supports  (cosine | spectral angle)")
    for name, key in [("Sage matched-ion support", "sage_support"),
                      ("independent full b/y    ", "independent")]:
        c, s = ho[key]["cosine"], ho[key]["spectral_angle"]
        print(f"  {name}: cos {c['pretrained']:.3f}->{c['finetuned']:.3f} ({c['gain']:+.3f})"
              f"   SA {s['pretrained']:.3f}->{s['finetuned']:.3f} ({s['gain']:+.3f})")
    cg_s = ho["sage_support"]["cosine"]["gain"]
    cg_i = ho["independent"]["cosine"]["gain"]
    print(f"  cosine gain retained on independent support: "
          f"{(cg_i / cg_s * 100) if cg_s else float('nan'):.0f}%")
    print(f"  predicted intensity on absent ions: "
          f"{result['pred_intensity_on_absent_ions_test']['pretrained']*100:.1f}% -> "
          f"{result['pred_intensity_on_absent_ions_test']['finetuned']*100:.1f}%")
    print("=" * 70)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
