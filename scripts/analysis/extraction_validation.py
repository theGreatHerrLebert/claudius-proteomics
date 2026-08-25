"""Validate the blob fragment extraction against Sage's reported intensities (v2).

v2: blob intensities are SUM-INTEGRATED within 20 ppm of each theoretical ion
(blob_intensity.integrate_peaks_for_psm) -- the correct readout for a raw
PASEF point cloud -- instead of single-closest-point sampling.

For ions Sage matched, two independent measurements of intensity:
  A = Sage matched_fragments intensities (the intensity fine-tune target)
  B = sum-integrated from the raw PASEF blob
High A-vs-B agreement => the blob extraction is accurate, not just coherent.

Dataset: PXD046777 held-out test split.
Output: notebook/analysis/intensity_exploration/cross_lab/extraction_validation_integrated.json + .png
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pyarrow.parquet as pq
from scipy.stats import pearsonr, spearmanr

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "analysis"))

from sagepy_rescore.sage_loader import load_psms_from_sage_parquets
from sagepy_rescore.predictors import apply_imspy_patches
from fragment_matching import FragmentMatcher, MatchConfig
from extract_fragment_peaks import read_blob, resolve_blob_path
from blob_intensity import integrate_peaks_for_psm

apply_imspy_patches()

DATASET   = "PXD046777"
SEED      = 0
TEST_FRAC = 0.20
TOL_PPM   = 20.0

DATA_ROOT = Path("/scratch/claudius-proteomics") / DATASET
CACHE     = PROJECT_ROOT / "notebook/analysis/intensity_exploration/cross_lab" / f"{DATASET}_cache"
STORE     = DATA_ROOT / "merged" / "precursor_store.parquet"
BLOB_DIR  = DATA_ROOT / "extracted"
OUT_JSON  = PROJECT_ROOT / "notebook/analysis/intensity_exploration/cross_lab" / "extraction_validation_integrated.json"
OUT_PNG   = PROJECT_ROOT / "notebook/analysis/intensity_exploration/cross_lab" / "extraction_validation_integrated.png"


def ion_id(x):
    s = str(x)
    return 0 if "B" in s else (1 if "Y" in s else int(x))


def cosine(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.sum() <= 0 or b.sum() <= 0:
        return np.nan
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def main():
    print(f"=== extraction validation (sum-integrated): {DATASET} ===\n", flush=True)

    psms = load_psms_from_sage_parquets(str(CACHE / "ft_results.parquet"),
                                        str(CACHE / "ft_fragments.parquet"),
                                        default_collision_energy=47.0, max_rank=1)
    rng = np.random.RandomState(SEED)
    peps = sorted({p.sequence for p in psms})
    rng.shuffle(peps)
    test_peps = set(peps[:int(TEST_FRAC * len(peps))])
    test = [p for p in psms if p.sequence in test_peps]
    print(f"held-out test PSMs: {len(test):,}", flush=True)

    t = pq.read_table(STORE, columns=["sage_psm_id", "raw_file", "blob_offset",
                                      "blob_size", "mono_mz"])
    df = t.filter(t["sage_psm_id"].is_valid()).to_pandas()
    df["sage_psm_id"] = df["sage_psm_id"].astype("int64")
    bidx = {r.sage_psm_id: (r.raw_file, int(r.blob_offset), int(r.blob_size), float(r.mono_mz))
            for r in df.itertuples()}

    matcher = FragmentMatcher(MatchConfig(mz_tolerance_ppm=TOL_PPM,
                                          ion_types=["b", "y"], max_fragment_charge=2))

    per_psm_pearson, per_psm_spearman, per_psm_cosine = [], [], []
    pool_sage, pool_blob = [], []
    all_ppm, blob_sage_ratio = [], []
    cover_b, cover_y, intens_expl = [], [], []
    presence_hits, presence_total = 0, 0
    n_done = 0

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
        n_done += 1
        cover_b.append(peaks["coverage_b"])
        cover_y.append(peaks["coverage_y"])
        intens_expl.append(peaks["intensity_explained"])

        blob_int = {}
        for it, num, chg, inten, err in zip(peaks["ion_type"], peaks["ion_number"],
                                            peaks["ion_charge"], peaks["intensity"],
                                            peaks["error_ppm"]):
            blob_int[(0 if it == "b" else 1, int(chg), int(num))] = float(inten)
            if inten > 0 and not np.isnan(err):
                all_ppm.append(float(err))

        fr = p.sage_feature.fragments
        s_int, b_int = [], []
        for it, ordn, chg, inten in zip(fr.ion_types, fr.fragment_ordinals,
                                        fr.charges, fr.intensities):
            key = (ion_id(it), int(chg), int(ordn))
            presence_total += 1
            bv = blob_int.get(key, 0.0)
            if bv > 0:
                presence_hits += 1
                s_int.append(float(inten))
                b_int.append(bv)
        if len(s_int) >= 5:
            s = np.array(s_int); b = np.array(b_int)
            per_psm_pearson.append(pearsonr(s, b)[0])
            per_psm_spearman.append(spearmanr(s, b)[0])
            per_psm_cosine.append(cosine(s, b))
            pool_sage.extend((s / s.sum()).tolist())
            pool_blob.extend((b / b.sum()).tolist())
            blob_sage_ratio.append(b.sum() / s.sum())

    per_psm_pearson = np.array(per_psm_pearson, float)
    per_psm_spearman = np.array(per_psm_spearman, float)
    per_psm_cosine = np.array(per_psm_cosine, float)
    pool_sage = np.array(pool_sage); pool_blob = np.array(pool_blob)
    all_ppm = np.array(all_ppm)

    result = {
        "dataset": DATASET, "n_test_psms": len(test), "n_extracted": n_done,
        "blob_readout": "sum-integrated within 20 ppm",
        "intensity_agreement_A_Sage_vs_B_blob": {
            "per_psm_pearson_median":  float(np.nanmedian(per_psm_pearson)),
            "per_psm_spearman_median": float(np.nanmedian(per_psm_spearman)),
            "per_psm_cosine_median":   float(np.nanmedian(per_psm_cosine)),
            "pooled_pearson":  float(pearsonr(pool_sage, pool_blob)[0]),
            "pooled_spearman": float(spearmanr(pool_sage, pool_blob)[0]),
            "n_psms_scored": int(len(per_psm_pearson)),
            "n_ion_pairs":   int(len(pool_sage)),
            "blob_to_sage_intensity_ratio_median": float(np.median(blob_sage_ratio)),
        },
        "presence_agreement": {
            "frac_sage_ions_blob_also_sees": presence_hits / max(presence_total, 1),
            "n_sage_matched_ions": presence_total,
        },
        "extraction_self_diagnostics": {
            "ppm_error_median": float(np.median(all_ppm)),
            "ppm_error_std":    float(np.std(all_ppm)),
            "ppm_abs_median":   float(np.median(np.abs(all_ppm))),
            "coverage_b_median": float(np.median(cover_b)),
            "coverage_y_median": float(np.median(cover_y)),
            "intensity_explained_median": float(np.median(intens_expl)),
        },
    }
    OUT_JSON.write_text(json.dumps(result, indent=2))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
        ax[0].hist(per_psm_pearson[~np.isnan(per_psm_pearson)], bins=50, color="#2563eb")
        ax[0].axvline(np.nanmedian(per_psm_pearson), color="red", ls="--",
                      label=f"median {np.nanmedian(per_psm_pearson):.3f}")
        ax[0].set_xlabel("per-PSM Pearson r (Sage vs blob intensity)")
        ax[0].set_ylabel("PSMs"); ax[0].set_title("A vs B agreement (sum-integrated)"); ax[0].legend()

        idx = np.random.RandomState(0).choice(len(pool_sage), min(8000, len(pool_sage)), replace=False)
        ax[1].scatter(pool_sage[idx], pool_blob[idx], s=3, alpha=0.25, color="#059669")
        lim = max(pool_sage[idx].max(), pool_blob[idx].max())
        ax[1].plot([0, lim], [0, lim], "k--", lw=1)
        ax[1].set_xlabel("Sage intensity (PSM-normalised)")
        ax[1].set_ylabel("blob intensity (PSM-normalised)")
        ax[1].set_title(f"pooled ions (pearson {result['intensity_agreement_A_Sage_vs_B_blob']['pooled_pearson']:.3f})")

        ax[2].hist(np.clip(all_ppm, -25, 25), bins=120, color="#7c3aed")
        ax[2].axvline(np.median(all_ppm), color="red", ls="--",
                      label=f"median {np.median(all_ppm):.2f} ppm")
        ax[2].set_xlabel("blob ion m/z error (ppm)")
        ax[2].set_ylabel("ions"); ax[2].set_title("extraction m/z self-diagnostic"); ax[2].legend()
        fig.tight_layout()
        fig.savefig(OUT_PNG, dpi=130)
        print(f"saved plot -> {OUT_PNG}", flush=True)
    except Exception as e:
        print(f"plot skipped: {e}", flush=True)

    a = result["intensity_agreement_A_Sage_vs_B_blob"]
    d = result["extraction_self_diagnostics"]
    print("\n" + "=" * 64)
    print(f"extracted {n_done:,}/{len(test):,} test PSMs")
    print("\nA (Sage) vs B (blob, sum-integrated) on Sage-matched ions:")
    print(f"  per-PSM Pearson  median {a['per_psm_pearson_median']:.3f}")
    print(f"  per-PSM Spearman median {a['per_psm_spearman_median']:.3f}")
    print(f"  per-PSM cosine   median {a['per_psm_cosine_median']:.3f}")
    print(f"  pooled Pearson {a['pooled_pearson']:.3f}  ({a['n_ion_pairs']:,} ion pairs)")
    print(f"  blob/Sage intensity ratio (median): {a['blob_to_sage_intensity_ratio_median']:.3f}")
    print(f"  Sage ions the blob also sees: "
          f"{result['presence_agreement']['frac_sage_ions_blob_also_sees']*100:.1f}%")
    print("\nextraction self-diagnostics:")
    print(f"  m/z error median {d['ppm_error_median']:+.2f} ppm (sd {d['ppm_error_std']:.2f}, "
          f"|med| {d['ppm_abs_median']:.2f})")
    print(f"  coverage b {d['coverage_b_median']:.2f}  y {d['coverage_y_median']:.2f}")
    print(f"  intensity explained median {d['intensity_explained_median']:.3f}")
    print("=" * 64)
    print(f"\nwrote {OUT_JSON}")


if __name__ == "__main__":
    main()
