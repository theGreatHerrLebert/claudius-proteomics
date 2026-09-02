"""Cluster-side independent intensity analysis -- one dataset, compact result.

Runs entirely on the cluster against /lustre, so the multi-100-GB blobs never
move. Per dataset it: calibrates NCE, fine-tunes the intensity predictor on
Sage matched-ion targets, and scores held-out test PSMs on two supports --
Sage matched ions and a search-engine-independent sum-integrated readout of
the raw PASEF blobs -- plus random-spectrum and ion-permutation chance floors.

Only a small JSON comes back: notebook/analysis ... no, written to
{data_root}/analysis/{accession}/intensity_independent.json  (~few KB).

Usage:  python cluster_intensity_analysis.py PXD050342 [--data-root DIR]

This consolidates, for batch use, the locally-validated logic of
scripts/analysis/{cross_lab_intensity_finetune,cross_lab_independent_rescore,
random_baseline}.py. See memory: blob-intensity-readout, intensity-finetune-cross-lab.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "analysis"))

from sagepy_rescore.sage_loader import load_psms_from_sage_parquets
from sagepy_rescore.predictors import apply_imspy_patches
from imspy_predictors.intensity.predictors import DeepPeptideIntensityPredictor, calibrate_nce
from fragment_matching import FragmentMatcher, MatchConfig
from extract_fragment_peaks import read_blob, resolve_blob_path
from blob_intensity import integrate_peaks_for_psm

apply_imspy_patches()

# --- fixed parameters (match the validated local scripts) -------------------
PEP_THRESHOLD     = 1e-3
MIN_MATCHED_PEAKS = 6
N_SAMPLE          = 3000
NCE_GRID          = list(range(12, 51))
FT_QVALUE         = 0.01
FT_MAX            = 20000
SEED              = 0
TEST_FRAC         = 0.20
TOL_PPM           = 20.0
GRID_DA           = 0.01
N_DRAWS           = 10
ION_AXIS          = {0: 1, 1: 0}        # grid[ordinal-1, ion_axis, charge-1]


def cosine(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size < 2 or a.sum() <= 0 or b.sum() <= 0:
        return np.nan
    return float(np.clip(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1))


def grid_value(grid, it, chg, ordn):
    if ordn < 1 or ordn > 29 or chg < 1 or chg > 3:
        return None
    return float(grid[ordn - 1, ION_AXIS[it], chg - 1])


def filter_results(res, *, pep=None, qval=None):
    mask = pc.greater_equal(res["matched_peaks"], MIN_MATCHED_PEAKS)
    mask = pc.and_(mask, pc.equal(res["is_decoy"], False))
    mask = pc.and_(mask, pc.less_equal(res["rank"], 1))
    if pep is not None:
        mask = pc.and_(mask, pc.less_equal(res["posterior_error"], float(np.log(pep))))
    if qval is not None:
        mask = pc.and_(mask, pc.less_equal(res["spectrum_q"], qval))
    return res.filter(mask)


def write_subset(res, frag, cache, tag):
    keep = pa.array(list(set(res.column("psm_id").to_pylist())))
    fsub = frag.filter(pc.is_in(frag["psm_id"], value_set=keep))
    rp, fp = cache / f"{tag}_results.parquet", cache / f"{tag}_fragments.parquet"
    pq.write_table(res, rp); pq.write_table(fsub, fp)
    return rp, fp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accession")
    ap.add_argument("--data-root", default="/lustre/project/ki-proanagi/dateschn/data")
    args = ap.parse_args()
    acc = args.accession
    root = Path(args.data_root)

    sage_dir = root / "processed" / acc / "sage"
    res_pq   = sage_dir / "results.sage.parquet"
    frag_pq  = sage_dir / "matched_fragments.sage.parquet"
    store    = root / "merged" / acc / "precursor_store.parquet"
    blob_dir = root / "extracted" / acc
    out_dir  = root / "analysis" / acc
    out_dir.mkdir(parents=True, exist_ok=True)
    out_json = out_dir / "intensity_independent.json"
    cache    = out_dir / "cache"
    cache.mkdir(exist_ok=True)

    t0 = time.time()
    result = {"accession": acc, "status": "running",
              "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

    missing = [str(p) for p in (res_pq, frag_pq, store) if not p.exists()]
    if missing or not blob_dir.exists():
        result["status"] = "missing_inputs"
        result["missing"] = missing + ([] if blob_dir.exists() else [str(blob_dir)])
        out_json.write_text(json.dumps(result, indent=2))
        print(f"[{acc}] missing inputs: {result['missing']}", flush=True)
        return

    try:
        matcher = FragmentMatcher(MatchConfig(mz_tolerance_ppm=TOL_PPM,
                                              ion_types=["b", "y"], max_fragment_charge=2))

        # --- calibrate NCE on the strict PEP<=1e-3 sample ----------------
        res = filter_results(pq.read_table(res_pq), pep=PEP_THRESHOLD)
        n_calib = res.num_rows
        if N_SAMPLE and n_calib > N_SAMPLE:
            idx = np.sort(np.random.RandomState(SEED).choice(n_calib, N_SAMPLE, replace=False))
            res = res.take(pa.array(idx))
        rp, fp = write_subset(res, pq.read_table(frag_pq), cache, "calib")
        calib_psms = load_psms_from_sage_parquets(str(rp), str(fp),
                                                  default_collision_energy=30.0, max_rank=1)
        model = DeepPeptideIntensityPredictor(verbose=False)
        cal = calibrate_nce(model, calib_psms, nce_grid=NCE_GRID, max_sample=N_SAMPLE)
        best_nce = int(cal["best_nce"])

        # --- fine-tune set: q<=0.01, peptide-level split -----------------
        res = filter_results(pq.read_table(res_pq), qval=FT_QVALUE)
        if FT_MAX and res.num_rows > FT_MAX:
            idx = np.sort(np.random.RandomState(SEED).choice(res.num_rows, FT_MAX, replace=False))
            res = res.take(pa.array(idx))
        rp, fp = write_subset(res, pq.read_table(frag_pq), cache, "ft")
        psms = load_psms_from_sage_parquets(str(rp), str(fp),
                                            default_collision_energy=float(best_nce), max_rank=1)
        for p in psms:
            p.collision_energy = p.collision_energy_calibrated = float(best_nce)
        rng = np.random.RandomState(SEED)
        peps = sorted({p.sequence for p in psms})
        rng.shuffle(peps)
        test_peps = set(peps[:int(TEST_FRAC * len(peps))])
        train = [p for p in psms if p.sequence not in test_peps]
        test  = [p for p in psms if p.sequence in test_peps]

        # --- blob index + independent integration of the test set -------
        t = pq.read_table(store, columns=["sage_psm_id", "raw_file", "blob_offset",
                                          "blob_size", "mono_mz"])
        df = t.filter(t["sage_psm_id"].is_valid()).to_pandas()
        df["sage_psm_id"] = df["sage_psm_id"].astype("int64")
        bidx = {r.sage_psm_id: (r.raw_file, int(r.blob_offset), int(r.blob_size), float(r.mono_mz))
                for r in df.itertuples()}

        indep = {}     # psm_id -> (keys, real_obs, blob_mz, blob_int)
        for p in test:
            loc = bidx.get(int(p.sage_feature.psm_id))
            if loc is None:
                continue
            raw_file, off, size, prec_mz = loc
            blob = read_blob(resolve_blob_path(blob_dir, raw_file), off, size)
            if blob is None or len(blob["frag_mz"]) == 0:
                continue
            theo = matcher.generate_theoretical_fragments(p.sequence_modified, p.charge)
            pk = integrate_peaks_for_psm(blob["frag_mz"], blob["frag_intensity"], theo,
                                         TOL_PPM, precursor_mz=prec_mz)
            keys = [(0 if it == "b" else 1, int(c), int(n))
                    for it, n, c in zip(pk["ion_type"], pk["ion_number"], pk["ion_charge"])]
            indep[int(p.sage_feature.psm_id)] = (theo, keys, np.array(pk["intensity"], float),
                                                 np.asarray(blob["frag_mz"], float),
                                                 np.asarray(blob["frag_intensity"], float))

        def score(m):
            seqs = [p.sequence_modified for p in test]
            chs  = [p.charge for p in test]
            grids = m.predict_intensities(seqs, chs, [float(best_nce)] * len(test),
                                          batch_size=1024, flatten=False)
            cs, ci = [], []
            for p, g in zip(test, grids):
                g = np.asarray(g)
                om = p.observed_fragments_map()
                o, pr = [], []
                for (it, chg, ordn), iv in om.items():
                    gv = grid_value(g, it, chg, ordn)
                    if gv is not None:
                        o.append(iv); pr.append(gv)
                cs.append(cosine(o, pr))
                rec = indep.get(int(p.sage_feature.psm_id))
                if rec is None:
                    ci.append(np.nan); continue
                _, keys, real, _, _ = rec
                pv = [grid_value(g, it, chg, ordn) or 0.0 for (it, chg, ordn) in keys]
                ci.append(cosine(real, pv))
            return float(np.nanmean(cs)), float(np.nanmean(ci))

        pre = DeepPeptideIntensityPredictor(verbose=False)
        pre_sage, pre_indep = score(pre)

        ft = DeepPeptideIntensityPredictor(verbose=False)
        ft.fine_tune_psms(train, batch_size=64, epochs=50, patience=5, verbose=False)
        ft_sage, ft_indep = score(ft)

        # --- random baselines on the fine-tuned predictions --------------
        rs = np.random.RandomState(SEED)
        seqs = [p.sequence_modified for p in test]
        chs  = [p.charge for p in test]
        ft_grids = [np.asarray(x) for x in
                    ft.predict_intensities(seqs, chs, [float(best_nce)] * len(test),
                                           batch_size=1024, flatten=False)]
        nullA, nullB = [], []
        for p, g in zip(test, ft_grids):
            rec = indep.get(int(p.sage_feature.psm_id))
            if rec is None:
                continue
            theo, keys, real, bmz, bint = rec
            pv = np.array([grid_value(g, it, chg, ordn) or 0.0 for (it, chg, ordn) in keys])
            lo, hi = bmz.min(), bmz.max()
            ng = max(int((hi - lo) / GRID_DA), 1)
            a = []
            for _ in range(N_DRAWS):
                rmz = lo + GRID_DA * rs.randint(0, ng + 1, size=bint.size)
                pk = integrate_peaks_for_psm(rmz, bint, theo, TOL_PPM)
                m = dict(zip([(0 if t2 == "b" else 1, int(c2), int(n2)) for t2, n2, c2 in
                              zip(pk["ion_type"], pk["ion_number"], pk["ion_charge"])],
                             pk["intensity"]))
                a.append(cosine(pv, np.array([m.get(k, 0.0) for k in keys])))
            nullA.append(np.nanmean(a))
            nullB.append(np.nanmean([cosine(pv, rs.permutation(real)) for _ in range(N_DRAWS)]))

        result.update({
            "status": "ok",
            "n_calib_pep1e-3": int(n_calib), "best_nce": best_nce,
            "n_finetune_psms": len(psms), "n_peptides": len(peps),
            "n_train": len(train), "n_test": len(test), "n_test_with_blob": len(indep),
            "held_out_test_cosine": {
                "sage_support":  {"pretrained": pre_sage,  "finetuned": ft_sage,
                                  "gain": ft_sage - pre_sage},
                "independent":   {"pretrained": pre_indep, "finetuned": ft_indep,
                                  "gain": ft_indep - pre_indep},
                "gain_retained_independent_frac":
                    ((ft_indep - pre_indep) / (ft_sage - pre_sage))
                    if (ft_sage - pre_sage) else None,
            },
            "chance_floor_finetuned": {
                "null_A_random_spectrum": float(np.nanmean(nullA)) if nullA else None,
                "null_B_ion_permutation": float(np.nanmean(nullB)) if nullB else None,
            },
            "runtime_seconds": round(time.time() - t0, 1),
        })
        print(f"[{acc}] ok | NCE {best_nce} | independent cos "
              f"{pre_indep:.3f}->{ft_indep:.3f} | nullA "
              f"{result['chance_floor_finetuned']['null_A_random_spectrum']}", flush=True)

    except Exception as e:
        import traceback
        result["status"] = "error"
        result["error"] = f"{type(e).__name__}: {e}"
        result["traceback"] = traceback.format_exc()
        print(f"[{acc}] ERROR: {e}", flush=True)

    out_json.write_text(json.dumps(result, indent=2))
    print(f"[{acc}] wrote {out_json}", flush=True)


if __name__ == "__main__":
    main()
