#!/usr/bin/env python3
"""build_hf_tier1.py — HF corpus Tier-1 builder (one accession = one Slurm task).

Implements TO_HF_CORPUS.md v5 §4.0/§4.1/§3.2/§3.4/§3.6:

  inner-join precursor_index + raw_features on (raw_file, precursor_id)
  (after a strict 1:1 alignment gate), apply the per-engine q<=0.01 rank-1
  target floor, then a PROJECTION blob pass (one blob at a time, MS1 + metadata
  only -- never the raw 4D cloud) to compute trace_snr_v1 and per-event
  collision_energies_v. Writes incremental Parquet row groups + a sidecar
  manifest with failure accounting and input checksums.

NOTE: this is the v5 first implementation. Two simplifications are flagged
inline (#SIMPLIFY) for review: peptidoform agreement uses sequence_normalized
(unmod) + charge; canonical mod-string normalization is a follow-up.

Run (compute node; venv needs its Python module):
  module load lang/Python/3.12.3-GCCcore-13.3.0
  .venv/bin/python scripts/build_hf_tier1.py \
      --accession PXD019086 --data-root /lustre/.../data --out hf_tier1
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

SCHEMA_VERSION = "tier1.v5"
SNR_VERSION_RT = "trace_snr_v1_rt"   # short-trace estimator for RT XICs
SNR_VERSION_IM = "trace_snr_v1"      # baseline-region estimator for IM mobilograms
ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"
GZIP_MAGIC = b"\x1f\x8b"
PROTON = 1.007276

# Explicit output schema (Codex §2: union null-fill must be schema-driven;
# also prevents per-row-group type drift -> the Mode-B ParquetWriter crash).
_F = pa.float64()
TIER1_SCHEMA = pa.schema([
    ("accession", pa.string()), ("raw_file", pa.string()),
    ("precursor_id", pa.int64()), ("sage_psm_id", _F),
    ("sequence", pa.string()), ("modified_sequence", pa.string()),
    ("sequence_normalized", pa.string()), ("charge", pa.int64()),
    ("protein", pa.string()), ("n_engines", pa.int64()),
    ("peptidoform_conflict", pa.bool_()),
    ("sage_assignment_passes_floor", pa.bool_()),
    ("fragpipe_assignment_passes_floor", pa.bool_()),
    ("sage_qvalue", _F), ("sage_pep", _F), ("sage_hyperscore", _F),
    ("sage_cosine", _F), ("fragpipe_qvalue", _F), ("fragpipe_pep", _F),
    ("fragpipe_probability", _F), ("fragpipe_hyperscore", _F),
    ("mz", _F), ("rt_seconds", _F), ("mobility", _F),
    ("sage_rt", _F), ("sage_mobility", _F), ("fragpipe_rt", _F),
    ("fragpipe_mobility", _F),
    ("rt_aligned", _F),   # Sage cross-run-aligned RT (~[0,1]); +3pp RT Pearson vs raw

    ("collision_energies_v", pa.list_(_F)), ("collision_energy_mean_v", _F),
    ("ms1_rt_apex", _F), ("ms1_rt_fwhm", _F), ("ms1_rt_sigma", _F),
    ("ms1_rt_skew", _F), ("ms1_rt_r2", _F), ("ms1_rt_snr", _F),
    ("ms1_rt_fit_status", pa.string()),
    ("ms1_im_apex", _F), ("ms1_im_fwhm", _F), ("ms1_im_sigma", _F),
    ("ms1_im_skew", _F), ("ms1_im_r2", _F), ("ms1_im_snr", _F),
    ("ms1_im_fit_status", pa.string()),
    ("isotope_cosim", _F), ("ms1_iso_0", _F), ("ms1_iso_1", _F),
    ("ms1_iso_2", _F), ("ms1_iso_3", _F), ("ms1_iso_4", _F),
    ("precursor_intensity", _F), ("ms1_total_intensity", _F),
    ("fragment_total_intensity", _F), ("n_peaks", pa.int64()),
    ("n_fragments_merged", pa.int64()),
    ("rt_width_reliable", pa.bool_()), ("im_width_reliable", pa.bool_()),
    ("width_reliable", pa.bool_()), ("strict", pa.bool_()),
])

# Floor (TO_HF_CORPUS §3.2)
Q_MAX = 0.01
# Width-reliability band (provisional, TO_HF_CORPUS §3.6)
R2_MIN = 0.8
SNR_MIN = 20.0


# --------------------------------------------------------------------------
# blob projection + SNR
# --------------------------------------------------------------------------
def _decompress(data: bytes) -> bytes:
    if data[:4] == ZSTD_MAGIC:
        from sanjose._compat import get_zstd_decompressor
        return get_zstd_decompressor().decompress(data)
    if data[:2] == GZIP_MAGIC:
        return gzip.decompress(data)
    raise ValueError(f"unknown compression {data[:4].hex()}")


def read_projection(fh, offset: int, size: int, file_size: int):
    """Decompress one blob; return (metadata, rt_c, rt_i, im_c, im_i).

    Touches ONLY metadata + the four MS1 marginal arrays -- the lazy NPZ
    means the raw 4D cloud / fragment arrays are never materialised.
    """
    if offset < 0 or size <= 0 or offset + size > file_size:
        raise ValueError(f"bad blob span offset={offset} size={size} file={file_size}")
    fh.seek(offset)
    combined = _decompress(fh.read(size))
    mlen = int.from_bytes(combined[:4], "little")
    md = json.loads(combined[4:4 + mlen].decode("utf-8"))
    npz = np.load(io.BytesIO(combined[4 + mlen:]))

    def g(name):
        return npz[name] if name in npz.files else np.array([], dtype=np.float32)

    return (md, g("ms1_rt_coords"), g("ms1_rt_intensities"),
            g("ms1_im_coords"), g("ms1_im_intensities"))


def trace_snr_v1_rt(coords, intens, apex_coord):
    """Short-trace SNR for timsTOF RT XICs (7-13 pts), TO_HF_CORPUS §3.6.

    A baseline-region SNR (trace_snr_v1) is undefined on such short traces
    (the peak fills most points). Instead: noise from the lowest floor(n/3)
    intensities (the baseline floor), with a successive-difference fallback
    when that floor is flat. Returns (snr|nan, status).
    """
    c = np.asarray(coords, dtype=float)
    y = np.asarray(intens, dtype=float)
    m = np.isfinite(c) & np.isfinite(y)
    c, y = c[m], y[m]
    if c.size < 7:
        return float("nan"), "too_few_points"
    y = np.clip(y, 0.0, None)
    order = np.argsort(c)
    c, y = c[order], y[order]
    if apex_coord is not None and np.isfinite(apex_coord):
        ai = int(np.argmin(np.abs(c - apex_coord)))
    else:
        ai = int(np.argmax(y))
    if ai == 0 or ai == c.size - 1:          # clear interior apex required
        return float("nan"), "apex_at_edge"
    k = max(3, c.size // 3)
    low = np.sort(y)[:k]                       # baseline-floor points
    base = float(np.median(low))
    noise = 1.4826 * float(np.median(np.abs(low - base)))
    if noise <= 0:                            # flat floor -> successive-diff noise
        d = np.diff(y)
        noise = 1.4826 * float(np.median(np.abs(d - np.median(d)))) / np.sqrt(2.0)
    if noise <= 0:
        return float("nan"), "zero_noise"
    return max(0.0, float(y[ai]) - base) / noise, "ok"


def trace_snr_v1(coords, intens, apex_coord, sigma):
    """Frozen baseline-region SNR (TO_HF_CORPUS §3.6); used for IM mobilograms
    (50-200 pts). Returns (snr|nan, status)."""
    if sigma is None or not np.isfinite(sigma) or sigma <= 0:
        return float("nan"), "no_sigma"
    if apex_coord is None or not np.isfinite(apex_coord):
        return float("nan"), "no_apex"
    c = np.asarray(coords, dtype=float)
    y = np.asarray(intens, dtype=float)
    m = np.isfinite(c) & np.isfinite(y)
    c, y = c[m], y[m]
    if c.size < 7:
        return float("nan"), "too_few_points"
    y = np.clip(y, 0.0, None)
    order = np.argsort(c)
    c, y = c[order], y[order]
    ai = int(np.argmin(np.abs(c - apex_coord)))
    if ai < 3 or ai > c.size - 4:          # >=3 points each side of apex
        return float("nan"), "apex_near_edge"
    lo, hi = apex_coord - 2.5 * sigma, apex_coord + 2.5 * sigma
    if lo <= c[0] or hi >= c[-1]:          # exclusion interval hits a boundary
        return float("nan"), "truncated"
    base_mask = (c < lo) | (c > hi)
    if int(base_mask.sum()) < 6:
        return float("nan"), "too_few_baseline"
    base = float(np.median(y[base_mask]))
    noise = 1.4826 * float(np.median(np.abs(y[base_mask] - base)))
    if noise <= 0:
        return float("nan"), "zero_noise"
    return max(0.0, float(y[ai]) - base) / noise, "ok"


# --------------------------------------------------------------------------
# join + 1:1 gate
# --------------------------------------------------------------------------
def _norm_raw(s):
    """Safe, deterministic join-key normalization: strip the trailing Bruker
    `.d` directory suffix (precursor_index drops it, raw_features keeps it).
    This is the ONLY rewriting applied — no %20/basename munging (§8)."""
    return s[:-2] if isinstance(s, str) and s.endswith(".d") else s


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def alignment_gate(pidx: pa.Table, rf: pa.Table) -> None:
    """Assert precursor_index <-> raw_features align 1:1 on (raw_file,
    precursor_id) (TO_HF_CORPUS §4.0/§8). Raises on PXD062839-class anomalies.
    """
    def keyset(tbl):
        rfn = [_norm_raw(x) for x in tbl["raw_file"].to_pylist()]
        pid = tbl["precursor_id"].to_pylist()
        if any(r is None or p is None for r, p in zip(rfn, pid)):
            raise AssertionError("null key in join input")
        keys = list(zip(rfn, pid))
        s = set(keys)
        if len(s) != len(keys):
            raise AssertionError("duplicate (raw_file, precursor_id) keys")
        return s

    ks_p, ks_r = keyset(pidx), keyset(rf)
    if ks_p != ks_r:
        raise AssertionError(
            f"key sets differ: pidx={len(ks_p)} rf={len(ks_r)} "
            f"only_pidx={len(ks_p - ks_r)} only_rf={len(ks_r - ks_p)}")


# --------------------------------------------------------------------------
# per-row assembly
# --------------------------------------------------------------------------
def first_present(row, *cols):
    for c in cols:
        if c in row and row[c] is not None:
            return row[c]
    return None


def build(acc: str, data_root: Path, out_dir: Path,
          q_max: float, fail_frac_max: float, limit: int = 0) -> dict:
    proc = data_root / "processed" / acc
    extr = data_root / "extracted" / acc
    # per-group/HLA datasets keep precursor_index + raw_features in subdirs and
    # may split across several files -> gather ALL and concat (not just [0]).
    pidx_hits = sorted({p for p in [proc / "precursor_index.parquet",
                                    *proc.rglob("precursor_index.parquet")]
                        if p.exists()})
    rf_hits = sorted(set(extr.rglob("raw_features.parquet")))
    if not pidx_hits or not rf_hits:
        raise SystemExit(f"missing inputs for {acc} "
                         f"(pidx={len(pidx_hits)} rf={len(rf_hits)})")

    def _read_concat(paths):
        tabs = [pq.read_table(p) for p in paths]
        return tabs[0] if len(tabs) == 1 else pa.concat_tables(
            tabs, promote_options="default")

    pidx = _read_concat(pidx_hits)
    rf = _read_concat(rf_hits)
    alignment_gate(pidx, rf)

    # Sage cross-run-aligned RT (validated +3pp RT Pearson vs raw rt_seconds):
    # {sage_psm_id: aligned_rt} from results.sage.parquet (rank-1 targets q<=q_max).
    aligned = {}
    for rp in sorted(proc.rglob("results.sage.parquet")):
        try:
            cols = [c for c in ["psm_id", "aligned_rt", "spectrum_q", "rank", "is_decoy"]
                    if c in pq.read_schema(rp).names]
            rt = pq.read_table(rp, columns=cols)
            if {"spectrum_q", "rank", "is_decoy"} <= set(cols):
                import pyarrow.compute as _pc
                rt = rt.filter(_pc.and_(_pc.and_(_pc.equal(rt["is_decoy"], False),
                                                 _pc.equal(rt["rank"], 1)),
                                        _pc.less_equal(rt["spectrum_q"], q_max)))
            for pid, a in zip(rt["psm_id"].to_pylist(), rt["aligned_rt"].to_pylist()):
                if pid is not None:
                    aligned[int(pid)] = a
        except Exception:
            pass

    pidx_cols = set(pidx.schema.names)
    have_sage = "sage_qvalue" in pidx_cols
    have_fp = "fragpipe_qvalue" in pidx_cols
    source_column_present = {c: (c in pidx_cols) for c in
                             ["sage_qvalue", "sage_pep", "fragpipe_qvalue",
                              "fragpipe_pep", "sequence_normalized", "n_engines"]}

    # index raw_features rows by key for the projection pass
    rf_d = rf.to_pylist()
    rf_by_key = {(_norm_raw(r["raw_file"]), r["precursor_id"]): r for r in rf_d}
    pidx_d = pidx.to_pylist()

    # group kept keys by raw_file (the .d form from raw_features, for blob
    # resolution), sorted by blob_offset (sequential reads)
    by_raw: dict[str, list] = {}
    n_floor = 0
    for prow in pidx_d:
        key = (_norm_raw(prow["raw_file"]), prow["precursor_id"])
        rrow = rf_by_key.get(key)
        if rrow is None:
            continue
        sage_q = prow.get("sage_qvalue") if have_sage else None
        fp_q = prow.get("fragpipe_qvalue") if have_fp else None
        sage_pass = sage_q is not None and sage_q <= q_max
        fp_pass = fp_q is not None and fp_q <= q_max
        if not (sage_pass or fp_pass):
            continue                       # PSM floor (§3.2)
        n_floor += 1
        by_raw.setdefault(rrow["raw_file"], []).append(
            (rrow["blob_offset"], prow, rrow, sage_pass, fp_pass))
        if limit and n_floor >= limit:
            break

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{acc}.parquet"
    writer = None
    fail = Counter()
    n_written = 0
    snr_ok = Counter()

    for raw_file in sorted(by_raw):
        items = sorted(by_raw[raw_file], key=lambda t: (t[0] if t[0] is not None else -1))
        blob_path = _resolve_blob(extr, raw_file)
        file_size = blob_path.stat().st_size if blob_path else 0
        fh = open(blob_path, "rb") if blob_path else None
        batch = []
        try:
            for offset, prow, rrow, sage_pass, fp_pass in items:
                rec = assemble(acc, raw_file, prow, rrow, sage_pass, fp_pass,
                               have_sage, have_fp, aligned)
                # projection blob read for SNR + per-event CE
                if fh is not None and offset is not None and rrow.get("blob_size"):
                    try:
                        md, rtc, rti, imc, imi = read_projection(
                            fh, offset, rrow["blob_size"], file_size)
                        if md.get("precursor_id") not in (None, prow["precursor_id"]):
                            raise ValueError("blob precursor_id mismatch")
                        ces = md.get("collision_energies") or []
                        rec["collision_energies_v"] = [float(x) for x in ces]
                        rec["collision_energy_mean_v"] = (
                            float(np.mean(ces)) if ces else rrow.get("collision_energy"))
                        snr_rt, st_rt = trace_snr_v1_rt(rtc, rti, rrow.get("ms1_rt_apex"))
                        snr_im, st_im = trace_snr_v1(imc, imi, rrow.get("ms1_im_apex"),
                                                     rrow.get("ms1_im_sigma"))
                        rec["ms1_rt_snr"] = None if np.isnan(snr_rt) else snr_rt
                        rec["ms1_im_snr"] = None if np.isnan(snr_im) else snr_im
                        rec["ms1_rt_fit_status"] = st_rt
                        rec["ms1_im_fit_status"] = st_im
                        snr_ok[st_rt == "ok"] += 1
                    except Exception as e:
                        fail[type(e).__name__] += 1
                        _null_blob_fields(rec, rrow, status="decode_failed")
                else:
                    fail["no_blob"] += 1
                    _null_blob_fields(rec, rrow, status="not_extracted")
                _derive_reliability(rec)
                batch.append(rec)
        finally:
            if fh is not None:
                fh.close()

        if batch:
            tbl = pa.Table.from_pylist(batch, schema=TIER1_SCHEMA)
            if writer is None:
                writer = pq.ParquetWriter(out_path, TIER1_SCHEMA, compression="zstd")
            writer.write_table(tbl)
            n_written += len(batch)

    if writer is not None:
        writer.close()

    n_attempt = sum(len(v) for v in by_raw.values())
    fail_frac = (sum(fail.values()) / n_attempt) if n_attempt else 0.0
    manifest = {
        "accession": acc, "schema_version": SCHEMA_VERSION,
        "snr_version_rt": SNR_VERSION_RT, "snr_version_im": SNR_VERSION_IM,
        "q_max": q_max,
        "input_shas": {str(p): sha256(p) for p in (*pidx_hits, *rf_hits)},
        "n_input_precursors": len(pidx_d), "n_passed_floor": n_floor,
        "n_written": n_written, "blob_failures": dict(fail),
        "blob_fail_frac": round(fail_frac, 5),
        "snr_ok_counts": {str(k): v for k, v in snr_ok.items()},
        "source_column_present": source_column_present,
        "built_unix": int(time.time()),
    }
    (out_dir / f"{acc}.manifest.json").write_text(json.dumps(manifest, indent=2))

    if fail_frac > fail_frac_max:
        raise SystemExit(
            f"QUARANTINE {acc}: blob_fail_frac={fail_frac:.3f} > {fail_frac_max}")
    return manifest


def _resolve_blob(extr: Path, raw_file: str):
    clean = raw_file.replace(".d", "")
    for p in (extr / f"{clean}.d" / "blobs.bin", extr / raw_file / "blobs.bin"):
        if p.exists():
            return p
    return None


def _null_blob_fields(rec, rrow, status):
    rec["collision_energies_v"] = None
    rec["collision_energy_mean_v"] = rrow.get("collision_energy")
    rec["ms1_rt_snr"] = None
    rec["ms1_im_snr"] = None
    rec["ms1_rt_fit_status"] = status
    rec["ms1_im_fit_status"] = status


def _derive_reliability(rec):
    def ok(r2, snr):
        return (r2 is not None and r2 >= R2_MIN and
                snr is not None and snr >= SNR_MIN)
    rec["rt_width_reliable"] = ok(rec.get("ms1_rt_r2"), rec.get("ms1_rt_snr"))
    rec["im_width_reliable"] = ok(rec.get("ms1_im_r2"), rec.get("ms1_im_snr"))
    rec["width_reliable"] = rec["rt_width_reliable"] and rec["im_width_reliable"]
    rec["strict"] = (rec.get("n_engines") == 2 and rec["width_reliable"]
                     and (rec.get("sage_qvalue") or 1) <= Q_MAX
                     and (rec.get("fragpipe_qvalue") or 1) <= Q_MAX)


def assemble(acc, raw_file, prow, rrow, sage_pass, fp_pass, have_sage, have_fp,
             aligned=None):
    # canonical assignment: Sage-preferred-else-FragPipe (§4.0)
    use_sage = sage_pass
    if use_sage:
        seq = prow.get("sage_peptide"); modseq = prow.get("sage_modified")
        prot = prow.get("sage_protein"); charge = prow.get("sage_charge")
    else:
        seq = prow.get("fragpipe_peptide"); modseq = prow.get("fragpipe_modified")
        prot = prow.get("fragpipe_protein"); charge = prow.get("fragpipe_charge")
    charge = int(charge) if charge is not None else prow.get("charge")
    # #SIMPLIFY: peptidoform agreement on sequence_normalized + charge (unmod)
    sn = prow.get("sequence_normalized")
    both_pass = sage_pass and fp_pass
    conflict = bool(both_pass and prow.get("sage_charge") != prow.get("fragpipe_charge"))
    n_eng = 1 + int(both_pass and not conflict)
    return {
        "accession": acc, "raw_file": raw_file,
        "precursor_id": prow["precursor_id"],
        "sage_psm_id": prow.get("sage_psm_id"),
        "sequence": seq, "modified_sequence": modseq,
        "sequence_normalized": sn, "charge": charge, "protein": prot,
        # engine status (§3.4)
        "n_engines": n_eng, "peptidoform_conflict": conflict,
        "sage_assignment_passes_floor": sage_pass,
        "fragpipe_assignment_passes_floor": fp_pass,
        "sage_qvalue": prow.get("sage_qvalue") if have_sage else None,
        "sage_pep": prow.get("sage_pep") if have_sage else None,
        "sage_hyperscore": prow.get("sage_hyperscore") if have_sage else None,
        "sage_cosine": rrow.get("sage_cosine"),
        "fragpipe_qvalue": prow.get("fragpipe_qvalue") if have_fp else None,
        "fragpipe_pep": prow.get("fragpipe_pep") if have_fp else None,
        "fragpipe_probability": prow.get("fragpipe_probability") if have_fp else None,
        "fragpipe_hyperscore": prow.get("fragpipe_hyperscore") if have_fp else None,
        # apex (per-engine kept; consensus convenience = canonical)
        "mz": rrow.get("mz"), "rt_seconds": rrow.get("rt_seconds"),
        "mobility": rrow.get("mobility"),
        "sage_rt": prow.get("sage_rt"), "sage_mobility": prow.get("sage_mobility"),
        "fragpipe_rt": prow.get("fragpipe_rt"),
        "fragpipe_mobility": prow.get("fragpipe_mobility"),
        "rt_aligned": ((aligned.get(int(prow["sage_psm_id"]))
                        if aligned and prow.get("sage_psm_id") is not None else None)),
        # width labels (RT s; IM 1/K0)
        "ms1_rt_apex": rrow.get("ms1_rt_apex"), "ms1_rt_fwhm": rrow.get("ms1_rt_fwhm"),
        "ms1_rt_sigma": rrow.get("ms1_rt_sigma"), "ms1_rt_skew": rrow.get("ms1_rt_skew"),
        "ms1_rt_r2": rrow.get("ms1_rt_r2"),
        "ms1_im_apex": rrow.get("ms1_im_apex"), "ms1_im_fwhm": rrow.get("ms1_im_fwhm"),
        "ms1_im_sigma": rrow.get("ms1_im_sigma"), "ms1_im_skew": rrow.get("ms1_im_skew"),
        "ms1_im_r2": rrow.get("ms1_im_r2"),
        # isotope + intensity
        "isotope_cosim": rrow.get("isotope_cosim"),
        "ms1_iso_0": rrow.get("ms1_iso_0"), "ms1_iso_1": rrow.get("ms1_iso_1"),
        "ms1_iso_2": rrow.get("ms1_iso_2"), "ms1_iso_3": rrow.get("ms1_iso_3"),
        "ms1_iso_4": rrow.get("ms1_iso_4"),
        "precursor_intensity": rrow.get("precursor_intensity"),
        "ms1_total_intensity": rrow.get("ms1_total_intensity"),
        "fragment_total_intensity": rrow.get("fragment_total_intensity"),
        "n_peaks": rrow.get("n_peaks"),
        "n_fragments_merged": rrow.get("n_fragments_merged"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True)
    ap.add_argument("--data-root", default="/lustre/project/ki-proanagi/dateschn/data")
    ap.add_argument("--out", default="hf_tier1")
    ap.add_argument("--repo-root", default="/lustre/project/ki-proanagi/dateschn/claudius-proteomics",
                    help="for importing sanjose")
    ap.add_argument("--q-max", type=float, default=Q_MAX)
    ap.add_argument("--fail-frac-max", type=float, default=0.02,
                    help="quarantine the dataset above this blob-failure fraction")
    ap.add_argument("--limit", type=int, default=0, help="smoke-test: cap kept precursors")
    args = ap.parse_args()
    sys.path.insert(0, args.repo_root)
    m = build(args.accession, Path(args.data_root), Path(args.out),
              args.q_max, args.fail_frac_max, args.limit)
    print(json.dumps(m, indent=2))


if __name__ == "__main__":
    main()
