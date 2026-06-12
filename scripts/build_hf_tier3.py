#!/usr/bin/env python3
"""build_hf_tier3.py — HF corpus Tier-3 builder (one accession = one Slurm task).

TO_HF_CORPUS.md v5 §3.3: the matcher-clean, permissive dual-engine fragment
table. Reworks the old extract_fragment_peaks.py path per Codex review #3:

  * NO quality gate / NO `--min-engines` filter — covers the full Tier-1 floor.
  * NO `consensus_peptide` fallback — match the canonical engine's peptide only.
  * DETERMINISTIC one-to-one matching: all ion<->peak candidates within
    tolerance, sorted (abs_ppm, -intensity, ion_type, ordinal, charge,
    peak_idx), each ion AND each peak consumed exactly once.
  * Unmatched ions -> not emitted (we emit one row per MATCHED fragment); a
    matched fragment carries real mz_exp/intensity (never the 0-sentinel).

Per precursor passing the floor: projection-read the blob fragment spectrum
(frag_mz/frag_intensity), pick the canonical peptidoform (Sage-else-FragPipe),
imspy-generate theoretical b/y, one-to-one match, emit matched-fragment rows
with engine-agreement context + Sage-native provenance.

#SIMPLIFY (v1, flagged for review): peptidoform agreement uses
sequence_normalized + charge (as in Tier 1). True mod-level conflict ->
separate per-engine dual-matching is deferred (assignment='consensus' covers
the agreed case; agreement_ppm between two different peptidoforms is undefined
and those are rare here because the merge keyed on sequence_normalized).
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

SCHEMA_VERSION = "tier3.v1"
MATCH_METHOD = "imspy_rematch"
ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"
GZIP_MAGIC = b"\x1f\x8b"
Q_MAX = 0.01
TOL_PPM = 20.0

_F = pa.float64()
TIER3_SCHEMA = pa.schema([
    ("accession", pa.string()), ("raw_file", pa.string()),
    ("precursor_id", pa.int64()), ("sage_psm_id", _F),
    ("assignment", pa.string()), ("canonical_engine", pa.string()),
    ("sequence", pa.string()), ("modified_sequence", pa.string()),
    ("sequence_normalized", pa.string()), ("precursor_charge", pa.int64()),
    ("n_engines", pa.int64()), ("peptidoform_conflict", pa.bool_()),
    ("fragment_type", pa.string()), ("fragment_ordinal", pa.int64()),
    ("fragment_charge", pa.int64()),
    ("fragment_mz_calculated", _F), ("fragment_mz_experimental", _F),
    ("fragment_intensity", _F), ("ppm_error", _F),
    ("match_method", pa.string()),
    ("sage_native_matched", pa.bool_()), ("sage_native_intensity", _F),
])


# --------------------------------------------------------------------------
# blob fragment-spectrum projection
# --------------------------------------------------------------------------
def _decompress(data: bytes) -> bytes:
    if data[:4] == ZSTD_MAGIC:
        from sanjose._compat import get_zstd_decompressor
        return get_zstd_decompressor().decompress(data)
    if data[:2] == GZIP_MAGIC:
        return gzip.decompress(data)
    raise ValueError(f"unknown compression {data[:4].hex()}")


def read_frag(fh, offset, size, file_size):
    """Projection: decompress one blob, return (frag_mz, frag_intensity) only."""
    if offset < 0 or size <= 0 or offset + size > file_size:
        raise ValueError(f"bad blob span offset={offset} size={size}")
    fh.seek(offset)
    combined = _decompress(fh.read(size))
    mlen = int.from_bytes(combined[:4], "little")
    md = json.loads(combined[4:4 + mlen].decode("utf-8"))
    npz = np.load(io.BytesIO(combined[4 + mlen:]))

    def g(name):
        return npz[name] if name in npz.files else np.array([], dtype=np.float64)

    return md, g("frag_mz"), g("frag_intensity")


# imspy's Rust core PANICS + core-dumps (uncatchable) on a bad sequence, so we
# MUST validate before calling. Strict: only valid AA + [UNIMOD:N] may remain
# after conversion; anything else (leftover [+mass], hyphens, B/J/O/X/Z) is
# skipped, never passed to imspy.
_VALID_SEQ = re.compile(r"^([ARNDCQEGHILKMFPSTWYVU]|\[UNIMOD:\d+\])+$")


def _clean_modseq(s):
    """Strip empty terminal-mod placeholders (`[]-PEP-[]` -> `PEP`) from the
    STORED sequence fields; real `[UNIMOD:N]` mods untouched."""
    if not isinstance(s, str):
        return s
    return s.replace("[]-", "").replace("-[]", "").replace("[]", "")


def clean_imspy_seq(seq):
    if not seq:
        return None
    s = seq.replace("[]", "")                     # empty terminal placeholders
    s = s.replace("]-", "]").replace("-[", "[")   # terminal-mod hyphens
    s = s.strip("-")
    return s if _VALID_SEQ.match(s) else None


# --------------------------------------------------------------------------
# deterministic one-to-one matcher (Codex #3)
# --------------------------------------------------------------------------
def match_one_to_one(theoretical, exp_mz, exp_int, tol_ppm=TOL_PPM):
    """theoretical: {ion_type: [(ordinal, charge, mz, seq), ...]}.
    Returns list of matched-fragment dicts; each ion and each peak used once."""
    exp_mz = np.asarray(exp_mz, dtype=float)
    exp_int = np.asarray(exp_int, dtype=float)
    if exp_mz.size == 0:
        return []
    order = np.argsort(exp_mz)
    exp_mz, exp_int = exp_mz[order], exp_int[order]
    cands = []  # (abs_ppm, -intensity, ion_type, ordinal, charge, peak_idx, tmz)
    for ion_type, ions in theoretical.items():
        for ordinal, fch, tmz, _seq in ions:
            if tmz <= 0:
                continue
            lo = tmz * (1 - tol_ppm * 1e-6)
            hi = tmz * (1 + tol_ppm * 1e-6)
            j0 = int(np.searchsorted(exp_mz, lo, "left"))
            j1 = int(np.searchsorted(exp_mz, hi, "right"))
            for pj in range(j0, j1):
                ppm = abs(exp_mz[pj] - tmz) / tmz * 1e6
                cands.append((ppm, -float(exp_int[pj]), ion_type, int(ordinal),
                              int(fch), pj, tmz))
    cands.sort(key=lambda c: (c[0], c[1], c[2], c[3], c[4], c[5]))
    used_ion, used_peak, out = set(), set(), []
    for ppm, negI, ion_type, ordinal, fch, pj, tmz in cands:
        ik = (ion_type, ordinal, fch)
        if ik in used_ion or pj in used_peak:
            continue
        used_ion.add(ik)
        used_peak.add(pj)
        out.append({
            "fragment_type": ion_type, "fragment_ordinal": ordinal,
            "fragment_charge": fch, "fragment_mz_calculated": tmz,
            "fragment_mz_experimental": float(exp_mz[pj]),
            "fragment_intensity": float(exp_int[pj]),
            "ppm_error": (float(exp_mz[pj]) - tmz) / tmz * 1e6,
        })
    return out


# --------------------------------------------------------------------------
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _norm_raw(s):
    return s[:-2] if isinstance(s, str) and s.endswith(".d") else s


def _resolve_blob(extr: Path, raw_file: str):
    clean = raw_file.replace(".d", "")
    for p in (extr / f"{clean}.d" / "blobs.bin", extr / raw_file / "blobs.bin"):
        if p.exists():
            return p
    return None


def _read_concat(paths):
    tabs = [pq.read_table(p) for p in paths]
    return tabs[0] if len(tabs) == 1 else pa.concat_tables(tabs, promote_options="default")


def build(acc, data_root, out_dir, matcher, mods, q_max, fail_frac_max, limit=0):
    from fragment_matching import (sage_to_imspy_sequence,
                                   fragpipe_to_imspy_sequence)
    proc = data_root / "processed" / acc
    extr = data_root / "extracted" / acc
    def _pick(top, base, name):
        # prefer the single top-level merge; else per-group subdir copies
        # (per-group datasets have both -> concatenating double-counts).
        return [top] if top.exists() else sorted(base.rglob(name))
    pidx_hits = _pick(proc / "precursor_index.parquet", proc, "precursor_index.parquet")
    rf_hits = _pick(extr / "raw_features.parquet", extr, "raw_features.parquet")
    if not pidx_hits or not rf_hits:
        raise SystemExit(f"missing inputs for {acc}")
    pidx = _read_concat(pidx_hits).to_pylist()
    rf = _read_concat(rf_hits)
    rf_by_key = {(_norm_raw(r["raw_file"]), r["precursor_id"]): r for r in rf.to_pylist()}
    have_sage = "sage_qvalue" in (pidx[0] if pidx else {})

    # Sage-native fragment provenance (best-effort): (sage_psm_id, type, ord, ch) -> intensity
    native = {}
    sage_frag = proc / "sage" / "matched_fragments.sage.parquet"
    if sage_frag.exists():
        nf = pq.read_table(sage_frag).to_pylist()
        for r in nf:
            ft = r.get("fragment_type")
            ft = chr(ft) if isinstance(ft, int) else (ft.decode() if isinstance(ft, bytes) else ft)
            native[(r.get("psm_id"), ft, r.get("fragment_ordinals"),
                    r.get("fragment_charge"))] = r.get("fragment_intensity")

    by_raw = {}
    n_floor = 0
    for prow in pidx:
        key = (_norm_raw(prow["raw_file"]), prow["precursor_id"])
        rrow = rf_by_key.get(key)
        if rrow is None:
            continue
        sage_pass = (prow.get("sage_qvalue") is not None and prow["sage_qvalue"] <= q_max)
        fp_pass = (prow.get("fragpipe_qvalue") is not None and prow["fragpipe_qvalue"] <= q_max)
        if not (sage_pass or fp_pass):
            continue
        n_floor += 1
        by_raw.setdefault(rrow["raw_file"], []).append(
            (rrow.get("blob_offset"), rrow.get("blob_size"), prow, sage_pass, fp_pass))
        if limit and n_floor >= limit:
            break

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{acc}.parquet"
    writer = None
    fail = Counter()      # real failures -> quarantine
    skip = Counter()      # coverage limits (unparseable mods, no match) -> reported only
    n_rows = 0
    n_matched_prec = 0

    for raw_file in sorted(by_raw):
        items = sorted(by_raw[raw_file], key=lambda t: (t[0] if t[0] is not None else -1))
        blob_path = _resolve_blob(extr, raw_file)
        file_size = blob_path.stat().st_size if blob_path else 0
        fh = open(blob_path, "rb") if blob_path else None
        batch = []
        try:
            for offset, bsize, prow, sage_pass, fp_pass in items:
                if fh is None or offset is None or not bsize:
                    fail["no_blob"] += 1
                    continue
                try:
                    md, fmz, fint = read_frag(fh, offset, bsize, file_size)
                except Exception as e:
                    fail[type(e).__name__] += 1
                    continue
                if fmz.size == 0:
                    fail["empty_frag"] += 1
                    continue
                use_sage = sage_pass
                canon = "sage" if use_sage else "fragpipe"
                assignment = "consensus" if (sage_pass and fp_pass) else canon
                if use_sage:
                    raw_seq = prow.get("sage_modified") or prow.get("sage_peptide")
                    conv = sage_to_imspy_sequence(raw_seq) if raw_seq else None
                    charge = prow.get("sage_charge")
                else:
                    raw_seq = prow.get("fragpipe_modified") or prow.get("fragpipe_peptide")
                    conv = fragpipe_to_imspy_sequence(raw_seq) if raw_seq else None
                    charge = prow.get("fragpipe_charge")
                charge = int(charge) if charge is not None else prow.get("charge")
                imspy_seq = clean_imspy_seq(conv)   # None if not imspy-safe
                if not imspy_seq or charge is None:
                    skip["unparseable_seq"] += 1     # converter limit, not a failure
                    continue
                try:
                    theo = matcher.generate_theoretical_fragments(imspy_seq, int(charge))
                except Exception as e:
                    skip["theo:" + type(e).__name__] += 1
                    continue
                matches = match_one_to_one(theo, fmz, fint)
                if not matches:
                    skip["no_match"] += 1
                    continue
                n_matched_prec += 1
                sage_pid = prow.get("sage_psm_id")
                pep_conf = bool(sage_pass and fp_pass
                                and prow.get("sage_charge") != prow.get("fragpipe_charge"))
                n_eng = 1 + int(sage_pass and fp_pass and not pep_conf)
                for m in matches:
                    nat = native.get((sage_pid, m["fragment_type"],
                                      m["fragment_ordinal"], m["fragment_charge"]))
                    batch.append({
                        "accession": acc, "raw_file": raw_file,
                        "precursor_id": prow["precursor_id"], "sage_psm_id": sage_pid,
                        "assignment": assignment, "canonical_engine": canon,
                        "sequence": _clean_modseq(prow.get(canon + "_peptide")),
                        "modified_sequence": _clean_modseq(raw_seq),
                        "sequence_normalized": _clean_modseq(prow.get("sequence_normalized")),
                        "precursor_charge": int(charge),
                        "n_engines": n_eng, "peptidoform_conflict": pep_conf,
                        "match_method": MATCH_METHOD,
                        "sage_native_matched": nat is not None,
                        "sage_native_intensity": (float(nat) if nat is not None else None),
                        **m,
                    })
        finally:
            if fh is not None:
                fh.close()
        if batch:
            tbl = pa.Table.from_pylist(batch, schema=TIER3_SCHEMA)
            if writer is None:
                writer = pq.ParquetWriter(out_path, TIER3_SCHEMA, compression="zstd")
            writer.write_table(tbl)
            n_rows += len(batch)

    if writer is not None:
        writer.close()

    n_attempt = sum(len(v) for v in by_raw.values())
    fail_frac = (sum(fail.values()) / n_attempt) if n_attempt else 0.0
    manifest = {
        "accession": acc, "schema_version": SCHEMA_VERSION,
        "match_method": MATCH_METHOD, "tol_ppm": TOL_PPM, "q_max": q_max,
        "n_passed_floor": n_floor, "n_matched_precursors": n_matched_prec,
        "n_fragment_rows": n_rows, "failures": dict(fail), "skipped": dict(skip),
        "fail_frac": round(fail_frac, 5),
        "input_shas": {str(p): sha256(p) for p in (*pidx_hits, *rf_hits)},
        "built_unix": int(time.time()),
    }
    (out_dir / f"{acc}.manifest.json").write_text(json.dumps(manifest, indent=2))
    if fail_frac > fail_frac_max:
        raise SystemExit(f"QUARANTINE {acc}: fail_frac={fail_frac:.3f} > {fail_frac_max}")
    return manifest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True)
    ap.add_argument("--data-root", default="/lustre/project/ki-proanagi/dateschn/data")
    ap.add_argument("--out", default="hf_tier3")
    ap.add_argument("--repo-root", default="/lustre/project/ki-proanagi/dateschn/claudius-proteomics")
    ap.add_argument("--q-max", type=float, default=Q_MAX)
    ap.add_argument("--tol-ppm", type=float, default=TOL_PPM)
    ap.add_argument("--fail-frac-max", type=float, default=0.05)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    sys.path.insert(0, args.repo_root)
    sys.path.insert(0, str(Path(args.repo_root) / "scripts"))
    from fragment_matching import FragmentMatcher, MatchConfig
    matcher = FragmentMatcher(MatchConfig(mz_tolerance_ppm=args.tol_ppm,
                                          ion_types=["b", "y"], max_fragment_charge=2))
    m = build(args.accession, Path(args.data_root), Path(args.out), matcher, None,
              args.q_max, args.fail_frac_max, args.limit)
    print(json.dumps(m, indent=2))


if __name__ == "__main__":
    main()
