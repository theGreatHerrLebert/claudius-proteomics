#!/usr/bin/env python3
"""aggregate_hf_corpus.py — assemble the HF corpus from per-dataset parquets.

TO_HF_CORPUS.md v5 §5/§6/§7. For each tier: cast every per-dataset parquet to
the canonical Arrow schema, add a peptide-hash `split` column, and stream rows
into per-split parquets. Emits a top-level manifest + corpus-wide build stats
(read from the per-dataset manifests).

Split (v0.1): peptide-hash on `sequence_normalized` -> 80/10/10, matching the
project's `peptide_split` convention. Guarantees NO peptide leakage across
splits. Dataset/group-level (sibling-PXD) splitting is deferred to v0.2 — it
needs group_id metadata we don't yet have; documented in the manifest.

Run on a COMPUTE node (heavy lustre I/O), not the login node.
"""
import glob
import hashlib
import json
import sys
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

BASE = Path("/lustre/project/ki-proanagi/dateschn")
sys.path.insert(0, str(BASE))
from build_hf_tier1 import TIER1_SCHEMA           # noqa: E402
from build_hf_tier3 import TIER3_SCHEMA           # noqa: E402

CORPUS_VERSION = "v0.1"
SPLIT_SEED = 0
VAL_FRAC, TEST_FRAC = 0.10, 0.10


def split_for(seq):
    if not seq:
        return "train"
    key = ("%d|%s" % (SPLIT_SEED, seq)).encode()
    h = (int(hashlib.sha1(key).hexdigest(), 16) % 10_000) / 10_000.0
    if h < TEST_FRAC:
        return "test"
    if h < TEST_FRAC + VAL_FRAC:
        return "val"
    return "train"


def _out_schema(base):
    return pa.schema(list(base) + [pa.field("split", pa.string())])


def aggregate(name, schema, in_glob, out_dir, batch_rows=200_000):
    """Stream each per-dataset parquet in row-group batches (memory-bounded,
    regardless of dataset size — a whole-table read OOM-killed the first run)."""
    osch = _out_schema(schema)
    names = [f.name for f in schema]
    out_dir.mkdir(parents=True, exist_ok=True)
    writers, counts, per_ds = {}, {"train": 0, "val": 0, "test": 0}, {}
    files = sorted(glob.glob(in_glob))
    for fp in files:
        acc = Path(fp).stem
        rows = 0
        for rb in pq.ParquetFile(fp).iter_batches(batch_size=batch_rows):
            t = pa.Table.from_batches([rb]).select(names).cast(schema)
            seqs = t.column("sequence_normalized").to_pylist()
            splits = [split_for(s) for s in seqs]
            t = t.append_column("split", pa.array(splits, pa.string()))
            rows += t.num_rows
            for sp in ("train", "val", "test"):
                sub = t.filter(pa.array([x == sp for x in splits]))
                if sub.num_rows == 0:
                    continue
                w = writers.get(sp)
                if w is None:
                    w = pq.ParquetWriter(out_dir / ("%s.parquet" % sp), osch, compression="zstd")
                    writers[sp] = w
                w.write_table(sub)
                counts[sp] += sub.num_rows
        per_ds[acc] = rows
    for w in writers.values():
        w.close()
    print("[%s] files=%d rows: %s" % (name, len(files), counts), flush=True)
    return {"n_datasets": len(files), "split_rows": counts,
            "total_rows": sum(counts.values()), "per_dataset_rows": per_ds}


def build_stats(man_glob, keys):
    out = {"n_manifests": 0, "quarantined": []}
    agg = {k: 0 for k in keys}
    skips, fails = {}, {}
    out["bad_manifests"] = []
    for m in sorted(glob.glob(man_glob)):
        try:
            d = json.load(open(m))
        except Exception:
            out["bad_manifests"].append(Path(m).name)   # empty/malformed -> skip, don't crash
            continue
        out["n_manifests"] += 1
        for k in keys:
            agg[k] += d.get(k, 0) or 0
        if d.get("fail_frac", 0) > 0.05 or d.get("blob_fail_frac", 0) > 0.02:
            out["quarantined"].append(d["accession"])
        for k, v in (d.get("skipped", {}) or {}).items():
            skips[k] = skips.get(k, 0) + v
        for k, v in (d.get("blob_failures", {}) or d.get("failures", {}) or {}).items():
            fails[k] = fails.get(k, 0) + v
    out.update(agg)
    out["skipped"] = skips
    out["failures"] = fails
    return out


def main():
    out_root = BASE / "hf_corpus" / CORPUS_VERSION
    t1 = aggregate("tier1", TIER1_SCHEMA, str(BASE / "hf_tier1" / "*.parquet"),
                   out_root / "tier1")
    t3 = aggregate("tier3", TIER3_SCHEMA, str(BASE / "hf_tier3" / "*.parquet"),
                   out_root / "tier3")
    s1 = build_stats(str(BASE / "hf_tier1" / "*.manifest.json"),
                     ["n_passed_floor", "n_written"])
    s3 = build_stats(str(BASE / "hf_tier3" / "*.manifest.json"),
                     ["n_passed_floor", "n_matched_precursors", "n_fragment_rows"])
    manifest = {
        "corpus_version": CORPUS_VERSION, "built_unix": int(time.time()),
        "split": {"method": "peptide_hash(sequence_normalized)",
                  "seed": SPLIT_SEED, "fracs": {"train": 0.8, "val": VAL_FRAC, "test": TEST_FRAC},
                  "note": "no peptide leakage; group/dataset-level split deferred to v0.2"},
        "tier1": t1, "tier3": t3,
        "tier1_build_stats": s1, "tier3_build_stats": s3,
        "filter": {"q_max": 0.01, "rank": 1, "target_only": True,
                   "floor_note": "union of per-engine reported q<=0.01; NOT a corpus-level FDR"},
        "schema_versions": {"tier1": "tier1.v5", "tier3": "tier3.v1"},
    }
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
