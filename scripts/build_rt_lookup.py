import pyarrow.parquet as pq, pyarrow as pa, pyarrow.compute as pc, glob
R = "/lustre/project/ki-proanagi/dateschn/data"
accs = open("/lustre/project/ki-proanagi/dateschn/val15_accs.txt").read().strip().split(",")
tabs = []
for acc in accs:
    hits = glob.glob(f"{R}/processed/{acc}/**/results.sage.parquet", recursive=True)
    if not hits:
        print("MISSING", acc); continue
    for h in hits:
        sch = pq.read_schema(h).names
        cols = [c for c in ["psm_id", "aligned_rt", "spectrum_q", "rank", "is_decoy"] if c in sch]
        t = pq.read_table(h, columns=cols)
        if {"spectrum_q", "rank", "is_decoy"} <= set(cols):
            m = pc.and_(pc.and_(pc.equal(t["is_decoy"], False), pc.equal(t["rank"], 1)),
                        pc.less_equal(t["spectrum_q"], 0.01))
            t = t.filter(m)
        if t.num_rows == 0:
            continue
        t = t.select(["psm_id", "aligned_rt"]).append_column(
            "accession", pa.array([acc] * t.num_rows, type=pa.string()))
        tabs.append(t)
out = pa.concat_tables(tabs)
pq.write_table(out, f"{R}/rt_aligned_lookup_val15.parquet", compression="zstd")
print("lookup rows:", out.num_rows, "datasets:", len(set(out['accession'].to_pylist())))
print("aligned_rt range:", round(pc.min(out['aligned_rt']).as_py(), 3), round(pc.max(out['aligned_rt']).as_py(), 3))
