#!/usr/bin/env python3
"""
A9 — provenance supplement generator (ECCB poster figure-2 reproducibility chain).

Emits one row per (accession, sample-group) linking, for each search engine:
  resolved mod_profile  ->  rendered config (hashed)  ->  engine versions  ->
  run (SLURM job + UTC)  ->  code (git commit + dirty flag)
plus the A2 conformance verdict, so the supplement *is* the figure-2 data table
with its provenance attached.

Sources (cluster layout; roots overridable for local reproduction):
  rendered:   <processed>/<ACC>[/<group>]/fragpipe_output/fragpipe.workflow
              <processed>/<ACC>[/<group>]/sage/sage_config.json
  provenance: <provenance>/<ACC>/run-*.json  (git_commit, working_tree_dirty,
              config sha256, sage version, job id, started_at_utc)
  profile:    config.mogon.yaml (dataset_metadata + mod_profiles)

FragPipe/MSFragger versions come from the workflow header ("# FragPipe version
24.0" / "# MSFragger version 4.4"); Sage version from the run json. Config sha is
compared against the current config to flag runs whose config has since changed
(a completeness signal, not an error). `.bak*`/`.failed*` dirs are excluded.

Honest limits recorded per row: working_tree_dirty (code not bit-pinned) and
config_matches_current (whether the current config still resolves the same
profile). The literal command line is not in the run json and is not claimed.

Usage (MOGON):
  python a9_provenance_supplement.py --config config/config.mogon.yaml \
      --processed <data>/processed --provenance <data>/provenance \
      --out A9_provenance_supplement.tsv
Local reproduction (separate rendered trees):
  ... --fp-root a2_workflows --sage-root sage_configs --provenance provenance
"""
import argparse, glob, json, re, hashlib, os, collections
import yaml

OX = 15.9949


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def load_expected(config_path):
    cfg = yaml.safe_load(open(config_path))
    dm, mp = cfg["dataset_metadata"], cfg["mod_profiles"]

    def expected_for(acc):
        meta = dm.get(acc)
        prof = meta.get("mod_profile") if meta else None
        p = mp.get(prof) if prof else None
        enzyme = (p or {}).get("enzyme_override") or "trypsin"
        return dict(profile=prof or "(none)", enzyme=enzyme,
                    termini=0 if enzyme == "nonspecific" else 2,
                    dmin=(p or {}).get("min_peptide_length", 7),
                    dmax=(p or {}).get("max_peptide_length", 50),
                    varmods=[(round(m["mass"], 4), tuple(m.get("residues", [m.get("site", "?")])))
                             for m in (p or {}).get("variable_modifications", [])])
    return expected_for


_STALE = re.compile(r"\.(?:bak|failed|oom|pre)(?:[-_.]|$)")  # delimiter boundary (#2): avoid over-matching e.g. .preprint


def acc_of(path, marker):
    rel = path.split(marker, 1)[1]
    # #5 fix: reject any stale path component (group-level .pre-resage etc.),
    # not just the accession-level dir.
    if any(_STALE.search(c) for c in rel.split("/")):
        return None, None
    a = rel.split("/", 1)[0]
    m = re.match(r"(PXD\d+)", a)
    return (m.group(1), a) if m else (None, None)


def group_of(path, marker):
    parts = path.split(marker, 1)[1].split("/")
    return parts[1] if parts[2] in ("fragpipe_output", "sage") else "(single)"


def parse_fp(path):
    kv, fp_ver, ms_ver = {}, "", ""
    for line in open(path, errors="ignore"):
        s = line.strip()
        m = re.match(r"#\s*FragPipe version\s+(\S+)", s)
        if m:
            fp_ver = m.group(1)
        m = re.match(r"#\s*MSFragger version\s+(\S+)", s)
        if m:
            ms_ver = m.group(1)
        if "=" in s and not s.startswith("#"):
            k, v = s.split("=", 1)
            kv[k.strip()] = v.strip()
    enabled = set()
    for e in kv.get("msfragger.table.var-mods", "").split(";"):
        p = [x.strip() for x in e.split(",")]
        if len(p) >= 3 and p[2].lower() == "true":
            try:
                enabled.add(round(float(p[0]), 4))
            except ValueError:
                pass
    return dict(fragpipe_ver=fp_ver, msfragger_ver=ms_ver,
                enzyme=kv.get("msfragger.search_enzyme_name_1", "?"),
                termini=kv.get("msfragger.num_enzyme_termini", "?"),
                dmin=kv.get("msfragger.digest_min_length", "?"),
                dmax=kv.get("msfragger.digest_max_length", "?"),
                enabled_masses=enabled)


def parse_sg(path):
    db = json.load(open(path)).get("database", {})
    enz = db.get("enzyme", {})
    cl = enz.get("cleave_at", "")
    ms = set()
    for _, mm in db.get("variable_mods", {}).items():
        for x in (mm if isinstance(mm, list) else [mm]):
            try:
                ms.add(round(float(x), 4))
            except (TypeError, ValueError):
                pass
    return dict(enzyme="nonspecific" if cl in ("", None) else "trypsin",
                termini=0 if cl in ("", None) else 2,
                dmin=enz.get("min_len", "?"), dmax=enz.get("max_len", "?"),
                enabled_masses=ms)


def conforms(r, exp):
    if r["enzyme"] != exp["enzyme"]:
        return False
    if str(r["termini"]) != str(exp["termini"]):
        return False
    if str(r["dmin"]) != str(exp["dmin"]) or str(r["dmax"]) != str(exp["dmax"]):
        return False
    for mass, _ in exp["varmods"]:
        if abs(mass - OX) < 0.01:
            continue
        if not any(abs(mass - em) < 0.02 for em in r["enabled_masses"]):
            return False
    return True


def load_provenance(prov_root):
    """Latest run-json per accession."""
    by = {}
    for f in glob.glob(f"{prov_root.rstrip('/')}/**/run-*.json", recursive=True):
        try:
            j = json.load(open(f))
        except Exception:
            continue
        acc = j.get("accession")
        if not acc:
            continue
        ts = j.get("started_at_utc", "")
        if acc not in by or ts > by[acc][0]:
            by[acc] = (ts, j)
    out = {}
    for acc, (_, j) in by.items():
        out[acc] = dict(
            git_commit=(j.get("code", {}).get("git_commit", "") or "")[:12],
            dirty=j.get("code", {}).get("working_tree_dirty"),
            run_config_sha=j.get("config", {}).get("sha256", ""),
            sage_ver=j.get("engines_versions", {}).get("sage", ""),
            job_id=j.get("slurm", {}).get("job_id", ""),
            started=j.get("started_at_utc", ""))
    return out


def collect(files, marker):
    by = collections.defaultdict(dict)  # acc -> {group -> path}
    for f in files:
        acc, _ = acc_of(f, marker)
        if acc:
            by[acc][group_of(f, marker)] = f
    return by


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.mogon.yaml")
    ap.add_argument("--processed", default="/lustre/project/ki-proanagi/dateschn/data/processed")
    ap.add_argument("--provenance", default="/lustre/project/ki-proanagi/dateschn/data/provenance")
    ap.add_argument("--fp-root", default=None)
    ap.add_argument("--sage-root", default=None)
    ap.add_argument("--out", default="A9_provenance_supplement.tsv")
    a = ap.parse_args()
    expected_for = load_expected(a.config)
    cur_cfg_sha = sha256_file(a.config)
    fp_root = a.fp_root or a.processed
    sg_root = a.sage_root or a.processed
    fp_marker = fp_root.rstrip("/").split("/")[-1] + "/"
    sg_marker = sg_root.rstrip("/").split("/")[-1] + "/"
    fp = collect(glob.glob(f"{fp_root.rstrip('/')}/**/fragpipe.workflow", recursive=True), fp_marker)
    sg = collect(glob.glob(f"{sg_root.rstrip('/')}/**/sage_config.json", recursive=True), sg_marker)
    prov = load_provenance(a.provenance)

    cols = ["accession", "group", "profile",
            "fragpipe_ver", "msfragger_ver", "fp_workflow_sha256", "fp_conforms",
            "sage_ver", "sage_config_sha256", "sage_conforms",
            "git_commit", "working_tree_dirty", "run_config_sha256",
            "config_matches_current", "job_id", "started_utc"]
    rows = []
    accs = sorted(set(fp) | set(sg))
    for acc in accs:
        exp = expected_for(acc)
        pv = prov.get(acc, {})
        groups = sorted(set(fp.get(acc, {})) | set(sg.get(acc, {})))
        for g in groups:
            row = {c: "" for c in cols}
            row.update(accession=acc, group=g, profile=exp["profile"])
            fpf = fp.get(acc, {}).get(g)
            if fpf:
                r = parse_fp(fpf)
                row.update(fragpipe_ver=r["fragpipe_ver"], msfragger_ver=r["msfragger_ver"],
                           fp_workflow_sha256=sha256_file(fpf)[:16],
                           fp_conforms=str(conforms(r, exp)))
            sgf = sg.get(acc, {}).get(g)
            if sgf:
                r = parse_sg(sgf)
                row.update(sage_config_sha256=sha256_file(sgf)[:16],
                           sage_conforms=str(conforms(r, exp)))
            row.update(sage_ver=pv.get("sage_ver", ""), git_commit=pv.get("git_commit", ""),
                       working_tree_dirty=str(pv.get("dirty", "")),
                       run_config_sha256=(pv.get("run_config_sha", "") or "")[:16],
                       config_matches_current=str(pv.get("run_config_sha", "") == cur_cfg_sha) if pv else "no_prov",
                       job_id=pv.get("job_id", ""), started_utc=pv.get("started", ""))
            rows.append(row)
    with open(a.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in cols) + "\n")

    # coverage / integrity summary
    n = len(rows)
    have_prov = sum(1 for r in rows if r["job_id"])
    dirty = sum(1 for r in rows if r["working_tree_dirty"] == "True")
    cfg_match = sum(1 for r in rows if r["config_matches_current"] == "True")
    fp_ver_seen = collections.Counter(r["fragpipe_ver"] for r in rows if r["fragpipe_ver"])
    ms_ver_seen = collections.Counter(r["msfragger_ver"] for r in rows if r["msfragger_ver"])
    print(f"A9 supplement: {n} (accession,group) rows -> {a.out}")
    print(f"  with provenance/run: {have_prov}/{n}")
    print(f"  working_tree_dirty=True: {dirty}/{n}  <-- code-pinning gap (the honest caveat)")
    print(f"  run config == current config: {cfg_match}/{n}")
    print(f"  FragPipe versions: {dict(fp_ver_seen)}   MSFragger: {dict(ms_ver_seen)}")
    print(f"  Sage versions: {dict(collections.Counter(r['sage_ver'] for r in rows if r['sage_ver']))}")


if __name__ == "__main__":
    main()
