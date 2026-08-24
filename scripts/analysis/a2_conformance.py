#!/usr/bin/env python3
"""
A2 configuration-conformance sweep (ECCB poster, figure 2).

For every rendered engine config under a processed-data root, compare the
*rendered* search configuration against the mod_profile that config.mogon.yaml
resolves for that accession. Reports FragPipe and Sage separately.

Rendered artifacts:
  FragPipe: <root>/<ACC>[/<group>]/fragpipe_output/fragpipe.workflow
            (variable mods in `msfragger.table.var-mods` = "mass,residues,enabled,max; ...")
  Sage:     <root>/<ACC>[/<group>]/sage/sage_config.json
            (`database.variable_mods` = {residue: [masses]}, `database.enzyme`)

Stale run dirs (`.bak*`, `.failed*`) are excluded. PTM matching is by mass
(±0.02) and skips ubiquitous Oxidation/M; a known limitation is that mass-only
matching conflates isobaric N-term vs side-chain sites (e.g. Acetyl 42.0106).

Usage (on MOGON):
  python a2_conformance.py \
      --config config/config.mogon.yaml \
      --root   /lustre/project/ki-proanagi/dateschn/data/processed
Or point --fragpipe-root / --sage-root at separate trees.
"""
import argparse, glob, json, re, collections
import yaml

OX = 15.9949  # Oxidation/M — ubiquitous, excluded from PTM matching


def load_expected(config_path):
    cfg = yaml.safe_load(open(config_path))
    dm, mp = cfg["dataset_metadata"], cfg["mod_profiles"]

    def expected_for(acc):
        meta = dm.get(acc)
        prof = meta.get("mod_profile") if meta else None
        p = mp.get(prof) if prof else None
        enzyme = (p or {}).get("enzyme_override") or "trypsin"
        return dict(
            profile=prof,
            enzyme=enzyme,
            termini=0 if enzyme == "nonspecific" else 2,
            dmin=(p or {}).get("min_peptide_length", 7),
            dmax=(p or {}).get("max_peptide_length", 50),
            varmods=[(round(m["mass"], 4), tuple(m.get("residues", [m.get("site", "?")])))
                     for m in (p or {}).get("variable_modifications", [])],
        )
    return expected_for


def acc_of(path, marker):
    accdir = path.split(marker, 1)[1].split("/", 1)[0]
    if ".bak" in accdir or ".failed" in accdir:
        return None
    m = re.match(r"(PXD\d+)", accdir)
    return m.group(1) if m else None


def parse_fragpipe(path):
    kv = {}
    for line in open(path, errors="ignore"):
        line = line.strip()
        if "=" in line and not line.startswith("#"):
            k, v = line.split("=", 1)
            kv[k.strip()] = v.strip()
    enabled = set()
    for entry in kv.get("msfragger.table.var-mods", "").split(";"):
        parts = [x.strip() for x in entry.split(",")]
        if len(parts) >= 3 and parts[2].lower() == "true":
            try:
                enabled.add(round(float(parts[0]), 4))
            except ValueError:
                pass
    return dict(enzyme=kv.get("msfragger.search_enzyme_name_1", "?"),
                termini=kv.get("msfragger.num_enzyme_termini", "?"),
                dmin=kv.get("msfragger.digest_min_length", "?"),
                dmax=kv.get("msfragger.digest_max_length", "?"),
                enabled_masses=enabled)


def parse_sage(path):
    sc = json.load(open(path))
    db = sc.get("database", {})
    enz = db.get("enzyme", {})
    cleave = enz.get("cleave_at", "")
    masses = set()
    for _, ms in db.get("variable_mods", {}).items():
        for mm in (ms if isinstance(ms, list) else [ms]):
            try:
                masses.add(round(float(mm), 4))
            except (TypeError, ValueError):
                pass
    return dict(enzyme="nonspecific" if cleave in ("", None) else "trypsin",
                termini=0 if cleave in ("", None) else 2,
                dmin=enz.get("min_len", "?"), dmax=enz.get("max_len", "?"),
                enabled_masses=masses)


def issues(rendered, exp):
    out = []
    if rendered["enzyme"] != exp["enzyme"]:
        out.append(f"enzyme:{rendered['enzyme']}!={exp['enzyme']}")
    if str(rendered["termini"]) != str(exp["termini"]):
        out.append(f"termini:{rendered['termini']}!={exp['termini']}")
    if str(rendered["dmin"]) != str(exp["dmin"]) or str(rendered["dmax"]) != str(exp["dmax"]):
        out.append(f"digest:{rendered['dmin']}-{rendered['dmax']}!={exp['dmin']}-{exp['dmax']}")
    for mass, res in exp["varmods"]:
        if abs(mass - OX) < 0.01:
            continue
        if not any(abs(mass - em) < 0.02 for em in rendered["enabled_masses"]):
            out.append(f"PTM:{mass}/{','.join(res)}")
    return out


def sweep(name, files, marker, parser, expected_for):
    byacc = collections.defaultdict(list)
    for f in files:
        acc = acc_of(f, marker)
        if acc:
            byacc[acc].append(f)
    rows, nonconf = [], []
    for acc in sorted(byacc):
        exp = expected_for(acc)
        for f in byacc[acc]:
            try:
                iss = issues(parser(f), exp)
            except Exception as e:  # noqa: BLE001
                iss = [f"parse-error:{e}"]
            rows.append((acc, exp["profile"], bool(iss)))
            if iss:
                nonconf.append((acc, exp["profile"], "; ".join(iss)))
    naccs = len(byacc)
    ok = sum(1 for *_, bad in rows if not bad)
    bad_accs = {a for a, _, _ in nonconf}
    print(f"\n=== {name} ===")
    print(f"datasets={naccs}  configs={len(rows)}")
    print(f"conforming configs:  {ok}/{len(rows)} = {100*ok/len(rows):.1f}%")
    print(f"conforming datasets: {naccs-len(bad_accs)}/{naccs} = {100*(naccs-len(bad_accs))/naccs:.1f}%")
    ptot = collections.Counter(expected_for(a)["profile"] for a in byacc)
    pbad = collections.Counter(expected_for(a)["profile"] for a in bad_accs)
    print("non-conformance by profile (datasets):")
    for prof, tot in sorted(ptot.items(), key=lambda x: str(x[0])):
        print(f"  {str(prof):16s} {pbad.get(prof,0):>2}/{tot:>2}")
    for acc, prof, why in nonconf:
        print(f"  - {acc} [{prof or 'none'}]: {why}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.mogon.yaml")
    ap.add_argument("--root", default="/lustre/project/ki-proanagi/dateschn/data/processed")
    ap.add_argument("--fragpipe-root", default=None)
    ap.add_argument("--sage-root", default=None)
    a = ap.parse_args()
    expected_for = load_expected(a.config)
    fp_root = a.fragpipe_root or a.root
    sg_root = a.sage_root or a.root
    fp = glob.glob(f"{fp_root.rstrip('/')}/**/fragpipe.workflow", recursive=True)
    sg = glob.glob(f"{sg_root.rstrip('/')}/**/sage_config.json", recursive=True)
    fp_marker = fp_root.rstrip("/").split("/")[-1] + "/"
    sg_marker = sg_root.rstrip("/").split("/")[-1] + "/"
    sweep("FragPipe conformance", fp, fp_marker, parse_fragpipe, expected_for)
    sweep("Sage conformance", sg, sg_marker, parse_sage, expected_for)


if __name__ == "__main__":
    main()
