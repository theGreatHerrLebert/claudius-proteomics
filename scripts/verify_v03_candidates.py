#!/usr/bin/env python3
"""verify_v03_candidates.py — protocol-level verification of v0.3 wave-1 candidates.

TO_V03_EXPANSION.md §5 hard prerequisite (Codex catch): title-only classification
has false positives. For each acyl/HLA candidate, pull the PRIDE
sampleProcessingProtocol + dataProcessingProtocol and classify the acyl subtype /
confirm HLA from the PROTOCOL TEXT, not the title. Emit a manifest for human review.

Output: corpus_expansion_verified_<date>.tsv  (+ stdout summary)
"""
import csv
import json
import re
import sys
import time
import urllib.request

API = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s"

# protocol-text signatures (word-boundary, case-insensitive)
ACYL_SIG = {
    "lactyl":   r"lactyl|lactylat|\bKla\b|L-lactyl|lactylome",
    "succinyl": r"succinyl|succinylat|\bKsucc\b|desuccinyl|succinylome",
    "crotonyl": r"crotonyl|crotonylat|\bKcr\b|crotonylome",
    "malonyl":  r"malonyl|malonylat|\bKmal\b|malonylome",
}
HLA_SIG = r"immunopeptidom|\bMHC\b|\bHLA\b|immunoprecipitat.*(MHC|HLA)|peptidome|class I|class II|W6/32|pan-HLA"
# enrichment evidence (antibody / affinity) — distinguishes a real PTM-proteome from incidental mention
ENRICH_SIG = r"enrich|antibod|immunoaffinit|immunoprecipit|affinity bead|pan-anti|PTMScan|agarose"


def fetch(acc):
    try:
        with urllib.request.urlopen(API % acc, timeout=30) as r:
            return json.load(r)
    except Exception as e:
        return {"_error": str(e)}


def classify(acc, title, expected_wf):
    d = fetch(acc)
    if "_error" in d:
        return {"accession": acc, "verdict": "FETCH_FAIL", "note": d["_error"]}
    spp = (d.get("sampleProcessingProtocol") or "")
    dpp = (d.get("dataProcessingProtocol") or "")
    text = (title + " " + spp + " " + dpp)
    instr = ",".join(i.get("name", "") for i in d.get("instruments", []))
    blob = text.lower()
    enrich = bool(re.search(ENRICH_SIG, blob, re.I))

    if expected_wf == "acyl":
        hits = {k: len(re.findall(v, text, re.I)) for k, v in ACYL_SIG.items()}
        present = {k: n for k, n in hits.items() if n > 0}
        if not present:
            verdict, sub = "EXCLUDE_no_acyl", ""
        elif len(present) == 1:
            sub = next(iter(present))
            verdict = "OK" if enrich else "REVIEW_no_enrich"
        else:
            sub = max(present, key=present.get)
            verdict = "REVIEW_multi_acyl"
        note = "hits=%s enrich=%s" % (present, enrich)
        profile = "acyl_%s" % sub if sub else ""
    else:  # HLA
        sub = ""
        is_hla = bool(re.search(HLA_SIG, text, re.I))
        verdict = ("OK" if (is_hla and enrich) else
                   "REVIEW_hla_no_enrich" if is_hla else "EXCLUDE_not_hla")
        note = "hla=%s enrich=%s" % (is_hla, enrich)
        profile = "nonspecific_hla"
    return {"accession": acc, "workflow": expected_wf, "subtype": sub,
            "mod_profile": profile, "verdict": verdict, "instrument": instr,
            "organisms": ",".join(o.get("name", "") for o in d.get("organisms", [])),
            "note": note, "title": title[:80], "protocol_snip": spp[:200].replace("\n", " ")}


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "corpus_expansion_cc0_2026-06-12.tsv"
    rows = list(csv.DictReader(open(src), delimiter="\t"))
    wave1 = [r for r in rows if r["workflow"] in ("acyl", "HLA")]
    print("verifying %d wave-1 candidates (acyl+HLA)..." % len(wave1), flush=True)
    out = []
    for r in wave1:
        res = classify(r["accession"], r.get("title", ""), r["workflow"])
        out.append(res)
        print("  %-10s %-7s %-22s %s" % (res["accession"], res.get("subtype", ""),
              res["verdict"], res["note"]), flush=True)
        time.sleep(0.4)
    cols = ["accession", "workflow", "subtype", "mod_profile", "verdict",
            "instrument", "organisms", "note", "title", "protocol_snip"]
    dst = "corpus_expansion_verified_2026-06-13.tsv"
    with open(dst, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for res in out:
            w.writerow({c: res.get(c, "") for c in cols})
    # summary
    from collections import Counter
    print("\n=== verdict summary ===")
    for k, n in Counter(r["verdict"] for r in out).most_common():
        print("  %-22s %d" % (k, n))
    print("\n=== acyl subtype counts (OK + REVIEW) ===")
    for k, n in Counter(r.get("subtype", "") for r in out
                        if r["workflow"] == "acyl" and r["verdict"].startswith(("OK", "REVIEW"))).most_common():
        print("  %-10s %d" % (k or "(none)", n))
    print("\nwrote %s" % dst)


if __name__ == "__main__":
    main()
