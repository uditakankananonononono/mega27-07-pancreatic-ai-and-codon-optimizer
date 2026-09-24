#!/usr/bin/env python3
"""Fetch gnomAD r4 population variants for the PDAC driver panel
(KRAS/TP53/CDKN2A/SMAD4) via the gnomAD GraphQL API (gnomad.broadinstitute.org,
no auth). Raw per-gene JSON -> data/gnomad/ (untracked); pass-filter screen +
flatten -> results/gnomad_panel_variants.csv (committed).
Question: are the somatic hotspot alleles a blood-based PDAC panel would
target actually absent from the healthy population (specificity), and what is
the germline LoF carrier background for the same genes?"""
import json, os, csv, urllib.request

API = "https://gnomad.broadinstitute.org/api"
GENES = {"KRAS": "ENSG00000133703", "TP53": "ENSG00000141510",
         "CDKN2A": "ENSG00000147889", "SMAD4": "ENSG00000141646"}
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data", "gnomad")
Q = """query($g:String!){ gene(gene_id:$g, reference_genome:GRCh38){ symbol
variants(dataset:gnomad_r4){ variant_id rsids hgvsp consequence
genome{af ac an filters} exome{af ac an filters} } } }"""


def fetch(gene, ensg):
    os.makedirs(RAW, exist_ok=True)
    dest = os.path.join(RAW, f"{gene}.json")
    if not os.path.exists(dest):
        req = urllib.request.Request(API, data=json.dumps({"query": Q, "variables": {"g": ensg}}).encode(),
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=120) as r:
            json.dump(json.load(r), open(dest, "w"))
    return json.load(open(dest))["data"]["gene"]["variants"]


def flatten():
    rows = []
    for gene in GENES:
        for v in fetch(gene, GENES[gene]):
            ex, ge = v.get("exome") or {}, v.get("genome") or {}
            pass_ex = ex.get("filters") == [] and ex.get("ac") is not None
            pass_ge = ge.get("filters") == [] and ge.get("ac") is not None
            rows.append({"gene": gene, "variant_id": v["variant_id"],
                         "rsids": ";".join(v.get("rsids") or []),
                         "hgvsp": v.get("hgvsp") or "", "consequence": v.get("consequence"),
                         "af_exome": ex.get("af") if pass_ex else "", "ac_exome": ex.get("ac") if pass_ex else "",
                         "an_exome": ex.get("an", ""), "af_genome": ge.get("af") if pass_ge else "",
                         "ac_genome": ge.get("ac") if pass_ge else "", "an_genome": ge.get("an", "")})
    out = os.path.join(HERE, "..", "results", "gnomad_panel_variants.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    return len(rows)


if __name__ == "__main__":
    n = flatten()
    print("variants flattened:", n)
