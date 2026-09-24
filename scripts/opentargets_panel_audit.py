#!/usr/bin/env python3
"""Open Targets association + tractability audit of the 4-gene PDAC panel.

Question: gnomAD (germline), GTEx (expression) and CIViC (curated clinical
evidence) all found no pancreatic specificity for KRAS/TP53/CDKN2A/SMAD4. Does
Open Targets' *systematic* cross-datasource aggregation agree, or does the
genetic evidence still concentrate the panel on pancreatic cancer? Controls:
BRCA1 must top breast cancer (MONDO_0007254) and JAK2 must top
myeloproliferative disorder (EFO_0004251); both must rank low for exocrine
pancreatic carcinoma (MONDO_0005192). If the controls fail, the rank metric is
broken and any panel result is void.

Reads data/opentargets/assoc_*.json + target_*.json (fetch:
scripts/opentargets_fetch.py). Writes results/opentargets_panel_audit.json
and results/opentargets_assoc_rows.csv (one row per target x disease
association record)."""
import csv, json, math
from collections import Counter

GENETIC = {"intogen", "cancer_gene_census", "eva", "eva_somatic",
           "genomics_england", "cancer_biomarkers", "crispr",
           "uniprot_variants", "gene_burden", "impc", "orphanet",
           "gene2phenotype", "clingen", "cancer_gene_census_somatic"}
LIT = "europepmc"
GENES = [("KRAS", "panel"), ("TP53", "panel"), ("CDKN2A", "panel"),
         ("SMAD4", "panel"), ("BRCA1", "breast_control"),
         ("JAK2", "mpn_control")]
DISEASES = {"MONDO_0005192": "exocrine pancreatic carcinoma",
            "MONDO_0007254": "breast cancer",
            "EFO_0004251": "myeloproliferative disorder"}

def wilson(k, n, z=1.959964):
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (p, c - h, c + h)

def load_assoc(dis):
    with open(f"data/opentargets/assoc_{dis}.json") as f:
        return json.load(f)

tables, ds_union = {}, set()
for dis in DISEASES:
    t = load_assoc(dis)
    assert t["count"] == len(t["rows"])
    rows = []
    for r in t["rows"]:
        ds = {d["id"]: d["score"] for d in r["datasourceScores"]}
        ds_union |= set(ds)
        gmax = max([ds.get(k, 0.0) for k in GENETIC] + [0.0])
        rows.append({"ensg": r["target"]["id"],
                     "symbol": r["target"]["approvedSymbol"],
                     "disease": dis, "score": r["score"],
                     "genetic_max": gmax, "europepmc": ds.get(LIT, 0.0),
                     "ds": ds})
    rows.sort(key=lambda r: -r["score"])
    for i, r in enumerate(rows):
        r["rank"] = i + 1
    tables[dis] = rows

# association rows CSV (one row per target x disease record)
with open("results/opentargets_assoc_rows.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["ensg", "symbol", "disease", "rank", "score",
                "genetic_max", "europepmc", "lit_only"])
    for dis, rows in tables.items():
        for r in rows:
            w.writerow([r["ensg"], r["symbol"], dis, r["rank"],
                        round(r["score"], 6), round(r["genetic_max"], 6),
                        round(r["europepmc"], 6),
                        int(r["genetic_max"] < 0.25 and r["europepmc"] >= 0.25)])

def rank_of(dis, sym):
    for r in tables[dis]:
        if r["symbol"] == sym:
            return r
    return None

per_gene = {}
for sym, group in GENES:
    entry = {"group": group, "diseases": {}}
    for dis in DISEASES:
        r = rank_of(dis, sym)
        if r is None:
            entry["diseases"][dis] = {"associated": False}
            continue
        n = len(tables[dis])
        entry["diseases"][dis] = {
            "associated": True, "rank": r["rank"], "n_targets": n,
            "score": round(r["score"], 6),
            "top_percentile": round(100 * r["rank"] / n, 4),
            "genetic_max": round(r["genetic_max"], 6),
            "europepmc": round(r["europepmc"], 6),
            "datasources": {k: round(v, 6) for k, v in
                            sorted(r["ds"].items(), key=lambda kv: -kv[1])}}
    with open(f"data/opentargets/target_{sym}.json") as f:
        t = json.load(f)
    cand = t["drugAndClinicalCandidates"]
    stages = Counter(row["maxClinicalStage"] for row in cand["rows"])
    panc_drugs = sorted({row["drug"]["name"] for row in cand["rows"]
                         if any(d.get("disease") and
                                "pancrea" in (d["disease"]["name"] or "").lower()
                                for d in row["diseases"])})
    tract = {}
    for b in t["tractability"]:
        if b["value"]:
            tract.setdefault(b["modality"], []).append(b["label"])
    entry["drug_candidates"] = {"count": cand["count"],
                                "by_stage": dict(stages),
                                "pancreatic_indication_drugs": panc_drugs}
    entry["tractability_true_buckets"] = tract
    per_gene[sym] = entry

# background literature-only fractions (pancreatic table)
prows = tables["MONDO_0005192"]
lit_only = [r for r in prows
            if r["genetic_max"] < 0.25 and r["europepmc"] >= 0.25]
genetic_backed = [r for r in prows if r["genetic_max"] >= 0.25]
top100 = prows[:100]
top100_lit_only = sum(1 for r in top100
                      if r["genetic_max"] < 0.25 and r["europepmc"] >= 0.25)
rest = prows[100:]
rest_lit_only = sum(1 for r in rest
                    if r["genetic_max"] < 0.25 and r["europepmc"] >= 0.25)
from scipy.stats import fisher_exact
fe = fisher_exact([[top100_lit_only, 100 - top100_lit_only],
                   [rest_lit_only, len(rest) - rest_lit_only]])

panel_ranks_panc = [per_gene[s]["diseases"]["MONDO_0005192"]["rank"]
                    for s, g in GENES if g == "panel"]
top_frac = max(panel_ranks_panc) / len(prows)
binom_all4 = top_frac ** 4  # chance all 4 panel genes land this high at random

out = {
    "source": "Open Targets Platform GraphQL API (api.platform.opentargets.org/api/v4/graphql), fetched 2026-09-25",
    "diseases": {d: {"name": DISEASES[d], "n_targets": len(tables[d])}
                 for d in DISEASES},
    "datasources_seen": sorted(ds_union),
    "per_gene": per_gene,
    "pancreatic_background": {
        "n": len(prows),
        "lit_only_n": len(lit_only),
        "lit_only_share": round(len(lit_only) / len(prows), 6),
        "lit_only_wilson95": [round(x, 6) for x in
                              wilson(len(lit_only), len(prows))],
        "genetic_backed_n": len(genetic_backed),
        "genetic_backed_share": round(len(genetic_backed) / len(prows), 6),
        "top100_lit_only": top100_lit_only,
        "rest_lit_only": rest_lit_only,
        "rest_n": len(rest),
        "fisher_top100_vs_rest_litonly_p": fe.pvalue,
        "fisher_oddsratio": fe.statistic},
    "panel_all4_top_rank_probability": binom_all4,
    "panel_ranks_pancreatic": panel_ranks_panc,
}
with open("results/opentargets_panel_audit.json", "w") as f:
    json.dump(out, f, indent=1)

print("panel pancreatic ranks:", panel_ranks_panc,
      "of", len(prows), "chance p<", f"{binom_all4:.2e}")
for sym, g in GENES:
    d = per_gene[sym]["diseases"]
    print(sym, g, {k: (v.get("rank"), v.get("top_percentile"))
                   for k, v in d.items()})
print("lit-only bg:", len(lit_only), "/", len(prows),
      "top100:", top100_lit_only, "rest:", rest_lit_only,
      "fisher p:", f"{fe.pvalue:.2e}")
for sym, g in GENES:
    print(sym, "drugs:", per_gene[sym]["drug_candidates"]["count"],
          per_gene[sym]["drug_candidates"]["by_stage"],
          "panc:", per_gene[sym]["drug_candidates"]["pancreatic_indication_drugs"])
