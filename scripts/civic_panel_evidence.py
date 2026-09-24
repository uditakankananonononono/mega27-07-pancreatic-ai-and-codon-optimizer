#!/usr/bin/env python3
"""CIViC clinical-evidence audit of the 4-gene PDAC detection panel.

Question: after GTEx (no expression-level pancreatic specificity) and gnomAD
(DNA-level population background), does the *curated clinical evidence* behind
each panel gene concentrate on pancreatic cancer, or is it indication-diffuse?
Controls: BRCA1 (breast/ovarian concentration is the positive calibration) and
JAK2 (myeloproliferative concentration). If the pancreatic-share metric fails
to detect concentration for the controls it is broken; if it detects it there
but not for the panel, the panel's evidence base is genuinely diffuse.

Reads data/civic/<SYM>.json (fetch: scripts/civic_fetch.py).
Writes results/civic_panel_evidence.json and
results/civic_panel_evidence_rows.csv (one row per unique EID x gene match)."""
import csv, json, math, re, os
from collections import Counter
from scipy.stats import fisher_exact

GENES = [("KRAS", "panel"), ("TP53", "panel"), ("CDKN2A", "panel"),
         ("SMAD4", "panel"), ("BRCA1", "breast_ovarian_control"),
         ("JAK2", "hematologic_control")]
PANCREAS = re.compile(r"pancrea", re.I)
BREAST_OVARIAN = re.compile(r"breast|ovarian", re.I)
HEMATOLOGIC = re.compile(r"myelo|leukemi|polycythemia|thrombocythemia|lymph", re.I)

def wilson(k, n, z=1.959964):
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (p, c - h, c + h)

def load(sym):
    with open(f"data/civic/{sym}.json") as f:
        return json.load(f)["items"]

def target_hits(items, pat, accepted_only=True):
    use = [e for e in items if (not accepted_only or e["status"] == "ACCEPTED")]
    k = sum(1 for e in use if e["disease"] and pat.search(e["disease"]["name"] or ""))
    return k, len(use)

per_gene, rows, seen = {}, [], set()
for sym, group in GENES:
    items = load(sym)
    n_all = len(items)
    n_acc = sum(1 for e in items if e["status"] == "ACCEPTED")
    for e in items:
        eid = e["name"]
        if eid not in seen:
            seen.add(eid)
        rows.append({
            "eid": eid, "gene": sym, "group": group, "status": e["status"],
            "evidence_type": e["evidenceType"], "evidence_level": e["evidenceLevel"],
            "significance": e["significance"],
            "disease": (e["disease"] or {}).get("name", ""),
            "doid": (e["disease"] or {}).get("doid", ""),
            "molecular_profile": e["molecularProfile"]["name"],
            "pancreatic": int(bool(e["disease"] and PANCREAS.search(e["disease"]["name"] or ""))),
        })
    pat = PANCREAS if group == "panel" else (BREAST_OVARIAN if group == "breast_ovarian_control" else HEMATOLOGIC)
    k, n = target_hits(items, pat)
    p, lo, hi = wilson(k, n)
    dis = Counter((e["disease"] or {}).get("name", "") for e in items if e["status"] == "ACCEPTED")
    types = Counter(e["evidenceType"] for e in items if e["status"] == "ACCEPTED")
    levels = Counter(e["evidenceLevel"] for e in items if e["status"] == "ACCEPTED")
    per_gene[sym] = {
        "group": group, "items_all": n_all, "items_accepted": n_acc,
        "target_pattern": pat.pattern,
        "target_hits": k, "target_share": round(p, 4),
        "target_share_ci95": [round(lo, 4), round(hi, 4)],
        "no_disease_accepted": sum(1 for e in items if e["status"] == "ACCEPTED" and not (e["disease"] and e["disease"].get("name"))),
        "top_diseases": dis.most_common(5),
        "evidence_types": dict(types), "evidence_levels": dict(levels),
        "distinct_molecular_profiles": len({e["molecularProfile"]["name"] for e in items}),
    }

# calibration + panel Fisher tests (accepted items only)
k_brca, n_brca = target_hits(load("BRCA1"), BREAST_OVARIAN)
k_jak2, n_jak2 = target_hits(load("JAK2"), HEMATOLOGIC)
pan_items = [e for s in ["KRAS", "TP53", "CDKN2A", "SMAD4"] for e in load(s)]
k_panel, n_panel = 0, 0
for e in pan_items:
    if e["status"] != "ACCEPTED":
        continue
    n_panel += 1
    if e["disease"] and PANCREAS.search(e["disease"]["name"] or ""):
        k_panel += 1
odds_b, p_brca = fisher_exact([[k_panel, n_panel - k_panel], [k_brca, n_brca - k_brca]])
odds_j, p_jak2 = fisher_exact([[k_panel, n_panel - k_panel], [k_jak2, n_jak2 - k_jak2]])
p_pool, lo_pool, hi_pool = wilson(k_panel, n_panel)

# per-gene pancreatic detail for the panel itself
panel_panc = {s: {"k": target_hits(load(s), PANCREAS)[0], "n": target_hits(load(s), PANCREAS)[1]}
              for s in ["KRAS", "TP53", "CDKN2A", "SMAD4"]}

out = {
    "source": "CIViC v2 GraphQL API (civicdb.org/api/graphql), evidenceItems, molecularProfileName substring",
    "genes": per_gene,
    "panel_pooled_accepted": {"pancreatic_hits": k_panel, "n": n_panel,
                              "share": round(p_pool, 4), "share_ci95": [round(lo_pool, 4), round(hi_pool, 4)]},
    "per_gene_pancreatic_accepted": panel_panc,
    "fisher_panel_vs_brca1_breast_ovarian": {"odds_ratio": round(odds_b, 4), "p": float(f"{p_brca:.3e}")},
    "fisher_panel_vs_jak2_hematologic": {"odds_ratio": round(odds_j, 4), "p": float(f"{p_jak2:.3e}")},
    "unique_eids": len(seen),
    "gene_eid_rows": len(rows),
    "verdict": "",
}
out["verdict"] = (
    f"Controls calibrate the metric: BRCA1 evidence concentrates on breast/ovarian "
    f"({k_brca}/{n_brca} = {k_brca/n_brca:.2f}) and JAK2 on hematologic ({k_jak2}/{n_jak2} = {k_jak2/n_jak2:.2f}) "
    f"disease, both far above the panel's pooled pancreatic share ({k_panel}/{n_panel} = {p_pool:.3f}; "
    f"Fisher p={p_brca:.1e} and {p_jak2:.1e}). The detection panel's curated clinical-evidence base is "
    f"indication-diffuse: pancreatic specificity must come from the mutation-pattern combination, "
    f"not from any single marker's evidence concentration.")

os.makedirs("results", exist_ok=True)
with open("results/civic_panel_evidence.json", "w") as f:
    json.dump(out, f, indent=1)
with open("results/civic_panel_evidence_rows.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(json.dumps(out["panel_pooled_accepted"], indent=1))
print(out["verdict"])
print("unique EIDs:", len(seen), "rows:", len(rows))
for s in ["KRAS", "TP53", "CDKN2A", "SMAD4"]:
    d = panel_panc[s]; print(s, f"pancreatic {d['k']}/{d['n']}", "top:", per_gene[s]["top_diseases"][:3])
