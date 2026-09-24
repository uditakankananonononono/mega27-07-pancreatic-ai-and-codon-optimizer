#!/usr/bin/env python3
"""Fetch Open Targets Platform association tables + per-target tractability.

Tables (full, all pages): exocrine pancreatic carcinoma (MONDO_0005192),
breast cancer (MONDO_0007254), myeloproliferative disorder (EFO_0004251).
Per-target: tractability + drugAndClinicalCandidates for the 4 panel genes
(KRAS, TP53, CDKN2A, SMAD4) and the 2 positive controls (BRCA1, JAK2).
Writes data/opentargets/*.json. Usage: python3 scripts/opentargets_fetch.py"""
import json, os, time, urllib.request

API = "https://api.platform.opentargets.org/api/v4/graphql"
OUT = "data/opentargets"
DISEASES = ["MONDO_0005192", "MONDO_0007254", "EFO_0004251"]
GENES = {"KRAS": "ENSG00000133703", "TP53": "ENSG00000141510",
         "CDKN2A": "ENSG00000147889", "SMAD4": "ENSG00000141646",
         "BRCA1": "ENSG00000012048", "JAK2": "ENSG00000096968"}

def gql(query):
    req = urllib.request.Request(
        API, data=json.dumps({"query": query}).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        d = json.load(r)
    if "errors" in d:
        raise RuntimeError(d["errors"][0]["message"])
    return d["data"]

def fetch_disease(dis):
    rows, idx, name, count = [], 0, None, None
    while True:
        q = ('{ disease(efoId:"%s") { name associatedTargets(page:{index:%d,size:3000}) '
             '{ count rows { target { id approvedSymbol } score '
             'datasourceScores { id score } } } } }') % (dis, idx)
        d = gql(q)["disease"]
        name, count = d["name"], d["associatedTargets"]["count"]
        batch = d["associatedTargets"]["rows"]
        rows.extend(batch)
        print(f"  {dis} page {idx}: +{len(batch)} (total {len(rows)}/{count})")
        if len(rows) >= count or not batch:
            break
        idx += 1
        time.sleep(0.3)
    assert len(rows) == count, f"{dis}: {len(rows)} != {count}"
    return {"disease_id": dis, "disease_name": name, "count": count, "rows": rows}

def fetch_target(sym, ensg):
    q = ('{ target(ensemblId:"%s") { approvedSymbol '
         'tractability { modality label value } '
         'drugAndClinicalCandidates { count rows { id maxClinicalStage '
         'drug { id name drugType } '
         'diseases { diseaseFromSource disease { id name } } } } } }') % ensg
    t = gql(q)["target"]
    assert t["approvedSymbol"] == sym, f"{ensg} -> {t['approvedSymbol']} != {sym}"
    return t

os.makedirs(OUT, exist_ok=True)
for dis in DISEASES:
    p = f"{OUT}/assoc_{dis}.json"
    if os.path.exists(p):
        print(f"  {dis} cached"); continue
    with open(p, "w") as f:
        json.dump(fetch_disease(dis), f)
for sym, ensg in GENES.items():
    p = f"{OUT}/target_{sym}.json"
    if os.path.exists(p):
        print(f"  {sym} cached"); continue
    with open(p, "w") as f:
        json.dump(fetch_target(sym, ensg), f)
    print(f"  {sym}: {t['drugAndClinicalCandidates']['count']} drug candidates")
    time.sleep(0.3)
print("done")
