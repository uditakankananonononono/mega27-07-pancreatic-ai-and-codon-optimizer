#!/usr/bin/env python3
"""Fetch CIViC v2 clinical-evidence items for the 4-gene PDAC detection panel plus controls.

Genes:
  panel:    KRAS, TP53, CDKN2A, SMAD4
  controls: BRCA1 (breast/ovarian tissue-specific evidence control)
            JAK2  (myeloproliferative hematologic-only evidence control)

Source: CIViC v2 GraphQL API (https://civicdb.org/api/graphql), evidenceItems
endpoint with molecularProfileName substring match (exactly as CIViC indexes
evidence; a gene's tag includes simple and complex molecular profiles).
All statuses are fetched and recorded; accepted vs non-accepted is tracked in
the analysis. Raw per-gene responses cached in data/civic/<SYM>.json
(untracked; refetch with this script if wiped)."""
import json, os, time, urllib.request

GENES = {
    "panel": ["KRAS", "TP53", "CDKN2A", "SMAD4"],
    "breast_ovarian_control": ["BRCA1"],
    "hematologic_control": ["JAK2"],
}
BASE = "https://civicdb.org/api/graphql"
OUT = "data/civic"
PAGE = 50

QUERY = """
query($name: String!, $after: String) {
  evidenceItems(molecularProfileName: $name, first: 50, after: $after) {
    totalCount
    pageInfo { hasNextPage endCursor }
    nodes {
      id name evidenceType evidenceLevel evidenceDirection significance
      status variantOrigin evidenceRating
      disease { name doid }
      molecularProfile { name }
      therapies { name }
    }
  }
}"""

def post(payload):
    req = urllib.request.Request(
        BASE, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "mega27/1.0"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())

os.makedirs(OUT, exist_ok=True)
for group, syms in GENES.items():
    for sym in syms:
        path = os.path.join(OUT, f"{sym}.json")
        if os.path.exists(path):
            print(sym, "cached"); continue
        nodes, after, total = [], None, None
        while True:
            d = post({"query": QUERY, "variables": {"name": sym, "after": after}})
            ei = d["data"]["evidenceItems"]
            total = ei["totalCount"]
            nodes.extend(ei["nodes"])
            if not ei["pageInfo"]["hasNextPage"]:
                break
            after = ei["pageInfo"]["endCursor"]
            time.sleep(0.3)
        payload = {"symbol": sym, "group": group, "totalCount": total,
                   "fetched": len(nodes), "items": nodes}
        with open(path, "w") as f:
            json.dump(payload, f)
        print(sym, "total", total, "fetched", len(nodes))
        time.sleep(0.5)
