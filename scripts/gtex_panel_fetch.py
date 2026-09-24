#!/usr/bin/env python3
"""Fetch GTEx v8 per-sample expression for the 4-gene PDAC detection panel plus controls.

Genes:
  panel:            KRAS, TP53, CDKN2A, SMAD4
  housekeeping:     GAPDH, ACTB
  pancreas-specific positive controls: PRSS1, INS

Source: GTEx Portal API v2 (https://gtexportal.org/api/v2), dataset gtex_v8.
Two calls per gene: reference/gene (symbol -> gencodeId) then
expression/geneExpression (per-sample TPM vectors keyed by tissueSiteDetailId).
Raw per-gene responses cached in data/gtex/<SYM>.json (untracked; refetch with
this script if wiped)."""
import json, os, time, urllib.request

GENES = {
    "panel": ["KRAS", "TP53", "CDKN2A", "SMAD4"],
    "housekeeping": ["GAPDH", "ACTB"],
    "pancreas_specific": ["PRSS1", "INS"],
}
BASE = "https://gtexportal.org/api/v2"
OUT = "data/gtex"

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "mega27/1.0"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())

os.makedirs(OUT, exist_ok=True)
for group, syms in GENES.items():
    for sym in syms:
        path = os.path.join(OUT, f"{sym}.json")
        if os.path.exists(path):
            print(sym, "cached"); continue
        ref = get(f"{BASE}/reference/gene?geneId={sym}")["data"]
        gid = ref[0]["gencodeId"]
        expr = get(f"{BASE}/expression/geneExpression?gencodeId={gid}&datasetId=gtex_v8")["data"]
        payload = {"symbol": sym, "group": group, "gencodeId": gid,
                   "datasetId": "gtex_v8",
                   "tissues": [{"tissueSiteDetailId": d["tissueSiteDetailId"],
                                "n_samples": len(d["data"]),
                                "tpm": d["data"]} for d in expr]}
        json.dump(payload, open(path, "w"))
        print(sym, gid, "tissues", len(expr), "samples", sum(len(d['data']) for d in expr))
        time.sleep(0.5)
print("done")
