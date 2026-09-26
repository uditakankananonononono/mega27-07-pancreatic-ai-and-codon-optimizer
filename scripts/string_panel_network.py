"""STRING network context for the PDAC driver panel + CDKN2A neighborhood.
Tool: STRING v12 REST API (string-db.org). Adds STRING to the strict tool ledger
with a real use: (1) pairwise evidence-channel decomposition for the four
drivers; (2) CDKN2A's top interactors as network context for the attribution
analysis (complements the Fisher-validated covariate story).
"""
import json, time, urllib.request, urllib.parse

API = "https://string-db.org/api/tsv"
def get(endpoint, params):
    url = f"{API}/{endpoint}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=60) as r:
        txt = r.read().decode()
    rows = [l.split("\t") for l in txt.strip().splitlines()]
    return rows[0], rows[1:]

DRIVERS = ["KRAS", "TP53", "CDKN2A", "SMAD4"]
hdr, net = get("network", {"identifiers": "%0d".join(DRIVERS), "species": 9606, "required_score": 400})
i = {h: k for k, h in enumerate(hdr)}
pairwise = [{"a": r[i["preferredName_A"]], "b": r[i["preferredName_B"]],
             "score": float(r[i["score"]]),
             "channels": {c: float(r[i[c]]) for c in ("nscore","fscore","pscore","ascore","escore","dscore","tscore")}}
            for r in net]

time.sleep(1)
hdr2, nb = get("interaction_partners", {"identifiers": "CDKN2A", "species": 9606, "limit": 20, "required_score": 700})
j = {h: k for k, h in enumerate(hdr2)}
partners = [{"partner": r[j["preferredName_B"]], "score": float(r[j["score"]])} for r in nb]

out = {"tool": "STRING v12 REST (string-db.org), accessed 2026-09-26",
       "driver_pairwise": pairwise, "cdkn2a_top_partners_score700": partners}
json.dump(out, open("results/string_panel_network.json", "w"), indent=1)
print(json.dumps({"n_pairwise": len(pairwise), "partners": [p["partner"] for p in partners[:10]]}, indent=1))
