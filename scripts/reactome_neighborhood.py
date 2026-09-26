"""Reactome pathway enrichment of CDKN2A's STRING neighborhood.
Tool: Reactome Analysis Service (reactome.org). Connects the network context
(string_panel_network.json) to pathway-level evidence for the attribution story.
"""
import json, time, urllib.request, urllib.parse

partners = [p["partner"] for p in json.load(open("results/string_panel_network.json"))["cdkn2a_top_partners_score700"]]
genes = sorted(set(partners + ["CDKN2A"]))
req = urllib.request.Request(
    "https://reactome.org/AnalysisService/identifiers/projection?pageSize=20&page=1&sortBy=ENTITIES_PVALUE&order=ASC",
    data="\n".join(genes).encode(), headers={"Content-Type": "text/plain", "User-Agent": "Mozilla/5.0 (research script)"})
with urllib.request.urlopen(req, timeout=90) as r:
    d = json.load(r)
pathways = [{"name": p["name"], "pValue": p["entities"]["pValue"], "fdr": p["entities"]["fdr"],
             "found": p["entities"]["found"], "total": p["entities"]["total"]}
            for p in d.get("pathways", [])]
out = {"tool": "Reactome Analysis Service (reactome.org), accessed 2026-09-26",
       "input_genes": genes, "pathways": pathways}
json.dump(out, open("results/reactome_neighborhood.json", "w"), indent=1)
print(json.dumps({"n_genes": len(genes), "top5": [(p["name"][:60], p["fdr"]) for p in pathways[:5]]}, indent=1))
