"""UniProt REST annotation check for the PDAC driver panel.
Tool: UniProt REST API (rest.uniprot.org). Verifies reviewed entry, gene name,
protein name and length for each driver - provenance for the panel's identities.
"""
import json, time, urllib.request, urllib.parse

GENES = ["KRAS", "TP53", "CDKN2A", "SMAD4"]
out = {"tool": "UniProt REST (rest.uniprot.org), accessed 2026-09-26", "entries": {}}
for g in GENES:
    q = urllib.parse.urlencode({"query": f"gene:{g} AND organism_id:9606 AND reviewed:true",
                                "fields": "accession,id,gene_names,protein_name,length", "format": "tsv"})
    with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/search?{q}", timeout=60) as r:
        rows = [l.split("\t") for l in r.read().decode().strip().splitlines()]
    hdr = rows[0]
    recs = [dict(zip(hdr, r)) for r in rows[1:]]
    out["entries"][g] = recs[:3]
    print(g, [(r["Entry"], r["Length"]) for r in recs[:3]])
    time.sleep(0.5)
json.dump(out, open("results/uniprot_panel_check.json", "w"), indent=1)
