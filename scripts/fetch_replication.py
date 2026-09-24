"""Fetch mutations + CNA for independent PDAC cohorts (QCMG 2016, MSK 2024)."""
import json, urllib.request, pathlib, sys

BASE = "https://www.cbioportal.org/api"

def get(url, data=None):
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.load(r)

for study in ["paad_qcmg_uq_2016", "pdac_msk_2024"]:
    out = pathlib.Path("data") / study
    out.mkdir(parents=True, exist_ok=True)
    try:
        profiles = get(f"{BASE}/studies/{study}/molecular-profiles")
        mut_prof = next((p for p in profiles if p["molecularAlterationType"] == "MUTATION_EXTENDED"), None)
        cna_prof = next((p for p in profiles if p["molecularAlterationType"] == "COPY_NUMBER_ALTERATION" and "DISCRETE" in p.get("datatype","")), None)
        slists = get(f"{BASE}/studies/{study}/sample-lists")
        all_list = next(s for s in slists if s["sampleListId"] == study + "_all")
        samples = get(f"{BASE}/sample-lists/{all_list['sampleListId']}/sample-ids")
        body = json.dumps({"sampleListId": all_list["sampleListId"]}).encode()
        if mut_prof:
            muts = get(f"{BASE}/molecular-profiles/{mut_prof['molecularProfileId']}/mutations/fetch?projection=DETAILED&pageSize=500000", data=body)
            slim = [{"sampleId": m["sampleId"], "hugo": m["gene"]["hugoGeneSymbol"],
                     "mutationType": m.get("mutationType", "")} for m in muts]
            (out / "mutations.json").write_text(json.dumps(slim))
            print(study, "muts", len(slim), "samples", len(samples), flush=True)
        if cna_prof:
            genes = ["CDKN2A", "ARID1A", "ARID1B", "ARID2", "SMARCA4", "SMARCB1", "PBRM1", "KRAS", "TP53"]
            ents = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in get(f"{BASE}/genes/fetch", data=json.dumps({"geneIds": genes, "geneIdType": "HUGO_GENE_SYMBOL"}).encode())}
            dbody = json.dumps({"sampleListId": all_list["sampleListId"], "entrezGeneIds": list(ents.values())}).encode()
            disc = get(f"{BASE}/molecular-profiles/{cna_prof['molecularProfileId']}/discrete-copy-number/fetch?discreteCopyNumberEventType=HOMDEL_AND_AMP&projection=DETAILED", data=dbody)
            slim = [{"sampleId": d["sampleId"], "entrez": d["entrezGeneId"], "alt": d["alteration"]} for d in disc]
            (out / "cna.json").write_text(json.dumps({"map": ents, "rows": slim}))
            print(study, "cna rows", len(slim), flush=True)
        (out / "samples.json").write_text(json.dumps(samples))
    except Exception as e:
        print(study, "FAILED", repr(e), flush=True)
