"""Fetch mutation data for ALL cBioPortal PDAC cohorts (multi-cohort replication)."""
import json, urllib.request, pathlib

BASE = "https://www.cbioportal.org/api"
STUDIES = ["paad_icgc", "paad_utsw_2015", "paad_cptac_2021", "paac_jhu_2014",
           "pancreas_msk_2024", "paad_msk_2025", "paad_tcga"]

def get(url, data=None):
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.load(r)

for study in STUDIES:
    out = pathlib.Path("data") / study
    out.mkdir(parents=True, exist_ok=True)
    try:
        profiles = get(f"{BASE}/studies/{study}/molecular-profiles")
        mut_prof = next((p for p in profiles if p["molecularAlterationType"] == "MUTATION_EXTENDED"), None)
        if not mut_prof:
            print(study, "no mutation profile", flush=True); continue
        slists = get(f"{BASE}/studies/{study}/sample-lists")
        all_list = next((s for s in slists if s["sampleListId"].endswith("_all") or "sequenced" in s["sampleListId"]), slists[0])
        body = json.dumps({"sampleListId": all_list["sampleListId"]}).encode()
        muts = get(f"{BASE}/molecular-profiles/{mut_prof['molecularProfileId']}/mutations/fetch?projection=DETAILED&pageSize=500000", data=body)
        slim = [{"sampleId": m["sampleId"], "hugo": m["gene"]["hugoGeneSymbol"],
                 "mutationType": m.get("mutationType", "")} for m in muts]
        (out / "mutations.json").write_text(json.dumps(slim))
        samples = get(f"{BASE}/sample-lists/{all_list['sampleListId']}/sample-ids")
        (out / "samples.json").write_text(json.dumps(samples))
        print(study, "muts", len(slim), "samples", len(samples), flush=True)
    except Exception as e:
        print(study, "FAIL", repr(e)[:120], flush=True)
