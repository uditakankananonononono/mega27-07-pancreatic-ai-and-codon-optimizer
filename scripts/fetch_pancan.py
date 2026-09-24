"""Pan-cancer CDKN2A/driver landscape: fetch a 4-gene mutation panel
(KRAS, TP53, CDKN2A, SMAD4) across cBioPortal studies. Each study = one
accession-level dataset; the analysis generalizes the PDAC CDKN2A finding
across cancer types. Output: data/pancan/<study>/{mutations.json,samples.json}
plus data/pancan/_manifest.json.
"""
import json, pathlib, urllib.request, time

BASE = "https://www.cbioportal.org/api"
PANEL = {"KRAS": 3845, "TP53": 7157, "CDKN2A": 1029, "SMAD4": 4089}
TARGET_N = 115
OUT = pathlib.Path("data/pancan"); OUT.mkdir(parents=True, exist_ok=True)

def get(url, data=None):
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.load(r)

studies = get(f"{BASE}/studies?projection=SUMMARY&pageSize=2000")
manifest = {}
done = {p.name for p in OUT.iterdir() if p.is_dir()}
n_ok = 0
# size-stratified spread: sort by sample count, take every k-th for a
# representative mix of small/medium/large cohorts (avoids tiny-cohort bias)
studies = sorted((s for s in studies if 25 <= s.get("allSampleCount", 0) <= 30000),
                 key=lambda s: s.get("allSampleCount", 0))
step = max(1, len(studies) // (TARGET_N * 2))
studies = studies[::step]
for st in studies:
    if n_ok >= TARGET_N:
        break
    sid = st["studyId"]
    if sid in done:
        continue
    try:
        profiles = get(f"{BASE}/studies/{sid}/molecular-profiles")
        mut_prof = next((p for p in profiles if p["molecularAlterationType"] == "MUTATION_EXTENDED"), None)
        if not mut_prof:
            continue
        slists = get(f"{BASE}/studies/{sid}/sample-lists")
        all_list = next((s for s in slists if s["sampleListId"].endswith("_all") or "sequenced" in s["sampleListId"]), None)
        if not all_list:
            continue
        body = json.dumps({"entrezGeneIds": list(PANEL.values()),
                           "sampleListId": all_list["sampleListId"]}).encode()
        muts = get(f"{BASE}/molecular-profiles/{mut_prof['molecularProfileId']}/mutations/fetch?projection=DETAILED&pageSize=500000", data=body)
        slim = [{"sampleId": m["sampleId"], "hugo": m["gene"]["hugoGeneSymbol"],
                 "mutationType": m.get("mutationType", "")} for m in muts]
        d = OUT / sid; d.mkdir(exist_ok=True)
        (d / "mutations.json").write_text(json.dumps(slim))
        (d / "samples.json").write_text(json.dumps(get(f"{BASE}/sample-lists/{all_list['sampleListId']}/sample-ids")))
        manifest[sid] = {"name": st.get("name", ""), "cancerType": st.get("cancerTypeId", ""),
                         "n_samples": len(json.loads((d / "samples.json").read_text())),
                         "n_panel_muts": len(slim)}
        (OUT / "_manifest.json").write_text(json.dumps(manifest, indent=1))
        n_ok += 1
        print(n_ok, sid, manifest[sid]["n_samples"], len(slim), flush=True)
    except Exception as e:
        print("FAIL", sid, repr(e)[:100], flush=True)
        time.sleep(2)
print("done:", n_ok, "studies")
