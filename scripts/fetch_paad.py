"""Fetch TCGA-PAAD mutation + clinical data from cBioPortal public API (no auth)."""
import json, urllib.request, pathlib

BASE = "https://www.cbioportal.org/api"
STUDY = "paad_tcga_pan_can_atlas_2018"
OUT = pathlib.Path(__file__).resolve().parent.parent / "data" / "paad"
OUT.mkdir(parents=True, exist_ok=True)

def get(url, data=None):
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)

# 1. mutation molecular profile
profiles = get(f"{BASE}/studies/{STUDY}/molecular-profiles")
mut_prof = next(p for p in profiles if p["molecularAlterationType"] == "MUTATION_EXTENDED")
print("mutation profile:", mut_prof["molecularProfileId"])

# 2. sample list
slists = get(f"{BASE}/studies/{STUDY}/sample-lists")
all_list = next(s for s in slists if s["sampleListId"] == STUDY + "_all")
samples = get(f"{BASE}/sample-lists/{all_list['sampleListId']}/sample-ids")
print("samples:", len(samples))

# 3. all mutations in the profile for these samples
body = json.dumps({"sampleListId": all_list["sampleListId"]}).encode()
muts = get(f"{BASE}/molecular-profiles/{mut_prof['molecularProfileId']}/mutations/fetch?projection=DETAILED&pageSize=200000", data=body)
print("mutation records:", len(muts))
slim = [{"sampleId": m["sampleId"], "hugo": m["gene"]["hugoGeneSymbol"],
         "proteinChange": m.get("proteinChange", ""), "mutationType": m.get("mutationType", ""),
         "ref": m.get("referenceAllele", ""), "alt": m.get("variantAllele", ""),
         "chr": m.get("chr", ""), "start": m.get("startPosition", 0),
         "ncbi": m["gene"].get("ncbiTaxaId", None)} for m in muts]
(OUT / "mutations.json").write_text(json.dumps(slim))

# 4. clinical data
clin = get(f"{BASE}/studies/{STUDY}/clinical-data?clinicalDataType=SAMPLE&projection=DETAILED&pageSize=20000")
(OUT / "clinical_sample.json").write_text(json.dumps(clin))
pat = get(f"{BASE}/studies/{STUDY}/clinical-data?clinicalDataType=PATIENT&projection=DETAILED&pageSize=20000")
(OUT / "clinical_patient.json").write_text(json.dumps(pat))
(OUT / "samples.json").write_text(json.dumps(samples))
print("clinical sample attrs:", len(clin), "| patient attrs:", len(pat))
