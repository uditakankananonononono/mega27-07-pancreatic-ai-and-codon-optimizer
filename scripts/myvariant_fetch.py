#!/usr/bin/env python3
"""Allele-level annotation of PDAC somatic variants via MyVariant.info.

Step 1: fetch allele-level somatic mutations (chr/start/ref/alt/proteinChange)
from cBioPortal (public API, no auth) for the 4 panel driver genes and 8
long, frequently-mutated passenger genes, across 7 non-overlapping hg19 PDAC
studies. Step 2: annotate every SNV with MyVariant.info (hg19 HGVS ids, batch
POST): CADD phred, dbNSFP AlphaMissense + REVEL, ClinVar RCV significance.
Raw JSON -> data/myvariant/ (untracked); flattened table ->
results/myvariant_variants.csv (committed). Resumable."""
import json, os, csv, time, urllib.request, urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data", "myvariant")
RES = os.path.join(HERE, "..", "results")
CBIO = "https://www.cbioportal.org/api"
MV = "https://myvariant.info/v1/variant"
# paad_tcga (legacy) duplicates the PanCan Atlas patients; pancreas_msk_2024 and
# paad_msk_2025 can share MSK patients with pdac_msk_2024 -> excluded.
STUDIES = ["paad_tcga_pan_can_atlas_2018", "paad_qcmg_uq_2016", "paad_utsw_2015",
           "paad_icgc", "paad_cptac_2021", "pdac_msk_2024", "paac_jhu_2014"]
GENES = {"KRAS": (3845, "panel"), "TP53": (7157, "panel"), "CDKN2A": (1029, "panel"),
         "SMAD4": (4089, "panel"), "TTN": (7273, "passenger"), "MUC16": (94025, "passenger"),
         "RYR1": (6261, "passenger"), "LRP1B": (53353, "passenger"), "CSMD3": (114788, "passenger"),
         "USH2A": (7399, "passenger"), "SYNE1": (23345, "passenger"), "FLG": (2312, "passenger")}
FIELDS = "cadd.phred,dbnsfp.alphamissense.score,dbnsfp.revel.score,clinvar.rcv.clinical_significance"


def req(url, data=None, ctype="application/json"):
    r = urllib.request.Request(url, data=data, headers={"Content-Type": ctype, "Accept": "application/json",
                                                        "User-Agent": "Mozilla/5.0 mega27-laneF"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(r, timeout=120) as f:
                return json.load(f)
        except Exception:
            if attempt == 3:
                raise
            time.sleep(3 * (attempt + 1))


def fetch_study(study):
    dest = os.path.join(RAW, f"cbio_{study}.json")
    if os.path.exists(dest):
        return json.load(open(dest))
    prof = next(p for p in req(f"{CBIO}/studies/{study}/molecular-profiles")
                if p["molecularAlterationType"] == "MUTATION_EXTENDED")["molecularProfileId"]
    body = json.dumps({"sampleListId": f"{study}_all",
                       "entrezGeneIds": [v[0] for v in GENES.values()]}).encode()
    muts = req(f"{CBIO}/molecular-profiles/{prof}/mutations/fetch?projection=DETAILED&pageSize=100000", body)
    slim = [{"study": study, "sampleId": m["sampleId"], "patientId": m.get("patientId", m["sampleId"]),
             "gene": m["gene"]["hugoGeneSymbol"], "proteinChange": m.get("proteinChange", ""),
             "mutationType": m.get("mutationType", ""), "chr": str(m.get("chr", "")),
             "start": m.get("startPosition"), "end": m.get("endPosition"),
             "ref": m.get("referenceAllele", ""), "alt": m.get("variantAllele", "")} for m in muts]
    json.dump(slim, open(dest, "w"))
    print(study, len(slim), flush=True)
    return slim


def hgvs(m):
    ref, alt = m["ref"], m["alt"]
    if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT" and m["start"]:
        c = m["chr"].replace("chr", "")
        c = {"23": "X", "24": "Y"}.get(c, c)
        return f"chr{c}:g.{m['start']}{ref}>{alt}"
    return None


def annotate(ids):
    dest = os.path.join(RAW, "myvariant_annotations.json")
    ann = json.load(open(dest)) if os.path.exists(dest) else {}
    todo = [i for i in ids if i not in ann]
    for k in range(0, len(todo), 500):
        chunk = todo[k:k + 500]
        data = urllib.parse.urlencode({"ids": ",".join(chunk), "fields": FIELDS, "assembly": "hg19"}).encode()
        out = req(MV, data, "application/x-www-form-urlencoded")
        for rec in out:
            ann[rec["query"]] = rec
        json.dump(ann, open(dest, "w"))
        print("annotated", len(ann), flush=True)
    return ann


def first(x):
    if isinstance(x, list):
        x = [v for v in x if v is not None]
        return max(x) if x else None
    return x


def clinvar(rec):
    rcv = (rec.get("clinvar") or {}).get("rcv")
    if rcv is None:
        return ""
    rcv = rcv if isinstance(rcv, list) else [rcv]
    sig = [str(r.get("clinical_significance", "")) for r in rcv]
    if any(s.lower() in ("pathogenic", "likely pathogenic", "pathogenic/likely pathogenic") for s in sig):
        return "P/LP"
    if any("benign" in s.lower() for s in sig):
        return "B/LB"
    return "other"


def main():
    os.makedirs(RAW, exist_ok=True)
    muts = [m for s in STUDIES for m in fetch_study(s)]
    ids = sorted({h for h in (hgvs(m) for m in muts) if h})
    ann = annotate(ids)
    rows = []
    for m in muts:
        h = hgvs(m)
        rec = ann.get(h, {}) if h else {}
        found = bool(rec) and not rec.get("notfound")
        db = rec.get("dbnsfp") or {}
        rows.append({**{k: m[k] for k in ("study", "sampleId", "patientId", "gene", "proteinChange",
                                         "mutationType", "chr", "start", "ref", "alt")},
                     "group": GENES[m["gene"]][1] if m["gene"] in GENES else "other",
                     "hgvs": h or "", "mv_found": int(found),
                     "cadd_phred": first((rec.get("cadd") or {}).get("phred")) if found else "",
                     "alphamissense": first((db.get("alphamissense") or {}).get("score")) if found else "",
                     "revel": first((db.get("revel") or {}).get("score")) if found else "",
                     "clinvar": clinvar(rec) if found else ""})
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "myvariant_variants.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("rows", len(rows), "snv ids", len(ids), "found", sum(r["mv_found"] for r in rows))


if __name__ == "__main__":
    main()
