"""Fisher-exact co-mutation enrichment: CDKN2A status vs pathway lesions
(non-KRAS RTK-RAS split from KRAS itself) in TCGA-PAAD. Falsifiable claim
for independent-cohort replication."""
import json, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from scipy.stats import fisher_exact
from pancreatic_ai import features

muts, samples, _ = features.load_raw()
status = features.build_targets(muts, samples)
y = np.array([status[s]["CDKN2A"] for s in samples])
per_sample = {s: {m["hugo"] for m in muts if m["sampleId"] == s
                  and m["mutationType"] in features.NONSYNONYMOUS} for s in samples}

def enrich(label, geneset):
    g = np.array([1.0 if per_sample[s] & set(geneset) else 0.0 for s in samples])
    tab = np.array([[int(((y == 1) & (g == 1)).sum()), int(((y == 1) & (g == 0)).sum())],
                    [int(((y == 0) & (g == 1)).sum()), int(((y == 0) & (g == 0)).sum())]])
    odds, p = fisher_exact(tab)
    return {"label": label, "table": tab.tolist(), "odds_ratio": round(float(odds), 3),
            "p": float(f"{p:.2e}"), "prev_in_cdkn2a": round(tab[0,0]/max(tab[0].sum(),1),3),
            "prev_in_wt": round(tab[1,0]/max(tab[1].sum(),1),3)}

ras_nonkras = [g for g in features.PATHWAYS["RTK_RAS"] if g != "KRAS"]
tests = [enrich("RTK_RAS_nonKRAS", ras_nonkras),
         enrich("KRAS", ["KRAS"]),
         enrich("SWI_SNF", features.PATHWAYS["SWI_SNF"]),
         enrich("DNA_REPAIR", features.PATHWAYS["DNA_REPAIR"]),
         enrich("TP53_pathway", features.PATHWAYS["TP53"]),
         enrich("TGF_BETA", features.PATHWAYS["TGF_BETA"])]
for t in tests:
    print(t, flush=True)
pathlib.Path("results/cdkn2a_comut_enrichment.json").write_text(json.dumps(tests, indent=2))
