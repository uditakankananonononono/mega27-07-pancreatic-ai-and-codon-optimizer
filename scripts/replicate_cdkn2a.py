"""Replicate CDKN2A co-mutation/exclusivity tests in independent cohorts
(mutation-defined status, matching the TCGA-PAAD analysis)."""
import json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from scipy.stats import fisher_exact
from pancreatic_ai import features

def run(study):
    d = pathlib.Path("data") / study
    muts = json.loads((d / "mutations.json").read_text())
    samples = json.loads((d / "samples.json").read_text())
    per = {s: {m["hugo"] for m in muts if m["sampleId"] == s
               and m["mutationType"] in features.NONSYNONYMOUS} for s in samples}
    y = np.array([1 if "CDKN2A" in per[s] else 0 for s in samples])
    out = {"study": study, "n": len(samples), "cdkn2a_mut_n": int(y.sum())}
    for label, gs in [("KRAS", ["KRAS"]), ("TP53_pathway", features.PATHWAYS["TP53"]),
                      ("SWI_SNF", features.PATHWAYS["SWI_SNF"]),
                      ("RTK_RAS_nonKRAS", [g for g in features.PATHWAYS["RTK_RAS"] if g != "KRAS"])]:
        g = np.array([1 if per[s] & set(gs) else 0 for s in samples])
        tab = np.array([[int(((y==1)&(g==1)).sum()), int(((y==1)&(g==0)).sum())],
                        [int(((y==0)&(g==1)).sum()), int(((y==0)&(g==0)).sum())]])
        odds, p = fisher_exact(tab)
        out[label] = {"table": tab.tolist(), "odds_ratio": round(float(odds),3) if np.isfinite(odds) else str(odds),
                      "p": float(f"{p:.2e}")}
    return out

res = [run(s) for s in ["paad_qcmg_uq_2016", "pdac_msk_2024"]]
print(json.dumps(res, indent=1))
pathlib.Path("results/cdkn2a_replication.json").write_text(json.dumps(res, indent=1))
