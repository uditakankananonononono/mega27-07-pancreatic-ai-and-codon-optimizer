"""Summary stats for the 5'-MFE audit from the crash-safe checkpoint
(results/mfe_audit_partial.json). Handles methods missing per gene by
aligning genes pairwise; commercial best is per-gene max over available.
Writes results/mfe_audit.json (same schema as mfe_audit.py intended).
"""
import json
from pathlib import Path
import numpy as np
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]
per_gene = json.load(open(ROOT / "results" / "mfe_audit_partial.json"))
METHODS = ["all_original", "icor", "BFC", "ERC", "HFC", "URC", "genscript", "ours"]
missing = [(g, m) for g, row in per_gene.items() for m in METHODS if row.get(m) is None]

def pair(method_a, method_b, key, genes):
    ga = [g for g in genes if per_gene[g].get(method_a) and per_gene[g].get(method_b)]
    a = np.array([per_gene[g][method_a][key] for g in ga])
    b = np.array([per_gene[g][method_b][key] for g in ga])
    d = a - b
    stat = wilcoxon(d)
    return {"mean_a": round(float(a.mean()), 2), "mean_b": round(float(b.mean()), 2),
            "mean_delta": round(float(d.mean()), 2), "wilcoxon_p": float(stat.pvalue), "n": int(len(d))}

genes = sorted(per_gene)
comm = ["BFC", "ERC", "HFC", "URC", "genscript"]
best_pairs = []
for g in genes:
    if per_gene[g].get("ours") is None: continue
    vals = [per_gene[g][m]["five"] for m in comm if per_gene[g].get(m)]
    if vals: best_pairs.append((per_gene[g]["ours"]["five"], max(vals)))
ba = np.array([p[0] for p in best_pairs]); bb = np.array([p[1] for p in best_pairs])
bd = ba - bb
bstat = wilcoxon(bd)
ours_vs_best = {"mean_a": round(float(ba.mean()), 2), "mean_b": round(float(bb.mean()), 2),
                "mean_delta": round(float(bd.mean()), 2), "wilcoxon_p": float(bstat.pvalue), "n": int(len(bd))}

summary = {
    "five_prime_window_nt": 42,
    "ours_vs_icor_five": pair("ours", "icor", "five", genes),
    "ours_vs_wt_five": pair("ours", "all_original", "five", genes),
    "icor_vs_wt_five": pair("icor", "all_original", "five", genes),
    "ours_vs_best_commercial_five": ours_vs_best,
    "ours_vs_icor_full": pair("ours", "icor", "full", genes),
    "method_mean_five": {m: round(float(np.mean([per_gene[g][m]["five"] for g in genes if per_gene[g].get(m)])), 2) for m in METHODS},
    "missing": missing,
    "note": "MFE in kcal/mol (ViennaRNA 2.7.2); less negative = less 5' structure = better initiation. mean_delta = ours - other, positive favors ours.",
}
json.dump({"summary": summary, "per_gene": per_gene}, open(ROOT / "results" / "mfe_audit.json", "w"), indent=1)
print(json.dumps(summary, indent=1))
