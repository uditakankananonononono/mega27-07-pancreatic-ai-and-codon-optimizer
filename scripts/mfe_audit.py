"""5'-structure audit (item 7): ViennaRNA MFE of the translation-initiation
region (first 42 nt of the CDS) across codon-optimization methods.

Strong 5'-CDS secondary structure is a known translation-initiation barrier,
so a good optimizer should not inflate 5' folding energy even while raising
tAI/CAI. We compute RNAfold MFE (ViennaRNA 2.7.2 Python API) for the first
42 nt and the full CDS for: wild-type (all_original), ICOR, five commercial
baselines (BFC, ERC, HFC, URC, GenScript), and our multiobjective CNN
designs (data/designs_mo). Paired Wilcoxon: ours vs ICOR and ours vs
wild-type on 5' MFE (positive delta = less structure = better initiation).
Output: results/mfe_audit.json
"""
import json
from pathlib import Path
import numpy as np
import RNA
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "data" / "icor" / "Lattice-Automation-icor-codon-optimization-c60b775" / "benchmark_sequences"
WIN = 42
METHODS = ["all_original", "icor", "BFC", "ERC", "HFC", "URC", "genscript"]

def read_fasta(p):
    lines = p.read_text().splitlines()
    return "".join(l.strip() for l in lines if not l.startswith(">")).upper().replace("T", "U")

def mfe(seq):
    return float(RNA.fold(seq)[1])

genes = sorted(p.name.replace("_mo.fasta", "") for p in (ROOT / "data" / "designs_mo").glob("*_mo.fasta"))
per_gene, by_method = {}, {m: {"five": [], "full": []} for m in METHODS + ["ours"]}
missing = []
for g in genes:
    row = {}
    for m in METHODS:
        cands = list((BENCH / m).glob(f"{g}*.fasta")) + list((BENCH / m).glob(f"{g}*.fa"))
        if not cands and m == "all_original":
            cands = list((BENCH / "dna").glob(f"{g}*_dna.fasta"))  # per-gene WT == all_original records (verified identical)
        if not cands:
            row[m] = None; missing.append((g, m)); continue
        s = read_fasta(cands[0])
        row[m] = {"five": round(mfe(s[:WIN]), 2), "full": round(mfe(s), 2)}
        by_method[m]["five"].append(row[m]["five"]); by_method[m]["full"].append(row[m]["full"])
    s = read_fasta(ROOT / "data" / "designs_mo" / f"{g}_mo.fasta")
    row["ours"] = {"five": round(mfe(s[:WIN]), 2), "full": round(mfe(s), 2)}
    by_method["ours"]["five"].append(row["ours"]["five"]); by_method["ours"]["full"].append(row["ours"]["full"])
    per_gene[g] = row
    json.dump(per_gene, open(ROOT / "results" / "mfe_audit_partial.json", "w"), indent=1)  # crash-safe checkpoint

def pair(a, b):
    a = np.array(a); b = np.array(b)
    d = a - b
    stat = wilcoxon(d)
    return {"mean_a": round(float(a.mean()), 2), "mean_b": round(float(b.mean()), 2),
            "mean_delta": round(float(d.mean()), 2), "wilcoxon_p": float(stat.pvalue),
            "n": int(len(d))}

five = by_method
summary = {
    "five_prime_window_nt": WIN,
    "ours_vs_icor_five": pair(five["ours"]["five"], five["icor"]["five"]),
    "ours_vs_wt_five": pair(five["ours"]["five"], five["all_original"]["five"]),
    "icor_vs_wt_five": pair(five["icor"]["five"], five["all_original"]["five"]),
    "ours_vs_best_commercial_five": pair(  # per-gene best of available commercial methods
        [per_gene[g]["ours"]["five"] for g in genes if per_gene[g].get("ours") and any(per_gene[g].get(m) for m in ("BFC","ERC","HFC","URC","genscript"))],
        [max(per_gene[g][m]["five"] for m in ("BFC","ERC","HFC","URC","genscript") if per_gene[g].get(m))
         for g in genes if per_gene[g].get("ours") and any(per_gene[g].get(m) for m in ("BFC","ERC","HFC","URC","genscript"))]),
    "ours_vs_icor_full": pair(five["ours"]["full"], five["icor"]["full"]),
    "method_mean_five": {m: round(float(np.mean(five[m]["five"])), 2) for m in five},
    "missing": missing,
}
json.dump({"summary": summary, "per_gene": per_gene}, open(ROOT / "results" / "mfe_audit.json", "w"), indent=1)
print(json.dumps(summary, indent=1)[:1600])
