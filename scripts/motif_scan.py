"""Regulatory-motif scan (verdict #11): codon optimization can inadvertently
create or destroy regulatory motifs. Scan all benchmark arms for:
  - canonical polyadenylation signals (AATAAA, ATTAAA)
  - cryptic splice-donor consensus (GT[AG]AGT core)
  - AU-rich elements (ATTTA pentamer; mRNA-destabilizing)
  - CpG dinucleotide frequency (immune-stimulatory / silencing risk)
Arms: all_original (WT), icor, BFC, ERC, HFC, URC, genscript, ours (designs_mo).
Output: results/motif_scan.json
"""
import json, re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "data" / "icor" / "Lattice-Automation-icor-codon-optimization-c60b775" / "benchmark_sequences"
METHODS = ["all_original", "icor", "BFC", "ERC", "HFC", "URC", "genscript"]

def read_fasta(p):
    return "".join(l.strip() for l in p.read_text().splitlines() if not l.startswith(">")).upper()

def scan(seq):
    n = len(seq)
    return {
        "polyA": len(re.findall(r"AATAAA|ATTAAA", seq)),
        "splice_donor": len(re.findall(r"GT[AG]AGT", seq)),
        "ARE_ATTTA": len(re.findall(r"ATTTA", seq)),
        "cpg_per_kb": 1000.0 * seq.count("CG") / max(n, 1),
        "len": n,
    }

genes = sorted(p.name.replace("_mo.fasta", "") for p in (ROOT / "data" / "designs_mo").glob("*_mo.fasta"))
rows = {}
for g in genes:
    row = {}
    for m in METHODS:
        cands = list((BENCH / m).glob(f"{g}*.fasta")) + list((BENCH / m).glob(f"{g}*.fa"))
        if not cands and m == "all_original":
            cands = list((BENCH / "dna").glob(f"{g}*_dna.fasta"))
        if cands:
            row[m] = scan(read_fasta(cands[0]))
    row["ours"] = scan(read_fasta(ROOT / "data" / "designs_mo" / f"{g}_mo.fasta"))
    rows[g] = row

def agg(m):
    vals = [r[m] for r in rows.values() if m in r]
    return {k: round(float(np.mean([v[k] for v in vals])), 3) for k in ("polyA", "splice_donor", "ARE_ATTTA", "cpg_per_kb")} | {"n": len(vals)}

summary = {m: agg(m) for m in METHODS + ["ours"]}
# paired deltas ours vs WT for the two risk-increasing motifs
d_polyA = [r["ours"]["polyA"] - r["all_original"]["polyA"] for r in rows.values() if "all_original" in r]
d_splice = [r["ours"]["splice_donor"] - r["all_original"]["splice_donor"] for r in rows.values() if "all_original" in r]
out = {"per_gene": rows, "summary_means": summary,
       "ours_minus_wt_polyA": {"mean": float(np.mean(d_polyA)), "n_up": int(sum(x > 0 for x in d_polyA)), "n_down": int(sum(x < 0 for x in d_polyA)), "n": len(d_polyA)},
       "ours_minus_wt_splice_donor": {"mean": float(np.mean(d_splice)), "n_up": int(sum(x > 0 for x in d_splice)), "n_down": int(sum(x < 0 for x in d_splice)), "n": len(d_splice)}}
json.dump(out, open(ROOT / "results" / "motif_scan.json", "w"), indent=1)
for m in METHODS + ["ours"]:
    print(m, summary[m])
print("ours-WT polyA:", out["ours_minus_wt_polyA"])
print("ours-WT splice:", out["ours_minus_wt_splice_donor"])
