"""Head-to-head vs ICOR (Jain et al. 2023, BMC Bioinformatics) on ITS OWN
40-gene benchmark: original vs ICOR vs GenSmart vs HFC/BFC/URC/ERC vs our
CNN-hillclimb, scored on CAI (E. coli high-expression weights), GC content,
and our PaxDb-trained ExpressionCNN predicted expression (z).
ICOR sequences: Zenodo 10.5281/zenodo.7487432 (v1.4 benchmark_sequences)."""
import json, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np
import torch
from codon_optimizer.codon import parse_cds_fasta, cai_weights, cai
from codon_optimizer.model import ExpressionCNN, encode
from codon_optimizer.optimize import cnn_hillclimb

ICOR = pathlib.Path("data/icor/Lattice-Automation-icor-codon-optimization-c60b775/benchmark_sequences")
METHODS = {"original": "dna", "icor": "icor", "gensmart": "genscript",
           "HFC": "HFC", "BFC": "BFC", "URC": "URC", "ERC": "ERC"}

def read_fasta_dir(d):
    out = {}
    for f in sorted(d.glob("*.fasta")):
        seq = "".join(l.strip() for l in f.read_text().splitlines() if not l.startswith(">"))
        out[f.stem.replace("_dna", "").replace("_aa", "")] = seq.upper().replace("U", "T")
    return out

def gc(s):
    return (s.count("G") + s.count("C")) / max(len(s), 1)

records = parse_cds_fasta()
w = cai_weights(records)
ck = torch.load("results_expr_model.pt", map_location="cpu", weights_only=False)
model = ExpressionCNN(); model.load_state_dict(ck["state"]); model.eval()
norm = ck["norm"]

def pred(s):
    with torch.no_grad():
        return float(model(torch.tensor(encode(s)).unsqueeze(0)))

aa_seqs = read_fasta_dir(ICOR / "aa")
method_seqs = {m: read_fasta_dir(ICOR / d) for m, d in METHODS.items()}
genes = sorted(aa_seqs)
print("benchmark genes:", len(genes), flush=True)
rows = {}
for g in genes:
    row = {}
    for m in METHODS:
        s = method_seqs[m].get(g)
        if s:
            row[m] = {"cai": round(cai(s, w), 3), "gc": round(gc(s), 3),
                      "pred_expr_z": round(pred(s), 2)}
    rows[g] = row
    print(g, {m: row[m]["cai"] for m in row}, flush=True)
json.dump({"genes": rows}, open("results/icor_headtohead_per_gene.json", "w"), indent=1)

# our optimizer on the same proteins
ours = {}
for g in genes:
    prot = aa_seqs[g].rstrip("*")
    if "X" in prot or len(prot) < 30:
        continue
    seq, _ = cnn_hillclimb(prot, w, model, norm, iters=120, batch=32, seed=1)
    ours[g] = {"cai": round(cai(seq, w), 3), "gc": round(gc(seq), 3),
               "pred_expr_z": round(pred(seq), 2), "dna": seq}
    rows.setdefault(g, {})["ours_cnn"] = {k: v for k, v in ours[g].items() if k != "dna"}
    print(g, "ours", ours[g]["cai"], ours[g]["pred_expr_z"], flush=True)

methods = list(METHODS) + ["ours_cnn"]
summary = {}
for m in methods:
    vals = [rows[g][m] for g in rows if m in rows[g]]
    if vals:
        summary[m] = {"n": len(vals),
                      "cai_mean": round(float(np.mean([v["cai"] for v in vals])), 3),
                      "gc_mean": round(float(np.mean([v["gc"] for v in vals])), 3),
                      "pred_expr_z_mean": round(float(np.mean([v["pred_expr_z"] for v in vals])), 2)}
print(json.dumps(summary, indent=1), flush=True)
json.dump({"metric_notes": "CAI weights from E. coli high-expression CDS (same reference as our codon_benchmark); pred_expr_z from PaxDb-trained ExpressionCNN (held-out spearman 0.61)",
           "summary": summary, "per_gene": rows},
          open("results/icor_headtohead.json", "w"), indent=1)
print("saved results/icor_headtohead.json")
