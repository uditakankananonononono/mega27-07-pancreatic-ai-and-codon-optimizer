"""Rescore the ICOR 40-gene benchmark under the source-pinned CLI checkpoint
(results_expr_model.pt, PaxDb 4.2 / MG1655 RefSeq, 5-epoch, holdout Spearman 0.564).

Why: the historical icor_headtohead.json was scored with an earlier checkpoint that
was never committed and is unrecoverable (overwritten by train_cli_checkpoint.py).
The historical ours_cnn z=+2.20 therefore stands as a historical claim only; this
script re-establishes the comparison end-to-end under the auditable model:
(a) all fixed comparator sequences rescored (ranking stability);
(b) fresh ours_cnn_v2 designs hill-climbed against the new model (same protocol:
    iters=120, batch=32, seed=1, gc_band default).
"""
import json, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np
import torch
from scipy import stats
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
assert "mean" in ck and "norm" not in ck, "expected source-pinned checkpoint"
model = ExpressionCNN(); model.load_state_dict(ck["state"]); model.eval()
norm = (ck["mean"], ck["sd"])

def pred(s):
    with torch.no_grad():
        return float(model(torch.tensor(encode(s)).unsqueeze(0)))

aa_seqs = read_fasta_dir(ICOR / "aa")
method_seqs = {m: read_fasta_dir(ICOR / d) for m, d in METHODS.items()}
genes = sorted(aa_seqs)
rows = {}
for g in genes:
    row = {}
    for m in METHODS:
        s = method_seqs[m].get(g)
        if s:
            row[m] = {"cai": round(cai(s, w), 3), "gc": round(gc(s), 3),
                      "pred_expr_z": round(pred(s), 3)}
    rows[g] = row
    print(g, {m: row[m]["pred_expr_z"] for m in row}, flush=True)

PROG = pathlib.Path("results/icor_rescore_progress.json")
done = json.load(open(PROG)) if PROG.exists() else {}
for g in genes:
    if g in done:
        rows.setdefault(g, {})["ours_cnn_v2"] = done[g]
        continue
    prot = aa_seqs[g].rstrip("*")
    if "X" in prot or len(prot) < 30:
        continue
    seq, _ = cnn_hillclimb(prot, w, model, norm, iters=120, batch=32, seed=1)
    rows.setdefault(g, {})["ours_cnn_v2"] = {"cai": round(cai(seq, w), 3), "gc": round(gc(seq), 3),
                                             "pred_expr_z": round(pred(seq), 3), "dna": seq}
    done[g] = rows[g]["ours_cnn_v2"]
    json.dump(done, open(PROG, "w"))
    print(g, "ours_v2", rows[g]["ours_cnn_v2"]["pred_expr_z"], flush=True)

methods = list(METHODS) + ["ours_cnn_v2"]
summary = {}
for m in methods:
    vals = [rows[g][m] for g in rows if m in rows[g]]
    if vals:
        summary[m] = {"n": len(vals),
                      "cai_mean": round(float(np.mean([v["cai"] for v in vals])), 3),
                      "gc_mean": round(float(np.mean([v["gc"] for v in vals])), 3),
                      "pred_expr_z_mean": round(float(np.mean([v["pred_expr_z"] for v in vals])), 3)}
# paired comparison ours vs icor on shared genes
shared = [g for g in rows if "ours_cnn_v2" in rows[g] and "icor" in rows[g]]
a = np.array([rows[g]["ours_cnn_v2"]["pred_expr_z"] for g in shared])
b = np.array([rows[g]["icor"]["pred_expr_z"] for g in shared])
wil = stats.wilcoxon(a, b)
out = {"model": "source-pinned CLI checkpoint (PaxDb 4.2/MG1655, holdout Spearman 0.564)",
       "historical_note": "prior icor_headtohead.json used an unrecoverable earlier checkpoint; its numbers are historical only",
       "summary": summary,
       "paired_ours_vs_icor": {"n": len(shared), "mean_diff": float((a - b).mean()),
                               "wilcoxon_p": float(wil.pvalue)},
       "per_gene": {g: {m: {k: v for k, v in r.items() if k != "dna"} for m, r in row.items()} for g, row in rows.items()}}
json.dump(out, open("results/icor_rescore.json", "w"), indent=1)
json.dump({g: rows[g]["ours_cnn_v2"]["dna"] for g in shared}, open("results/icor_rescore_ours_dna.json", "w"), indent=1)
print(json.dumps({"summary": summary, "paired": out["paired_ours_vs_icor"]}, indent=1), flush=True)
print("saved results/icor_rescore.json")
