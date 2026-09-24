"""Multi-objective codon optimization: maximize CNN-predicted expression (z)
subject to tAI >= per-gene ICOR floor. Closes the tAI deficit documented in
results/tai_benchmark.json (ours 0.3342 vs ICOR 0.3451) while keeping the
z-score win (+2.20 vs +0.68). Designs -> data/designs_mo/, metrics ->
results/multiobjective_benchmark.json.
"""
import json, sys, pathlib, random
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np, torch
from codon_optimizer.codon import parse_cds_fasta, cai_weights, cai, SYNONYMS
from codon_optimizer.model import ExpressionCNN, encode
from codon_optimizer.optimize import cai_greedy
from codon_optimizer.tai import copy_numbers, relative_weights, TaiScorer

ICOR = pathlib.Path("data/icor/Lattice-Automation-icor-codon-optimization-c60b775/benchmark_sequences")

def read_fasta_dir(d):
    out = {}
    for f in sorted(d.glob("*.fasta")):
        seq = "".join(l.strip() for l in f.read_text().splitlines() if not l.startswith(">"))
        out[f.stem.replace("_dna", "").replace("_aa", "")] = seq.upper().replace("U", "T")
    return out

copy = copy_numbers("data/gtrnadb/eschColi_K_12_MG1655-tRNAs.out")
scorer = TaiScorer(relative_weights(copy))
records = parse_cds_fasta()
wcai = cai_weights(records)
ck = torch.load("results_expr_model.pt", map_location="cpu", weights_only=False)
model = ExpressionCNN(); model.load_state_dict(ck["state"]); model.eval()
norm = ck["norm"]
mu, sd = norm

aa_seqs = read_fasta_dir(ICOR / "aa")
icor = read_fasta_dir(ICOR / "icor")

def gc(s):
    return (s.count("G") + s.count("C")) / max(len(s), 1)

def mo_hillclimb(protein, tai_floor, iters=150, batch=32, seed=7, gc_band=(0.35, 0.70)):
    """Hill-climb CNN score; only accept mutants with tAI >= floor."""
    rng = random.Random(seed)
    seq = cai_greedy(protein, wcai)
    if scorer.tai(seq) < tai_floor:  # nudge seed up to floor via max-tAI codons
        best = {}
        for aa, cods in SYNONYMS.items():
            best[aa] = max(cods, key=lambda c: scorer.logw.get(c) or -1e9)
        seq = "".join(best[aa] for aa in protein)
    positions = [i for i, aa in enumerate(protein) if len(SYNONYMS[aa]) > 1]
    def cnn(s):
        with torch.no_grad():
            return float(model(torch.tensor(encode(s)).unsqueeze(0)))
    cur = cnn(seq)
    for _ in range(iters):
        cands = []
        for _ in range(batch):
            pos = rng.choice(positions)
            aa = protein[pos]
            alt = rng.choice([c for c in SYNONYMS[aa] if c != seq[pos*3:pos*3+3]])
            cand = seq[:pos*3] + alt + seq[pos*3+3:]
            if gc_band[0] <= gc(cand) <= gc_band[1] and scorer.tai(cand) >= tai_floor:
                cands.append(cand)
        if not cands:
            continue
        with torch.no_grad():
            scores = model(torch.tensor(np.stack([encode(c) for c in cands]))).numpy()
        k = int(np.argmax(scores))
        if scores[k] > cur:
            cur = float(scores[k]); seq = cands[k]
    return seq, cur

desdir = pathlib.Path("data/designs_mo"); desdir.mkdir(exist_ok=True)
rows = {}
for g in sorted(aa_seqs):
    prot = aa_seqs[g].rstrip("*")
    if "X" in prot or len(prot) < 30 or g not in icor:
        continue
    floor = scorer.tai(icor[g])
    f = desdir / f"{g}_mo.fasta"
    if f.exists():
        seq = "".join(l.strip() for l in f.read_text().splitlines() if not l.startswith(">"))
    else:
        seq, _ = mo_hillclimb(prot, floor)
        f.write_text(f">{g}_mo\n{seq}\n")
    z = float((float(model(torch.tensor(encode(seq)).unsqueeze(0))) - float(mu)) / float(sd))
    with torch.no_grad():
        pass
    icor_z = None
    with torch.no_grad():
        icor_z = float((float(model(torch.tensor(encode(icor[g])).unsqueeze(0))) - float(mu)) / float(sd))
    rows[g] = {"tai": round(scorer.tai(seq), 4), "tai_icor": round(floor, 4),
               "cai": round(cai(seq, wcai), 4), "z": round(z, 3), "z_icor": round(icor_z, 3),
               "gc": round(gc(seq), 4)}
    print(g, rows[g], flush=True)

n = len(rows)
summ = {
    "n_genes": n,
    "tai_mean": round(float(np.mean([r["tai"] for r in rows.values()])), 4),
    "tai_icor_mean": round(float(np.mean([r["tai_icor"] for r in rows.values()])), 4),
    "cai_mean": round(float(np.mean([r["cai"] for r in rows.values()])), 4),
    "z_mean": round(float(np.mean([r["z"] for r in rows.values()])), 3),
    "z_icor_mean": round(float(np.mean([r["z_icor"] for r in rows.values()])), 3),
    "genes_tai_ge_icor": sum(1 for r in rows.values() if r["tai"] >= r["tai_icor"] - 1e-9),
    "genes_z_gt_icor": sum(1 for r in rows.values() if r["z"] > r["z_icor"]),
}
out = {"objective": "maximize CNN z s.t. tAI >= per-gene ICOR tAI floor",
       "summary": summ, "per_gene": rows}
json.dump(out, open("results/multiobjective_benchmark.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
