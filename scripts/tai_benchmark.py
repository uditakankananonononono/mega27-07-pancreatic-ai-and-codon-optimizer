"""tRNA Adaptation Index (tAI, dos Reis et al. 2004 NAR) on the ICOR 40-gene
benchmark for all methods + our regenerated designs (saved to data/designs/).
tRNA gene copy numbers: GtRNAdb E. coli K-12 MG1655 (tRNAscan-SE). Wobble
s-values from dos Reis 2004: WC 0.0, G:U 0.41, I:U/C/A 0.28, U:G 0.68 (rare),
I:G 0.9999. Also saves our design fastas (reproducibility fix: previously unsaved).
"""
import json, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np, torch
from codon_optimizer.codon import parse_cds_fasta, cai_weights, cai, STANDARD_CODE
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

# --- tRNA copy numbers from GtRNAdb tRNAscan-SE output
copy = {}
for line in open("data/gtrnadb/eschColi_K_12_MG1655-tRNAs.out"):
    p = line.split()
    if len(p) >= 6 and p[0] == "chr" and p[4] != "Pseudo" and len(p[5]) == 3:
        ac = p[5].replace("T", "U")
        copy[ac] = copy.get(ac, 0) + 1
print("anticodons:", len(copy), "total tRNAs:", sum(copy.values()), flush=True)

# wobble pairing: codon position 1 pairs anticodon position 3 (written 5'->3')
def revcomp(s):
    return s.translate(str.maketrans("ACGU", "UGCA"))[::-1]
S_WC, S_GU, S_IU, S_IG = 0.0, 0.41, 0.28, 0.9999
codons = [c for c in STANDARD_CODE if STANDARD_CODE[c] != "*"]
W = {}
for cod in codons:
    c = cod.replace("T", "U")
    acc = 0.0
    for ac, n in copy.items():
        anti = revcomp(ac)  # align anticodon 5'->3' to codon
        if anti[1:] != c[1:]:
            continue
        a1, c1 = ac[2], c[0]  # anticodon wobble base (pos 34) vs codon first base
        pair = a1 + c1
        if pair in ("AU", "UA", "GC", "CG"): s = S_WC
        elif pair == "GU": s = S_GU
        elif pair in ("IU", "IC", "IA"): s = S_IU
        elif pair == "IG": s = S_IG
        else: continue
        acc += (1 - s) * n
    W[cod] = acc
wmax = max(v for v in W.values() if v > 0)
wrel = {c: (v / wmax if v > 0 else 0.0) for c, v in W.items()}

def tai(seq):
    vals = [wrel.get(seq[i:i+3], 0.0) for i in range(0, len(seq) - 2, 3)]
    vals = [v for v in vals if v > 0]
    return float(np.exp(np.mean(np.log(vals)))) if vals else 0.0

records = parse_cds_fasta()
wcai = cai_weights(records)
ck = torch.load("results_expr_model.pt", map_location="cpu", weights_only=False)
model = ExpressionCNN(); model.load_state_dict(ck["state"]); model.eval()
norm = ck["norm"]
aa_seqs = read_fasta_dir(ICOR / "aa")
method_seqs = {m: read_fasta_dir(ICOR / d) for m, d in METHODS.items()}

desdir = pathlib.Path("data/designs"); desdir.mkdir(exist_ok=True)
ours = {}
genes = sorted(aa_seqs)
for g in genes:
    f = desdir / f"{g}_codonopt.fasta"
    if f.exists():
        seq = "".join(l.strip() for l in f.read_text().splitlines() if not l.startswith(">"))
    else:
        prot = aa_seqs[g].rstrip("*")
        if "X" in prot or len(prot) < 30: continue
        seq, _ = cnn_hillclimb(prot, wcai, model, norm, iters=120, batch=32, seed=1)
        f.write_text(f">{g}_codonopt\n{seq}\n")
    ours[g] = seq
print("our designs:", len(ours), flush=True)

rows = {}
for g in genes:
    row = {}
    for m in METHODS:
        s = method_seqs[m].get(g)
        if s: row[m] = round(tai(s), 4)
    if g in ours: row["ours_cnn"] = round(tai(ours[g]), 4)
    rows[g] = row
summ = {m: round(float(np.mean([r[m] for r in rows.values() if m in r])), 4)
        for m in list(METHODS) + ["ours_cnn"]}
out = {"metric": "tAI (dos Reis 2004), GtRNAdb E.coli K-12 tRNA copy numbers",
       "summary_tai": summ, "per_gene": rows, "anticodon_copy": copy}
json.dump(out, open("results/tai_benchmark.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
