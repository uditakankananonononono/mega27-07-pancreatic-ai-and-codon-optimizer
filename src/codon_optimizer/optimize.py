"""Synonymous recoding optimizer.

Strategies compared on identical protein sequences:
- wildtype: native CDS
- cai_greedy: highest-frequency synonymous codon per position
- cnn_hillclimb: CAI-max seed + CNN-guided hill climbing (batch-scored mutants)
Metrics: CAI and CNN-predicted expression (z-scored log10 ppm).
"""
import math
import random

import numpy as np
import torch

from .codon import STANDARD_CODE, SYNONYMS, cai, cai_weights
from .model import encode


def cai_greedy(protein, w):
    best = {aa: max(cods, key=lambda c: w.get(c, 0.0)) for aa, cods in SYNONYMS.items()}
    return "".join(best[aa] for aa in protein)


def to_protein(cds):
    prot = "".join(STANDARD_CODE.get(cds[i:i+3], "X") for i in range(0, len(cds) - 2, 3))
    return prot.rstrip("*")


def cnn_hillclimb(protein, w, model, norm, iters=250, batch=48, seed=0, gc_band=(0.35, 0.70)):
    """Hill-climb synonymous codons to maximize CNN-predicted expression."""
    rng = random.Random(seed)
    seq = cai_greedy(protein, w)
    mu, sd = norm
    positions = [i for i, aa in enumerate(protein) if len(SYNONYMS[aa]) > 1]

    def score(s):
        with torch.no_grad():
            return float(model(torch.tensor(encode(s)).unsqueeze(0)))

    def gc(s):
        return (s.count("G") + s.count("C")) / max(len(s), 1)

    cur = score(seq)
    for _ in range(iters):
        cands = []
        for _ in range(batch):
            pos = rng.choice(positions)
            aa = protein[pos]
            alt = rng.choice([c for c in SYNONYMS[aa] if c != seq[pos*3:pos*3+3]])
            cand = seq[:pos*3] + alt + seq[pos*3+3:]
            if gc_band[0] <= gc(cand) <= gc_band[1]:
                cands.append((pos, cand))
        if not cands:
            continue
        with torch.no_grad():
            scores = model(torch.tensor(np.stack([encode(c) for _, c in cands]))).numpy()
        k = int(np.argmax(scores))
        if scores[k] > cur:
            cur = float(scores[k]); seq = cands[k][1]
    return seq, cur


def benchmark(proteins_cds, model, norm, w, n=25, seed=0):
    m = model
    rng = random.Random(seed)
    picks = rng.sample(proteins_cds, min(n, len(proteins_cds)))
    rows = []
    for lt, gene, cds, abun in picks:
        prot = to_protein(cds)
        if "X" in prot or "*" in prot or len(prot) < 50:
            continue
        opt_seq, _ = cnn_hillclimb(prot, w, model, norm)
        cai_seq = cai_greedy(prot, w)
        def pred(s):
            with torch.no_grad():
                return float(m(torch.tensor(encode(s)).unsqueeze(0)))
        rows.append({
            "locus": lt, "gene": gene, "len_aa": len(prot),
            "wt_cai": round(cai(cds, w), 3), "cai_greedy_cai": round(cai(cai_seq, w), 3),
            "cnn_opt_cai": round(cai(opt_seq, w), 3),
            "wt_pred": round(pred(cds), 3), "cai_greedy_pred": round(pred(cai_seq), 3),
            "cnn_opt_pred": round(pred(opt_seq), 3),
        })
    return rows
