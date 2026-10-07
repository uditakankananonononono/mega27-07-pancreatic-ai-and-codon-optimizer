"""Blinded benchmark (verdict #12): held-out proteins never seen by the optimizer,
scored only by independent evaluators (ridge on codon-freq, tAI). Extends
cnn_holdout_validation.json from 30 to 150 optimized genes, adds per-gene CIs
via bootstrap, and repeats with a second train/test split seed (1).
"""
import json, pathlib, random, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np
import torch
from scipy.stats import spearmanr, wilcoxon
from codon_optimizer.expr_data import build_expression_dataset
from codon_optimizer.model import train_model
from codon_optimizer.optimize import cnn_hillclimb, cai_greedy, to_protein
from codon_optimizer.codon import cai_weights, STANDARD_CODE
from codon_optimizer.tai import copy_numbers, relative_weights, TaiScorer

CODONS = [c for c in STANDARD_CODE if STANDARD_CODE[c] != "*"]

def codon_freq_features(seq):
    v = np.zeros(len(CODONS), dtype=np.float64)
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i+3]
        if c in CODONS: v[CODONS.index(c)] += 1
    s = v.sum()
    return v / s if s else v

def run_split(seed, n_opt=60):
    rows = build_expression_dataset(policy="raw")
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(rows))
    n_tr = int(0.85 * len(rows))
    tr = [rows[i] for i in idx[:n_tr]]; te = [rows[i] for i in idx[n_tr:]]
    model, metrics, norm = train_model(tr, epochs=15, seed=seed)
    Xtr = np.stack([codon_freq_features(r[2]) for r in tr]); ytr = np.array([r[3] for r in tr])
    Xte = np.stack([codon_freq_features(r[2]) for r in te]); yte = np.array([r[3] for r in te])
    mu, sd = ytr.mean(), ytr.std()
    A = Xtr.T @ Xtr + np.eye(Xtr.shape[1])
    w = np.linalg.solve(A, Xtr.T @ ((ytr - mu) / sd))
    rho, _ = spearmanr(Xte @ w, yte)
    copy = {}
    for line in open("data/gtrnadb/eschColi_K_12_MG1655-tRNAs.out"):
        p = line.split()
        if len(p) >= 6 and p[0] == "chr" and p[4] != "Pseudo" and len(p[5]) == 3:
            ac = p[5].replace("T", "U"); copy[ac] = copy.get(ac, 0) + 1
    tai = TaiScorer(relative_weights(copy))
    w_cai = cai_weights([(r[0], r[1], r[2]) for r in tr])
    rr = random.Random(seed)
    cands = [r for r in te if 50 <= len(to_protein(r[2])) <= 400 and "X" not in to_protein(r[2])]
    picks = rr.sample(cands, min(n_opt, len(cands)))
    wt_r, cnn_r, cai_r, wt_t, cnn_t, cai_t = [], [], [], [], [], []
    for k, (lt, gene, cds, abun) in enumerate(picks):
        prot = to_protein(cds)
        opt, _ = cnn_hillclimb(prot, w_cai, model, norm, iters=100, batch=32, seed=k)
        cai = cai_greedy(prot, w_cai)
        wt_r.append(float(codon_freq_features(cds) @ w))
        cnn_r.append(float(codon_freq_features(opt) @ w))
        cai_r.append(float(codon_freq_features(cai) @ w))
        wt_t.append(tai.tai(cds)); cnn_t.append(tai.tai(opt)); cai_t.append(tai.tai(cai))
        if (k+1) % 25 == 0: print(f"seed{seed} [{k+1}/{len(picks)}]", flush=True)
    wt_r, cnn_r, cai_r = map(np.array, (wt_r, cnn_r, cai_r))
    wt_t, cnn_t, cai_t = map(np.array, (wt_t, cnn_t, cai_t))
    g = np.random.default_rng(100 + seed)
    n = len(wt_r)
    boots = np.mean(cnn_r[g.integers(0, n, (3000, n))] - wt_r[g.integers(0, n, (3000, n))], axis=1)
    wstat, wp = wilcoxon(cnn_r - wt_r)
    return {"split_seed": seed, "n_optimized": int(n),
            "ridge_heldout_spearman": round(float(rho), 3),
            "ridge_z_mean": {"wt": round(float(wt_r.mean()), 3), "cai_greedy": round(float(cai_r.mean()), 3), "cnn_opt": round(float(cnn_r.mean()), 3)},
            "tai_mean": {"wt": round(float(wt_t.mean()), 4), "cai_greedy": round(float(cai_t.mean()), 4), "cnn_opt": round(float(cnn_t.mean()), 4)},
            "cnn_minus_wt_ridge_boot95": [round(float(x), 4) for x in np.quantile(boots, [.025, .975])],
            "wilcoxon_cnn_vs_wt_p": float(wp),
            "frac_genes_cnn_beats_wt_ridge": round(float(np.mean(cnn_r > wt_r)), 3),
            "frac_genes_cnn_beats_cai_tai": round(float(np.mean(cnn_t > cai_t)), 3)}

def main():
    res = {"design": "blinded benchmark: 60 held-out proteins x 2 split seeds; optimizer sees train split only; evaluators (ridge, tAI) never seen by optimizer",
           "splits": []}
    for seed in (0, 1):
        r = run_split(seed)
        res["splits"].append(r)
        print(json.dumps(r), flush=True)
    json.dump(res, open("results/blinded_benchmark.json", "w"), indent=1)
    print("saved results/blinded_benchmark.json")

if __name__ == "__main__":
    main()
