"""Independent-evaluator hold-out validation of CNN-guided codon optimization.

Circularity problem: cnn_hillclimb maximizes ExpressionCNN-predicted expression,
and codon_benchmark.json scores designs with the SAME network. Here we:
  1. split PaxDb E. coli expression genes into train/test (seed 0, 85/15)
  2. train ExpressionCNN on train only (optimizer sees only this)
  3. train an INDEPENDENT ridge evaluator on 64-dim codon-frequency features,
     train genes only - never used by the optimizer
  4. on held-out test genes, run wildtype / cai_greedy / cnn_hillclimb and
     score all designs with the ridge evaluator AND tAI (two metrics the
     optimizer never saw)
  5. Wilcoxon signed-rank wt vs cnn-opt under the ridge evaluator.
Output: results/cnn_holdout_validation.json
"""
import json, math, pathlib, random, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np
import torch
from codon_optimizer.expr_data import build_expression_dataset
from codon_optimizer.model import train_model, encode
from codon_optimizer.optimize import cnn_hillclimb, cai_greedy, to_protein
from codon_optimizer.codon import cai_weights, STANDARD_CODE, DATA
from codon_optimizer.tai import copy_numbers, relative_weights, TaiScorer

CODONS = [c for c in STANDARD_CODE if STANDARD_CODE[c] != "*"]

def codon_freq_features(seq):
    v = np.zeros(len(CODONS), dtype=np.float64)
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i+3]
        if c in CODONS:
            v[CODONS.index(c)] += 1
    s = v.sum()
    return v / s if s else v

def main():
    rows = build_expression_dataset()
    print(f"dataset: {len(rows)} genes", flush=True)
    rng = np.random.default_rng(0)
    idx = rng.permutation(len(rows))
    n_tr = int(0.85 * len(rows))
    tr_rows = [rows[i] for i in idx[:n_tr]]
    te_rows = [rows[i] for i in idx[n_tr:]]
    print(f"train {len(tr_rows)} / test {len(te_rows)}", flush=True)

    # 1) optimizer's model: CNN on train only
    model, cnn_metrics, norm = train_model(tr_rows, epochs=15, seed=0)
    print(f"CNN (train-only) internal holdout: {cnn_metrics}", flush=True)

    # 2) independent ridge evaluator on codon-frequency features, train only
    Xtr = np.stack([codon_freq_features(r[2]) for r in tr_rows])
    ytr = np.array([r[3] for r in tr_rows])
    Xte = np.stack([codon_freq_features(r[2]) for r in te_rows])
    yte = np.array([r[3] for r in te_rows])
    mu, sd = ytr.mean(), ytr.std()
    lam = 1.0
    A = Xtr.T @ Xtr + lam * np.eye(Xtr.shape[1])
    w_ridge = np.linalg.solve(A, Xtr.T @ ((ytr - mu) / sd))
    from scipy.stats import spearmanr
    ridge_rho, _ = spearmanr(Xte @ w_ridge, yte)
    print(f"ridge independent evaluator held-out spearman: {ridge_rho:.3f}", flush=True)

    # 3) tAI scorer (second independent metric)
    copy = {}
    for line in open("data/gtrnadb/eschColi_K_12_MG1655-tRNAs.out"):
        p = line.split()
        if len(p) >= 6 and p[0] == "chr" and p[4] != "Pseudo" and len(p[5]) == 3:
            ac = p[5].replace("T", "U")
            copy[ac] = copy.get(ac, 0) + 1
    tai_scorer = TaiScorer(relative_weights(copy))

    # 4) optimize held-out genes, score with both independent metrics
    w_cai = cai_weights([(r[0], r[1], r[2]) for r in tr_rows])
    rr = random.Random(0)
    cands = [r for r in te_rows if 50 <= len(to_protein(r[2])) <= 400 and "X" not in to_protein(r[2])]
    picks = rr.sample(cands, min(30, len(cands)))
    out = []
    for k, (lt, gene, cds, abun) in enumerate(picks):
        prot = to_protein(cds)
        opt_seq, _ = cnn_hillclimb(prot, w_cai, model, norm, iters=100, batch=32, seed=k)
        cai_seq = cai_greedy(prot, w_cai)
        row = {"locus": lt, "gene": gene}
        for name, s in (("wt", cds), ("cai_greedy", cai_seq), ("cnn_opt", opt_seq)):
            row[f"{name}_ridge_z"] = round(float((codon_freq_features(s) @ w_ridge)), 3)
            row[f"{name}_tai"] = round(tai_scorer.tai(s), 4)
        out.append(row)
        print(f"[{k+1}/{len(picks)}] {gene}: ridge_z wt {row['wt_ridge_z']} cnn {row['cnn_opt_ridge_z']} | tAI wt {row['wt_tai']} cnn {row['cnn_opt_tai']}", flush=True)

    from scipy.stats import wilcoxon
    wt_ridge = np.array([r["wt_ridge_z"] for r in out])
    cnn_ridge = np.array([r["cnn_opt_ridge_z"] for r in out])
    cai_ridge = np.array([r["cai_greedy_ridge_z"] for r in out])
    wt_tai = np.array([r["wt_tai"] for r in out])
    cnn_tai = np.array([r["cnn_opt_tai"] for r in out])
    try:
        wstat, wp = wilcoxon(cnn_ridge - wt_ridge)
    except ValueError:
        wstat, wp = None, None
    res = {
        "design": "independent-evaluator hold-out: optimizer=CNN(train genes only); evaluators=codon-freq ridge (train only) + tAI; test=held-out genes",
        "n_train": len(tr_rows), "n_test_pool": len(te_rows), "n_optimized": len(out),
        "cnn_internal_holdout": cnn_metrics,
        "ridge_heldout_spearman": round(float(ridge_rho), 3),
        "ridge_z": {"wt": round(float(wt_ridge.mean()), 3),
                     "cai_greedy": round(float(cai_ridge.mean()), 3),
                     "cnn_opt": round(float(cnn_ridge.mean()), 3)},
        "tai": {"wt": round(float(wt_tai.mean()), 4),
                "cnn_opt": round(float(cnn_tai.mean()), 4)},
        "wilcoxon_cnn_vs_wt_ridge": {"stat": (float(wstat) if wstat is not None else None),
                                      "p": (float(wp) if wp is not None else None)},
        "cnn_wins_ridge": int((cnn_ridge > wt_ridge).sum()),
        "cnn_wins_tai": int((cnn_tai > wt_tai).sum()),
        "per_gene": out,
    }
    json.dump(res, open("results/cnn_holdout_validation.json", "w"), indent=2)
    print("saved results/cnn_holdout_validation.json", flush=True)

if __name__ == "__main__":
    main()
