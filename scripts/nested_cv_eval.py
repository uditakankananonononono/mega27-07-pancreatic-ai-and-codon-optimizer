"""Nested cross-validation statistics (verdict #5): single-split held-out scores
carry split noise; nested CV gives bias-corrected performance with a defensible
uncertainty estimate. 5 outer folds x 3 inner folds (alpha selection), on:
  (a) 61 codon-frequency features + ridge  (matches the RF evaluator's inputs)
  (b) frozen NT v2-50M embeddings + ridge  (matches the #18 evaluator)
Reports per-outer-fold Pearson r, mean +/- SE, and bootstrap 95% CI of the mean.
"""
import json, pathlib, sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from codon_optimizer.expr_data import build_expression_dataset
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
CODONS = [c for c in CODONS if c not in ("TAA", "TAG", "TGA")]

def codon_freq(seq):
    seq = seq.upper().replace("U", "T")
    n = len(seq) // 3
    counts = dict.fromkeys(CODONS, 0)
    for i in range(0, 3 * n, 3):
        c = seq[i:i + 3]
        if c in counts: counts[c] += 1
    return np.array([counts[c] / n for c in CODONS]) if n else np.zeros(len(CODONS))

def nested_cv(X, y, tag):
    outer = KFold(5, shuffle=True, random_state=7)
    rs = []
    for k, (tr, te) in enumerate(outer.split(X)):
        m = RidgeCV(alphas=np.logspace(-3, 3, 13), cv=KFold(3, shuffle=True, random_state=k)).fit(X[tr], y[tr])
        r = float(np.corrcoef(m.predict(X[te]), y[te])[0, 1])
        rs.append(r)
        print(f"{tag} outer fold {k+1}: r={r:.3f} alpha={m.alpha_}", flush=True)
    rs = np.array(rs)
    rng = np.random.default_rng(0)
    boots = rng.choice(rs, size=(10000, len(rs)), replace=True).mean(axis=1)
    return {"per_fold_r": [round(float(x), 4) for x in rs],
            "mean_r": float(rs.mean()), "se": float(rs.std(ddof=1) / np.sqrt(len(rs))),
            "boot95_mean": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]}

def main():
    rows = [(lt, seq, logab) for lt, gene, seq, logab in build_expression_dataset() if len(seq) >= 90]
    y = np.array([r[2] for r in rows])
    out = {"design": "5x3 nested CV, RidgeCV inner alpha selection, Pearson r on outer folds",
           "n_genes": len(rows), "views": {}}
    Xc = np.array([codon_freq(r[1]) for r in rows])
    out["views"]["codon_freq_ridge"] = nested_cv(Xc, y, "codon")
    z = np.load(ROOT / "results" / "nt_embeddings.npz", allow_pickle=True)
    emb = dict(zip(z["keys"].tolist(), z["vals"]))
    Xn = np.array([emb[r[0]] for r in rows])
    out["views"]["nt_frozen_ridge"] = nested_cv(Xn, y, "nt")
    json.dump(out, open(ROOT / "results" / "nested_cv_eval.json", "w"), indent=1)
    print(json.dumps(out["views"], indent=1))

if __name__ == "__main__":
    main()
