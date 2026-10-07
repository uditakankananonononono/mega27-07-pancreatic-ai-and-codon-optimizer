"""Sequence-feature analysis (verdict #19): which sequence features carry the
expression signal? Permutation importance of the 61 codon frequencies under the
ridge evaluator (5-fold CV, 10 permutations per feature), plus a 16-dinucleotide
view to separate codon usage from local sequence context.
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
DINUC = [a + b for a in "ACGT" for b in "ACGT"]

def codon_freq(seq):
    seq = seq.upper().replace("U", "T")
    n = len(seq) // 3
    counts = dict.fromkeys(CODONS, 0)
    for i in range(0, 3 * n, 3):
        c = seq[i:i + 3]
        if c in counts: counts[c] += 1
    return np.array([counts[c] / n for c in CODONS]) if n else np.zeros(len(CODONS))

def dinuc_freq(seq):
    seq = seq.upper().replace("U", "T")
    n = len(seq) - 1
    counts = dict.fromkeys(DINUC, 0)
    for i in range(n):
        d = seq[i:i+2]
        if d in counts: counts[d] += 1
    return np.array([counts[d] / n for d in DINUC]) if n > 0 else np.zeros(len(DINUC))

def perm_importance(X, y, names, tag):
    kf = KFold(5, shuffle=True, random_state=3)
    base, imps = [], np.zeros((5, X.shape[1]))
    rng = np.random.default_rng(0)
    for k, (tr, te) in enumerate(kf.split(X)):
        m = RidgeCV(alphas=np.logspace(-3, 3, 13)).fit(X[tr], y[tr])
        r0 = np.corrcoef(m.predict(X[te]), y[te])[0, 1]
        base.append(r0)
        for j in range(X.shape[1]):
            drops = []
            for _ in range(10):
                Xp = X[te].copy(); rng.shuffle(Xp[:, j])
                drops.append(r0 - np.corrcoef(m.predict(Xp), y[te])[0, 1])
            imps[k, j] = np.mean(drops)
    imp = imps.mean(0)
    order = np.argsort(-imp)
    print(tag, "base r:", round(np.mean(base), 3), flush=True)
    return {"base_r_mean": float(np.mean(base)),
            "top_features": [{"feature": names[i], "importance": round(float(imp[i]), 4)} for i in order[:15]]}

def main():
    rows = [(lt, seq, logab) for lt, gene, seq, logab in build_expression_dataset(policy="raw") if len(seq) >= 90]
    y = np.array([r[2] for r in rows])
    out = {"design": "permutation importance, ridge, 5-fold CV, 10 permutations/feature; importance = mean r drop"}
    out["codon_view"] = perm_importance(np.array([codon_freq(r[1]) for r in rows]), y, CODONS, "codon")
    out["dinuc_view"] = perm_importance(np.array([dinuc_freq(r[1]) for r in rows]), y, DINUC, "dinuc")
    json.dump(out, open(ROOT / "results" / "sequence_features.json", "w"), indent=1)
    print(json.dumps(out["codon_view"]["top_features"][:10], indent=1))
    print(json.dumps(out["dinuc_view"]["top_features"][:8], indent=1))

if __name__ == "__main__":
    main()
