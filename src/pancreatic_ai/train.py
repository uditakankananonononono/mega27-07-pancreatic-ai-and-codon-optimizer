"""Stratified 5-fold CV benchmark: CNN vs GNN vs logistic baseline."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "src"))
import numpy as np
import torch
import torch.nn as nn
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from pancreatic_ai import features
from pancreatic_ai.models import MutationCNN, CoMutationGNN, comutation_adjacency, baseline_logreg

torch.manual_seed(0); np.random.seed(0)
EPOCHS = 60

def train_eval(model_fn, X, C, y, use_adj=False):
    skf = StratifiedKFold(5, shuffle=True, random_state=42)
    aucs = []
    Xt = torch.tensor(X, dtype=torch.float32)
    Ct = torch.tensor(C, dtype=torch.float32)
    for tr, te in skf.split(X, y):
        model = model_fn(X[tr])
        pos = max(int(y[tr].sum()), 1); neg = len(tr) - pos
        crit = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(neg / pos))
        opt = torch.optim.Adam(model.parameters(), lr=3e-3, weight_decay=1e-3)
        for _ in range(EPOCHS):
            model.train(); opt.zero_grad()
            loss = crit(model(Xt[tr], Ct[tr]), torch.tensor(y[tr], dtype=torch.float32))
            loss.backward(); opt.step()
        model.eval()
        with torch.no_grad():
            p = torch.sigmoid(model(Xt[te], Ct[te])).numpy()
        aucs.append(roc_auc_score(y[te], p))
    return float(np.mean(aucs)), float(np.std(aucs))

def run():
    results = {}
    for gene in features.DRIVERS:
        X, C, y, _ = features.build_dataset(gene)
        if y.sum() < 10 or (len(y) - y.sum()) < 10:
            results[gene] = {"skipped": "too few positives/negatives"}; continue
        cnn = train_eval(lambda xtr: MutationCNN(n_genes=X.shape[1]), X, C, y)
        X100 = X[:, :100, :].copy()
        gnn = train_eval(lambda xtr: CoMutationGNN(comutation_adjacency(xtr)), X100, C, y)
        flat = X.reshape(len(X), -1)
        skf = StratifiedKFold(5, shuffle=True, random_state=42)
        b_aucs = []
        for tr, te in skf.split(flat, y):
            clf = baseline_logreg().fit(np.hstack([flat[tr], C[tr]]), y[tr])
            b_aucs.append(roc_auc_score(y[te], clf.predict_proba(np.hstack([flat[te], C[te]]))[:, 1]))
        results[gene] = {
            "prevalence": round(float(y.mean()), 3),
            "cnn_auc": [round(a, 3) for a in cnn],
            "gnn_auc": [round(a, 3) for a in gnn],
            "logreg_auc": [round(float(np.mean(b_aucs)), 3), round(float(np.std(b_aucs)), 3)],
        }
        print(gene, results[gene], flush=True)
    import json
    pathlib.Path("results").mkdir(exist_ok=True)
    pathlib.Path("results/pancreatic_ai_benchmark.json").write_text(json.dumps(results, indent=2))
    return results

if __name__ == "__main__":
    run()
