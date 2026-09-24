"""CDKN2A-loss signal attribution: which co-mutated genes carry the signal the
GNN exploits (AUC 0.728 vs logreg 0.679)? Permutation importance per gene on a
held-out fold, GNN trained on the rest. Discovery output: ranked candidate
genes, to be tested on an independent PDAC cohort (ICGC/CPTAC via cBioPortal).
"""
import json, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np
import torch
import torch.nn as nn
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import roc_auc_score
from pancreatic_ai import features
from pancreatic_ai.models import CoMutationGNN, comutation_adjacency

torch.manual_seed(0); np.random.seed(0)
REPEATS = 4

X, C, y, genes = features.build_dataset("CDKN2A")
X100 = X[:, :100, :].copy()
genes100 = genes[:100]
sss = StratifiedShuffleSplit(1, test_size=0.25, random_state=7)
(tr, te), = sss.split(X100, y)
adj = comutation_adjacency(X100[tr])
model = CoMutationGNN(adj, clin_dim=C.shape[1])
pos = int(y[tr].sum()); neg = len(tr) - pos
crit = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(neg / pos))
opt = torch.optim.Adam(model.parameters(), lr=3e-3, weight_decay=1e-3)
Xt = torch.tensor(X100); Ct = torch.tensor(C); yt = torch.tensor(y, dtype=torch.float32)
for ep in range(120):
    model.train(); opt.zero_grad()
    loss = crit(model(Xt[tr], Ct[tr]), yt[tr]); loss.backward(); opt.step()
model.eval()
with torch.no_grad():
    base = roc_auc_score(y[te], torch.sigmoid(model(Xt[te], Ct[te])).numpy())
print(f"held-out AUC {base:.3f}  (n_test={len(te)}, pos={int(y[te].sum())})", flush=True)

rng = np.random.default_rng(0)
drops = np.zeros(len(genes100))
for g in range(len(genes100)):
    d = []
    for r in range(REPEATS):
        Xp = X100[te].copy()
        Xp[:, g, :] = Xp[rng.permutation(len(te)), g, :]
        with torch.no_grad():
            p = torch.sigmoid(model(torch.tensor(Xp), Ct[te])).numpy()
        d.append(base - roc_auc_score(y[te], p))
    drops[g] = np.mean(d)
rank = np.argsort(-drops)
out = {"target": "CDKN2A", "held_out_auc": round(float(base), 3),
       "n_test": int(len(te)), "method": "permutation importance, GNN, 4 repeats",
       "ranked_genes": [{"gene": genes100[i], "auc_drop": round(float(drops[i]), 4)}
                        for i in rank[:25]]}
for e in out["ranked_genes"][:12]:
    print(e, flush=True)
pathlib.Path("results/cdkn2a_attribution.json").write_text(json.dumps(out, indent=2))
print("saved results/cdkn2a_attribution.json")
