"""Follow-up: permute the 11 covariate features (clinical + burden + pathway)
on the same held-out fold to locate the CDKN2A signal the GNN uses."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np, torch, torch.nn as nn
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import roc_auc_score
from pancreatic_ai import features
from pancreatic_ai.models import CoMutationGNN, comutation_adjacency

torch.manual_seed(0); np.random.seed(0)
X, C, y, genes = features.build_dataset("CDKN2A")
X100 = X[:, :100, :].copy()
(tr, te), = StratifiedShuffleSplit(1, test_size=0.25, random_state=7).split(X100, y)
model = CoMutationGNN(comutation_adjacency(X100[tr]), clin_dim=C.shape[1])
pos = int(y[tr].sum()); neg = len(tr) - pos
crit = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(neg / pos))
opt = torch.optim.Adam(model.parameters(), lr=3e-3, weight_decay=1e-3)
Xt = torch.tensor(X100); Ct = torch.tensor(C); yt = torch.tensor(y, dtype=torch.float32)
for _ in range(120):
    opt.zero_grad(); loss = crit(model(Xt[tr], Ct[tr]), yt[tr]); loss.backward(); opt.step()
model.eval()
with torch.no_grad():
    base = roc_auc_score(y[te], torch.sigmoid(model(Xt[te], Ct[te])).numpy())
names = ["age", "smoker", "log_burden", "frac_truncating", "log_n_genes"] + list(features.PATHWAYS)
rng = np.random.default_rng(1)
print(f"base {base:.3f}")
for j in range(C.shape[1]):
    d = []
    for r in range(6):
        Cp = C[te].copy()
        Cp[:, j] = Cp[rng.permutation(len(te)), j]
        with torch.no_grad():
            p = torch.sigmoid(model(Xt[te], torch.tensor(Cp))).numpy()
        d.append(base - roc_auc_score(y[te], p))
    print(f"{names[j]:16s} drop {np.mean(d):+.4f}", flush=True)
