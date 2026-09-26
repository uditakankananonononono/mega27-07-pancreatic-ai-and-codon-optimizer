"""Seed-ensemble attribution: single-split permutation importance was shown to be
seed noise (cdkn2a_attribution_stability.json: pairwise Spearman -0.013/-0.016/0.042
across splits, held-out AUC 0.565-0.844). This script tests whether ENSEMBLING over
many split seeds recovers a stable gene ranking.

Protocol: 16 split seeds; per seed, identical GNN training + 100-gene permutation
importance (2 repeats) as the stability audit. Then split-half reliability: mean
ranking from seeds 1-8 vs seeds 9-16, Spearman rho. If rho(split-half ensemble) is
high while single-split rho is ~0, the ensemble ranking is a stable object and the
paper's attribution claim can be re-stated at ensemble level only.
"""
import json, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np
import torch
import torch.nn as nn
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
from pancreatic_ai import features
from pancreatic_ai.models import CoMutationGNN, comutation_adjacency

X, C, y, genes = features.build_dataset("CDKN2A")
X100 = X[:, :100, :].copy()
genes100 = list(genes[:100])

def attribute(split_seed, repeats=2):
    torch.manual_seed(split_seed); np.random.seed(split_seed)
    sss = StratifiedShuffleSplit(1, test_size=0.25, random_state=split_seed)
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
    rng = np.random.default_rng(split_seed)
    drops = np.zeros(len(genes100))
    for g in range(len(genes100)):
        d = []
        for r in range(repeats):
            Xp = X100[te].copy()
            Xp[:, g, :] = Xp[rng.permutation(len(te)), g, :]
            with torch.no_grad():
                p = torch.sigmoid(model(torch.tensor(Xp), Ct[te])).numpy()
            d.append(base - roc_auc_score(y[te], p))
        drops[g] = np.mean(d)
    return base, drops

SEEDS = list(range(101, 117))  # 16 fresh seeds, disjoint from prior audits
bases, D = [], []
for sd in SEEDS:
    b, dr = attribute(sd)
    bases.append(float(b)); D.append(dr)
    print(f"seed {sd}: AUC {b:.3f}", flush=True)
D = np.stack(D)
np.save("/tmp/ensemble_drops.npy", D)

mean_ranking = D.mean(0)
h1, h2 = D[:8].mean(0), D[8:].mean(0)
rho_half = float(spearmanr(h1, h2).statistic)
# single-split baseline rho for contrast (adjacent seed pairs)
pair_rhos = [float(spearmanr(D[i], D[i+1]).statistic) for i in range(0, 16, 2)]
top10 = sorted(zip(genes100, mean_ranking), key=lambda kv: -kv[1])[:10]
out = {
 "target": "CDKN2A", "seeds": SEEDS, "held_out_aucs": [round(b,3) for b in bases],
 "held_out_auc_median": round(float(np.median(bases)),3),
 "split_half_spearman_8v8": round(rho_half, 4),
 "single_split_pairwise_spearman": [round(r,3) for r in pair_rhos],
 "single_split_rho_median": round(float(np.median(pair_rhos)), 4),
 "ensemble_top10": [{"gene": g, "mean_auc_drop": round(float(v),5)} for g, v in top10],
 "reading": "if split-half rho >> single-split rho, ensemble ranking is stable though single splits are noise",
}
json.dump(out, open("results/cdkn2a_attribution_ensemble.json", "w"), indent=1)
print(json.dumps(out, indent=1))
