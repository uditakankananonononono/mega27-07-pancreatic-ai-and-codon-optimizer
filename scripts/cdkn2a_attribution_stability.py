"""Attribution stability audit: does the CDKN2A permutation-importance ranking
replicate across train/test split seeds? The committed ranking
(results/cdkn2a_attribution.json, split seed 7) shows tiny drops (max 0.0038)
led by passenger-looking genes. If the ranking is seed-noise, the attribution
claim must be downgraded. Spearman rank correlation of per-gene AUC drops
across 3 fresh splits (seeds 11, 23, 42), 2 permutation repeats each.
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
genes100 = genes[:100]

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

bases, all_drops = [], []
for sd in (11, 23, 42):
    b, dr = attribute(sd)
    bases.append(b); all_drops.append(dr)
    print(f"seed {sd}: held-out AUC {b:.3f}", flush=True)

all_drops = np.stack(all_drops)
rhos, _ = spearmanr(all_drops.T, axis=0) if False else (None, None)
from itertools import combinations
pair_rhos = []
for a, b in combinations(range(3), 2):
    r = spearmanr(all_drops[a], all_drops[b]).statistic
    pair_rhos.append(float(r))
mean_drop = all_drops.mean(0)
rank = np.argsort(-mean_drop)
out = {"target": "CDKN2A", "split_seeds": [11, 23, 42], "held_out_aucs": [round(float(b), 3) for b in bases],
       "pairwise_spearman_of_gene_drops": [round(r, 3) for r in pair_rhos],
       "mean_drop_top10": [{"gene": genes100[i], "mean_auc_drop": round(float(mean_drop[i]), 4)} for i in rank[:10]],
       "verdict": "stable" if min(pair_rhos) > 0.5 else "seed-noise: attribution ranking does NOT replicate across splits"}
json_path = pathlib.Path("results/cdkn2a_attribution_stability.json")
json_path.write_text(json.dumps(out, indent=2))
print(json.dumps({k: out[k] for k in ("held_out_aucs", "pairwise_spearman_of_gene_drops", "verdict")}, indent=1))
print("top10:", [e["gene"] for e in out["mean_drop_top10"]])
