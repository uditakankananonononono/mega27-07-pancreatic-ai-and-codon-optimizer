"""CNN regressor: CDS one-hot sequence -> log10 protein abundance (PaxDb)."""
import numpy as np
import torch
import torch.nn as nn

BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def encode(seq, max_len=1500):
    x = np.zeros((4, max_len), dtype=np.float32)
    for i, b in enumerate(seq[:max_len]):
        j = BASE_IDX.get(b)
        if j is not None:
            x[j, i] = 1.0
    return x


class ExpressionCNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv1d(4, 48, kernel_size=9, padding=4), nn.ReLU(),
            nn.Conv1d(48, 64, kernel_size=9, padding=4, stride=3), nn.ReLU(),
            nn.Conv1d(64, 96, kernel_size=7, padding=3, stride=3), nn.ReLU(),
        )
        self.head = nn.Sequential(nn.Linear(192, 64), nn.ReLU(), nn.Linear(64, 1))

    def forward(self, x):
        h = self.net(x)
        h = torch.cat([h.mean(dim=2), h.amax(dim=2)], dim=1)
        return self.head(h).squeeze(-1)


def train_model(rows, epochs=15, seed=0, max_len=1500):
    from scipy.stats import spearmanr, pearsonr
    torch.manual_seed(seed); np.random.seed(seed)
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(rows))
    n_tr = int(0.85 * len(rows))
    tr, te = idx[:n_tr], idx[n_tr:]
    X = np.stack([encode(r[2], max_len) for r in rows])
    y = np.array([r[3] for r in rows], dtype=np.float32)
    mu, sd = y[tr].mean(), y[tr].std()
    yz = (y - mu) / sd
    model = ExpressionCNN()
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    crit = nn.MSELoss()
    bs = 64
    for ep in range(epochs):
        model.train()
        order = rng.permutation(tr)
        tot = 0.0
        for i in range(0, len(order), bs):
            b = order[i:i + bs]
            opt.zero_grad()
            loss = crit(model(torch.tensor(X[b])), torch.tensor(yz[b]))
            loss.backward(); opt.step()
            tot += loss.item() * len(b)
    model.eval()
    with torch.no_grad():
        preds = np.concatenate([
            model(torch.tensor(X[te[i:i + 256]])).numpy() for i in range(0, len(te), 256)])
    rho, _ = spearmanr(preds, yz[te]); r, _ = pearsonr(preds, yz[te])
    return model, {"spearman": round(float(rho), 3), "pearson": round(float(r), 3),
                   "n_train": len(tr), "n_test": len(te)}, (mu, sd)
