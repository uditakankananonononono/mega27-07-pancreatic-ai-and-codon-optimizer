"""Models: 1D-CNN over the gene x mutation-type tensor, GNN over the
co-mutation graph, and a logistic-regression baseline."""
import numpy as np
import torch
import torch.nn as nn
from sklearn.linear_model import LogisticRegression


class MutationCNN(nn.Module):
    """1D CNN across the gene axis of the (genes, 3) mutation-type tensor,
    with clinical covariates fused at the readout."""

    def __init__(self, n_genes=300, in_ch=3, clin_dim=2, hid=32):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv1d(in_ch, hid, kernel_size=5, padding=2), nn.ReLU(),
            nn.Conv1d(hid, hid, kernel_size=5, padding=2, stride=2), nn.ReLU(),
            nn.AdaptiveAvgPool1d(4),
        )
        self.head = nn.Sequential(
            nn.Linear(hid * 4 + clin_dim, 32), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(32, 1),
        )

    def forward(self, x, clin):
        h = self.conv(x.transpose(1, 2)).flatten(1)
        return self.head(torch.cat([h, clin], dim=1)).squeeze(-1)


class CoMutationGNN(nn.Module):
    """GCN over a fixed gene co-mutation graph; per-sample node features are
    the 3 mutation-type counts; clinical covariates fused at readout."""

    def __init__(self, adj, in_dim=3, clin_dim=2, hid=24):
        super().__init__()
        self.register_buffer("adj", torch.tensor(adj, dtype=torch.float32))
        self.lin1 = nn.Linear(in_dim, hid)
        self.lin2 = nn.Linear(hid, hid)
        self.head = nn.Sequential(
            nn.Linear(hid + clin_dim, 32), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(32, 1),
        )

    def forward(self, x, clin):
        h = torch.relu(self.lin1(torch.einsum("ij,bjk->bik", self.adj, x)))
        h = torch.relu(self.lin2(torch.einsum("ij,bjk->bik", self.adj, h)))
        h = h.mean(dim=1)
        return self.head(torch.cat([h, clin], dim=1)).squeeze(-1)


def comutation_adjacency(X, top_k=100, min_jaccard=0.05):
    """Binary adjacency over the top_k genes from co-mutation Jaccard overlap."""
    binary = (X[:, :top_k, :].sum(axis=2) > 0).astype(np.float32)
    inter = binary.T @ binary
    union = binary.sum(0)[:, None] + binary.sum(0)[None, :] - inter
    j = np.divide(inter, np.maximum(union, 1e-9))
    adj = (j >= min_jaccard).astype(np.float32)
    np.fill_diagonal(adj, 1.0)
    deg = adj.sum(1)
    d_inv_sqrt = np.power(np.maximum(deg, 1e-9), -0.5)
    return d_inv_sqrt[:, None] * adj * d_inv_sqrt[None, :]  # symmetric norm


def baseline_logreg():
    return LogisticRegression(max_iter=2000, C=0.5, class_weight="balanced")
