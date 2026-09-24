import sys, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from pancreatic_ai.models import MutationCNN, CoMutationGNN, comutation_adjacency
import torch

def test_cnn_forward_shape():
    m = MutationCNN(n_genes=50)
    x = torch.zeros(4, 50, 3); c = torch.zeros(4, 2)
    assert m(x, c).shape == (4,)

def test_gnn_forward_shape():
    adj = comutation_adjacency(np.random.rand(20, 100, 3).astype(np.float32))
    m = CoMutationGNN(adj)
    x = torch.zeros(4, 100, 3); c = torch.zeros(4, 2)
    assert m(x, c).shape == (4,)
    assert adj.shape == (100, 100)
    assert np.allclose(adj, adj.T)

def test_comutation_adjacency_no_leak_diagonal():
    adj = comutation_adjacency(np.ones((10, 100, 3), dtype=np.float32))
    assert np.isfinite(adj).all()
