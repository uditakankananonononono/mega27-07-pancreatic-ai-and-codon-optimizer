import pytest
from codon_optimizer.expr_data import build_expression_dataset
from collections import Counter
def test_required_policy():
 with pytest.raises(TypeError):build_expression_dataset()
def test_invalid_policy():
 with pytest.raises(ValueError):build_expression_dataset(policy='longest')
def test_unique_quarantines_all_repeats():
 raw=build_expression_dataset(policy='raw');unique=build_expression_dataset(policy='unique');c=Counter(r[0] for r in raw)
 assert len(raw)==3746 and len(unique)==3729
 assert unique==[r for r in raw if c[r[0]]==1]
 assert len({r[0] for r in unique})==len(unique)
