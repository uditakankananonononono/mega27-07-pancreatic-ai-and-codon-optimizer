import json,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from external_joint_audit import compute

def test_joint_selection_and_covariate_reproduction():
 j=json.loads((ROOT/'results/external_joint_audit.json').read_text());assert compute()==j
 a=json.loads((ROOT/'results/external_similarity_top10_audit.json').read_text())
 for k,q in [('retained','retained'),('top10_flagged','flagged'),('top10_unflagged','unflagged')]:
  assert j['rows'][k]['n']==a[q]['n']
  assert j['rows'][k]['category_counts']==a[q]['category_counts']
  assert j['rows'][k]['raw_prediction_spearman']==pytest.approx(a[q]['rho'])
 assert j['rows']['top10_flagged']['n']+j['rows']['top10_unflagged']['n']==j['rows']['retained']['n']
 old=json.loads((ROOT/'results/external_covariate_audit.json').read_text())
 for k in ['gc_spearman','length_spearman','gc_length_partial_rank_correlation']:assert j['rows']['retained'][k]==pytest.approx(old[k])
