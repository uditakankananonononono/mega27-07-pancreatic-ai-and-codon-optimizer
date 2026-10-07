import json,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from pancreatic_probability_audit import compute

def test_proper_scores_and_rank_invariance():
 j=json.loads((ROOT/'results/pancreatic_probability_audit.json').read_text());assert compute()==j
 old=json.loads((ROOT/'results/pancreatic_control_audit.json').read_text());assert len(j['rows'])==9
 for r,p in zip(j['rows'],old['rows']):
  assert r['seed']==p['seed'] and r['view']==p['view']
  for metric in ['auc','average_precision','brier','logloss']:assert r['scores']['balanced'][metric]==pytest.approx(p[metric])
  for metric in ['auc','average_precision']:assert r['scores']['train_prior_odds_shift'][metric]==pytest.approx(r['scores']['balanced'][metric])
  assert r['scores']['train_prevalence_constant']['auc']==.5
  assert r['scores']['balanced']['brier']>r['scores']['train_prevalence_constant']['brier']
  assert r['scores']['balanced']['logloss']>r['scores']['train_prevalence_constant']['logloss']
  assert r['train_prevalence']==26/138 and r['test_prevalence']==9/46
