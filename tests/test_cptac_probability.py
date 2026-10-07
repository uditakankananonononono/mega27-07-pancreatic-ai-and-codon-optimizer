import hashlib,json
from pathlib import Path
import numpy as np
from sklearn.metrics import brier_score_loss,log_loss,roc_auc_score
ROOT=Path(__file__).resolve().parents[1]
def test_external_scores_and_controls():
 j=json.loads((ROOT/'results/cptac_probability_audit.json').read_text())
 assert j['plan_sha256']==hashlib.sha256((ROOT/'scripts/cptac_probability_plan.json').read_bytes()).hexdigest()
 assert len(j['records'])==560
 for r in j['rows']:
  records=[x for x in j['records'] if x['target']==r['target']];y=[x['label'] for x in records]
  for key,col in [('balanced','balanced'),('train_prior_odds_shift','shift')]:
   p=[x[col] for x in records]
   assert np.isclose(brier_score_loss(y,p),r['scores'][key]['brier'])
   assert np.isclose(log_loss(y,p),r['scores'][key]['logloss'])
  assert r['scores']['balanced']['auc']==r['scores']['train_prior_odds_shift']['auc']
  wins=r['target']=='TP53'
  assert (r['scores']['balanced']['brier']<r['scores']['train_prior_constant']['brier'])==wins
  assert (r['scores']['train_prior_odds_shift']['logloss']<r['scores']['train_prior_constant']['logloss'])==wins
