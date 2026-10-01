import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def test_saved_oof_ledger_and_conditional_comparisons():
 j=json.loads((ROOT/'results/oof_evaluator_audit.json').read_text());raw=(ROOT/'results/oof_evaluator_predictions.jsonl').read_bytes();rows=[json.loads(x) for x in raw.splitlines()]
 assert len(rows)==j['n_genes']==3745 and len({r['locus_tag'] for r in rows})==3736
 assert hashlib.sha256(raw).hexdigest()==j['prediction_ledger_sha256']
 assert hashlib.sha256((ROOT/'scripts/oof_evaluator_plan.json').read_bytes()).hexdigest()==j['plan_sha256']
 y=np.array([r['log10_ppm'] for r in rows])
 for view,s in j['summary'].items():
  p=np.array([r['oof_predictions'][view] for r in rows]);assert np.isclose(np.corrcoef(y,p)[0,1],s['r']);assert np.isclose(np.sqrt(np.mean((y-p)**2)),s['rmse'])
 assert len(j['folds'])==15
 old=json.loads((ROOT/'results/nested_cv_eval.json').read_text())['views']['codon_freq_ridge']['per_fold_r']
 new=[r['r'] for r in j['folds'] if r['view']=='codon'];assert np.allclose(old,new,atol=.000051)
 assert j['summary']['length_gc']['r']<.16
 assert 0<j['differences']['codon_length_gc']['paired_conditional_bootstrap_r_delta95'][0]
