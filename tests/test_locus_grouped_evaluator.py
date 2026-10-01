import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from locus_grouped_evaluator_audit import grouped_folds

def test_grouped_inner_and_outer_example():
 groups=np.array(['a','a','b','c','d','e','f','g'])
 for tr,te in grouped_folds(groups,7,3):
  assert not set(groups[tr])&set(groups[te])
  for a,b in grouped_folds(groups[tr],0,2):assert not set(groups[tr][a])&set(groups[tr][b])

def test_actual_grouped_ledger_and_duplicate_audit():
 j=json.loads((ROOT/'results/locus_grouped_evaluator_audit.json').read_text());raw=(ROOT/'results/locus_grouped_predictions.jsonl').read_bytes();rows=[json.loads(x) for x in raw.splitlines()]
 assert j['n_rows']==3745 and j['n_loci']==len(rows)==3736
 assert len(j['duplicate_loci'])==8 and sum(r['rows']-1 for r in j['duplicate_loci'])==9 and j['leaked_locus_count_old_outer']==6
 assert all(r['shared_target'] and r['distinct_sequences']==r['rows'] for r in j['duplicate_loci'])
 assert len({r['locus_tag'] for r in rows})==3736 and sum(r['cds_row_count'] for r in rows)==3745
 assert hashlib.sha256(raw).hexdigest()==j['prediction_ledger_sha256']
 y=np.array([r['log10_ppm'] for r in rows])
 for k,v in j['summary'].items():
  pred=np.array([r['oof_prediction_mean'][k] for r in rows]);assert np.isclose(np.corrcoef(y,pred)[0,1],v['r'])
 assert j['delta']['conditional_paired_delta95'][0]>0
