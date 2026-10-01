import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def test_complete_refit_grid_and_ledger():
 j=json.loads((ROOT/'results/nested_partition_audit.json').read_text());p=(ROOT/'scripts/nested_partition_plan.json').read_bytes();assert j['plan_sha256']==hashlib.sha256(p).hexdigest()
 assert [r['seed'] for r in j['records']]==[17,29,43] and j['n_loci']==3736 and j['n_rows']==3745
 raw=(ROOT/'results/nested_partition_predictions.jsonl').read_bytes();assert hashlib.sha256(raw).hexdigest()==j['prediction_ledger_sha256']
 rows=[json.loads(x) for x in raw.splitlines()];assert len(rows)==3*3736
 for result in j['records']:
  rr=[r for r in rows if r['seed']==result['seed']];assert len({r['locus_tag'] for r in rr})==3736
  assert len(result['folds'])==15
  for view,m in result['summary'].items():
   y=np.array([r['log10_ppm'] for r in rr]);v=np.array([r['predictions'][view] for r in rr]);assert np.isclose(np.corrcoef(y,v)[0,1],m['r']);assert np.isclose(np.sqrt(np.mean((y-v)**2)),m['rmse'])
  assert np.isclose(result['r_delta_combined_minus_codon'],result['summary']['codon_length_gc']['r']-result['summary']['codon']['r'])
