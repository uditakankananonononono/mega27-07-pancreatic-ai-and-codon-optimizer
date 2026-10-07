import json,hashlib
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr,pearsonr
ROOT=Path(__file__).resolve().parents[1]
def test_transfer_counts_and_metrics():
 j=json.loads((ROOT/'results/gse63789_transfer_audit.json').read_text());rr=j['records']
 assert j['plan_sha256']==hashlib.sha256((ROOT/'scripts/gse63789_transfer_plan.json').read_bytes()).hexdigest()
 assert j['training_genes']==5068 and len(rr)==5110
 assert sum(r['source_overlap'] for r in rr)==4755
 assert j['excluded']=={'length_mismatch':32,'unmatched':83}
 for name in j['models']:
  for label in ['all','overlap','target_only']:
   sub=[r for r in rr if label=='all' or r['source_overlap']==(label=='overlap')]
   y=[r['target'] for r in sub];p=[r['predictions'][name] for r in sub];s=j['models'][name][label]
   assert s['n']==len(sub)
   assert np.isclose(s['spearman'],spearmanr(y,p).statistic)
   assert np.isclose(s['pearson'],pearsonr(y,p).statistic)
