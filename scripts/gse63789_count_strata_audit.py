"""Frozen descriptive RNA count bins of archived source-only predictions."""
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1]
def compute():
 p=ROOT/'scripts/gse63789_count_strata_plan.json';plan=json.loads(p.read_text());f=ROOT/plan['input'];j=json.loads(f.read_text());rows=sorted(j['records'],key=lambda r:(r['RNA'],r['gene']));assert len(rows)==5110;names=list(rows[0]['predictions']);groups=np.array_split(np.arange(len(rows)),4);bins=[]
 for n,idx in enumerate(groups,1):
  rs=[rows[i] for i in idx];bins.append({'quartile':n,'n':len(rs),'RNA_min':rs[0]['RNA'],'RNA_max':rs[-1]['RNA'],'source_overlap_n':sum(r['source_overlap'] for r in rs),'models':{name:float(spearmanr([r['target'] for r in rs],[r['predictions'][name] for r in rs]).statistic) for name in names}})
 return {'plan_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bins':bins,'component_correlations':{name:{key:float(spearmanr([r[key] for r in rows],[r['predictions'][name] for r in rows]).statistic) for key in ['RNA','FP']} for name in names},'limits':plan['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/gse63789_count_strata_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j)
