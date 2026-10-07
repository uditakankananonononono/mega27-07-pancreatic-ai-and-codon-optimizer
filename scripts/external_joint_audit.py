"""Joint descriptive selection and covariate sensitivity; no prediction refit."""
import hashlib,json,sys
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from external_covariate_audit import aggregate

def compute():
 raw=(ROOT/'scripts/external_joint_plan.json').read_bytes();p=json.loads(raw)
 paths=[ROOT/'results/external_covariates.jsonl',ROOT/'results/external_similarity_top10.jsonl']
 a=[json.loads(s) for s in paths[0].read_text().splitlines()];b=[json.loads(s) for s in paths[1].read_text().splitlines()]
 if len({r['name'] for r in a})!=len(a) or len({r['name'] for r in b})!=len(b):raise ValueError('duplicate identity')
 lookup={r['name']:r for r in b}
 if set(lookup)!={r['name'] for r in a}:raise ValueError('identity set')
 for r in a:
  s=lookup[r['name']]
  if any(r[k]!=s[q] for k,q in [('expression_category','E'),('natural_abundance_ridge_prediction','prediction'),('protein_containment_overlap','protein_containment_overlap')]):raise ValueError('ledger mismatch')
  r['similarity_flag']=s['similarity_flag']
 retained=[r for r in a if not r['protein_containment_overlap']]
 subsets={'retained':retained,'top10_flagged':[r for r in retained if r['similarity_flag']],'top10_unflagged':[r for r in retained if not r['similarity_flag']]}
 rows={k:aggregate(v)|{'category_counts':dict(sorted(Counter(str(r['expression_category']) for r in v).items()))} for k,v in subsets.items()}
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'input_sha256':{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},'rows':rows,'limits':p['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/external_joint_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j,indent=2))
