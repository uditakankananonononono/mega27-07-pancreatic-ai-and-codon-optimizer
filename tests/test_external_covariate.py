import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from external_covariate_audit import aggregate,compute
def test_covariate_ledger_and_rank_aggregate():
 raw=(ROOT/'results/external_covariates.jsonl').read_bytes();j=json.loads((ROOT/'results/external_covariate_audit.json').read_text());assert hashlib.sha256(raw).hexdigest()==j['ledger_sha256'];rows=[json.loads(r) for r in raw.splitlines()];a=aggregate(rows)
 assert len(rows)==6348 and a['n']==6305
 for k,v in a.items():assert np.isclose(v,j[k])
 assert j['plan_sha256']==hashlib.sha256((ROOT/'scripts/external_covariate_plan.json').read_bytes()).hexdigest()
 assert all(0<=r['gc_fraction']<=1 and r['nt_length']%3==0 for r in rows)
 if (ROOT/'data/external/nature16509_data2.xlsx').exists():assert compute()[0]==j
