import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from external_similarity_audit import aggregate
def test_similarity_screen_ledger_frozen_metrics():
 raw=(ROOT/'results/external_similarity.jsonl').read_bytes();j=json.loads((ROOT/'results/external_similarity_audit.json').read_text());assert hashlib.sha256(raw).hexdigest()==j['ledger_sha256'];rows=[json.loads(x) for x in raw.splitlines()];assert len(rows)==6348 and aggregate(rows)=={k:j[k] for k in aggregate(rows)}
 assert j['plan_sha256']==hashlib.sha256((ROOT/'scripts/external_similarity_plan.json').read_bytes()).hexdigest()
 assert j['retained']['n']==6305 and j['flagged']['n']==506 and j['unflagged']['n']==5799
 for r in rows:
  assert 0<=r['max_construct_5mer_fraction']<=1 and r['candidate_count']<=3
  b=r['best'];assert r['similarity_flag']==bool(b and b['identity']>=.4 and min(b['construct_coverage'],b['training_coverage'])>=.7)
  if b:assert all(0<=b[k]<=1 for k in ['identity','construct_coverage','training_coverage'])
