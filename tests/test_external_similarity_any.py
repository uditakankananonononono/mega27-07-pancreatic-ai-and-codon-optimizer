import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from external_similarity_any_audit import aggregate
def test_all_candidate_threshold_ledger():
 raw=(ROOT/'results/external_similarity_any.jsonl').read_bytes();j=json.loads((ROOT/'results/external_similarity_any_audit.json').read_text());rows=[json.loads(x) for x in raw.splitlines()];assert hashlib.sha256(raw).hexdigest()==j['ledger_sha256'];assert len(rows)==6348
 assert aggregate(rows)=={k:j[k] for k in aggregate(rows)}
 prior=[json.loads(x) for x in (ROOT/'results/external_similarity.jsonl').read_text().splitlines()]
 for r,p in zip(rows,prior):
  assert r['name']==p['name'] and r['best']==p['best'] and r['prior_best_flag']==p['similarity_flag']
  assert r['candidate_count']==len(r['aligned_candidates'])<=3
  assert r['similarity_flag']==any(v['identity']>=.4 and min(v['construct_coverage'],v['training_coverage'])>=.7 for v in r['aligned_candidates'])
 assert j['flag_disagreement_names']==[] and j['flagged']['n']==506
 assert hashlib.sha256((ROOT/'scripts/external_similarity_any_plan.json').read_bytes()).hexdigest()==j['plan_sha256']
def test_highest_product_need_not_meet_both_thresholds():
 values=[(.9,.5),(.4,.7)];assert max(values,key=lambda v:v[0]*v[1])==(.9,.5)
 assert any(i>=.4 and c>=.7 for i,c in values)
