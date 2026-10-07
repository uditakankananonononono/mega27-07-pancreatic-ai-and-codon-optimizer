import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from external_similarity_top10_audit import aggregate
def test_top10_ledger_and_top3_comparison():
 raw=(ROOT/'results/external_similarity_top10.jsonl').read_bytes();j=json.loads((ROOT/'results/external_similarity_top10_audit.json').read_text());rows=[json.loads(x) for x in raw.splitlines()];old={r['name']:r for r in map(json.loads,(ROOT/'results/external_similarity_any.jsonl').read_text().splitlines())}
 assert len(rows)==6348 and hashlib.sha256(raw).hexdigest()==j['ledger_sha256']
 assert aggregate(rows)=={k:j[k] for k in aggregate(rows)}
 changed=[]
 for r in rows:
  assert 1<=r['candidate_count']<=10 and r['candidate_count']==len(r['aligned_candidates'])
  assert r['similarity_flag']==any(v['identity']>=.4 and min(v['construct_coverage'],v['training_coverage'])>=.7 for v in r['aligned_candidates'])
  assert r['aligned_candidates'][:old[r['name']]['candidate_count']]==old[r['name']]['aligned_candidates']
  if old[r['name']]['similarity_flag']:assert r['similarity_flag']
  if r['similarity_flag']!=old[r['name']]['similarity_flag']:changed.append(r['name'])
 assert changed==j['top3_vs_top10_flag_changes'] and len(changed)==18
 assert j['flagged']['n']==524 and j['unflagged']['n']==5781 and j['flag_disagreement_names']==['CtR128']
 assert hashlib.sha256((ROOT/'scripts/external_similarity_top10_plan.json').read_bytes()).hexdigest()==j['plan_sha256']
