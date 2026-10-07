"""Candidate-bounded protein similarity and selection-sensitivity screen."""
import hashlib,json,sys
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np,openpyxl
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner,substitution_matrices
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from codon_optimizer.expr_data import build_expression_dataset
def kmers(p):return {p[i:i+5] for i in range(max(0,len(p)-4)) if '*' not in p[i:i+5]}
def aggregate(rows):
 r=[x for x in rows if not x['protein_containment_overlap']];flag=[x for x in r if x['similarity_flag']];other=[x for x in r if not x['similarity_flag']]
 def summary(rr):return {'n':len(rr),'rho':float(spearmanr([x['prediction'] for x in rr],[x['E'] for x in rr]).statistic) if len(rr)>2 else None,'category_counts':dict(sorted(Counter(str(x['E']) for x in rr).items()))}
 return {'retained':summary(r),'flagged':summary(flag),'unflagged':summary(other),'all_constructs_flagged':sum(x['similarity_flag'] for x in rows),'maximum_construct_5mer_overlap_expression_rho':float(spearmanr([x['max_construct_5mer_fraction'] for x in r],[x['E'] for x in r]).statistic)}
def compute():
 raw=(ROOT/'scripts/external_similarity_top10_plan.json').read_bytes();prior=[json.loads(x) for x in (ROOT/'results/external_transfer_predictions.jsonl').read_text().splitlines()];f=ROOT/'data/external/nature16509_data2.xlsx';assert hashlib.sha256(f.read_bytes()).hexdigest()=='7e32a22549f4565a08284e64f4187c5562181b909bfb382849cd9227fe4ef978';sheet=list(openpyxl.load_workbook(f,read_only=True,data_only=True)['E6348.csv'].values);seqs={r[0]:str(Seq(r[-1]).translate()).rstrip('*') for r in sheet[1:]};train=[(lt,str(Seq(s).translate()).rstrip('*')) for lt,g,s,y in build_expression_dataset() if len(s)>=90];inv=defaultdict(list)
 for i,(lt,p) in enumerate(train):
  for k in kmers(p):inv[k].append(i)
 align=PairwiseAligner();align.mode='local';align.substitution_matrix=substitution_matrices.load('BLOSUM62');align.open_gap_score=-10;align.extend_gap_score=-.5;rows=[]
 for ix,r in enumerate(prior):
  p=seqs[r['name']];p=p.replace('*','X');counts=Counter(j for k in kmers(p) for j in inv.get(k,[]));candidates=sorted(counts,key=lambda j:(-counts[j],j))[:10];best=None;aligned=[]
  for j in candidates:
   q=train[j][1].replace('*','X');a=align.align(p,q)[0];pairs=list(zip(a.aligned[0],a.aligned[1]));nongap=sum(int(e-s) for (s,e),_ in pairs);matched=sum(sum(x==y for x,y in zip(p[s:e],q[t:u])) for (s,e),(t,u) in pairs);ident=matched/nongap if nongap else 0;covp=nongap/len(p);covq=nongap/len(q);v={'training_locus':train[j][0],'identity':ident,'construct_coverage':covp,'training_coverage':covq,'shared_unique_5mers':counts[j]}
   aligned.append(v)
   if best is None or ident*min(covp,covq)>best['identity']*min(best['construct_coverage'],best['training_coverage']):best=v
  rows.append({'name':r['name'],'E':r['expression_category'],'prediction':r['natural_abundance_ridge_prediction'],'protein_containment_overlap':r['protein_containment_overlap'],'max_construct_5mer_fraction':max(counts.values(),default=0)/max(1,len(kmers(p))),'candidate_count':len(candidates),'best':best,'similarity_flag':any(v['identity']>=.4 and min(v['construct_coverage'],v['training_coverage'])>=.7 for v in aligned),'prior_best_flag':bool(best and best['identity']>=.4 and min(best['construct_coverage'],best['training_coverage'])>=.7),'aligned_candidates':aligned})
  if ix%1000==0:print('processed',ix,flush=True)
 old={r['name']:r for r in map(json.loads,(ROOT/'results/external_similarity_any.jsonl').read_text().splitlines())}
 text=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in rows);return aggregate(rows)|{'top3_vs_top10_flag_changes':[r['name'] for r in rows if r['similarity_flag']!=old[r['name']]['similarity_flag']],'prior_best_flag_note':'best-versus-any disagreement within top10, not top3 comparison','flag_disagreement_names':[r['name'] for r in rows if r['prior_best_flag']!=r['similarity_flag']],'plan_sha256':hashlib.sha256(raw).hexdigest(),'ledger_sha256':hashlib.sha256(text.encode()).hexdigest(),'limits':json.loads(raw)['limits'],'selection':'Top10 candidate expansion; same frozen predictions, no E-dependent fitting.'},text
if __name__=='__main__':
 j,text=compute();(ROOT/'results/external_similarity_top10_audit.json').write_text(json.dumps(j,indent=2)+'\n');(ROOT/'results/external_similarity_top10.jsonl').write_text(text);print(json.dumps(j,indent=2))
