"""Separate alternate CDS rows by locus in both outer and inner validation."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
from sklearn.linear_model import RidgeCV
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
from oof_evaluator_audit import freq,metrics
from codon_optimizer.expr_data import build_expression_dataset

def grouped_folds(groups,seed,k):
 unique=np.unique(groups);unique=np.random.default_rng(seed).permutation(unique)
 for test in np.array_split(unique,k):
  mask=np.isin(groups,test);yield np.flatnonzero(~mask),np.flatnonzero(mask)

def main():
 raw=(ROOT/'scripts/nested_partition_plan.json').read_bytes();plan=json.loads(raw)
 rows=[(lt,seq,logab) for lt,gene,seq,logab in build_expression_dataset() if len(seq)>=90];groups=np.array([r[0] for r in rows]);y=np.array([r[2] for r in rows]);c=np.array([freq(r[1]) for r in rows]);simple=np.array([[np.log10(len(r[1])),sum(b in 'GC' for b in r[1])/len(r[1])] for r in rows]);views={'codon':c,'length_gc':simple,'codon_length_gc':np.c_[c,simple]};unique=np.unique(groups);records=[];ledger=[]
 for seed in plan['outer_seeds']:
  pred={k:np.zeros(len(y)) for k in views};folds=[];row_outer=np.zeros(len(y),dtype=int)
  for fold,(tr,te) in enumerate(grouped_folds(groups,seed,5)):
   assert not set(groups[tr])&set(groups[te]);inner=list(grouped_folds(groups[tr],seed*100+fold,3));row_outer[te]=fold
   for a,b in inner:assert not set(groups[tr][a])&set(groups[tr][b])
   for name,X in views.items():
    m=RidgeCV(alphas=np.array(plan['alphas']),cv=inner).fit(X[tr],y[tr]);pred[name][te]=m.predict(X[te]);folds.append({'fold':fold,'view':name,'alpha':float(m.alpha_),**metrics(y[te],pred[name][te])})
  ly=[];lp={k:[] for k in views}
  for tag in unique:
   ids=np.flatnonzero(groups==tag);assert len(set(y[ids]))==1 and len(set(row_outer[ids]))==1
   ly.append(float(y[ids[0]]));ps={k:float(v[ids].mean()) for k,v in pred.items()}
   for k,v in ps.items():lp[k].append(v)
   ledger.append({'seed':seed,'locus_tag':str(tag),'outer_fold':int(row_outer[ids[0]]),'log10_ppm':ly[-1],'predictions':ps})
  summary={k:metrics(np.array(ly),np.array(v)) for k,v in lp.items()};result={'seed':seed,'summary':summary,'r_delta_combined_minus_codon':summary['codon_length_gc']['r']-summary['codon']['r'],'folds':folds};records.append(result);print(result,flush=True)
 text=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in ledger);(ROOT/'results/nested_partition_predictions.jsonl').write_text(text)
 out={'plan_sha256':hashlib.sha256(raw).hexdigest(),'n_rows':len(rows),'n_loci':len(unique),'records':records,'prediction_ledger_sha256':hashlib.sha256(text.encode()).hexdigest(),'limits':plan['limits']};(ROOT/'results/nested_partition_audit.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
