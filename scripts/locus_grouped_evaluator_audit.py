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
 raw=(ROOT/'scripts/locus_grouped_plan.json').read_bytes();rows=[(lt,seq,logab) for lt,gene,seq,logab in build_expression_dataset() if len(seq)>=90];groups=np.array([r[0] for r in rows]);y=np.array([r[2] for r in rows]);c=np.array([freq(r[1]) for r in rows]);simple=np.array([[np.log10(len(r[1])),sum(b in 'GC' for b in r[1])/len(r[1])] for r in rows]);views={'codon':c,'length_gc':simple,'codon_length_gc':np.c_[c,simple]};pred={k:np.zeros(len(y)) for k in views};folds=[];row_outer=np.zeros(len(y),dtype=int)
 old_folds=np.array_split(np.random.default_rng(7).permutation(len(y)),5);old_id=np.zeros(len(y),int)
 for i,te in enumerate(old_folds):old_id[te]=i
 duplicated=[]
 for tag in np.unique(groups):
  ids=np.flatnonzero(groups==tag)
  if len(ids)>1:duplicated.append({'locus':str(tag),'rows':len(ids),'lengths':[len(rows[i][1]) for i in ids],'distinct_sequences':len(set(rows[i][1] for i in ids)),'shared_target':len(set(y[ids]))==1,'historical_outer_folds':old_id[ids].tolist(),'old_locus_leakage':len(set(old_id[ids]))>1})
 for fold,(tr,te) in enumerate(grouped_folds(groups,7,5)):
  assert not set(groups[tr])&set(groups[te]);inner=list(grouped_folds(groups[tr],fold,3));row_outer[te]=fold
  for a,b in inner:assert not set(groups[tr][a])&set(groups[tr][b])
  for name,X in views.items():
   m=RidgeCV(alphas=np.logspace(-3,3,13),cv=inner).fit(X[tr],y[tr]);pred[name][te]=m.predict(X[te]);folds.append({'fold':fold,'view':name,'n_train_rows':len(tr),'n_test_rows':len(te),'n_train_loci':len(set(groups[tr])),'n_test_loci':len(set(groups[te])),'alpha':float(m.alpha_),**metrics(y[te],pred[name][te])});print(fold,name,folds[-1],flush=True)
 unique=np.unique(groups);ly=[];lp={k:[] for k in views};ledger=[]
 for tag in unique:
  ids=np.flatnonzero(groups==tag);assert len(set(y[ids]))==1 and len(set(row_outer[ids]))==1
  ly.append(float(y[ids[0]]));ps={k:float(v[ids].mean()) for k,v in pred.items()}
  for k,v in ps.items():lp[k].append(v)
  ledger.append({'locus_tag':str(tag),'cds_row_count':len(ids),'outer_fold':int(row_outer[ids[0]]),'log10_ppm':ly[-1],'oof_prediction_mean':ps})
 ly=np.array(ly);lp={k:np.array(v) for k,v in lp.items()};rng=np.random.default_rng(1001593);boot={k:[] for k in lp}
 for _ in range(5000):
  ids=rng.integers(0,len(unique),len(unique))
  for k,v in lp.items():boot[k].append(metrics(ly[ids],v[ids])['r'])
 summary={k:{**metrics(ly,v),'conditional_locus_bootstrap_r95':np.percentile(boot[k],[2.5,97.5]).tolist()} for k,v in lp.items()};delta={'r_delta_combined_minus_codon':summary['codon_length_gc']['r']-summary['codon']['r'],'conditional_paired_delta95':np.percentile(np.array(boot['codon_length_gc'])-np.array(boot['codon']),[2.5,97.5]).tolist()}
 text=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in ledger);(ROOT/'results/locus_grouped_predictions.jsonl').write_text(text)
 out={'plan_sha256':hashlib.sha256(raw).hexdigest(),'n_rows':len(rows),'n_loci':len(unique),'duplicate_loci':duplicated,'leaked_locus_count_old_outer':sum(r['old_locus_leakage'] for r in duplicated),'summary':summary,'delta':delta,'folds':folds,'prediction_ledger_sha256':hashlib.sha256(text.encode()).hexdigest(),'limits':json.loads(raw)['limits']};(ROOT/'results/locus_grouped_evaluator_audit.json').write_text(json.dumps(out,indent=2)+'\n');print('GROUPED',summary,delta,'old leaks',out['leaked_locus_count_old_outer'])
if __name__=='__main__':main()
