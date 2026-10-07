"""Fixed nested-CV OOF controls: separate composition signal from simple length/GC."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from codon_optimizer.expr_data import build_expression_dataset
CODONS=[a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG' if a+b+c not in ['TAA','TAG','TGA']]
def freq(seq):
 n=len(seq)//3;return np.array([sum(seq[i:i+3]==c for i in range(0,n*3,3))/n for c in CODONS])
def metrics(y,p):return {'r':float(np.corrcoef(y,p)[0,1]),'rmse':float(np.sqrt(np.mean((p-y)**2)))}
def main():
 raw=(ROOT/'scripts/oof_evaluator_plan.json').read_bytes();rows=[(lt,seq,logab) for lt,gene,seq,logab in build_expression_dataset(policy="raw") if len(seq)>=90];y=np.array([r[2] for r in rows]);c=np.array([freq(r[1]) for r in rows]);simple=np.array([[np.log10(len(r[1])),sum(b in 'GC' for b in r[1])/len(r[1])] for r in rows]);views={'codon':c,'length_gc':simple,'codon_length_gc':np.c_[c,simple]};pred={k:np.zeros(len(y)) for k in views};folds=[]
 for fold,(tr,te) in enumerate(KFold(5,shuffle=True,random_state=7).split(y)):
  for name,X in views.items():
   m=RidgeCV(alphas=np.logspace(-3,3,13),cv=KFold(3,shuffle=True,random_state=fold)).fit(X[tr],y[tr]);pred[name][te]=m.predict(X[te]);folds.append({'fold':fold,'view':name,'n_train':len(tr),'n_test':len(te),'alpha':float(m.alpha_),**metrics(y[te],pred[name][te])});print(fold,name,folds[-1],flush=True)
 rng=np.random.default_rng(1001571);boot={k:[] for k in pred}
 for _ in range(5000):
  ids=rng.integers(0,len(y),len(y))
  for k,p in pred.items():boot[k].append(metrics(y[ids],p[ids])['r'])
 summary={k:{**metrics(y,p),'conditional_gene_bootstrap_r95':np.percentile(boot[k],[2.5,97.5]).tolist()} for k,p in pred.items()}
 deltas={k:{'pooled_r_delta_vs_codon':summary[k]['r']-summary['codon']['r'],'paired_conditional_bootstrap_r_delta95':np.percentile(np.array(boot[k])-np.array(boot['codon']),[2.5,97.5]).tolist()} for k in ['length_gc','codon_length_gc']}
 records=[{'locus_tag':row[0],'log10_ppm':row[2],'oof_predictions':{k:float(v[i]) for k,v in pred.items()}} for i,row in enumerate(rows)]
 ledger=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in records);(ROOT/'results/oof_evaluator_predictions.jsonl').write_text(ledger)
 out={'plan_sha256':hashlib.sha256(raw).hexdigest(),'n_genes':len(rows),'input_sha256':hashlib.sha256(json.dumps(rows).encode()).hexdigest(),'prediction_ledger_sha256':hashlib.sha256(ledger.encode()).hexdigest(),'summary':summary,'differences':deltas,'folds':folds,'limits':json.loads(raw)['limits']};(ROOT/'results/oof_evaluator_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(summary,deltas,flush=True)
if __name__=='__main__':main()
