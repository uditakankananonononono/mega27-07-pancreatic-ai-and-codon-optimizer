"""Saved-probability scoring with a training-prior-only odds heuristic."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.special import expit,logit
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss,log_loss
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from pancreatic_ai import features as F

def compute():
 raw=(ROOT/'scripts/pancreatic_probability_plan.json').read_bytes();p=json.loads(raw);old=json.loads((ROOT/'results/pancreatic_control_audit.json').read_text());path=ROOT/'results/pancreatic_control_predictions.jsonl';b=path.read_bytes()
 if hashlib.sha256(b).hexdigest()!=old['prediction_sha256'] or hashlib.sha256((F.DATA/'mutations.json').read_bytes()).hexdigest()!=old['mutation_input_sha256']:raise ValueError('source hash')
 muts,samples,_=F.load_raw();target=F.build_targets(muts,samples);labels={s:int(target[s]['CDKN2A']) for s in samples};patients={s:'-'.join(s.split('-')[:3]) for s in samples};prior=[json.loads(s) for s in b.decode().splitlines()];rows=[]
 def score(y,v):return {'auc':float(roc_auc_score(y,v)),'average_precision':float(average_precision_score(y,v)),'brier':float(brier_score_loss(y,v)),'logloss':float(log_loss(y,v)),'mean_probability':float(np.mean(v))}
 for split in old['splits']:
  seed=split['seed'];tr=set(split['train_patients']);te=set(split['test_patients']);assert not tr&te
  train=[s for s in samples if patients[s] in tr];test=[s for s in samples if patients[s] in te]
  if len(train)!=len(tr) or len(test)!=len(te) or len(train)+len(test)!=len(samples):raise ValueError('patient/sample geometry')
  prevalence=sum(labels[s] for s in train)/len(train)
  for view in ['clinical_age_sex','burden_pathways','combined']:
   rr=[r for r in prior if r['seed']==seed and r['view']==view]
   if {r['sample'] for r in rr}!=set(test) or len(rr)!=len(test):raise ValueError('test identities')
   if any(r['patient']!=patients[r['sample']] or r['target']!=labels[r['sample']] or not 0<r['prediction']<1 for r in rr):raise ValueError('prediction identity/range')
   y=np.array([r['target'] for r in rr]);pred=np.array([r['prediction'] for r in rr]);shift=expit(logit(pred)+logit(prevalence));scores={'balanced':score(y,pred),'train_prevalence_constant':score(y,np.full(len(y),prevalence)),'train_prior_odds_shift':score(y,shift)}
   rows.append({'seed':seed,'view':view,'n_train':len(train),'n_test':len(test),'train_prevalence':prevalence,'test_prevalence':float(np.mean(y)),'scores':scores})
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'prediction_sha256':hashlib.sha256(b).hexdigest(),'rows':rows,'limits':p['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/pancreatic_probability_audit.json').write_text(json.dumps(j,indent=2)+'\n')
 for r in j['rows']:print(r['seed'],r['view'],r['train_prevalence'],{k:(round(v['brier'],6),round(v['logloss'],6)) for k,v in r['scores'].items()})
