"""External proper scores with source-prior controls, no target calibration."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.special import expit,logit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss,log_loss
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from pancreatic_ai import features as F
def score(y,p):return {'auc':float(roc_auc_score(y,p)),'ap':float(average_precision_score(y,p)),'brier':float(brier_score_loss(y,p)),'logloss':float(log_loss(y,p)),'mean_probability':float(p.mean())}
def compute():
 raw=(ROOT/'scripts/cptac_probability_plan.json').read_bytes();old=json.loads((ROOT/'results/cptac_external_validation.json').read_text())
 for p,d in old['inputs_sha256'].items():
  if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=d:raise ValueError('input digest')
 a=json.loads((ROOT/'data/paad/mutations.json').read_text());b=json.loads((ROOT/'data/paad_cptac_2021/mutations.json').read_text());sa=json.loads((ROOT/'data/paad/samples.json').read_text());sb=json.loads((ROOT/'data/paad_cptac_2021/samples.json').read_text());rows=[];records=[]
 for target in F.DRIVERS:
  ya=np.array([F.build_targets(a,sa)[s][target] for s in sa]);yb=np.array([F.build_targets(b,sb)[s][target] for s in sb]);xa=F.extra_features(a,sa,exclude=[target]);xb=F.extra_features(b,sb,exclude=[target]);m=make_pipeline(StandardScaler(),LogisticRegression(C=1,class_weight='balanced',max_iter=2000,random_state=741));m.fit(xa,ya);p=m.predict_proba(xb)[:,1];prior=float(ya.mean());shift=expit(logit(p)+logit(prior));scores={'balanced':score(yb,p),'train_prior_constant':score(yb,np.full(len(yb),prior)),'train_prior_odds_shift':score(yb,shift)};priorrow=next(r for r in old['rows'] if r['target']==target)
  for k in ['auc','ap','brier']:
   if abs(scores['balanced'][k]-priorrow[k])>1e-10:raise ValueError('old replay '+k)
  if abs(scores['balanced']['auc']-scores['train_prior_odds_shift']['auc'])>1e-12:raise ValueError('rank invariance')
  rows.append({'target':target,'n_train':len(sa),'n_test':len(sb),'train_prevalence':prior,'test_prevalence_descriptive':float(yb.mean()),'scores':scores})
  records.extend({'sample':s,'target':target,'label':int(y),'balanced':float(v),'shift':float(w)} for s,y,v,w in zip(sb,yb,p,shift))
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'input_sha256':old['inputs_sha256'],'rows':rows,'records':records,'limits':json.loads(raw)['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/cptac_probability_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j['rows'])
