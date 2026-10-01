"""Patient-grouped simple CDKN2A controls with training-only age preprocessing."""
import hashlib,json,sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss,log_loss
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from pancreatic_ai import features as F

def main():
 raw=(ROOT/'scripts/pancreatic_control_plan.json').read_bytes();muts,samples,_=F.load_raw();status=F.build_targets(muts,samples);y=np.array([status[s]['CDKN2A'] for s in samples]);patient=np.array(['-'.join(s.split('-')[:3]) for s in samples]);unique=np.unique(patient);py=[]
 for pid in unique:
  ids=np.flatnonzero(patient==pid);assert len(set(y[ids]))==1;py.append(y[ids[0]])
 clin=F.clinical_features();clinical=np.array([clin.get(pid,(np.nan,0.)) for pid in patient]);extra=F.extra_features(muts,samples,exclude=['CDKN2A']);rows=[];predictions=[];splits=[]
 for seed in [17,29,43]:
  (a,b),=StratifiedShuffleSplit(1,test_size=.25,random_state=seed).split(unique,py);tr=np.flatnonzero(np.isin(patient,unique[a]));te=np.flatnonzero(np.isin(patient,unique[b]));assert not set(patient[tr])&set(patient[te]);splits.append({'seed':seed,'train_patients':unique[a].tolist(),'test_patients':unique[b].tolist()});C=clinical.copy();median=float(np.nanmedian(C[tr,0]));C[:,0]=np.where(np.isnan(C[:,0]),median,C[:,0]);views={'clinical_age_sex':C,'burden_pathways':extra,'combined':np.c_[C,extra]}
  for view,X in views.items():
   scaler=StandardScaler().fit(X[tr]);model=LogisticRegression(C=1,class_weight='balanced',max_iter=2000,random_state=seed).fit(scaler.transform(X[tr]),y[tr]);p=model.predict_proba(scaler.transform(X[te]))[:,1]
   row={'seed':seed,'view':view,'n_train':len(tr),'n_test':len(te),'n_train_patients':len(a),'n_test_patients':len(b),'test_positives':int(y[te].sum()),'train_age_median':median,'auc':float(roc_auc_score(y[te],p)),'average_precision':float(average_precision_score(y[te],p)),'brier':float(brier_score_loss(y[te],p)),'logloss':float(log_loss(y[te],p))};rows.append(row);print(row,flush=True)
   for i,v in zip(te,p):predictions.append({'seed':seed,'view':view,'sample':samples[i],'patient':patient[i],'target':int(y[i]),'prediction':float(v)})
 ledger=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in predictions);(ROOT/'results/pancreatic_control_predictions.jsonl').write_text(ledger)
 out={'plan_sha256':hashlib.sha256(raw).hexdigest(),'n_samples':len(samples),'n_patients':len(unique),'positive_samples':int(y.sum()),'missing_age':int(np.isnan(clinical[:,0]).sum()),'clinical_second_feature':'is_male, not smoking','rows':rows,'splits':splits,'prediction_sha256':hashlib.sha256(ledger.encode()).hexdigest(),'mutation_input_sha256':hashlib.sha256((F.DATA/'mutations.json').read_bytes()).hexdigest(),'limits':json.loads(raw)['limits']};(ROOT/'results/pancreatic_control_audit.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
