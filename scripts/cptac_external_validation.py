"""External CPTAC mutation-inference audit. No external-label tuning.
Transparent logistic pathway/burden baseline, not the existing GNN.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from pancreatic_ai import features as F
def main():
 train=ROOT/'data/paad';test=ROOT/'data/paad_cptac_2021';a=json.loads((train/'mutations.json').read_text());b=json.loads((test/'mutations.json').read_text());sa=json.loads((train/'samples.json').read_text());sb=json.loads((test/'samples.json').read_text());assert not set(sa)&set(sb)
 out={'source':'https://www.cbioportal.org/study/summary?id=paad_cptac_2021','train':'TCGA-PAAD','test':'CPTAC-PDAC','n_train':len(sa),'n_test':len(sb),'model':'fixed StandardScaler + balanced LogisticRegression C=1, max_iter=2000; pathway/burden only','target_definition':'nonsynonymous mutation only, does not include deletion or methylation','sample_identifier_overlap':0,'rows':[],'limits':['Distinct study identifiers are not a patient-identity overlap audit','Variant calling, tumor purity and platforms differ across cohorts','Not validation of the previously fitted GNN/CNN','Features contain burden, so association can reflect mutation load','This is genomic inference, not diagnosis or treatment selection']}
 rng=np.random.default_rng(741)
 for target in F.DRIVERS:
  ya=np.array([F.build_targets(a,sa)[s][target] for s in sa]);yb=np.array([F.build_targets(b,sb)[s][target] for s in sb]);xa=F.extra_features(a,sa,exclude=[target]);xb=F.extra_features(b,sb,exclude=[target]);m=make_pipeline(StandardScaler(),LogisticRegression(C=1,class_weight='balanced',max_iter=2000,random_state=741));m.fit(xa,ya);p=m.predict_proba(xb)[:,1];boot=[]
  for _ in range(2000):
   ix=np.r_[rng.choice(np.where(yb==0)[0],sum(yb==0),replace=True),rng.choice(np.where(yb==1)[0],sum(yb==1),replace=True)];boot.append(roc_auc_score(yb[ix],p[ix]))
  row={'target':target,'train_positive':int(ya.sum()),'test_positive':int(yb.sum()),'auc':float(roc_auc_score(yb,p)),'auc95_conditional_stratified_bootstrap':np.quantile(boot,[.025,.975]).tolist(),'ap':float(average_precision_score(yb,p)),'prevalence_baseline_ap':float(yb.mean()),'brier':float(brier_score_loss(yb,p))};out['rows'].append(row);print(row,flush=True)
 out['inputs_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [train/'samples.json',train/'mutations.json',test/'samples.json',test/'mutations.json']}
 (ROOT/'results/cptac_external_validation.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
