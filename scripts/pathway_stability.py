"""Repeated train/test CDKN2A pathway-group permutation audit.
Training-only gene selection and age imputation; group permutations preserve
within-pathway covariance. Noncausal, curated membership, short-model audit.
"""
import sys,json,itertools
from pathlib import Path
import numpy as np
import torch
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from pancreatic_ai import features as f
from pancreatic_ai.models import CoMutationGNN,comutation_adjacency
torch.set_num_threads(1)
def main():
 muts,samples,_=f.load_raw();status=f.build_targets(muts,samples);y=np.array([status[s]['CDKN2A'] for s in samples]);clin=f.clinical_features();sp=f.sample_patient_map()
 C=np.array([clin.get(sp[s],(np.nan,0)) for s in samples],dtype=np.float32)
 extras=f.extra_features(muts,samples,exclude=['CDKN2A']);names=list(f.PATHWAYS);rows=[]
 for seed in [11,23,42,71,89,101]:
  torch.manual_seed(seed);rng=np.random.default_rng(seed);(tr,te),=StratifiedShuffleSplit(1,test_size=.25,random_state=seed).split(np.zeros(len(y)),y)
  train_samples={samples[i] for i in tr};training_muts=[m for m in muts if m['sampleId'] in train_samples]
  _,genes=f.build_matrix(training_muts,[samples[i] for i in tr],n_genes=100,exclude=['CDKN2A'])
  # Build the selected training vocabulary without selecting genes from test data.
  X=np.zeros((len(samples),100,3),np.float32);idx={g:i for i,g in enumerate(genes)};si={s:i for i,s in enumerate(samples)}
  for m in muts:
   g=m['hugo'];t=m['mutationType']
   if g in idx and t in f.NONSYNONYMOUS:X[si[m['sampleId']],idx[g],f.TYPE_INDEX['missense' if t in f.MISSENSE else 'truncating' if t in f.TRUNCATING else 'other']]+=1
  X=np.log1p(X);Cc=C.copy();age=np.nanmean(Cc[tr,0]);Cc[:,0]=np.where(np.isnan(Cc[:,0]),age,Cc[:,0])/100.;Cc=np.hstack([Cc,extras])
  Xt=torch.tensor(X);Ct=torch.tensor(Cc);yt=torch.tensor(y,dtype=torch.float32)
  m=CoMutationGNN(comutation_adjacency(X[tr]),clin_dim=Cc.shape[1]);opt=torch.optim.Adam(m.parameters(),lr=.003,weight_decay=.001)
  crit=torch.nn.BCEWithLogitsLoss(pos_weight=torch.tensor((len(tr)-y[tr].sum())/y[tr].sum()))
  for _ in range(120):m.train();opt.zero_grad();loss=crit(m(Xt[tr],Ct[tr]),yt[tr]);loss.backward();opt.step()
  m.eval()
  with torch.no_grad():base=roc_auc_score(y[te],torch.sigmoid(m(Xt[te],Ct[te])).numpy())
  drops=[]
  for pi,name in enumerate(names):
   selected=[idx[g] for g in f.PATHWAYS[name] if g in idx];ds=[]
   for _ in range(10):
    order=rng.permutation(len(te));xp=X[te].copy();cp=Cc[te].copy();xp[:,selected,:]=xp[order][:,selected,:];cp[:,5+pi]=cp[order,5+pi]
    with torch.no_grad():p=torch.sigmoid(m(torch.tensor(xp),torch.tensor(cp))).numpy()
    ds.append(base-roc_auc_score(y[te],p))
   drops.append({'pathway':name,'selected_genes':[genes[k] for k in selected],'mean_auc_drop':float(np.mean(ds)),'repeat_sd':float(np.std(ds,ddof=1))})
  rows.append({'split_seed':seed,'test_auc':float(base),'drops':drops});print('seed',seed,'AUC',base,flush=True)
 values=np.array([[d['mean_auc_drop'] for d in r['drops']] for r in rows]);cor=[]
 for a,b in itertools.combinations(range(len(rows)),2):
  rho=float(spearmanr(values[a],values[b]).statistic);cor.append({'seeds':[rows[a]['split_seed'],rows[b]['split_seed']],'rho':rho if np.isfinite(rho) else None})
 out={'target':'CDKN2A','design':'6 stratified 75/25 splits, 120 epochs, 10 grouped permutations; train-only vocabulary, adjacency and age imputation','pathways':names,'rows':rows,'rank_correlations':cor,'mean_auc_drop':dict(zip(names,values.mean(0).tolist())),'limits':['Overlapping test sets: split replicates are not independent patient cohorts','Permutation is association, not causal explanation','Hand-curated pathways and fixed model hyperparameters; no new pathway discovery','Model excludes target gene but global burden remains a confounder']}
 (ROOT/'results/pathway_stability.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['mean_auc_drop'],indent=2))
if __name__=='__main__':main()
