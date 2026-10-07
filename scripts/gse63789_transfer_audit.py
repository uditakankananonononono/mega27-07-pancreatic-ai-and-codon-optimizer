"""Frozen source-only fit and raw-count-ratio rank transfer, not calibration."""
import argparse,gzip,hashlib,itertools,json,re
from pathlib import Path
from collections import Counter
import numpy as np
from Bio import SeqIO
from scipy.stats import pearsonr,spearmanr
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import RidgeCV
ROOT=Path(__file__).resolve().parents[1]
def values(p):return {g:float(v) for g,v in (s.split() for s in gzip.open(p,'rt'))}
def stats(y,p):
 if len(y)<3:return {'n':len(y),'pearson':None,'spearman':None}
 return {'n':len(y),'pearson':float(pearsonr(y,p).statistic),'spearman':float(spearmanr(y,p).statistic)}
def compute(source):
 planpath=ROOT/'scripts/gse63789_transfer_plan.json';plan=json.loads(planpath.read_bytes());source=Path(source)
 if hashlib.sha256(source.read_bytes()).hexdigest()!=plan['source_sha256']:raise ValueError('source digest')
 D=ROOT/'data/yeast';rpf=values(D/'GSE75897_RPF_RPKMs.txt.gz');rna=values(D/'GSE75897_RiboZero_RPKMs.txt.gz');cds={}
 for r in SeqIO.parse(gzip.open(D/'GCF_000146045.2_R64_cds_from_genomic.fna.gz','rt'),'fasta'):
  m=re.search(r'\[locus_tag=([^\]]+)\]',r.description)
  if m:cds.setdefault(m[1],[]).append(str(r.seq).upper())
 codons=[''.join(x) for x in itertools.product('ACGT',repeat=3)];feat={}
 for g,ss in cds.items():
  if len(ss)!=1:continue
  s=ss[0]
  if len(s)%3 or set(s)-set('ACGT'):continue
  c=Counter(s[i:i+3] for i in range(0,len(s),3));freq=[c[k]/(len(s)/3) for k in codons];simple=[np.log(len(s)),(s.count('G')+s.count('C'))/len(s)]
  feat[g]={'length_gc':simple,'codon_freq':freq,'codon_length_gc':freq+simple}
 train=[g for g in feat if g in rpf and g in rna and min(rpf[g],rna[g])>0];y=np.array([np.log2(rpf[g]/rna[g]) for g in train]);test=[];excluded=Counter()
 with gzip.open(source,'rt') as f:
  next(f)
  for line in f:
   g,L,R,F,_=line.rstrip().split('\t');L,R,F=map(int,[L,R,F])
   if g not in cds:excluded['unmatched']+=1;continue
   if len(cds[g])!=1:excluded['ambiguous']+=1;continue
   if 3*L!=len(cds[g][0]):excluded['length_mismatch']+=1;continue
   if g not in feat:excluded['invalid_sequence']+=1;continue
   if min(R,F)<=0:excluded['nonpositive_counts']+=1;continue
   test.append({'gene':g,'source_overlap':g in train,'target':float(np.log2(F/R)),'RNA':R,'FP':F})
 target=np.array([r['target'] for r in test]);models={};preds={}
 for name in plan['models']:
  X=np.array([feat[g][name] for g in train]);Z=np.array([feat[r['gene']][name] for r in test]);m=make_pipeline(StandardScaler(),RidgeCV(alphas=[.1,1,10,100,1000]));m.fit(X,y);p=m.predict(Z);preds[name]=p
  models[name]={'alpha':float(m[-1].alpha_),'all':stats(target,p)}
  for label,flag in [('overlap',True),('target_only',False)]:
   ii=np.array([r['source_overlap']==flag for r in test]);models[name][label]=stats(target[ii],p[ii])
 for i,r in enumerate(test):r['predictions']={k:float(p[i]) for k,p in preds.items()}
 overlap=[r for r in test if r['source_overlap']];baseline=stats([r['target'] for r in overlap],[np.log2(rpf[r['gene']]/rna[r['gene']]) for r in overlap])
 return {'plan_sha256':hashlib.sha256(planpath.read_bytes()).hexdigest(),'source_sha256':plan['source_sha256'],'training_genes':len(train),'target_genes':len(test),'excluded':dict(excluded),'models':models,'old_endpoint_persistence_overlap':baseline,'records':test,'limits':plan['limits']}
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--source',required=True);x=a.parse_args();j=compute(x.source);(ROOT/'results/gse63789_transfer_audit.json').write_text(json.dumps(j,indent=2)+'\n');print({k:v for k,v in j.items() if k!='records'})
