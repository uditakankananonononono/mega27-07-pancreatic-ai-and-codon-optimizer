"""Measured yeast RPF/RNA endpoint; native CDS sequence models, not optimized-design validation."""
import json,gzip,sys,re,itertools,hashlib
from pathlib import Path
import numpy as np
from Bio import SeqIO
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold
from scipy.stats import pearsonr,spearmanr
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data/yeast'
def read(p):return dict((g,float(v)) for g,v in (line.split() for line in gzip.open(p,'rt')))
def main():
 rpf=read(D/'GSE75897_RPF_RPKMs.txt.gz');rna=read(D/'GSE75897_RiboZero_RPKMs.txt.gz');codons=[''.join(x) for x in itertools.product('ACGT',repeat=3)];rows=[];seen=set()
 for r in SeqIO.parse(gzip.open(D/'GCF_000146045.2_R64_cds_from_genomic.fna.gz','rt'),'fasta'):
  m=re.search(r'\[locus_tag=([^\]]+)\]',r.description)
  if not m:continue
  g=m.group(1);s=str(r.seq).upper()
  if g in seen or g not in rpf or g not in rna or rpf[g]<=0 or rna[g]<=0 or len(s)%3 or set(s)-set('ACGT'):continue
  seen.add(g);freq=np.array([sum(s[i:i+3]==c for i in range(0,len(s),3)) for c in codons])/ (len(s)/3)
  rows.append((g,freq,[np.log(len(s)),(s.count('G')+s.count('C'))/len(s)],np.log2(rpf[g]/rna[g])))
 X=np.array([x[1] for x in rows]);simple=np.array([x[2] for x in rows]);y=np.array([x[3] for x in rows]);splits=list(KFold(5,shuffle=True,random_state=307).split(y));out={}
 for name,xx in [('length_gc',simple),('codon_freq',X),('codon_length_gc',np.c_[X,simple])]:
  pred=np.zeros(len(y));folds=[]
  for tr,te in splits:
   model=make_pipeline(StandardScaler(),RidgeCV(alphas=[.1,1,10,100,1000]));model.fit(xx[tr],y[tr]);pred[te]=model.predict(xx[te]);folds.append(float(pearsonr(y[te],pred[te]).statistic))
  out[name]={'fold_pearson':folds,'oof_pearson':float(pearsonr(y,pred).statistic),'oof_spearman':float(spearmanr(y,pred).statistic),'oof_mse':float(np.mean((y-pred)**2))}
 result={'dataset':'GSE75897 RPF / RiboZero RPKM ratio','source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75897','n_genes':len(rows),'gene_ids':[x[0] for x in rows],'models':out,'target':'log2 RPF/RNA density ratio, not protein yield or initiation efficiency','seed':307,'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in D.glob('*.gz')},'limits':['Observed native yeast CDS only, not intervention or validation of E. coli optimized designs','One measurement pair, no replicate uncertainty','Homologous genes may cross folds; this is not family-held-out transfer','Ratio shares measurement noise and does not correct elongation-rate differences','No new biological mechanism or expression improvement established']}
 (ROOT/'results/yeast_translation_audit.json').write_text(json.dumps(result,indent=2)+'\n');print('genes',len(rows));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
