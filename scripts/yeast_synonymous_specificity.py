"""Chromosome-held-out yeast audit separates amino-acid from synonymous signal.
Association only. Chromosomes are not gene-family holdouts.
"""
import gzip,json,re,hashlib,itertools
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.Data import CodonTable
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import RidgeCV
from scipy.stats import pearsonr
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data/yeast'
def abund(p):return dict((g,float(v)) for g,v in (s.split() for s in gzip.open(p,'rt')))
def main():
 codons=sorted(CodonTable.unambiguous_dna_by_id[1].forward_table);table=CodonTable.unambiguous_dna_by_id[1].forward_table;aas=sorted(set(table.values()));rpf=abund(D/'GSE75897_RPF_RPKMs.txt.gz');rna=abund(D/'GSE75897_RiboZero_RPKMs.txt.gz');rows=[];seen=set()
 for r in SeqIO.parse(gzip.open(D/'GCF_000146045.2_R64_cds_from_genomic.fna.gz','rt'),'fasta'):
  g=re.search(r'\[locus_tag=([^\]]+)\]',r.description)
  if not g:continue
  g=g.group(1);s=str(r.seq).upper();chr_=re.search(r'(NC_\d+\.\d+)',r.id)
  if g in seen or g not in rpf or g not in rna or rpf[g]<=0 or rna[g]<=0 or len(s)%3 or set(s)-set('ACGT') or not chr_:continue
  seen.add(g);counts=np.array([sum(s[k:k+3]==c for k in range(0,len(s),3)) for c in codons],float);total=counts.sum();ac=np.array([sum(counts[i] for i,c in enumerate(codons) if table[c]==a) for a in aas]);den={a:ac[i] for i,a in enumerate(aas)}
  aa=ac/total;syn=np.array([counts[i]/den[table[c]] if den[table[c]] else 0. for i,c in enumerate(codons)]);basic=[np.log(len(s)),(s.count('G')+s.count('C'))/len(s)]
  rows.append({'gene':g,'chromosome':chr_.group(1),'aa':aa,'syn':syn,'basic':basic,'y':np.log2(rpf[g]/rna[g])})
 y=np.array([r['y'] for r in rows]);groups=np.array([r['chromosome'] for r in rows]);basic=np.array([r['basic'] for r in rows]);aa=np.array([r['aa'] for r in rows]);syn=np.array([r['syn'] for r in rows]);views={'length_gc':basic,'aa_length_gc':np.c_[aa,basic],'aa_syn_length_gc':np.c_[aa,syn,basic]};splits=list(GroupKFold(5).split(y,y,groups));results={};predictions={}
 for name,x in views.items():
  pred=np.zeros(len(y));folds=[]
  for tr,te in splits:
   m=make_pipeline(StandardScaler(),RidgeCV(alphas=[.1,1,10,100,1000]));m.fit(x[tr],y[tr]);pred[te]=m.predict(x[te]);folds.append({'pearson':float(pearsonr(y[te],pred[te]).statistic),'mse':float(np.mean((y[te]-pred[te])**2)),'test_chromosomes':sorted(set(groups[te])),'n_test':len(te)})
  predictions[name]=pred;results[name]={'folds':folds,'oof_pearson':float(pearsonr(y,pred).statistic),'oof_mse':float(np.mean((y-pred)**2))}
 # Paired chromosome-cluster resampling, conditional on fitted OOF models.
 base=(y-predictions['aa_length_gc'])**2;full=(y-predictions['aa_syn_length_gc'])**2;rng=np.random.default_rng(778);unique=np.unique(groups);diff=[]
 for _ in range(3000):
  ix=np.concatenate([np.where(groups==g)[0] for g in rng.choice(unique,len(unique),replace=True)]);diff.append(float(np.mean(base[ix]-full[ix])))
 out={'n_genes':len(y),'chromosomes':unique.tolist(),'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75897','target':'log2 RPF/RNA','design':'five chromosome-disjoint GroupKFold folds; inner RidgeCV on training folds; standard genetic code excludes stops from composition features','models':results,'synonymous_incremental_mse_reduction':float(np.mean(base-full)),'conditional_chromosome_bootstrap95':np.quantile(diff,[.025,.975]).tolist(),'limits':['Chromosome separation does not remove paralog/gene-family leakage','Conditional cluster bootstrap does not refit models','Amino-acid and within-AA codon fractions correlate with unmeasured regulation','Association only; no redesigned-sequence validation','S288C reference may differ from experimental strain','No replicate biological uncertainty']}
 (ROOT/'results/yeast_synonymous_specificity.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
