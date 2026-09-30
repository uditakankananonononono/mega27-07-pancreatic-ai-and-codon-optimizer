"""Protein-similarity cluster holdouts with nested group-disjoint ridge tuning.
Association only. Representative clusters do not exclude all homologs.
"""
import gzip,json,re,hashlib,itertools,subprocess,argparse
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV
from scipy.stats import pearsonr
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data/yeast'
def abund(p):return dict((g,float(v)) for g,v in (s.split() for s in gzip.open(p,'rt')))
def main(cdhit):
 codons=sorted(CodonTable.unambiguous_dna_by_id[1].forward_table);table=CodonTable.unambiguous_dna_by_id[1].forward_table;aas=sorted(set(table.values()));rpf=abund(D/'GSE75897_RPF_RPKMs.txt.gz');rna=abund(D/'GSE75897_RiboZero_RPKMs.txt.gz');rows=[];seen=set()
 for r in SeqIO.parse(gzip.open(D/'GCF_000146045.2_R64_cds_from_genomic.fna.gz','rt'),'fasta'):
  g=re.search(r'\[locus_tag=([^\]]+)\]',r.description)
  if not g:continue
  g=g.group(1);s=str(r.seq).upper();chr_=re.search(r'(NC_\d+\.\d+)',r.id)
  if g in seen or g not in rpf or g not in rna or rpf[g]<=0 or rna[g]<=0 or len(s)%3 or set(s)-set('ACGT') or not chr_:continue
  seen.add(g);counts=np.array([sum(s[k:k+3]==c for k in range(0,len(s),3)) for c in codons],float);total=counts.sum();ac=np.array([sum(counts[i] for i,c in enumerate(codons) if table[c]==a) for a in aas]);den={a:ac[i] for i,a in enumerate(aas)}
  aa=ac/total;syn=np.array([counts[i]/den[table[c]] if den[table[c]] else 0. for i,c in enumerate(codons)]);basic=[np.log(len(s)),(s.count('G')+s.count('C'))/len(s)]
  rows.append({'protein':str(Seq(s).translate()).rstrip('*'),'gene':g,'chromosome':chr_.group(1),'aa':aa,'syn':syn,'basic':basic,'y':np.log2(rpf[g]/rna[g])})
 protein_file=D/'matched_proteins.fasta';protein_file.write_text(''.join('>'+r['gene']+'\n'+r['protein']+'\n' for r in rows))
 cluster_file=D/'protein_clusters50';cmd=[cdhit,'-i',str(protein_file),'-o',str(cluster_file),'-c','.5','-n','3','-aS','.8','-aL','.8','-g','1','-d','0','-T','2','-M','1200']
 proc=subprocess.run(cmd,capture_output=True,text=True,timeout=90);assert proc.returncode==0,proc.stderr
 (ROOT/'results/yeast_cdhit50.log').write_text(proc.stdout+proc.stderr)
 assignment={};cluster=None
 for line in Path(str(cluster_file)+'.clstr').read_text().splitlines():
  if line.startswith('>Cluster'):cluster=line.split()[-1]
  else:
   hit=re.search(r'>([^.]+)\.\.\.',line);assert hit,line;assignment[hit.group(1)]=cluster
 assert all(r['gene'] in assignment for r in rows)
 y=np.array([r['y'] for r in rows]);groups=np.array([assignment[r['gene']] for r in rows]);basic=np.array([r['basic'] for r in rows]);aa=np.array([r['aa'] for r in rows]);syn=np.array([r['syn'] for r in rows]);views={'length_gc':basic,'aa_length_gc':np.c_[aa,basic],'aa_syn_length_gc':np.c_[aa,syn,basic]};splits=list(GroupKFold(5).split(y,y,groups));results={};predictions={}
 for name,x in views.items():
  pred=np.zeros(len(y));folds=[]
  for tr,te in splits:
   m=GridSearchCV(make_pipeline(StandardScaler(),Ridge()),{'ridge__alpha':[.1,1,10,100,1000]},cv=list(GroupKFold(4).split(x[tr],y[tr],groups[tr])),scoring='neg_mean_squared_error');m.fit(x[tr],y[tr]);pred[te]=m.predict(x[te]);folds.append({'pearson':float(pearsonr(y[te],pred[te]).statistic),'mse':float(np.mean((y[te]-pred[te])**2)),'alpha':m.best_params_['ridge__alpha'],'test_clusters':sorted(set(groups[te])),'n_test':len(te)})
  predictions[name]=pred;results[name]={'folds':folds,'oof_pearson':float(pearsonr(y,pred).statistic),'oof_mse':float(np.mean((y-pred)**2))}
 # Paired protein-cluster resampling, conditional on fitted OOF models.
 base=(y-predictions['aa_length_gc'])**2;full=(y-predictions['aa_syn_length_gc'])**2;rng=np.random.default_rng(778);unique=np.unique(groups);diff=[]
 cluster_delta=np.array([np.sum((base-full)[groups==g]) for g in unique]);cluster_n=np.array([np.sum(groups==g) for g in unique])
 for _ in range(3000):
  ix=rng.integers(len(unique),size=len(unique));diff.append(float(cluster_delta[ix].sum()/cluster_n[ix].sum()))
 out={'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in list(D.glob('*.gz'))+[protein_file,cluster_file,Path(str(cluster_file)+'.clstr')]},'cdhit_source':'https://github.com/weizhongli/cdhit','cdhit_source_commit':subprocess.check_output(['git','-C',str(Path(cdhit).parent),'rev-parse','HEAD']).decode().strip(),'cdhit_command':cmd,'cluster_assignment':assignment,'n_clusters':len(unique),'n_genes':len(y),'cluster_sizes':{g:int(np.sum(groups==g)) for g in unique},'oof_predictions':{name:pred.tolist() for name,pred in predictions.items()},'gene_ids':[r['gene'] for r in rows],'target_values':y.tolist(),'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75897','target':'log2 RPF/RNA','design':'five protein-cluster-disjoint GroupKFold folds; four inner group-disjoint folds tune ridge with standardization fitted separately in each inner train split; standard genetic code excludes stops from composition features','models':results,'synonymous_incremental_mse_reduction':float(np.mean(base-full)),'conditional_protein_cluster_bootstrap95':np.quantile(diff,[.025,.975]).tolist(),'limits':['CD-HIT50 clusters at80% coverage are representative-based similarity groups, not exhaustive homology families; distant homologs and cross-cluster similarity remain','Conditional cluster bootstrap does not refit models','Amino-acid and within-AA codon fractions correlate with unmeasured regulation','Association only; no redesigned-sequence validation','S288C reference may differ from experimental strain','No replicate biological uncertainty']}
 (ROOT/'results/yeast_protein_cluster_specificity.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--cdhit',required=True);x=a.parse_args();main(x.cdhit)
