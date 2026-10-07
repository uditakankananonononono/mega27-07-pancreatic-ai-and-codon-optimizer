"""Train a source-pinned exploratory CLI checkpoint from public MG1655/PaxDb data.

Not the historical validation model; no held-out genes enter training.
"""
import hashlib, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import torch
from scipy.stats import spearmanr
from codon_optimizer.expr_data import build_expression_dataset
from codon_optimizer.model import ExpressionCNN, encode
ROOT=pathlib.Path(__file__).resolve().parents[1]
FILES=['data/gtrnadb/eschColi_K_12_MG1655-tRNAs.out','data/ecoli/mg1655_cds.fna.gz','data/ecoli/paxdb_abundance.tsv']
from collections import Counter
prior=json.loads((ROOT/'results/cli_checkpoint_provenance.json').read_text());assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h for f,h in prior['input_sha256'].items());raw=build_expression_dataset(policy="raw");counts=Counter(r[0] for r in raw);rows=build_expression_dataset(policy="unique");assert len({r[0] for r in rows})==len(rows); rng=np.random.default_rng(0); ids=rng.permutation(len(rows)); n=int(.85*len(ids)); tr,te=ids[:n],ids[n:]
X=np.stack([encode(r[2]) for r in rows]); y=np.array([r[3] for r in rows],dtype=np.float32); mu,sd=y[tr].mean(),y[tr].std(); yz=(y-mu)/sd
print('dataset',len(rows),'train',n,'test',len(te),flush=True)
torch.manual_seed(0); torch.set_num_threads(2); model=ExpressionCNN(); opt=torch.optim.Adam(model.parameters(),lr=1e-3); lossfn=torch.nn.MSELoss()
for ep in range(5):
 model.train(); order=rng.permutation(tr)
 for i in range(0,len(order),64):
  b=order[i:i+64]; opt.zero_grad(); loss=lossfn(model(torch.tensor(X[b])),torch.tensor(yz[b])); loss.backward(); opt.step()
 print('epoch',ep+1,'last_loss',round(float(loss.detach()),4),flush=True)
model.eval()
with torch.no_grad(): pred=np.concatenate([model(torch.tensor(X[te[i:i+64]])).numpy() for i in range(0,len(te),64)])
res={'plan_sha256':hashlib.sha256((ROOT/'scripts/cli_unique_locus_plan.json').read_bytes()).hexdigest(),'excluded_ambiguous_records':len(raw)-len(rows),'train_locus_ids':[rows[i][0] for i in tr],'test_locus_ids':[rows[i][0] for i in te],'test_predictions_normalized':pred.tolist(),'scope':'unique-locus exploratory rerun on same admitted4.2snapshot;not family-independent/wetlabvalidation', 'source_urls':['https://gtrnadb.ucsc.edu/genomes/bacteria/Esch_coli_K_12_MG1655/eschColi_K_12_MG1655-tRNAs.tar.gz','https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz','https://pax-db.org/downloads/4.2/datasets/paxdb-abundance-files-v4.2/511145/511145-WHOLE_ORGANISM-integrated.txt'], 'input_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in FILES},'matched_rows':len(rows),'train_rows':len(tr),'test_rows':len(te),'seed':0,'epochs':5,'holdout_spearman':float(spearmanr(pred,yz[te]).statistic),'holdout_pearson':float(np.corrcoef(pred,yz[te])[0,1]),'test_targets_normalized':yz[te].tolist(),'label':'log10 protein abundance, normalized using train only'}
torch.save({'state':model.state_dict(),'mean':float(mu),'sd':float(sd),'provenance':res},ROOT/'results_expr_unique_locus_model.pt'); (ROOT/'results/cli_unique_locus_audit.json').write_text(json.dumps(res,indent=2)+'\n'); print('saved',{k:v for k,v in res.items() if not isinstance(v,list)},flush=True)
