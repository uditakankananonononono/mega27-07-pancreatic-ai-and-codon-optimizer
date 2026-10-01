"""Exploratory frozen covariate audit of retained ordinal transfer predictions."""
import hashlib,json
from pathlib import Path
import numpy as np,openpyxl
from scipy.stats import rankdata,spearmanr
ROOT=Path(__file__).resolve().parents[1]
def aggregate(rows):
 r=[x for x in rows if not x['protein_containment_overlap']];y=np.array([x['expression_category'] for x in r]);pred=np.array([x['natural_abundance_ridge_prediction'] for x in r]);gc=np.array([x['gc_fraction'] for x in r]);length=np.log([x['nt_length'] for x in r]);design=np.column_stack([np.ones(len(r)),rankdata(gc),rankdata(length)])
 def residual(v):
  ranks=rankdata(v);return ranks-design@np.linalg.lstsq(design,ranks,rcond=None)[0]
 return {'n':len(r),'raw_prediction_spearman':float(spearmanr(pred,y).statistic),'gc_spearman':float(spearmanr(gc,y).statistic),'length_spearman':float(spearmanr(length,y).statistic),'prediction_gc_spearman':float(spearmanr(pred,gc).statistic),'prediction_length_spearman':float(spearmanr(pred,length).statistic),'gc_length_partial_rank_correlation':float(np.corrcoef(residual(pred),residual(y))[0,1])}
def compute():
 raw=(ROOT/'scripts/external_covariate_plan.json').read_bytes();p=json.loads(raw);source=ROOT/'data/external/nature16509_data2.xlsx';assert hashlib.sha256(source.read_bytes()).hexdigest()==p['source_sha256'];sheet=list(openpyxl.load_workbook(source,read_only=True,data_only=True)['E6348.csv'].values);seq={r[0]:r[-1] for r in sheet[1:]};prior=[json.loads(x) for x in (ROOT/'results/external_transfer_predictions.jsonl').read_text().splitlines()];rows=[]
 for r in prior:
  s=seq[r['name']];assert hashlib.sha256(s.encode()).hexdigest()==r['sequence_sha256'];rows.append({k:r[k] for k in ['name','expression_category','natural_abundance_ridge_prediction','protein_containment_overlap']}|{'nt_length':len(s),'gc_fraction':(s.count('G')+s.count('C'))/len(s)})
 text=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in rows);out=aggregate(rows)|{'plan_sha256':hashlib.sha256(raw).hexdigest(),'ledger_sha256':hashlib.sha256(text.encode()).hexdigest(),'limits':p['limits']};return out,text
if __name__=='__main__':
 out,text=compute();(ROOT/'results/external_covariate_audit.json').write_text(json.dumps(out,indent=2)+'\n');(ROOT/'results/external_covariates.jsonl').write_text(text);print(json.dumps(out,indent=2))
