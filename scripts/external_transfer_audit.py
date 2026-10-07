"""Frozen natural-abundance ridge transferred to an ordinal T7 assay; no external tuning."""
import hashlib,json,sys
from pathlib import Path
import numpy as np,openpyxl
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from Bio.Seq import Seq
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from codon_optimizer.expr_data import build_expression_dataset
from oof_evaluator_audit import freq
from locus_grouped_evaluator_audit import grouped_folds
def compute():
 raw=(ROOT/'scripts/external_transfer_plan.json').read_bytes();plan=json.loads(raw);file=ROOT/'data/external/nature16509_data2.xlsx';assert hashlib.sha256(file.read_bytes()).hexdigest()==plan['file_sha256']
 rr=list(openpyxl.load_workbook(file,read_only=True,data_only=True)['E6348.csv'].values);rows=[dict(zip(rr[0],r)) for r in rr[1:]];train=[(lt,seq,logab) for lt,gene,seq,logab in build_expression_dataset(policy="raw") if len(seq)>=90];groups=np.array([r[0] for r in train]);X=np.array([freq(r[1]) for r in train]);y=np.array([r[2] for r in train]);proteins=sorted(set(str(Seq(r[1]).translate()).rstrip('*') for r in train));assert all(len(p)>=30 for p in proteins)
 overlaps=[]
 for r in rows:
  prot=str(Seq(r['seq']).translate()).rstrip('*');hits=[hashlib.sha256(p.encode()).hexdigest() for p in proteins if p in prot];overlaps.append(hits)
 model=RidgeCV(alphas=np.logspace(-3,3,13),cv=list(grouped_folds(groups,501,3))).fit(X,y);pred=model.predict(np.array([freq(r['seq']) for r in rows]));target=np.array([r['E'] for r in rows]);keep=np.array([not bool(x) for x in overlaps]);ledger=[{'name':r['name'],'expression_category':r['E'],'sequence_sha256':hashlib.sha256(r['seq'].encode()).hexdigest(),'natural_abundance_ridge_prediction':float(pred[i]),'protein_containment_overlap':bool(overlaps[i]),'overlapping_training_protein_hashes':overlaps[i]} for i,r in enumerate(rows)];text=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in ledger)
 out={'plan_sha256':hashlib.sha256(raw).hexdigest(),'assay_file_sha256':plan['file_sha256'],'selected_training_alpha':float(model.alpha_),'n_train_cds':len(train),'n_assay_constructs':len(rows),'protein_containment_overlap_constructs':int(np.sum(~keep)),'all_construct_spearman':float(spearmanr(pred,target).statistic),'remaining_constructs':int(keep.sum()),'protein_containment_excluded_spearman':float(spearmanr(pred[keep],target[keep]).statistic),'prediction_ledger_sha256':hashlib.sha256(text.encode()).hexdigest(),'limits':plan['limits']};return out,text
def main():
 out,text=compute();(ROOT/'results/external_transfer_audit.json').write_text(json.dumps(out,indent=2)+'\n');(ROOT/'results/external_transfer_predictions.jsonl').write_text(text);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
