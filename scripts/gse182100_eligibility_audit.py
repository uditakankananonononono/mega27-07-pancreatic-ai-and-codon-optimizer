"""Pre-admission eligibility only. No ratios or validation scores."""
import json,hashlib
from pathlib import Path
import numpy as np,openpyxl
ROOT=Path(__file__).resolve().parents[1]
def compute(folder):
 p=ROOT/'scripts/gse182100_eligibility_plan.json';plan=json.loads(p.read_text());tables=[]
 for name,h in plan['hashes'].items():
  f=Path(folder)/name;assert hashlib.sha256(f.read_bytes()).hexdigest()==h;tables.append(list(openpyxl.load_workbook(f,read_only=True,data_only=True).active.values))
 a,b=tables;assert a[0]==b[0];assert [r[:3] for r in a[1:]]==[r[:3] for r in b[1:]];keys=[r[:3] for r in a[1:]];assert len(set(keys))==len(keys)==4321
 RNA=np.array([r[3:] for r in a[1:]],float);FP=np.array([r[3:] for r in b[1:]],float);assert np.isfinite(RNA).all() and np.isfinite(FP).all() and (RNA>=0).all() and (FP>=0).all();eligible=(RNA>0)&(FP>0)
 return {'plan_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'source_sha256':plan['hashes'],'rows':len(keys),'unique_coordinate_keys':len(set(keys)),'duplicate_gene_name_extra_rows':len(keys)-len(set(r[0] for r in keys)),'RNA_zero_cells':int((RNA==0).sum()),'FP_zero_cells':int((FP==0).sum()),'all36_nominal_columns_positive_both_rows':int(eligible.all(1).sum()),'columns':[{'label':name,'RNA_zero':int((RNA[:,i]==0).sum()),'FP_zero':int((FP[:,i]==0).sum()),'positive_both':int(eligible[:,i].sum())} for i,name in enumerate(a[0][3:])],'ratios_constructed':False,'validation_performed':False,'limits':plan['limits']}
if __name__=='__main__':
 import argparse
 a=argparse.ArgumentParser();a.add_argument('--folder',required=True);v=a.parse_args();j=compute(v.folder);(ROOT/'results/gse182100_eligibility_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j)
