"""Version and input identity only, no scoring/fitting."""
import csv,json,hashlib,sys
from pathlib import Path
from collections import Counter
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from codon_optimizer.expr_data import build_expression_dataset
ROOT=Path(__file__).resolve().parents[1]
def compute(folder):
 folder=Path(folder);q=json.loads((folder/'qualification.json').read_text());newhash={}
 for r in q['files']:
  f=folder/r['file'];h=hashlib.sha256(f.read_bytes()).hexdigest();assert h==r['sha256'];newhash[r['file']]=h
 prov=json.loads((ROOT/'results/cli_checkpoint_provenance.json').read_text())
 for f,h in prov['input_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
 rows=build_expression_dataset(policy="raw");assert len(rows)==prov['matched_rows'];counts=Counter(r[0] for r in rows);ambiguous={k for k,v in counts.items() if v>1};old={r[0]:r[2] for r in rows if r[0] not in ambiguous};perm=np.random.default_rng(0).permutation(len(rows));train={rows[i][0] for i in perm[:int(.85*len(rows))]};test={rows[i][0] for i in perm[int(.85*len(rows)):]} 
 seq={};key=None
 for line in (folder/'reference-selected-CDS.fasta').read_text().splitlines():
  if line.startswith('>'):key=line[1:].split()[0];assert key not in seq;seq[key]=''
  else:seq[key]+=line.strip()
 join=list(csv.DictReader((folder/'identity-join.csv').open()));eligible=[r for r in join if r['exact_protein_sequence_match']=='True'];overlap=[r for r in eligible if r['locus_tag'] in old];shared=[r['locus_tag'] for r in overlap];changed=[k for k in shared if old[k]!=seq[k]]
 return {'plan_sha256':hashlib.sha256((ROOT/'scripts/paxdb_v6_geometry_plan.json').read_bytes()).hexdigest(),'qualified_source_sha256':newhash,'old_input_sha256':prov['input_sha256'],'old_dataset_rows':len(rows),'old_unique_loci':len(counts),'old_ambiguous_loci':sorted(ambiguous),'old_ambiguous_record_count':sum(counts[k] for k in ambiguous),'old_train_rows':int(.85*len(rows)),'old_test_rows':len(rows)-int(.85*len(rows)),'old_train_test_shared_loci':sorted(train&test),'new_exact_protein_rows':len(eligible),'protein_mismatch_quarantined':len(join)-len(eligible),'new_exact_protein_old_locus_overlap':len(shared),'overlap_old_train':sum(k in train for k in shared),'overlap_old_test':sum(k in test for k in shared),'new_exact_protein_not_old_input':len(eligible)-len(shared),'full_CDS_DNA_changed_overlap':len(changed),'encoded1500prefix_changed_overlap':sum(old[k][:1500]!=seq[k][:1500] for k in shared),'new_exact_CDS_over1500':sum(len(seq[r['locus_tag']])>1500 for r in eligible),'old_CDS_over1500':sum(len(r[2])>1500 for r in rows),'changed_locus_tags':changed,'no_scoring_or_fitting':True,'limits':json.loads((ROOT/'scripts/paxdb_v6_geometry_plan.json').read_text())['limits']}
if __name__=='__main__':
 import argparse
 a=argparse.ArgumentParser();a.add_argument('--folder',required=True);v=a.parse_args();j=compute(v.folder);(ROOT/'results/paxdb_v6_geometry_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j)
