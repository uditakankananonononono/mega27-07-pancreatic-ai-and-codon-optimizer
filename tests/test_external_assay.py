import json,sys,hashlib
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from external_assay_audit import compute
def test_landed_provenance_and_exact_sequence_scope():
 j=json.loads((ROOT/'results/external_assay_audit.json').read_text());source=ROOT/'data/external/nature16509_data2.xlsx'
 if source.exists():assert compute()==j
 assert j['file_sha256']=='7e32a22549f4565a08284e64f4187c5562181b909bfb382849cd9227fe4ef978'
 assert j['n_rows']==6348 and j['unique_exact_sequences']==6348 and j['non_acgt_rows']==0 and j['not_triplet_length_rows']==0
 assert j['exact_full_construct_overlap_rows_with_current_paxdb_cds']==0
 assert sum(j['expression_categories'].values())==6348
def test_frozen_transfer_ledger_metrics_and_overlap():
 j=json.loads((ROOT/'results/external_transfer_audit.json').read_text());raw=(ROOT/'results/external_transfer_predictions.jsonl').read_bytes();assert hashlib.sha256(raw).hexdigest()==j['prediction_ledger_sha256']
 rows=[json.loads(s) for s in raw.splitlines()];assert len(rows)==6348
 pred=np.array([r['natural_abundance_ridge_prediction'] for r in rows]);target=np.array([r['expression_category'] for r in rows]);keep=np.array([not r['protein_containment_overlap'] for r in rows]);assert int((~keep).sum())==43 and keep.sum()==6305
 assert np.isclose(spearmanr(pred,target).statistic,j['all_construct_spearman']) and np.isclose(spearmanr(pred[keep],target[keep]).statistic,j['protein_containment_excluded_spearman'])
 assert j['plan_sha256']==hashlib.sha256((ROOT/'scripts/external_transfer_plan.json').read_bytes()).hexdigest()
