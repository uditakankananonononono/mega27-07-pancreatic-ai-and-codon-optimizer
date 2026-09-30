"""Integrity of archived similarity-separated OOF result, no fitting at test time."""
import json
from pathlib import Path
import numpy as np

def test_cluster_fold_disjoint_and_loss_recomputable():
 j=json.loads((Path(__file__).resolve().parents[1]/'results/yeast_protein_cluster_specificity.json').read_text())
 assert j['n_genes']==len(j['gene_ids'])==len(set(j['gene_ids']))
 assert j['n_clusters']==len(set(j['cluster_assignment'].values()))
 seen=set()
 for f in j['models']['aa_length_gc']['folds']:
  groups=set(f['test_clusters']);assert not seen & groups;seen|=groups
 assert seen==set(j['cluster_assignment'].values())
 y=np.asarray(j['target_values']);p=j['oof_predictions']
 a=np.mean((y-np.asarray(p['aa_length_gc']))**2);b=np.mean((y-np.asarray(p['aa_syn_length_gc']))**2)
 assert np.isclose(a-b,j['synonymous_incremental_mse_reduction'])
 assert np.isclose(a,j['models']['aa_length_gc']['oof_mse'])
