import hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from codon_optimizer.expr_data import build_expression_dataset
from collections import Counter
def test_unique_locus_split_and_metrics():
 j=json.loads((ROOT/'results/cli_unique_locus_audit.json').read_text());assert j['plan_sha256']==hashlib.sha256((ROOT/'scripts/cli_unique_locus_plan.json').read_bytes()).hexdigest()
 raw=build_expression_dataset(policy="raw");c=Counter(r[0] for r in raw);rows=[r for r in raw if c[r[0]]==1];assert len(rows)==3729;assert len(raw)-len(rows)==17
 perm=np.random.default_rng(0).permutation(len(rows));assert [rows[i][0] for i in perm[:3169]]==j['train_locus_ids'];assert [rows[i][0] for i in perm[3169:]]==j['test_locus_ids'];assert not set(j['train_locus_ids'])&set(j['test_locus_ids'])
 assert np.isclose(spearmanr(j['test_predictions_normalized'],j['test_targets_normalized']).statistic,j['holdout_spearman'])
