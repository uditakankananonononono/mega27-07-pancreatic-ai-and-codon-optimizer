import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def test_scope_corrections_preserve_archived_values():
 text=(ROOT/'paper/main.tex').read_text();j=json.loads((ROOT/'results/manuscript_claim_consistency.json').read_text())
 assert len(j)>=20 and all(r['before']!=r['after'] for r in j)
 assert '3.2$\\times$' not in text and '$4.0\\times$' not in text
 assert 'above the reported bootstrap interval' in text
 assert 'ridge uses codon-frequency features' in text
 result=json.loads((ROOT/'results/nested_cv_eval.json').read_text())
 assert result['views']['nt_frozen_ridge']['boot95_mean'][1]<.627
