import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def test_no_validation_eligibility_scope():
 j=json.loads((ROOT/'results/gse182100_eligibility_audit.json').read_text());assert j['plan_sha256']==hashlib.sha256((ROOT/'scripts/gse182100_eligibility_plan.json').read_bytes()).hexdigest();assert not j['validation_performed'] and not j['ratios_constructed'];assert len(j['columns'])==36;assert j['unique_coordinate_keys']==4321
