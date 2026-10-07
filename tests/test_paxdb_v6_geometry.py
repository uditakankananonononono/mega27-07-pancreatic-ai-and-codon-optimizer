import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def test_identity_geometry_record():
 j=json.loads((ROOT/'results/paxdb_v6_geometry_audit.json').read_text());assert j['plan_sha256']==hashlib.sha256((ROOT/'scripts/paxdb_v6_geometry_plan.json').read_bytes()).hexdigest();assert j['old_unique_loci']==3737;assert j['old_train_test_shared_loci']==['b2592'];assert j['full_CDS_DNA_changed_overlap']==0;assert j['no_scoring_or_fitting']
