import importlib.util,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];s=importlib.util.spec_from_file_location('strata',ROOT/'scripts/gse63789_count_strata_audit.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def test_stratum_replay():
 j=m.compute();assert j==json.loads((ROOT/'results/gse63789_count_strata_audit.json').read_text());assert sum(r['n'] for r in j['bins'])==5110;assert sum(r['source_overlap_n'] for r in j['bins'])==4755
