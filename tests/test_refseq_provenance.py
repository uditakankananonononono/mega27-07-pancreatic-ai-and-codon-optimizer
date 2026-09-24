import importlib.util, os
spec = importlib.util.spec_from_file_location('rp', os.path.join(os.path.dirname(__file__), '..', 'scripts', 'refseq_provenance.py'))
rp = importlib.util.module_from_spec(spec)
import sys
sys.modules['rp'] = rp


def test_accession_regex():
    import re
    acc = re.compile(r'^>?([A-Z]{1,3}_?\d+\.\d+)\s')
    assert acc.match('>NM_139071.3 Homo sapiens').group(1) == 'NM_139071.3'
    assert acc.match('>U65407.1 Plasmodium').group(1) == 'U65407.1'
    assert acc.match('> Human peptide deformylase') is None


def test_summary_counts_from_committed_json():
    import json
    d = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'results', 'refseq_provenance.json')))
    s = d['summary']
    assert s['n_genes'] == 40
    assert s['n_exact_cds_match'] + s['n_mismatch'] + len(s['n_fetch_failed']) == s['n_with_accession']
