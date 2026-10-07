import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def test_qualified_source_counts_and_units():
 j=json.loads((ROOT/'results/gse63789_qualification_audit.json').read_text())
 assert j['plan_sha256']==hashlib.sha256((ROOT/'scripts/gse63789_qualification_plan.json').read_bytes()).hexdigest()
 assert j['rows']==j['unique_genes']==5225 and j['duplicates']==0
 assert j['position_length_mismatches']==j['position_sum_mismatches']==0
 assert j['zero_rna']==j['zero_fp']==0
 assert j['unique_cds_joined']==5142 and len(j['unmatched_genes'])==83
 assert sum(j['three_times_reported_minus_cds_nt_histogram'].values())==5142
 assert j['three_times_reported_minus_cds_nt_histogram']['0']==5110
 assert j['old_gse75897_gene_overlap']==4785
 assert j['qualification'].startswith('conditional:')
 assert len(j['length_mismatch_genes'])==32
 assert 'reported_minus_cds_length_histogram' not in j
 assert j['source_sha256']=='1177ee0766cdc664787b1b32044fcbd4a24a55c921ec6844cd74080371c26e93'
