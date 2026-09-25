"""PaxDb 4.2 abundance column mapping stays pinned to actual file schema."""
from codon_optimizer.expr_data import load_abundance

def test_three_column_paxdb_join(tmp_path):
    p = tmp_path / 'paxdb.tsv'
    p.write_text('#internal_id\tstring_external_id\tabundance\n1520\t511145.b0001\t12.5\n1521\t511145.b0002\t0\n')
    assert load_abundance(p) == {'b0001': 12.5, 'b0002': 0.0}
