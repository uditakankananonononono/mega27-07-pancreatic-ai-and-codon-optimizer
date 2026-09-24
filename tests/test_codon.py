import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from codon_optimizer import codon, optimize

def test_standard_code():
    assert codon.STANDARD_CODE["ATG"] == "M"
    assert codon.STANDARD_CODE["TAA"] == "*"
    assert len(codon.SYNONYMS["L"]) == 6

def test_cai_perfect_sequence():
    w = {c: 1.0 for c in codon.CODONS}
    seq = "ATGGCAGCATAA"
    assert abs(codon.cai(seq, w) - 1.0) < 1e-9

def test_cai_greedy_synonymous():
    w = codon.cai_weights  # noqa - just ensure import
    weights = {c: (1.0 if c.endswith("G") else 0.1) for c in codon.CODONS}
    prot = "MASL"
    opt = optimize.cai_greedy(prot, weights)
    assert optimize.to_protein(opt + "TAA") == prot  # preserves protein
    # each chosen codon should be the max-weight synonym
    from codon_optimizer.codon import STANDARD_CODE, SYNONYMS
    for i, aa in enumerate(prot):
        chosen = opt[i*3:i*3+3]
        assert weights[chosen] == max(weights[c] for c in SYNONYMS[aa])

def test_to_protein_strips_stop():
    assert optimize.to_protein("ATGAAATAA") == "MK"
