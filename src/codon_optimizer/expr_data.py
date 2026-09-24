"""Join PaxDb protein abundances to MG1655 CDS via locus_tag (511145.bXXXX)."""
import math
import pathlib
from .codon import parse_cds_fasta, DATA

def load_abundance(path=None):
    path = path or (DATA / "paxdb_abundance.tsv")
    ab = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            sid = parts[0].split(".")[-1]  # 511145.b0001 -> b0001
            try:
                ab[sid] = float(parts[1])  # ppm abundance
            except ValueError:
                pass
    return ab

def build_expression_dataset():
    records = parse_cds_fasta()
    ab = load_abundance()
    rows = []
    for lt, gene, seq in records:
        if lt in ab and ab[lt] > 0:
            rows.append((lt, gene, seq, math.log10(ab[lt])))
    return rows

if __name__ == "__main__":
    rows = build_expression_dataset()
    import numpy as np
    vals = np.array([r[3] for r in rows])
    print(f"matched {len(rows)} CDS with abundance; log10 ppm range {vals.min():.1f}..{vals.max():.1f}, median {np.median(vals):.2f}")
