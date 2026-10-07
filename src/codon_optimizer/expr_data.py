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
        if len(parts) >= 3:
            sid = parts[1].split(".")[-1]  # STRING ID 511145.b0001 -> b0001
            try:
                ab[sid] = float(parts[2])  # ppm abundance
            except ValueError:
                pass
    return ab

def build_expression_dataset(*, policy):
    """Explicit policy: raw preserves historical rows; unique quarantines every repeated locus.

    No first/longest sequence selection and no silent default. Counts are applied
    after the historical positive-abundance join, before any caller length filter.
    """
    if policy not in {"raw", "unique"}:
        raise ValueError("dataset policy must be raw or unique")
    records = parse_cds_fasta()
    ab = load_abundance()
    rows = []
    for lt, gene, seq in records:
        if lt in ab and ab[lt] > 0:
            rows.append((lt, gene, seq, math.log10(ab[lt])))
    if policy == "unique":
        from collections import Counter
        counts = Counter(row[0] for row in rows)
        rows = [row for row in rows if counts[row[0]] == 1]
    return rows

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--policy", choices=("raw", "unique"), required=True)
    rows = build_expression_dataset(policy=p.parse_args().policy)
    import numpy as np
    vals = np.array([r[3] for r in rows])
    print(f"matched {len(rows)} CDS with abundance; log10 ppm range {vals.min():.1f}..{vals.max():.1f}, median {np.median(vals):.2f}")
