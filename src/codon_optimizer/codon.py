"""Codon utilities: standard code, usage tables, CAI (Sharp & Li 1987)."""
import gzip
import re
import math
import pathlib
from collections import Counter, defaultdict

DATA = pathlib.Path(__file__).resolve().parents[2] / "data" / "ecoli"

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
STANDARD_CODE = {}
_i = 0
for aa in "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG":
    STANDARD_CODE[CODONS[_i]] = aa
    _i += 1
SYNONYMS = defaultdict(list)
for cod, aa in STANDARD_CODE.items():
    if aa != "*":
        SYNONYMS[aa].append(cod)


def parse_cds_fasta(path=None):
    """Return [(locus_tag, gene_name, sequence)] for complete CDS records."""
    path = path or (DATA / "mg1655_cds.fna.gz")
    records, header, seq = [], None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header, seq = line, []
            else:
                seq.append(line)
        if header is not None:
            records.append((header, "".join(seq)))
    out = []
    for header, s in records:
        lt_m = re.search(r"\[locus_tag=([^\]]+)\]", header)
        g_m = re.search(r"\[gene=([^\]]+)\]", header)
        lt = lt_m.group(1).strip() if lt_m else None
        gene = g_m.group(1).strip() if g_m else None
        s = s.upper()
        if len(s) % 3 == 0 and lt:
            out.append((lt, gene or lt, s))
    return out


def codon_counts(records):
    c = Counter()
    for _, _, s in records:
        for i in range(0, len(s) - 2, 3):
            cod = s[i:i + 3]
            if cod in STANDARD_CODE:
                c[cod] += 1
    return c


def cai_weights(records):
    """Relative synonymous codon usage w_i = f_i / max(f_syn) per amino acid."""
    c = codon_counts(records)
    w = {}
    for aa, cods in SYNONYMS.items():
        fmax = max((c[x] for x in cods), default=0)
        for x in cods:
            w[x] = (c[x] / fmax) if fmax > 0 else 0.0
    return w


def cai(seq, w):
    logs = []
    for i in range(0, len(seq) - 2, 3):
        cod = seq[i:i + 3]
        aa = STANDARD_CODE.get(cod)
        if aa in (None, "*", "M", "W"):
            continue
        wi = w.get(cod, 0.0)
        if wi > 0:
            logs.append(math.log(wi))
    return math.exp(sum(logs) / len(logs)) if logs else 0.0
