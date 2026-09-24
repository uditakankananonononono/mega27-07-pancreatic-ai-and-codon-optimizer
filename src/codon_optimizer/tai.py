"""tRNA Adaptation Index (tAI, dos Reis et al. 2004 NAR) shared module.
tRNA gene copy numbers: GtRNAdb tRNAscan-SE output. Wobble s-values from
dos Reis 2004: WC 0.0, G:U 0.41, I:U/C/A 0.28, I:G 0.9999.
"""
import math
from .codon import STANDARD_CODE

S_WC, S_GU, S_IU, S_IG = 0.0, 0.41, 0.28, 0.9999


def revcomp(s):
    return s.translate(str.maketrans("ACGU", "UGCA"))[::-1]


def copy_numbers(trnascan_out_path):
    copy = {}
    for line in open(trnascan_out_path):
        p = line.split()
        if len(p) >= 6 and p[0] == "chr" and p[4] != "Pseudo" and len(p[5]) == 3:
            ac = p[5].replace("T", "U")
            copy[ac] = copy.get(ac, 0) + 1
    return copy


def relative_weights(copy):
    codons = [c for c in STANDARD_CODE if STANDARD_CODE[c] != "*"]
    W = {}
    for cod in codons:
        c = cod.replace("T", "U")
        acc = 0.0
        for ac, n in copy.items():
            anti = revcomp(ac)
            if anti[1:] != c[1:]:
                continue
            pair = ac[2] + c[0]
            if pair in ("AU", "UA", "GC", "CG"): s = S_WC
            elif pair == "GU": s = S_GU
            elif pair in ("IU", "IC", "IA"): s = S_IU
            elif pair == "IG": s = S_IG
            else: continue
            acc += (1 - s) * n
        W[cod] = acc
    wmax = max(v for v in W.values() if v > 0)
    return {c: (v / wmax if v > 0 else 0.0) for c, v in W.items()}


class TaiScorer:
    """Vectorized-ish tAI via log-sum over a codon list."""

    def __init__(self, wrel):
        self.logw = {c: (math.log(w) if w > 0 else None) for c, w in wrel.items()}

    def tai(self, seq):
        tot, n = 0.0, 0
        lw = self.logw
        for i in range(0, len(seq) - 2, 3):
            v = lw.get(seq[i:i+3])
            if v is not None:
                tot += v; n += 1
        return math.exp(tot / n) if n else 0.0
