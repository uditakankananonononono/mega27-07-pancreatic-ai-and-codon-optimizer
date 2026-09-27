#!/usr/bin/env python3
"""Independent random-forest evaluator (verdict #9): does an RF trained on
codon-usage features + PaxDb abundance rank CodonOpt above ICOR too?

Independent of the CNN expression model: features are 61 codon frequencies
only, target is log PaxDb abundance, held-out gene split. Scores the published
benchmark sequences (wild type, ICOR, GenScript, HFC/BFC/URC/ERC) plus
CodonOpt's archived optimized DNA (results/icor_rescore_ours_dna.json).
"""
import gzip
import json
import pathlib

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

ROOT = pathlib.Path(__file__).resolve().parent.parent
BENCH = ROOT / "data/icor/Lattice-Automation-icor-codon-optimization-c60b775/benchmark_sequences"

CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
CODONS = [c for c in CODONS if c not in ("TAA", "TAG", "TGA")]


def codon_freq(seq):
    seq = seq.upper().replace("U", "T")
    n = (len(seq) // 3)
    if n == 0:
        return np.zeros(len(CODONS))
    counts = {c: 0 for c in CODONS}
    for i in range(0, 3 * n, 3):
        c = seq[i:i + 3]
        if c in counts:
            counts[c] += 1
    return np.array([counts[c] / n for c in CODONS])


def read_fasta(path):
    recs, name, buf = {}, None, []
    fh = path if hasattr(path, 'read') else open(path)
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            if name:
                recs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        recs[name] = "".join(buf)
    return recs


def main():
    # training data via the repo's own locus_tag join
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from codon_optimizer.expr_data import build_expression_dataset
    rows = build_expression_dataset()
    X, y = [], []
    for lt, gene, seq, logab in rows:
        if len(seq) >= 90:
            X.append(codon_freq(seq))
            y.append(logab)
    X, y = np.array(X), np.array(y)
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=0)
    rf = RandomForestRegressor(n_estimators=300, n_jobs=-1, random_state=0)
    rf.fit(Xtr, ytr)
    r_holdout = float(np.corrcoef(rf.predict(Xte), yte)[0, 1])

    # score benchmark
    ours = json.load(open(ROOT / "results/icor_rescore_ours_dna.json"))
    arms = {}
    for method, sub in [("wild_type", "all_original"), ("icor", "icor"), ("genscript", "genscript"),
                        ("HFC", "HFC"), ("BFC", "BFC"), ("URC", "URC"), ("ERC", "ERC")]:
        d = BENCH / sub
        if not d.exists():
            continue
        seqs = {}
        for f in sorted(d.glob("*.faa")) + sorted(d.glob("*.fna")) + sorted(d.glob("*.fasta")) + sorted(d.glob("*.txt")):
            seqs.update(read_fasta(f))
        if seqs:
            arms[method] = seqs
    arms["codonopt"] = ours

    # restrict to the 40 shared genes where possible
    out = {"design": "RF on 61 codon frequencies -> log10 PaxDb abundance; 80/20 gene split; independent of the CNN evaluator",
           "n_train": len(Xtr), "n_test": len(Xte), "holdout_pearson_r": r_holdout, "arms": {}}
    for method, seqs in arms.items():
        vals = [float(rf.predict(codon_freq(s).reshape(1, -1))[0]) for s in seqs.values() if len(s) >= 90]
        if vals:
            out["arms"][method] = {"n": len(vals), "mean_log10_abund": float(np.mean(vals)),
                                   "sd": float(np.std(vals))}
    (ROOT / "results/rf_evaluator.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1)[:1200])


if __name__ == "__main__":
    main()
