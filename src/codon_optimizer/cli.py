"""codonopt - multi-objective codon optimization CLI.

Maximizes CNN-predicted expression (z) subject to a tAI floor, under the explicitly selected checkpoint. Historical benchmark results do not
certify a different checkpoint or wet-lab yield. Input: protein FASTA. Output: optimized DNA FASTA + metrics JSON.

Usage:
  codonopt optimize input.fasta -o out_dir [--tai-floor F] [--iters N]
  codonopt metrics designs.fasta
"""
import argparse, json, pathlib, random, sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]


def gc(s):
    return (s.count("G") + s.count("C")) / max(len(s), 1)


def read_protein_fasta(path):
    out, name, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(buf).upper().rstrip("*")
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name is not None:
        out[name] = "".join(buf).upper().rstrip("*")
    return out


def read_dna_fasta(path):
    seqs = read_protein_fasta(path)
    return {k: v.replace("U", "T") for k, v in seqs.items()}


def load_scorer(trnascan=None):
    from .tai import copy_numbers, relative_weights, TaiScorer
    trnascan = trnascan or (REPO / "data" / "gtrnadb" / "eschColi_K_12_MG1655-tRNAs.out")
    return TaiScorer(relative_weights(copy_numbers(str(trnascan))))


def load_cai_weights(cds=None):
    from .codon import parse_cds_fasta, cai_weights
    return cai_weights(parse_cds_fasta(str(cds)) if cds else parse_cds_fasta())


def load_model(ckpt, provenance):
    import torch, hashlib
    from .model import ExpressionCNN
    if not ckpt or not provenance:
        raise ValueError("Explicit checkpoint and provenance paths are required")
    ckpt = pathlib.Path(ckpt)
    provenance = pathlib.Path(provenance)
    expected = json.loads(provenance.read_text())
    ck = torch.load(str(ckpt), map_location="cpu", weights_only=True)
    if ck.get("provenance") != expected:
        raise ValueError("Checkpoint embedded provenance does not match selected ledger")
    model = ExpressionCNN()
    model.load_state_dict(ck["state"])
    model.eval()
    model.selection_metadata = {"checkpoint_sha256": hashlib.sha256(ckpt.read_bytes()).hexdigest(),
        "provenance_sha256": hashlib.sha256(provenance.read_bytes()).hexdigest(),
        "scope": expected.get("scope"), "mean": ck["mean"], "sd": ck["sd"],
        "score_units": "checkpoint-specific normalized log10 abundance prediction, not fold change or wet-lab yield"}
    return model


def cnn_z(model, seq):
    import torch
    from .model import encode
    with torch.no_grad():
        return float(model(torch.tensor(encode(seq)).unsqueeze(0)))


def mo_hillclimb(protein, scorer, wcai, model, tai_floor, iters=150, batch=32,
                 seed=7, gc_band=(0.35, 0.70)):
    """Hill-climb CNN score; accept a mutation only if tAI >= floor and GC in band."""
    from .codon import SYNONYMS
    from .optimize import cai_greedy
    rng = random.Random(seed)
    seq = cai_greedy(protein, wcai)
    if scorer.tai(seq) < tai_floor:
        best = {aa: max(cods, key=lambda c: scorer.logw.get(c) or -1e9)
                for aa, cods in SYNONYMS.items()}
        seq = "".join(best[aa] for aa in protein)
    positions = [i for i, aa in enumerate(protein) if len(SYNONYMS[aa]) > 1]
    if not positions:
        return seq, cnn_z(model, seq)
    cur = cnn_z(model, seq)
    for _ in range(iters):
        cands = []
        for _ in range(batch):
            pos = rng.choice(positions)
            aa = protein[pos]
            alt = rng.choice([c for c in SYNONYMS[aa] if c != seq[pos*3:pos*3+3]])
            cand = seq[:pos*3] + alt + seq[pos*3+3:]
            if gc_band[0] <= gc(cand) <= gc_band[1] and scorer.tai(cand) >= tai_floor:
                cands.append(cand)
        if not cands:
            continue
        import torch
        from .model import encode
        with torch.no_grad():
            scores = model(torch.tensor(np.stack([encode(c) for c in cands]))).numpy()
        k = int(np.argmax(scores))
        if scores[k] > cur:
            cur = float(scores[k]); seq = cands[k]
    return seq, cur


def cmd_optimize(a):
    from .codon import cai
    scorer = load_scorer(a.trnascan)
    wcai = load_cai_weights(a.cds)
    model = load_model(a.ckpt, a.provenance)
    proteins = read_protein_fasta(a.input)
    outdir = pathlib.Path(a.out); outdir.mkdir(parents=True, exist_ok=True)
    rows = {}
    for name, prot in proteins.items():
        if "X" in prot or len(prot) < 5:
            print(f"skip {name}: too short or ambiguous", file=sys.stderr)
            continue
        floor = a.tai_floor
        if floor is None:
            from .optimize import cai_greedy
            floor = scorer.tai(cai_greedy(prot, wcai))
        seq, z = mo_hillclimb(prot, scorer, wcai, model, floor,
                              iters=a.iters, batch=a.batch, seed=a.seed)
        (outdir / f"{name}_mo.fasta").write_text(f">{name}_mo\n{seq}\n")
        rows[name] = {"tai": round(scorer.tai(seq), 4), "tai_floor": round(floor, 4),
                      "cai": round(cai(seq, wcai), 4), "z": round(z, 3),
                      "gc": round(gc(seq), 4), "len_nt": len(seq)}
        print(name, rows[name], flush=True)
    metrics_path = outdir / "codonopt_metrics.json"
    json.dump({"tool": "codonopt", "objective": "maximize CNN z s.t. tAI >= floor",
               "model_selection": model.selection_metadata, "genes": rows}, open(metrics_path, "w"), indent=1)
    print(f"wrote {metrics_path}")


def cmd_metrics(a):
    scorer = load_scorer(a.trnascan)
    wcai = load_cai_weights(a.cds)
    model = load_model(a.ckpt, a.provenance)
    from .codon import cai
    rows = {}
    for name, seq in read_dna_fasta(a.input).items():
        rows[name] = {"tai": round(scorer.tai(seq), 4), "cai": round(cai(seq, wcai), 4),
                      "z": round(cnn_z(model, seq), 3), "gc": round(gc(seq), 4),
                      "len_nt": len(seq)}
    print(json.dumps({"model_selection": model.selection_metadata, "genes": rows}, indent=1))


def main(argv=None):
    p = argparse.ArgumentParser(prog="codonopt",
        description="Multi-objective codon optimization (expression z subject to tAI floor).")
    sub = p.add_subparsers(dest="cmd", required=True)
    for name, fn in (("optimize", cmd_optimize), ("metrics", cmd_metrics)):
        sp = sub.add_parser(name)
        sp.add_argument("input", help="protein FASTA (optimize) or DNA FASTA (metrics)")
        sp.add_argument("--trnascan", help="GtRNAdb tRNA scan output (default: bundled E. coli K-12)")
        sp.add_argument("--cds", help="reference CDS FASTA for CAI weights (default: bundled MG1655)")
        sp.add_argument("--ckpt", required=True, help="Explicit trusted ExpressionCNN checkpoint path")
        sp.add_argument("--provenance", required=True, help="JSON ledger matching checkpoint embedded provenance")
        if name == "optimize":
            sp.add_argument("-o", "--out", default="codonopt_out")
            sp.add_argument("--tai-floor", type=float, default=None,
                            help="minimum tAI; default = tAI of CAI-greedy design")
            sp.add_argument("--iters", type=int, default=150)
            sp.add_argument("--batch", type=int, default=32)
            sp.add_argument("--seed", type=int, default=7)
        sp.set_defaults(fn=fn)
    a = p.parse_args(argv)
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
