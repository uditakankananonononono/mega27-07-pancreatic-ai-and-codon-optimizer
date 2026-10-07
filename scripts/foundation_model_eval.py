"""Foundation-model comparison (verdict #18): frozen Nucleotide Transformer
v2-50M (multi-species) embeddings + ridge head vs our CNN evaluator and the
codon-usage RF, on the SAME expression dataset and SAME 80/20 gene split,
then scored on the benchmark arms (WT / ICOR / GenScript / ours).

Frozen embeddings, ridge head only - no fine-tuning. This tests whether a
general DNA foundation model's sequence representations already contain the
expression signal our task-specific CNN learned, and how it ranks the
codon-optimization arms.
Checkpoints embeddings to results/nt_embeddings.npz so partial progress persists.
"""
import json, os, pathlib, sys, time
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from codon_optimizer.expr_data import build_expression_dataset

MAXLEN = 512  # nt tokens kept per sequence (uniform cap, documented)
BATCH = 16
CKPT = ROOT / "results/nt_embeddings.npz"

def main():
    import torch
    from transformers import AutoTokenizer, AutoModelForMaskedLM
    torch.set_num_threads(os.cpu_count())
    t0 = time.time()
    rows = [(lt, seq, logab) for lt, gene, seq, logab in build_expression_dataset(policy="raw") if len(seq) >= 90]
    print(f"dataset: {len(rows)} genes", flush=True)
    tok = AutoTokenizer.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-50m-multi-species", trust_remote_code=True)
    model = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-50m-multi-species", trust_remote_code=True)
    model.eval()
    emb = {}
    if CKPT.exists():
        z = np.load(CKPT, allow_pickle=True)
        emb = dict(zip(z["keys"].tolist(), z["vals"]))
        print(f"resumed: {len(emb)} embeddings", flush=True)
    todo = [(lt, s, y) for lt, s, y in rows if lt not in emb]
    for i in range(0, len(todo), BATCH):
        chunk = todo[i:i+BATCH]
        enc = tok([s[:MAXLEN] for _, s, _ in chunk], return_tensors="pt",
                  padding=True, truncation=True, max_length=MAXLEN//6+8)
        with torch.no_grad():
            out = model(**enc, output_hidden_states=True).hidden_states[-1]  # B,T,H
            mask = enc["attention_mask"].unsqueeze(-1).float()
            pooled = (out * mask).sum(1) / mask.sum(1).clamp(min=1)
        for (lt, _, _), e in zip(chunk, pooled.numpy()):
            emb[lt] = e
        if (i // BATCH) % 10 == 0:
            np.savez(CKPT, keys=np.array(list(emb.keys())), vals=np.array(list(emb.values())))
            print(f"embedded {len(emb)}/{len(rows)} ({time.time()-t0:.0f}s)", flush=True)
    np.savez(CKPT, keys=np.array(list(emb.keys())), vals=np.array(list(emb.values())))
    embmap = {k: v for k, v in zip(np.load(CKPT, allow_pickle=True)["keys"].tolist(),
                                   np.load(CKPT, allow_pickle=True)["vals"])}
    X = np.array([embmap[lt] for lt, _, _ in rows])
    y = np.array([logab for _, _, logab in rows])
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import train_test_split
    from scipy.stats import spearmanr
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=0)
    ridge = RidgeCV(alphas=np.logspace(-3, 3, 13)).fit(Xtr, ytr)
    pred = ridge.predict(Xte)
    res = {"model": "nucleotide-transformer-v2-50m-multi-species, frozen, mean-pooled, ridge head",
           "maxlen_nt": MAXLEN, "n_train": len(Xtr), "n_test": len(Xte),
           "holdout_pearson_r": float(np.corrcoef(pred, yte)[0, 1]),
           "holdout_spearman": float(spearmanr(pred, yte).statistic),
           "alpha": float(ridge.alpha_)}
    print(json.dumps(res, indent=1), flush=True)
    # score benchmark arms with the NT+ridge evaluator
    from rf_evaluator import read_fasta, BENCH  # reuse loader
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
    arms["codonopt"] = json.load(open(ROOT / "results/icor_rescore_ours_dna.json"))
    def embed_seqs(seqs):
        vals = {}
        items = list(seqs.items())
        for i in range(0, len(items), BATCH):
            chunk = items[i:i+BATCH]
            enc = tok([s[:MAXLEN] for _, s in chunk], return_tensors="pt",
                      padding=True, truncation=True, max_length=MAXLEN//6+8)
            with torch.no_grad():
                out = model(**enc, output_hidden_states=True).hidden_states[-1]
                mask = enc["attention_mask"].unsqueeze(-1).float()
                pooled = (out * mask).sum(1) / mask.sum(1).clamp(min=1)
            for (k, _), e in zip(chunk, pooled.numpy()):
                vals[k] = ridge.predict(e.reshape(1, -1))[0]
        return vals
    res["arms"] = {}
    for method, seqs in arms.items():
        v = embed_seqs({k: s for k, s in seqs.items() if len(s) >= 90})
        if v:
            arr = np.array(list(v.values()))
            res["arms"][method] = {"n": len(arr), "mean_pred_log10_abund": float(arr.mean()), "sd": float(arr.std())}
            print(method, res["arms"][method], flush=True)
    json.dump(res, open(ROOT / "results/foundation_model_eval.json", "w"), indent=1)
    print("saved results/foundation_model_eval.json", flush=True)

if __name__ == "__main__":
    main()
