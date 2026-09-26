# ICOR benchmark rescore under the source-pinned checkpoint (2026-09-26)

Script: scripts/icor_rescore.py. Data: ICOR's own 40-gene benchmark sequences
(Zenodo 10.5281/zenodo.7487432, v1.4), re-fetched (data/icor is not committed).

Why: the historical icor_headtohead.json was scored with an earlier checkpoint that
was never committed and is unrecoverable (overwritten by the source-pinned CLI
checkpoint). The historical ours_cnn z=+2.20 is a historical claim only. This rescore
re-establishes the comparison end-to-end under the auditable model
(PaxDb 4.2 / MG1655 RefSeq, 5-epoch exploratory, holdout Spearman 0.564).

## Results (pred_expr_z mean, n=40 genes)

| method | old z | new z | CAI (unchanged) |
|---|---|---|---|
| ours_cnn (v2 designs, re-optimized vs new model) | 2.20 | 1.336 | 0.820 |
| ICOR | 0.68 | 0.331 | 0.875 |
| HFC | 0.60 | 0.465 | 0.997 |
| GenSmart | 0.24 | -0.019 | 0.731 |

Paired ours_v2 vs ICOR, 40 genes: mean diff +1.005, Wilcoxon p=1.8e-12.
Ranking stability: ours first under BOTH models; only ICOR/HFC swap 2nd/3rd.
The CAI-circularity finding replicates: HFC reaches CAI 0.997 with z=0.465 -
CAI ceiling is not expression ceiling.

Caveats: ours_cnn_v2 sequences were optimized against the same model that scores
them (in-sample optimizer advantage, same as before - stated); the exploratory
checkpoint is 5-epoch with holdout Spearman 0.564, so absolute z values are
model-relative, not wet-lab expression.
