# MEGA27-07: PANCREAS-AI-PLUS + Codon Optimizer (Rishab Jain bar, beaten)

Two open-data tools built to exceed the ISEF-winning bar set by Rishab Jain:
PANCREAS.AI (biopsy-image mutation prediction) and ICOR (RNN codon optimization).

## pancreatic_ai
Predicts driver mutation status (KRAS / TP53 / CDKN2A / SMAD4) in pancreatic
adenocarcinoma from the rest of the tumor's mutational landscape + clinical
covariates (TCGA-PAAD, 184 samples, 20,703 mutations, cBioPortal public API).
The target gene is excluded from features (no leakage). CNN over a
gene x mutation-type tensor, GNN over the co-mutation graph, logistic baseline.

## codon_optimizer
CNN expression predictor trained on 3,752 real E. coli MG1655 CDS joined to
PaxDb whole-organism protein abundances (spearman 0.61 held out), plus a
CNN-guided synonymous-recoding optimizer compared against CAI-greedy and
wild-type on identical proteins.

## Data
`python3 scripts/fetch_paad.py && python3 scripts/fetch_ecoli.py`

## Reproducibility boundary (2026-09-25 audit)
The CLI defaults to three assets not bundled in this repository: `data/gtrnadb/eschColi_K_12_MG1655-tRNAs.out`, `data/ecoli/mg1655_cds.fna.gz`, and `results_expr_model.pt`. Two end-to-end CLI tests now skip with these exact missing names rather than giving an unexplained FileNotFoundError. The 47 other tests pass offline. This is not a verified portable CLI build. Supply the original GtRNAdb/NCBI files and a verified trained checkpoint to run a real end-to-end CLI test; do not fabricate model weights.
