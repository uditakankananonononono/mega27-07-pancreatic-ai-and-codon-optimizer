# MEGA27-07: PANCREAS-AI-PLUS + Codon Optimizer

Two open-data tools modelled on the problems in Rishab Jain's PANCREAS.AI (biopsy-image mutation prediction) and ICOR (RNN codon optimization). No superiority over either is established: see the status notes below.

## pancreatic_ai
Predicts driver mutation status (KRAS / TP53 / CDKN2A / SMAD4) in pancreatic
adenocarcinoma from the rest of the tumor's mutational landscape + clinical
covariates (TCGA-PAAD, 184 samples, 20,703 mutations, cBioPortal public API).
The target gene is excluded from features (no leakage). CNN over a
gene x mutation-type tensor, GNN over the co-mutation graph, logistic baseline.

## codon_optimizer
CNN expression predictor trained on 3,752 real E. coli MG1655 CDS joined to
PaxDb whole-organism protein abundances (spearman 0.61 held out in the historical run; the current auditable checkpoint gets 0.564), plus a
CNN-guided synonymous-recoding optimizer compared against CAI-greedy and
wild-type on identical proteins.

## Status (2026-10-07 audit)
- Pancreas: within-repo TCGA-PAAD results are KRAS AUC 0.85-0.86, TP53 0.79-0.80, CDKN2A 0.68-0.73, SMAD4 0.51-0.55 (chance). On the external CPTAC-PDAC cohort (`results/cptac_external_validation.json`, n=140) external AUCs are KRAS 0.384 (95% CI 0.157-0.621; 96% prevalence), TP53 0.741 (0.643-0.830), CDKN2A 0.497 (0.376-0.612), SMAD4 0.347 (0.253-0.449), i.e. only TP53 transfers; no published like-for-like co-mutation baseline was located (`results/PANCREAS_BENCHMARK_CONTEXT.md`).
- Codon: under the auditable 5-epoch checkpoint the optimizer's predicted-expression z is 1.336 vs ICOR 0.331 (n=40, paired Wilcoxon p=1.8e-12), but the designs were optimized against the same model that scores them, and predicted z is model-relative, not wet-lab expression (`results/RESCORE.md`). No wet-lab or independent-model test.

## Data
`python3 scripts/fetch_paad.py && python3 scripts/fetch_ecoli.py`

## Reproducibility boundary (2026-09-25 audit)
The CLI now includes a source-pinned public GtRNAdb tRNAscan output, NCBI MG1655 RefSeq CDS, versioned PaxDb 4.2 abundance file, and a five-epoch *exploratory* checkpoint trained locally from those inputs (`scripts/train_cli_checkpoint.py`; `results/cli_checkpoint_provenance.json`). The current source snapshot yields 3,746 CDS-abundance pairs, 3,184 train and 562 held-out genes, and held-out Spearman 0.564. This new checkpoint is **not** the previously reported 3,752-pair model or the separately archived 0.627 hold-out experiment, and has not been used to rescore the historical ICOR comparison. Both end-to-end CLI tests pass with these bundled assets, and 68 test functions exist (static count, suite not rerun in this audit; the earlier README said 50 pass). Source rights and software dependencies still matter for use outside this repo. The strict 40-tool gate remains open. The old PaxDb `latest/abundances` URL now returns 404; `scripts/fetch_ecoli.py` uses the pinned versioned URL.

## Install and CLI
```
pip install -e .        # installs the codonopt command (setup.py shim enables PEP 660 editable installs)
codonopt --help         # subcommands: optimize, metrics
codonopt optimize --gene BIRC5 ...   # multi-objective optimization (expression z subject to tAI floor)
codonopt metrics ...                 # CAI/tAI/predicted-expression scoring of a sequence
```
The package pins its public data sources (GtRNAdb tRNAscan output, NCBI MG1655 RefSeq CDS, versioned PaxDb 4.2) and ships an exploratory checkpoint for evaluation; see the reproducibility boundary section below for exactly what the bundled checkpoint is and is not.

### Admitted external ratio source: fixed RNA-count strata
No newfit:5110GSE63789genes/fourRNAcountbins. CodonfrequencySpearman.359/.462/.510/.417,combined.355/.412/.459/.363,lengthGC.186/.223/.283/.179. PredictionalsoassociateswithRNA/FPcomponents;stratificationnotcausaladjustment ornoiseestimate. Samegene/sourceandratioendpointlimitationsremain. See `results/gse63789_count_strata_audit.json`.

### CLI input identity correction
The archived3746CLIinputrowsare3737uniqueloci,8ambiguousloci/17records;originalrow-split sharesb2592acrosstrain/test. Historical.564score isnot certifiedlocus-disjoint. NewPaxDb6.0geometry(notrefit/validation):3708unambiguousoverlapCDSDNA/prefixidentical;29newproteinmismatchesquarantined. See `results/paxdb_v6_geometry_audit.json`.
Unique-locusexploratory4.2rerun excludesall17ambiguousrecords:3729loci,3169train/560test,0locusoverlap,5epochs. Spearman.594738vsold.564187,not causalestimateofleakageeffect(changedsplit/dataset). Separatecheckpoint/historicalretained,no6.0externalvalidationscore. See `results/cli_unique_locus_audit.json`.

### Explicit CLI model selection (breaking safety change)
Both `optimize` and `metrics` now require `--ckpt` and `--provenance`; there is no historical-model default. Example:
```
codonopt metrics designs.fasta --ckpt results_expr_unique_locus_model.pt --provenance results/cli_unique_locus_audit.json
```
The selected JSON must equal the checkpoint's embedded provenance or loading fails. Checkpoints use restricted weights-only loading; use trusted files. Output records checkpoint/provenance SHA256, training scope and normalization. `metrics` now returns `{model_selection, genes}` rather than a bare gene map. Scores remain checkpoint-specific normalized predictions, not fold changes or demonstrated protein yield. Historical checkpoint selection is possible only by explicitly naming its matching historical ledger, whose locus-overlap caveat remains in the geometry audit.

The expression dataset API also requires `build_expression_dataset(policy="raw"|"unique")`. `raw` is historical replay, retaining ambiguous locus rows; it is not a locus-disjoint validation dataset. `unique` removes all records of every repeated locus after the abundance join, before downstream filters, preserving survivor order. Existing historical audit/training scripts explicitly name `raw`; the repaired trainer obtains `unique` from this API. Direct module invocation requires `--policy`. No source file or historical result is changed by this migration.
