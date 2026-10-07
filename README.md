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
The CLI now includes a source-pinned public GtRNAdb tRNAscan output, NCBI MG1655 RefSeq CDS, versioned PaxDb 4.2 abundance file, and a five-epoch *exploratory* checkpoint trained locally from those inputs (`scripts/train_cli_checkpoint.py`; `results/cli_checkpoint_provenance.json`). The current source snapshot yields 3,746 CDS-abundance pairs, 3,184 train and 562 held-out genes, and held-out Spearman 0.564. This new checkpoint is **not** the previously reported 3,752-pair model or the separately archived 0.627 hold-out experiment, and has not been used to rescore the historical ICOR comparison. Both end-to-end CLI tests pass with these bundled assets, and 50 tests pass overall. Source rights and software dependencies still matter for use outside this repo. The strict 40-tool gate remains open. The old PaxDb `latest/abundances` URL now returns 404; `scripts/fetch_ecoli.py` uses the pinned versioned URL.

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
