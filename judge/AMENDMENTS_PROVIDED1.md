# Lane-07 AMENDMENT QUEUE - LOCKED BEFORE EXECUTION
Verdict: PROVIDED round 1 (WhatsApp 10:48:50 IST, wamid...RUVCOTk4MQA=, verbatim in judge/round_provided1_verdict_whatsapp.txt).
Locked: 2026-09-27 10:50 IST, against 52pp build 9ef60fc. No execution began before this lock.

## Weakness items (#1-#19) mapped vs current paper
- #1 ONE flagship (CodonOpt primary OR PANCREAS.AI++ primary, not equal): FRAMING - execute immediately. DECISION: CodonOpt primary (it carries the clean benchmark win: 3.2x ICOR on its own benchmark; the cancer arm becomes the replication-rigor supporting study). Title/abstract restructure.
- #2 novelty = replication-aware cancer genomics INTERPRETATION framework, not another classifier: FRAMING (cross-cohort replication already the core) + title shift
- #3 "research-grade genomic inference" explicit; no diagnosis/treatment claims: EDITORIAL
- #4 add ICGC/CPTAC/GENIE validation: PARTIAL (QCMG-UQ n=456 + MSK n=2336 already replicate; add CPTAC-PAAD if accessible via cBioPortal)
- #5 repeated nested CV + CIs + statistical comparison of CNN/GNN: NEW
- #6 SHAP/integrated-gradients/stability-selection at pathway level: PARTIAL (permutation importance + 16-split ensemble stability already; add one model-agnostic method at pathway level)
- #7 driver-discovery experiment beyond the 4 known drivers: NEW (score all genes, report novel candidates)
- #8 expression-independent metrics (ribosome occupancy/proteomics/translation efficiency): PARTIAL (MFE audit exists) -> add Ribo-seq/proteomics correlation vs public datasets
- #9 multiple independent evaluators (RF/transformer evaluator): PARTIAL (independent ridge validation exists) -> add RF evaluator
- #10 expand beyond E. coli (yeast or mammalian): NEW (yeast PaxDb arm cheapest)
- #11 more biological constraints (mRNA structure/regulatory motifs/cryptic splice/RNA stability): PARTIAL (ViennaRNA 5' MFE audit exists) -> add motif scan
- #12 blinded benchmark (held-out proteins never seen by optimization): PARTIAL (held-out gene split exists for the expression model; add blinded optimization set)
- #13 compare vs published expression measurements: PARTIAL (PaxDb abundances are exactly that; state explicitly)
- #14 unifying theme "Reliable sequence-to-function modeling in biology" if keeping both: FRAMING (adopt with CodonOpt flagship)
- #15 2-page executive summary + one central figure + one main claim: NEW - cheap, high value
- #16 move audits to supporting, lead with biological insight: EDITORIAL restructure
- #17 package both as reproducible tools: PARTIAL (CLI exists for codon pipeline - cli_checkpoint_provenance.json; extend + docs)
- #18 foundation-model comparison (ESM embeddings/DNA transformers): NEW - biggest compute item
- #19 mechanism: pathways for cancer (exists), sequence-features-learned analysis for codons: PARTIAL/NEW

## NOT-to-do (implied by #1/#3): no equal-footing two-domain presentation; no clinical-deployment implication.

## First deliverables (cheap, existing data)
1. Flagship restructure: CodonOpt primary, new title/abstract, executive summary (#1/#14/#15)
2. RF independent evaluator + blinded held-out-protein optimization benchmark (#9/#12)
3. Driver-discovery scan (#7) from existing mutation matrix
