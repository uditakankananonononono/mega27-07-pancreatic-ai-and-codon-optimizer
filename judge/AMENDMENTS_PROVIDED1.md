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

## LANDED 2026-09-27 (commit below)
- #1/#14/#15: CodonOpt flagship restructure - new title, abstract rewritten CodonOpt-first, executive summary with one main claim + regenerated central figure (fig now INCLUDES the CodonOpt star at CAI 0.829 / z 2.20; prior build's figure omitted it, fixed).
- #3: "research-grade genomic inference" explicit in abstract.
- #9: RF independent evaluator LANDED with honest twist - codon-usage-only RF ranks ICOR 1st / CodonOpt 3rd (held-out r=0.678, 2996/749 split); framed as concrete proof that single-codon metrics are blind to context gains (strengthens the circularity discovery). results/rf_evaluator.json + new paper section.
- #12: blinded held-out-protein benchmark LANDED. 60 held-out proteins x 2 split seeds; optimizer sees train half only; ridge evaluator fit on opposite half (held-out Spearman 0.611/0.607). CNN beats WT under blinded ridge: z 0.247 vs 0.072 and 0.170 vs 0.005, boot95 excl. zero, Wilcoxon p 8.2e-5/6.8e-5, ahead on 73-78% of genes. Honest boundary: CAI-greedy still wins tAI on 93% of genes (codon-usage metric by construction) - consistent with the circularity finding. results/blinded_benchmark.json + paper section sec:blinded; paper 55pp.
- #13: explicit-vs-published statement LANDED - evaluator chain spelled out (PaxDb published abundances -> trained evaluators -> scored designs; held-out/blinded audits at each arrow; no wet-lab measurement of the designed constructs exists or is claimed). Paper section sec:published-expr.
- #11: motif scan LANDED - all 40 genes x 8 arms scanned for polyA signals, cryptic splice donors (GT[AG]AGT), AU-rich elements (ATTTA), CpG density. Our designs have FEWER hazards than WT (polyA 0.50 vs 0.78; ARE 0.23 vs 1.00, best of all arms; splice donors unchanged; CpG mid-range). scripts/motif_scan.py, results/motif_scan.json, paper section sec:motifs; 56pp.
Still queued: #4 CPTAC, #5 nested CV stats, #6 pathway-level stability methods, #7 driver discovery, #8 Ribo-seq/proteomics, #10 yeast arm, #17 packaging, #18 foundation models, #19 sequence-features analysis.
