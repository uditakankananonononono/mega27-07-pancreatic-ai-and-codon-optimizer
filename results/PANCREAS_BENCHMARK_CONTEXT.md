# PANCREAS.AI published-benchmark context (2026-09-26, live web search)

Question: what published baselines exist for predicting PDAC driver status
(KRAS/TP53/CDKN2A/SMAD4) so our mutation-profile classifiers can be compared?

Findings (live search, 2026-09-26):
- Published driver-status classifiers for PDAC use DIFFERENT input modalities:
  deep learning on H&E pathology images predicting KRAS status
  (MDPI Journal of Personalized Medicine 14(7):249), pathomics-transcriptomics
  multimodal frameworks, epigenomic prognostic signatures. None predicts driver
  status FROM the tumor mutational profile (co-mutation patterns) as
  PANCREAS.AI does.
- Our current numbers (results/pancreatic_ai_benchmark.json): KRAS 0.85-0.86,
  TP53 0.79-0.80, CDKN2A 0.68-0.73, SMAD4 0.51-0.55 (at chance - documented
  boundary).

Consequence: a like-for-like published baseline for co-mutation-based driver
inference has not been located; the honest benchmark frame is currently
within-repo (CNN vs GNN vs logreg on identical splits) plus modality-adjacent
published AUCs cited with their input difference stated. A ChatGPT literature
triangulation consult is queued (judge-round route) to either locate a true
published comparator or confirm that none exists, in which case the SMAD4
chance-level boundary and the modality gap become the reported state.
