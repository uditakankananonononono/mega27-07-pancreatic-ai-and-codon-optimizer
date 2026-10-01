# External-validation source audit, 1 October 2026

This is source scouting, not a completed independent validation. The current codon target is PaxDb protein abundance (log10 ppm), not mRNA abundance or expression measured after synonymous edits. The current pancreatic outcome is the cached mutation-status definition from TCGA-PAAD; a new source must preserve or explicitly change that definition.

## Pancreatic cohort candidate and a label mismatch
The Pancreas Genome Phenome Atlas ICGC Australian-cohort page exposes clinical information and columns explicitly titled Missense Mutation - KRAS/TP53/CDKN2A/SMAD4. It includes DO-prefixed patient identifiers, ages and missing ages. The retrieved rendering stops partway through a row, so no complete cohort denominator or mutation prevalence is asserted. Hidden columns and vertical bars within cell values make naive Markdown-delimiter parsing unsafe. A missense-only CDKN2A column is not the full mutation-status label, and these four genes do not supply the original burden/pathway feature panel. Treating every No as absence of all CDKN2A alteration would silently change the estimand. No model was run against this portal and no patient-independent external AUC is claimed.

A 2017 primary pancreatic-genomics article explicitly discusses highly prevalent CDKN2A homozygous deletions. This reinforces that missense status is narrower than overall genomic alteration, but the article does not prove the cached training label includes deletions. The newer portal paper describes multiple curated pancreatic sources. Neither study alone proves the accessible Australian records are nonoverlapping with every training/selection source. Full mutation and clinical files, sample-to-patient mapping, source-specific assay definitions and provenance/overlap checks remain necessary before external scoring.

Sources inspected:
- https://pancreasexpression.org/analytics/cohort/icgc_paca_au/ (operator cohort page, title/column names/record format)
- https://www.nature.com/articles/bjc2017209 (primary article, Clinical study of genomic drivers in pancreatic ductal adenocarcinoma, 18 July 2017; deletion scope)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12523802/ (primary portal article; multiple source/cohort context)

## Codon-source candidates have different endpoints
The 2016 Nature article reports 6,348 T7-polymerase expression experiments in E. coli and redesigned synonymous constructs, with protein and mRNA analyses. This is a promising design/assay-transfer source, but not an automatic held-out natural-MG1655 protein-abundance benchmark. Strain, promoter, expression scoring, construct sequences, supplements and training/selection overlap must be checked. No supplement contents or completed join is claimed.

The 2022 PLOS article links RNA-seq and ribosome profiling to GSE182100 and studies nutrient/growth conditions. Its relative translation efficiency is footprint density divided by RNA density, not absolute protein abundance or synthetic-design gain. The article reports 2,914 genes after its own mRNA cutoff and twelve conditions. Those are the paper's analysis counts, not a join count for this repo. A separate attempt to fetch the official GEO record failed with SOURCE_NOT_AVAILABLE; no official series contents were verified. Repeated measurements of the same loci across conditions are not new held-out genes; transfer across assays/conditions must be distinguished from gene-level generalization. No external score or independence claim is made.

Sources inspected (plus one failed official-record access):
- https://www.nature.com/articles/nature16509 (primary article, Codon influence on protein expression in E. coli correlates with mRNA levels, 13 January 2016)
- https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1010641 (primary article, Global and gene-specific translational regulation in Escherichia coli across different conditions, 20 October 2022)
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182100 (accession URL observed in the PLOS article; attempted fetch failed, so not an inspected series record)

## Decision
Keep these as source candidates. Do not turn their labels into the current target or call them independent validation until full files, identity/overlap, assay definitions and frozen transfer protocol have been verified. The nested partition rerun in this patch remains internal robustness on the same protein-abundance join.
