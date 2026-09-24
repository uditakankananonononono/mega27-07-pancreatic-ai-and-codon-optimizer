#!/usr/bin/env python3
"""Generate paper/myvariant_sec.tex from results/myvariant_panel_audit.json
(token replacement from committed JSON; no numbers typed by hand)."""
import json

a = json.load(open("results/myvariant_panel_audit.json"))
q1, q2, q3 = a["q1_positive_control"], a["q2_recurrence"], a["q3_patient_support"]
rep, cpg, ceil, pg = q2["cohort_split_replication"], q2["cpg_control"]["panel"], q2["am_ceiling"], q2["per_gene_am"]
f2 = lambda x: f"{x:.2f}"
f3 = lambda x: f"{x:.3f}"
def pv(x):
    return f"{x:.2g}".replace("e-0", "e-").replace("e-", r"\times10^{-") + ("}" if "e-" in f"{x:.2g}" else "")
pct = lambda x: f"{100*x:.1f}"
npat = sum(a["n_patients_by_study"].values())

rows = []
for sc, lab in (("am", "AlphaMissense"), ("revel", "REVEL"), ("cadd", "CADD phred")):
    rows.append(f"{lab} & {f3(q1[sc]['panel_median'])} & {f3(q1[sc]['passenger_median'])} & {f3(q1[sc]['auroc_panel_vs_passenger'])} "
                f"& {f3(q2['panel'][sc]['auroc_recurrent_vs_singleton'])} & {f3(rep['msk_impact'][sc]['auroc_recurrent_vs_singleton'])} "
                f"($p={pv(rep['msk_impact'][sc]['p'])}$) & {f3(rep['wes_cohorts'][sc]['auroc_recurrent_vs_singleton'])} ($p={pv(rep['wes_cohorts'][sc]['p'])}$) "
                f"& {f3(cpg[sc]['noncpg_only_rho'])} ($p={pv(cpg[sc]['noncpg_only_p'])}$) \\\\")
top = ", ".join(f"{t['gene']} {t['protein']} ({t['n_pat']})" for t in q2["top_panel_alleles"][:6])

tex = r"""\section{Allele-level audit: AlphaMissense saturates on panel alleles (MyVariant.info)}
The earlier audits worked at gene level. A mutation-based detection panel calls
\emph{alleles}, so we asked three allele-level questions. We fetched every
somatic mutation in the four panel genes and in eight long, frequently mutated
passenger genes (TTN, MUC16, RYR1, LRP1B, CSMD3, USH2A, SYNE1, FLG) from seven
hg19 PDAC cohorts on cBioPortal (NPAT patients carrying at least one mutation in these genes; legacy \texttt{paad\_tcga} and
the two smaller MSK cohorts were excluded to avoid double-counting patients),
NROWS mutation records in total. Each of the NSNV distinct SNVs was queried by
hg19 HGVS identifier against MyVariant.info (batch POST), which returned CADD
phred, dbNSFP AlphaMissense and REVEL scores, and ClinVar RCV significance
(FOUND\% of SNV records resolved; \texttt{results/myvariant\_variants.csv}).
Trinucleotide context for every missense allele came from the UCSC Genome
Browser REST API (hg19); the reference base matched the cBioPortal allele for
all NCTX of NCTX alleles (\texttt{results/myvariant\_allele\_context.csv}).

\paragraph{Positive control (passes).} Panel missense alleles score far above
passenger missense on every predictor (Table~\ref{tab:myvariant}; AlphaMissense
AUROC AMAUC, $p=AMP$). This calibrates the annotations; had it failed, the
rest would be void.

\paragraph{Within-panel recurrence: AlphaMissense has no room left.} We then
asked whether a predictor ranks the panel's alleles by how many patients carry
them (recurrent $\ge3$ patients vs.\ singletons; the most recurrent are TOP).
AlphaMissense does not: pooled AUROC AMREC ($p=AMRECP$), and it fails in both
cohort splits. The reason is a ceiling: CEILLP\% of panel \emph{singleton}
missense alleles already exceed the published likely-pathogenic threshold
(0.564), and CEIL9\% exceed 0.9. REVEL and CADD do track recurrence in the
pooled data (AUROC REVREC and CADDREC), and the signal survives a mutability
control --- CpG transitions are only CPGF\% of panel alleles but are
enriched among recurrent ones ($\rho=CPGRHO$, $p=CPGP$), and restricting to
non-CpG alleles leaves REVEL $\rho=REVNON$ and CADD $\rho=CADDNON$. But the
signal does not replicate across the cohort split: it holds in MSK-IMPACT
(CADD AUROC CADDMSK) and is null in the pooled WES/WGS cohorts (CADD AUROC
CADDWES, $p=CADDWESP$, recurrent $=\ge2$ patients). REVEL is also trained on
disease-mutation databases that contain the TP53 hotspots, so its advantage is
partly circular. \textbf{Verdict:} no predictor gives a cohort-robust ranking of
panel alleles; the robust finding is the negative one, that AlphaMissense
cannot rank them at all. We name this the \emph{AlphaMissense panel ceiling} and state it falsifiably: in any new PDAC cohort with $\ge100$ recurrent panel missense alleles, AlphaMissense AUROC for recurrent vs singleton will stay below 0.60; an AUROC $\ge0.60$ with $p<0.01$ refutes it.

\begin{table}[h]\centering\footnotesize
\caption{Allele-level predictor audit. Columns: median score (panel / passenger
missense alleles); AUROC panel vs passenger; AUROC recurrent vs singleton
within panel (pooled, MSK-IMPACT, WES/WGS); Spearman $\rho$ recurrence vs score
on non-CpG panel alleles.}\label{tab:myvariant}
\begin{tabular}{p{2.0cm}p{0.9cm}p{0.9cm}p{1.0cm}p{1.0cm}p{2.0cm}p{2.0cm}p{2.0cm}}
\hline
Predictor & panel & pass. & drv/pass & pooled & MSK & WES/WGS & non-CpG $\rho$ \\ \hline
ROWS
\hline\end{tabular}\end{table}

\paragraph{CDKN2A is the exception.} Patient-weighted, KRAS KRASW\%, TP53 TP53W\% and
SMAD4 SMAD4W\% of missense carriers hold an AlphaMissense likely-pathogenic allele,
but only CDKNW\% of CDKN2A missense carriers do. One candidate explanation (tested and not supported in the next section) is that
the scores are defined on a single canonical protein, while CDKN2A also
encodes p14$^{ARF}$ from an alternate reading frame, so an allele's effect on
that product is not scored. We treat the low scores as a limit of the
predictor, not as evidence that those alleles are passengers.

\paragraph{How many panel calls rest on unsupported alleles?} Calling an
alteration \emph{supported} if it is truncating, ClinVar pathogenic/likely
pathogenic, or AlphaMissense likely-pathogenic, SUPP\% of the NPOS
panel-positive patients carry at least one supported alteration (ClinVar P/LP
alone: CLIN\%); NUNS patients carry only unsupported alterations (mostly
in-frame indels and low-scoring missense). The denominator is panel-positive
patients only, so this is a bound on annotation support, not on sensitivity.

\begin{figure}[h]\centering
\includegraphics[width=\linewidth]{figs/fig_myvariant.pdf}
\caption{MyVariant allele audit. A: AlphaMissense separates driver from
passenger missense. B: no predictor ranks panel alleles by recurrence
robustly across cohorts. C: CDKN2A missense carriers fall below the
AlphaMissense threshold far more often than the other panel genes.}
\label{fig:myvariant}\end{figure}
"""
rep_tok = {
    "NPAT": f"{npat:,}", "NROWS": f"{a['n_rows']:,}", "NSNV": f"{a['n_snv_ids']:,}",
    "FOUND": pct(a["myvariant_found_fraction_snv"]),
    "AMAUC": f3(q1["am"]["auroc_panel_vs_passenger"]), "AMP": pv(q1["am"]["mwu_p"]),
    "TOP": top, "AMRECP": pv(q2["panel"]["am"]["p_auroc"]), "AMREC": f3(q2["panel"]["am"]["auroc_recurrent_vs_singleton"]),
    "CEILLP": pct(ceil["panel_singletons_frac_am_lp"]), "CEIL9": pct(ceil["panel_singletons_frac_am_ge_0_9"]),
    "REVREC": f3(q2["panel"]["revel"]["auroc_recurrent_vs_singleton"]), "CADDREC": f3(q2["panel"]["cadd"]["auroc_recurrent_vs_singleton"]),
    "CPGF": pct(cpg["frac_cpg_transition"]), "CPGRHO": f2(cpg["spearman_recurrence_vs_cpg_ts"]), "CPGP": pv(cpg["p"]),
    "REVNON": f2(cpg["revel"]["noncpg_only_rho"]), "CADDNON": f2(cpg["cadd"]["noncpg_only_rho"]),
    "CADDMSK": f3(rep["msk_impact"]["cadd"]["auroc_recurrent_vs_singleton"]),
    "CADDWESP": pv(rep["wes_cohorts"]["cadd"]["p"]), "CADDWES": f3(rep["wes_cohorts"]["cadd"]["auroc_recurrent_vs_singleton"]),
    "KRASW": pct(pg["KRAS"]["frac_patient_weighted_am_lp"]), "TP53W": pct(pg["TP53"]["frac_patient_weighted_am_lp"]),
    "SMAD4W": pct(pg["SMAD4"]["frac_patient_weighted_am_lp"]), "CDKNW": pct(pg["CDKN2A"]["frac_patient_weighted_am_lp"]),
    "SUPP": pct(q3["frac_panel_positive_with_supported_allele"]), "NPOS": f"{q3['n_panel_positive_patients']:,}",
    "CLIN": pct(q3["frac_panel_positive_with_clinvar_plp"]), "NUNS": str(q3["n_panel_positive_unsupported_only"]),
    "ROWS": "\n".join(rows),
}
import csv
nctx = sum(1 for _ in open("results/myvariant_allele_context.csv")) - 1
nm = sum(int(r["ref_matches"]) for r in csv.DictReader(open("results/myvariant_allele_context.csv")))
assert nm == nctx
rep_tok["NCTX"] = f"{nctx:,}"
for k in sorted(rep_tok, key=len, reverse=True):
    tex = tex.replace(k, rep_tok[k])
open("paper/myvariant_sec.tex", "w").write(tex)
print("wrote paper/myvariant_sec.tex")
