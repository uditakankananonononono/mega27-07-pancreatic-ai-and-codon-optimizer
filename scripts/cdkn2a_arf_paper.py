#!/usr/bin/env python3
"""Generate paper/cdkn2a_arf_sec.tex from results/cdkn2a_arf_frame.json (token replace)."""
import json
a = json.load(open("results/cdkn2a_arf_frame.json"))
c, h, tc = a["codon_position_control_shared_region"], a["cancerhotspots"], a["transcript_checks"]
def pv(x):
    s = f"{x:.2g}"
    return s.replace("e-0", "e-").replace("e-", r"\times10^{-") + "}" if "e-" in s else s
t, sr, st = a["table_low_am_x_arf_altered"], a["shared_region_only"]["table"], c["strata_tables"]
cp = c["low_am_by_codon_pos_table_rows_pos1_pos2_cols_high_low"]
h83 = next(r for r in a["low_am_alleles"] if r["protein_cbio"] == "H83Y")
tok = {
 "P16LEN": str(tc["p16_len_aa"]), "ARFLEN": str(tc["arf_len_aa"]), "SHBP": str(tc["shared_coding_bp"]),
 "NALL": str(a["n_alleles"]), "NEXC": str(a["n_excluded_not_p16_missense"]),
 "LOWALT": f"{t[1][1]}/{t[1][0]+t[1][1]}", "HIALT": f"{t[0][1]}/{t[0][0]+t[0][1]}",
 "FOR": f"{a['fisher_or']:.1f}", "FPV": pv(a["fisher_p"]),
 "SRLOW": f"{sr[1][1]}/{sr[1][0]+sr[1][1]}", "SRHI": f"{sr[0][1]}/{sr[0][0]+sr[0][1]}",
 "SROR": f"{a['shared_region_only']['fisher_or']:.1f}", "SRP": pv(a["shared_region_only"]["fisher_p"]),
 "POS1ARF": f"{100*c['frac_arf_altered_by_codon_pos_shared_region']['1']:.0f}",
 "POS2ARF": f"{100*c['frac_arf_altered_by_codon_pos_shared_region']['2']:.0f}",
 "CP1LOW": f"{cp[0][1]}/{cp[0][0]+cp[0][1]}", "CP2LOW": f"{cp[1][1]}/{cp[1][0]+cp[1][1]}", "CPP": pv(c["low_am_pos1_vs_pos2_fisher_p"]),
 "S1HI": f"{st['1'][0][1]}/{sum(st['1'][0])}", "S1LO": f"{st['1'][1][1]}/{sum(st['1'][1])}", "MHP": f"{c['p_mh']:.2f}",
 "H83ARF": h83["arf_change"], "H83N": str(h83["n_pat"]), "H83AM": f"{h83['am']:.3f}",
 "NHRES": str(sum(h["n_hotspot_residues"].values())), "NHAL": str(h["n_hotspot_alleles"]), "NPAN": str(h["n_panel_alleles"]),
 "HLOWAL": str(h["hotspot_alleles_below_am_lp"]), "HLOWPAT": str(h["hotspot_patients_below_am_lp"]), "HPAT": f"{h['hotspot_patients_total']:,}",
 "HAMA": f"{h['auroc_hotspot_vs_not_am']:.3f}", "HAMP": pv(h["p_am"]),
 "HREA": f"{h['auroc_hotspot_vs_not_revel']:.3f}", "HREP": pv(h["p_revel"]),
 "HCAA": f"{h['auroc_hotspot_vs_not_cadd']:.3f}", "HCAP": pv(h["p_cadd"]),
}
tex = r"""\section{Testing the p14\textsuperscript{ARF} explanation of the CDKN2A exception (negative)}
The allele-level section left one explanation untested: low AlphaMissense
scores on CDKN2A might mark alleles that act through p14$^{ARF}$, the second
protein encoded in an alternate reading frame. We tested it directly. From UCSC
hg19 refGene exon coordinates and the hg19 sequence we rebuilt the NM\_000077
(p16$^{INK4a}$) and NM\_058195 (p14$^{ARF}$) coding sequences; both translate
to full-length proteins (P16LEN and ARFLEN aa, Met start, one terminal stop)
and share SHBP coding bp in exon 2. Each of the NALL CDKN2A missense SNVs that
match their cBioPortal p16 call (NEXC excluded as annotated on other
transcripts) was translated in both frames (\texttt{results/cdkn2a\_arf\_alleles.csv}).

\paragraph{The raw association is strong.} Low-AlphaMissense alleles alter
p14$^{ARF}$ in LOWALT cases versus HIALT for the rest (Fisher OR FOR,
$p=FPV$), and inside the shared exon-2 region SRLOW versus SRHI (OR SROR,
$p=SRP$). The most common low-scoring allele, p16 H83Y (H83N patients,
AlphaMissense H83AM, ClinVar P/LP), is ARF H83ARF.

\paragraph{It is explained by reading-frame geometry.} In the shared region
the ARF frame is offset by one base, so a change at p16 codon position 1 alters
ARF in POS1ARF\% of alleles and a change at position 2 in POS2ARF\%.
Independently, low AlphaMissense scores concentrate at p16 codon position 1
(CP1LOW alleles vs CP2LOW at position 2, $p=CPP$). Stratified by codon
position the association vanishes: in the only informative stratum (position
1) ARF is altered in S1HI high-scoring and S1LO low-scoring alleles
(Mantel--Haenszel $p=MHP$). \textbf{Verdict:} the data cannot separate an ARF
mechanism from frame geometry; the explanation is not supported, and the
CDKN2A exception is best described as AlphaMissense rating many
codon-position-1 substitutions (e.g.\ H83Y) as benign.

\paragraph{Independent recurrence definition.} We repeated the recurrence test
with residue hotspots from cancerhotspots.org (NHRES single-residue hotspots in
the four panel genes; statistically defined on large pan-cancer cohorts that
overlap ours, so this is a check of definition, not of data). NHAL of NPAN
panel missense alleles sit at a hotspot residue. AlphaMissense does not
separate hotspot from non-hotspot residues (AUROC HAMA, $p=HAMP$), while
REVEL (HREA, $p=HREP$) and CADD (HCAA, $p=HCAP$) do, reproducing the
AlphaMissense panel ceiling. HLOWAL hotspot alleles (HLOWPAT of HPAT hotspot
carriers) score below the likely-pathogenic threshold, led by CDKN2A H83Y.

\begin{figure}[h]\centering
\includegraphics[width=0.85\linewidth]{figs/fig_cdkn2a_arf.pdf}
\caption{A: fraction of CDKN2A missense alleles that alter p14$^{ARF}$, split by
p16 codon position and AlphaMissense class (shared exon-2 region). B:
AlphaMissense does not separate cancerhotspots.org residues from the rest of
the panel; REVEL and CADD do.}\label{fig:cdkn2aarf}\end{figure}
"""
for k in sorted(tok, key=len, reverse=True):
    tex = tex.replace(k, tok[k])
open("paper/cdkn2a_arf_sec.tex", "w").write(tex); print("ok")
