import json, csv, os

R = "/home/sandbox/repos/mega27-07-pancreatic-ai-and-codon-optimizer/results"
OUT = "/home/sandbox/repos/mega27-07-pancreatic-ai-and-codon-optimizer/paper/appendix.tex"

def esc(s):
    s = str(s)
    for a, b in [("\\","\\textbackslash "),("_","\\_"),("%","\\%"),("&","\\&"),("#","\\#"),("$","\\$"),("{","\\{"),("}","\\}"),("~","\\textasciitilde "),("^","\\textasciicircum ")]:
        s = s.replace(a, b)
    return s

def num(x, nd=3):
    try:
        f = float(x)
        if abs(f) < 1e-4 and f != 0:
            return f"{f:.2e}"
        return f"{f:.{nd}f}"
    except Exception:
        return esc(x)

L = []
A = L.append

A(r"\appendix")
A(r"\section{Numbered formula compendium}")
A(r"Every quantitative claim in the body traces to one of the following definitions. Equation numbering continues from the body.")
A(r"\subsection{Codon adaptation index}")
A(r"For codon $c$ encoding amino acid $a$, the relative adaptiveness is")
A(r"\begin{equation} w_c = \frac{f_c}{\max_{c' \in \mathrm{syn}(a)} f_{c'}}, \end{equation}")
A(r"where $f_c$ is the codon's frequency in the reference set of highly expressed genes. CAI of a coding sequence of length $L$ codons is the geometric mean")
A(r"\begin{equation} \mathrm{CAI} = \left( \prod_{i=1}^{L} w_{c_i} \right)^{1/L}. \end{equation}")
A(r"The high-frequency-codon (HFC) baseline sets every $w_{c_i}=1$ by construction, giving")
A(r"\begin{equation} \mathrm{CAI}_{\mathrm{HFC}} = 1 \quad \text{by definition}, \end{equation}")
A(r"which is the circularity the body reports empirically (HFC reaches CAI 0.997 on the benchmark while its predicted expression trails CodonOpt).")
A(r"\subsection{Predicted-expression z-score}")
A(r"The CNN expression model outputs a score calibrated against PaxDb whole-proteome abundances; we standardize per gene as")
A(r"\begin{equation} z_g = \frac{s_g - \mu_{\mathrm{wt}}}{\sigma_{\mathrm{wt}}}, \end{equation}")
A(r"with $\mu_{\mathrm{wt}}, \sigma_{\mathrm{wt}}$ the mean and standard deviation of wild-type scores over the benchmark panel. The head-to-head metric is the mean $z_g$ over the 40 benchmark genes.")
A(r"\subsection{Co-mutation enrichment}")
A(r"For target alteration $T$ and covariate set $S$, the $2\times2$ contingency table over $n$ tumours is $[[a,b],[c,d]]$ with $a$ the count of tumours altered in both. The odds ratio and Fisher exact $p$-value are")
A(r"\begin{equation} \mathrm{OR} = \frac{a\,d}{b\,c}, \qquad p = \sum_{k:\,P(k)\le P(a)} \frac{\binom{a+b}{k}\binom{c+d}{a+c-k}}{\binom{n}{a+c}}. \end{equation}")
A(r"Replication requires sign-consistency of $\log \mathrm{OR}$ and $p < 0.05$ in an independent cohort with no shared samples.")
A(r"\subsection{Permutation importance}")
A(r"For covariate $j$ and held-out metric $M$ (AUC), with $R=4$ permuted replicates,")
A(r"\begin{equation} \Delta_j = M_{\mathrm{base}} - \frac{1}{R}\sum_{r=1}^{R} M_{j}^{(r)}, \end{equation}")
A(r"ranked descending. Stability is the Spearman rank correlation between rankings from disjoint splits:")
A(r"\begin{equation} \rho = 1 - \frac{6 \sum_k d_k^2}{m(m^2-1)}, \end{equation}")
A(r"where $d_k$ is the rank difference of gene $k$ between the two half-ensembles and $m$ the number of ranked genes.")
A(r"\subsection{5$'$ mRNA structure}")
A(r"ViennaRNA fold reports the minimum free energy of the 5$'$ window of length $w$:")
A(r"\begin{equation} \mathrm{MFE}_g = \min_{s \in \mathcal{S}} E(s, \mathrm{seq}_g[1..w]), \end{equation}")
A(r"over the secondary-structure space $\mathcal{S}$; less negative MFE near the start codon is associated with higher initiation efficiency.")
A(r"\subsection{Bootstrap confidence intervals}")
A(r"Paired bootstrap over $B=10{,}000$ resamples of the gene axis gives the percentile interval")
A(r"\begin{equation} \mathrm{CI}_{95} = \left[ Q_{0.025}(\Delta_b),\; Q_{0.975}(\Delta_b) \right],\quad \Delta_b = \bar z^{\mathrm{CodonOpt}}_b - \bar z^{\mathrm{ICOR}}_b. \end{equation}")

# Appendix B: 40-gene benchmark
A(r"\section{Per-gene benchmark table (ICOR 40-gene panel)}")
A(r"All six methods, all 40 genes. pred\_z is the predicted-expression z-score under our CNN; original = wild type. Full machine-readable record: \texttt{results/icor\_headtohead\_per\_gene.json}.")
d = json.load(open(f"{R}/icor_headtohead_per_gene.json"))["genes"]
methods = ["original", "icor", "gensmart", "HFC", "BFC", "URC", "ERC", "codonopt"]
avail = None
rows = []
for g, v in d.items():
    if avail is None:
        avail = [m for m in methods if any(m in gg for gg in d.values())]
    row = [esc(g)]
    for m in avail:
        s = v.get(m)
        row.append((num(s.get("cai")) + "/" + num(s.get("pred_expr_z"), 2)) if s else "--")
    rows.append(row)
hdr = "Gene & " + " & ".join(esc(m) + " CAI/$z$" for m in avail) + r" \\"
ncol = 1 + len(avail)
A(r"\begin{footnotesize}")
A(r"\begin{longtable}{l" + "l" * len(avail) + "}")
A(r"\caption{Per-gene codon-optimization benchmark, all methods. Cells show CAI/predicted-expression $z$.}\\")
A(r"\toprule " + hdr + r" \midrule \endfirsthead")
A(r"\toprule " + hdr + r" \midrule \endhead")
A(r"\midrule \multicolumn{" + str(ncol) + r"}{r}{continued}\\ \endfoot \bottomrule \endlastfoot")
for r_ in rows:
    A(" & ".join(r_) + r" \\")
A(r"\end{longtable}")
A(r"\end{footnotesize}")

# Appendix C: ensemble
A(r"\section{Seed-ensemble attribution audit (16 splits)}")
e = json.load(open(f"{R}/cdkn2a_attribution_ensemble.json"))
A(r"Per-seed held-out AUCs and the stability statistics reported in the body. Median held-out AUC = " + num(e["held_out_auc_median"]) + r"; split-half Spearman (8v8) = " + num(e["split_half_spearman_8v8"]) + r"; median single-split pairwise Spearman = " + num(e["single_split_rho_median"]) + r".")
A(r"\begin{center}\begin{tabular}{lr}")
A(r"\toprule Seed & held-out AUC \\ \midrule")
for s, a in zip(e["seeds"], e["held_out_aucs"]):
    A(f"{s} & {num(a)} \\\\")
A(r"\bottomrule\end{tabular}\end{center}")
A(r"Single-split pairwise Spearman values: " + ", ".join(num(x) for x in e["single_split_pairwise_spearman"]) + r".")
A(r"Ensemble top-10 genes by mean permutation importance: " + ", ".join(esc(g if isinstance(g, str) else g.get("gene", g)) for g in e["ensemble_top10"]) + r". The top genes are known long-gene passengers (TTN, MACF1, SYNE2 class); the rescue is rejected, as reported in the body.")

# Appendix D: replication tables
A(r"\section{Replication contingency tables}")
for f, label in [("cdkn2a_replication.json", "per-study"), ("cdkn2a_meta_replication.json", "meta")]:
    try:
        rep = json.load(open(f"{R}/{f}"))
    except Exception:
        continue
    studies = rep["per_study"] if isinstance(rep, dict) and "per_study" in rep else (rep if isinstance(rep, list) else [])
    for st in studies:
        A(r"\subsection{" + esc(st["study"]) + r" ($n=" + str(st["n"]) + r"$, CDKN2A-altered $=" + str(st["cdkn2a_mut_n"]) + r"$)}")
        A(r"\begin{center}\begin{tabular}{lrrrr}")
        A(r"\toprule Covariate & OR & $p$ & $a$ & table \\ \midrule")
        for cov, v in st.items():
            if isinstance(v, dict) and "table" in v:
                t = v["table"]
                A(esc(cov) + " & " + num(v["odds_ratio"], 2) + " & " + num(v["p"], 4) + " & " + str(t[0][0]) + " & " + esc(str(t)) + r" \\")
        A(r"\bottomrule\end{tabular}\end{center}")

A(r"\end{document-marker}")
open(OUT, "w").write("\n".join(L) + "\n")
print("wrote", OUT, len(L), "lines")
