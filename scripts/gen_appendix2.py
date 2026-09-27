import json, csv

R = "/home/sandbox/repos/mega27-07-pancreatic-ai-and-codon-optimizer/results"
OUT = "/home/sandbox/repos/mega27-07-pancreatic-ai-and-codon-optimizer/paper/appendix_tables.tex"

def esc(s):
    s = str(s)
    for a, b in [("\\","\\textbackslash "),("_","\\_"),("%","\\%"),("&","\\&"),("#","\\#"),("$","\\$"),("{","\\{"),("}","\\}")]:
        s = s.replace(a, b)
    return s[:60]

def num(x, nd=3):
    try:
        f = float(x)
        if f == 0: return "0"
        if abs(f) < 1e-4 or abs(f) > 1e5: return f"{f:.2e}"
        return f"{f:.{nd}f}"
    except Exception:
        return esc(x)

L = []
A = L.append

# gnomAD top 150 by total AC
A(r"\section{gnomAD r4 panel variants (top 150 by allele count)}")
A(r"Population-frequency audit of the detection panel. Full 9{,}410-row record: \texttt{results/gnomad\_panel\_variants.csv}.")
rows = list(csv.DictReader(open(f"{R}/gnomad_panel_variants.csv")))
def ac(r):
    try: return int(r["ac_exome"] or 0) + int(r["ac_genome"] or 0)
    except Exception: return 0
rows.sort(key=ac, reverse=True)
A(r"\begin{footnotesize}\begin{longtable}{llllrrr}")
A(r"\caption{Top 150 panel variants by combined gnomAD allele count.}\\ \toprule")
A(r"Gene & Variant & Consequence & rsID & AF (exome) & AC & AN \\ \midrule \endfirsthead")
A(r"\toprule Gene & Variant & Consequence & rsID & AF (exome) & AC & AN \\ \midrule \endhead")
A(r"\midrule \multicolumn{7}{r}{continued}\\ \endfoot \bottomrule \endlastfoot")
for r in rows[:150]:
    A(" & ".join([esc(r["gene"]), esc(r["variant_id"]), esc(r["consequence"]), esc(r["rsids"] or "--"), num(r["af_exome"], 6) if r["af_exome"] else "--", r["ac_exome"] or "--", r["an_exome"] or "--"]) + r" \\")
A(r"\end{longtable}\end{footnotesize}")

# CIViC: pancreatic-first, 150
A(r"\section{CIViC evidence rows for the driver panel}")
A(r"Curated clinical evidence supporting panel genes. Full 1{,}153-row record: \texttt{results/civic\_panel\_evidence\_rows.csv}. Pancreatic-relevant rows first.")
rows = list(csv.DictReader(open(f"{R}/civic_panel_evidence_rows.csv")))
rows.sort(key=lambda r: -int(r["pancreatic"] or 0))
npanc = sum(1 for r in rows if r["pancreatic"] == "1")
A(f"Of {{len(rows)}} evidence rows, {npanc} are pancreatic-cancer relevant.".replace("{len(rows)}", str(len(rows))))
A(r"\begin{footnotesize}\begin{longtable}{llllllp{4.2cm}}")
A(r"\caption{CIViC evidence rows (pancreatic-relevant first, top 150 shown).}\\ \toprule")
A(r"EID & Gene & Level & Type & Significance & Panc. & Molecular profile \\ \midrule \endfirsthead")
A(r"\toprule EID & Gene & Level & Type & Significance & Panc. & Molecular profile \\ \midrule \endhead")
A(r"\midrule \multicolumn{7}{r}{continued}\\ \endfoot \bottomrule \endlastfoot")
for r in rows[:150]:
    A(" & ".join([esc(r["eid"]), esc(r["gene"]), esc(r["evidence_level"]), esc(r["evidence_type"][:12]), esc(r["significance"][:18]), r["pancreatic"], esc(r["molecular_profile"][:40])]) + r" \\")
A(r"\end{longtable}\end{footnotesize}")

# DepMap sample 150
A(r"\section{DepMap 24Q4 panel gene-effect sample}")
A(r"Chronos gene-effect scores for the seven panel genes across screened lines (sample of 150 of 1{,}178). Full record: \texttt{results/depmap\_panel\_subset.csv}.")
rows = list(csv.DictReader(open(f"{R}/depmap_panel_subset.csv")))
genes = ["CDKN2A","KRAS","PCNA","PSMA1","RPS3","SMAD4","TP53"]
A(r"\begin{footnotesize}\begin{longtable}{llrrrrrrr}")
A(r"\caption{DepMap Chronos scores, sample.}\\ \toprule")
A(r"Model & Lineage & " + " & ".join(genes) + r" \\ \midrule \endfirsthead")
A(r"\toprule Model & Lineage & " + " & ".join(genes) + r" \\ \midrule \endhead")
A(r"\midrule \multicolumn{9}{r}{continued}\\ \endfoot \bottomrule \endlastfoot")
for r in rows[:150]:
    A(" & ".join([esc(r["ModelID"]), esc(r["OncotreeLineage"][:18])] + [num(r[g], 2) for g in genes]) + r" \\")
A(r"\end{longtable}\end{footnotesize}")

# GTEx all 433
A(r"\section{GTEx v8 tissue medians for the panel (complete)}")
A(r"Median TPM by tissue for panel genes; the complete 432-row record. Background for the expression-specificity discussion in the body.")
rows = list(csv.DictReader(open(f"{R}/gtex_panel_tissue_medians.csv")))
A(r"\begin{footnotesize}\begin{longtable}{llrr}")
A(r"\caption{GTEx median TPM by tissue, panel genes.}\\ \toprule")
A(r"Gene & Tissue & Median TPM & $n$ \\ \midrule \endfirsthead")
A(r"\toprule Gene & Tissue & Median TPM & $n$ \\ \midrule \endhead")
A(r"\midrule \multicolumn{4}{r}{continued}\\ \endfoot \bottomrule \endlastfoot")
for r in rows:
    A(" & ".join([esc(r["gene"]), esc(r["tissueSiteDetailId"]), num(r["median_tpm"], 2), r["n_samples"]]) + r" \\")
A(r"\end{longtable}\end{footnotesize}")

open(OUT, "w").write("\n".join(L) + "\n")
print("wrote", len(L), "lines")
