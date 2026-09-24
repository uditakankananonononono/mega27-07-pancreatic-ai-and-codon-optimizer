#!/usr/bin/env python3
"""Generate paper/structure_sec.tex from results/structure_panel_audit.json (token replace)."""
import json
a = json.load(open("results/structure_panel_audit.json"))
def pv(x):
    s = f"{x:.2g}"
    return s.replace("e-0", "e-").replace("e-", r"\times10^{-") + "}" if "e-" in s else s
S = a["structures"]
c, pc, bu, lg = a["low_am_x_contact"], a["low_am_x_partner_contact"], a["low_am_x_buried"], a["logit"]
t, tp, tb = c["table_rows_highAM_lowAM_cols_no_yes"], pc["table_rows_highAM_lowAM_cols_no_yes"], bu["table_rows_highAM_lowAM_cols_no_yes"]
h83 = next(r for r in a["low_am_alleles_mapped"] if r["protein"] == "H83Y")
srow = "\n".join(f"{g} & {S[g]['pdb']} ({S[g]['chain']}) & {S[g]['n_residues']} & {100*S[g]['uniprot_identity']:.1f} & {S[g]['n_partner_contact']} & {S[g]['n_ligand_contact']} & {a['mapped_by_gene'].get(g,0)} \\\\"
                 for g in ("KRAS", "TP53", "CDKN2A", "SMAD4"))
tok = {"NMAP": str(a["n_alleles_mapped"]), "NPAN": str(a["n_panel_alleles"]),
       "RAM": f"{a['spearman_rsa_vs_am']:.3f}", "PAM": pv(a["p_rsa_am"]),
       "RRE": f"{a['spearman_rsa_vs_revel']:.3f}", "PRE": pv(a["p_rsa_revel"]),
       "RCA": f"{a['spearman_rsa_vs_cadd']:.3f}", "PCA": pv(a["p_rsa_cadd"]),
       "CLO": f"{t[1][1]}/{sum(t[1])}", "CHI": f"{t[0][1]}/{sum(t[0])}", "COR": f"{c['or']:.2f}", "CP": pv(c["p"]),
       "PLO": f"{tp[1][1]}/{sum(tp[1])}", "PHI": f"{tp[0][1]}/{sum(tp[0])}", "PPV": pv(pc["p"]),
       "BLO": f"{tb[1][1]}/{sum(tb[1])}", "BHI": f"{tb[0][1]}/{sum(tb[0])}", "BOR": f"{bu['or']:.2f}", "BP": pv(bu["p"]),
       "LRSA": f"{lg['rsa']['coef']:.2f}", "LRSAP": pv(lg["rsa"]["p"]), "LCON": f"{lg['contact']['coef']:.2f}", "LCONP": pv(lg["contact"]["p"]),
       "H83RSA": f"{h83['rsa']:.2f}", "SROWS": srow}
tex = r"""\section{Structural context of low-scoring panel alleles (RCSB PDB, negative)}
A natural reason for AlphaMissense to miss damaging alleles is that it models
a single chain, so damage at a partner interface might be under-weighted. We
stated this as a falsifiable hypothesis --- alleles below the likely-pathogenic
threshold are \emph{enriched} at partner or ligand contact residues --- and
tested it on four experimental complexes from the RCSB PDB
(Table~\ref{tab:structure}). For each residue we computed relative solvent
accessibility (RSA) of the isolated chain with the Biopython Shrake--Rupley
implementation (maximum ASA from Tien et al.\ 2013) and contact with any
partner chain, DNA, or bound ligand (any heavy atom within 5\,\AA, Biopython
NeighborSearch). Author numbering was checked residue by residue against the
UniProt canonical sequence, and an allele was mapped only when its reference
amino acid matched the structure: NMAP of NPAN panel missense alleles mapped
(\texttt{results/structure\_panel\_alleles.csv}).

\begin{table}[h]\centering\footnotesize
\caption{Structures used. Identity = share of modelled residues matching UniProt.}\label{tab:structure}
\begin{tabular}{p{1.4cm}p{1.6cm}p{1.2cm}p{1.4cm}p{1.4cm}p{1.4cm}p{1.4cm}}
\hline
Gene & PDB (chain) & residues & identity \% & partner contacts & ligand contacts & alleles mapped \\ \hline
SROWS
\hline\end{tabular}\end{table}

\paragraph{Control.} Buried residues score higher, as they should:
Spearman $\rho$(RSA, AlphaMissense) $=RAM$ ($p=PAM$); REVEL $\rho=RRE$
($p=PRE$); CADD, which has no structural input, shows none ($\rho=RCA$,
$p=PCA$).

\paragraph{Result: the hypothesis is falsified in the opposite direction.}
Low-scoring alleles are \emph{depleted} at contact residues: CLO versus CHI
for the rest (Fisher OR COR, $p=CP$; protein/DNA partner contacts alone PLO
versus PHI, $p=PPV$). They are also less often buried (BLO versus BHI, OR
BOR, $p=BP$). Both effects hold jointly with gene fixed effects in a logistic
model (RSA coefficient LRSA, $p=LRSAP$; contact coefficient LCON, $p=LCONP$).
\textbf{Verdict:} AlphaMissense's low scores fall mostly on exposed,
non-contact residues, where low damage is structurally plausible;
interface blindness does not explain them. The main exception is CDKN2A H83Y,
which is fully buried (RSA H83RSA) in the p16 ankyrin fold yet scored benign.
The CDK6--p16 structure is 3.4\,\AA, and the p53 interface here is DNA only
(other p53 partners are not modelled), so contact calls are a lower bound.

\begin{figure}[h]\centering
\includegraphics[width=0.85\linewidth]{figs/fig_structure.pdf}
\caption{A: AlphaMissense versus relative solvent accessibility for mapped
panel alleles. B: low-scoring alleles are less often at partner/ligand
contacts and less often buried.}\label{fig:structure}\end{figure}
"""
for k in sorted(tok, key=len, reverse=True):
    tex = tex.replace(k, tok[k])
open("paper/structure_sec.tex", "w").write(tex); print("ok")
