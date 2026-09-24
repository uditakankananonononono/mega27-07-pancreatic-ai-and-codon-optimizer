#!/usr/bin/env python3
"""Structural context of PDAC panel missense alleles vs AlphaMissense scores.

Hypothesis (falsifiable): alleles AlphaMissense rates below its likely-
pathogenic threshold (0.564) sit disproportionately at partner-contact
residues (protein/DNA interfaces or bound ligand), where a single-chain
evolutionary+structure model may under-weight damage.
Structures (RCSB PDB, mmCIF): KRAS-RAF1 RBD 6VJJ (KRAS chain A; partner B;
ligands GNP/MG), p53 core-DNA 1TSR (p53 chain B; DNA chains E,F; ZN),
CDK6-p16 1BI7 (p16 chain B; CDK6 chain A), SMAD3/SMAD4 MH2 1U7F (SMAD4 chain
B; SMAD3 chains A,C). Numbering verified residue-by-residue against the
UniProt canonical sequence (UniProt REST); domains from UniProt features.
Per residue: relative solvent accessibility (Biopython Shrake-Rupley on the
isolated target chain, Tien et al. 2013 max ASA) and partner contact (any
heavy atom within 5 A of a partner chain / ligand, Biopython NeighborSearch).
Raw files -> data/structure/ (untracked). Writes results/structure_panel_audit.json
and results/structure_panel_residues.csv."""
import json, os, time, urllib.request
import numpy as np, pandas as pd
from scipy import stats
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import three_to_index, index_to_one

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data", "structure")
RES = os.path.join(HERE, "..", "results")
AM_LP = 0.564
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225, "G": 104, "H": 224, "I": 197,
          "L": 201, "K": 236, "M": 224, "F": 240, "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}
CFG = {"KRAS": ("6VJJ", "A", ["B"], ["GNP", "MG"], "P01116"),
       "TP53": ("1TSR", "B", ["E", "F"], ["ZN"], "P04637"),
       "CDKN2A": ("1BI7", "B", ["A"], [], "P42771"),
       "SMAD4": ("1U7F", "B", ["A", "C"], [], "Q13485")}


def fetch(url, dest, js=False):
    dest = os.path.join(RAW, dest)
    if not os.path.exists(dest):
        for a in range(4):
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=120) as r:
                    open(dest, "wb").write(r.read()); break
            except Exception:
                if a == 3: raise
                time.sleep(4)
    return json.load(open(dest)) if js else dest


def residue_table(gene):
    pdb, ch, partners, ligs, acc = CFG[gene]
    path = fetch(f"https://files.rcsb.org/download/{pdb}.cif", f"{pdb}.cif")
    up = fetch(f"https://rest.uniprot.org/uniprotkb/{acc}.json", f"uniprot_{acc}.json", js=True)
    seq = up["sequence"]["value"]
    doms = [(f["type"], f.get("description", ""), f["location"]["start"]["value"], f["location"]["end"]["value"])
            for f in up.get("features", []) if f["type"] in ("Domain", "Repeat", "DNA binding", "Region", "Zinc finger")]
    model = next(iter(MMCIFParser(QUIET=True).get_structure(pdb, path)))
    chain = model[ch]
    # isolated-chain SASA
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    s2 = Structure("x"); m2 = Model(0); s2.add(m2); m2.add(chain.copy())
    ShrakeRupley().compute(s2[0], level="R")
    iso = s2[0][ch]
    partner_atoms = [a for p in partners for a in model[p].get_atoms() if a.element != "H"]
    lig_atoms = [a for c in model for r in c for a in r if r.id[0] != " " and r.get_resname() in ligs]
    ns_p = NeighborSearch(partner_atoms) if partner_atoms else None
    ns_l = NeighborSearch(lig_atoms) if lig_atoms else None
    rows, match, tot = [], 0, 0
    for r in chain:
        if r.id[0] != " ":
            continue
        try:
            aa = index_to_one(three_to_index(r.get_resname()))
        except Exception:
            continue
        n = r.id[1]; tot += 1
        ok = 1 <= n <= len(seq) and seq[n - 1] == aa
        match += ok
        heavy = [a for a in r if a.element != "H"]
        cp = any(ns_p.search(a.coord, 5.0) for a in heavy) if ns_p else False
        cl = any(ns_l.search(a.coord, 5.0) for a in heavy) if ns_l else False
        rsa = min(1.0, iso[r.id].sasa / MAXASA[aa])
        dom = ";".join(d[1] or d[0] for d in doms if d[2] <= n <= d[3])
        rows.append({"gene": gene, "pdb": pdb, "resnum": n, "aa": aa, "uniprot_match": int(ok), "rsa": rsa,
                     "partner_contact": int(cp), "ligand_contact": int(cl), "domain": dom})
    return rows, {"pdb": pdb, "chain": ch, "n_residues": tot, "uniprot_identity": match / tot,
                  "n_partner_contact": sum(r["partner_contact"] for r in rows),
                  "n_ligand_contact": sum(r["ligand_contact"] for r in rows)}


def main():
    os.makedirs(RAW, exist_ok=True)
    allrows, meta = [], {}
    for g in CFG:
        rows, m = residue_table(g); allrows += rows; meta[g] = m
    rt = pd.DataFrame(allrows)
    rt.to_csv(os.path.join(RES, "structure_panel_residues.csv"), index=False)
    al = pd.read_csv(os.path.join(RES, "myvariant_missense_alleles.csv"))
    al = al[al.group == "panel"].dropna(subset=["am"]).copy()
    al["pos"] = al.protein.astype(str).str.extract(r"^[A-Z](\d+)[A-Z]$")[0].astype(float)
    al["refaa"] = al.protein.astype(str).str[0]
    m = al.merge(rt, left_on=["gene", "pos"], right_on=["gene", "resnum"], how="inner")
    m = m[m.refaa == m.aa].copy()          # reference residue must match structure
    m["low_am"] = (m.am < AM_LP).astype(int)
    m["contact"] = ((m.partner_contact == 1) | (m.ligand_contact == 1)).astype(int)
    m["buried"] = (m.rsa < 0.2).astype(int)
    out = {"structures": meta, "n_alleles_mapped": int(len(m)), "n_panel_alleles": int(len(al)),
           "mapped_by_gene": m.gene.value_counts().to_dict()}
    # positive control: AM and burial (buried residues should score higher)
    for sc in ("am", "revel", "cadd"):
        rho, p = stats.spearmanr(m.rsa, m[sc], nan_policy="omit")
        out[f"spearman_rsa_vs_{sc}"] = float(rho); out[f"p_rsa_{sc}"] = float(p)
    def fisher(col):
        t = pd.crosstab(m.low_am, m[col]).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
        o, p = stats.fisher_exact(t.values)
        return {"table_rows_highAM_lowAM_cols_no_yes": t.values.tolist(), "or": float(o), "p": float(p)}
    out["low_am_x_contact"] = fisher("contact")
    out["low_am_x_partner_contact"] = fisher("partner_contact")
    out["low_am_x_buried"] = fisher("buried")
    # per gene contact test (where both classes exist)
    pg = {}
    for g, s in m.groupby("gene"):
        t = pd.crosstab(s.low_am, s.contact).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
        o, p = stats.fisher_exact(t.values)
        pg[g] = {"n": int(len(s)), "n_low_am": int(s.low_am.sum()), "table": t.values.tolist(), "or": float(o), "p": float(p),
                 "median_rsa_low_am": float(s[s.low_am == 1].rsa.median()) if s.low_am.sum() else None,
                 "median_rsa_high_am": float(s[s.low_am == 0].rsa.median())}
    out["per_gene"] = pg
    # logistic: low_am ~ rsa + contact + gene
    import statsmodels.formula.api as smf
    try:
        fit = smf.logit("low_am ~ rsa + contact + C(gene)", data=m).fit(disp=0)
        out["logit"] = {k: {"coef": float(fit.params[k]), "p": float(fit.pvalues[k])} for k in ("rsa", "contact")}
    except Exception as e:
        out["logit"] = {"error": repr(e)[:200]}
    out["low_am_alleles_mapped"] = m[m.low_am == 1].sort_values("n_pat", ascending=False)[
        ["gene", "protein", "n_pat", "am", "rsa", "partner_contact", "ligand_contact", "domain"]].to_dict("records")
    json.dump(out, open(os.path.join(RES, "structure_panel_audit.json"), "w"), indent=1, default=str)
    m.to_csv(os.path.join(RES, "structure_panel_alleles.csv"), index=False)
    print(json.dumps({k: v for k, v in out.items() if k != "low_am_alleles_mapped"}, indent=1, default=str))
    print(pd.DataFrame(out["low_am_alleles_mapped"]).to_string())


if __name__ == "__main__":
    main()
