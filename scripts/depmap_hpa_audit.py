"""PDAC-driver feasibility audit: DepMap 24Q4 dependency + HPA cell-line expression.

Questions: (1) Are KRAS/TP53/CDKN2A/SMAD4 dependencies pancreatic-selective?
Tumor-suppressor loss (TP53/CDKN2A/SMAD4) should be tolerated (effect >= 0),
which is what makes them sustainable, detectable tumor states; KRAS is the
expected pancreatic dependency. (2) Are the four driver transcripts measurable
in pancreatic cell lines (HPA nTPM), and does CDKN2A zero-expression in lines
mirror its deletion prevalence? Reads the committed subsets only.
"""
import json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

GENES = ["KRAS", "TP53", "CDKN2A", "SMAD4"]
CTRL = ["RPS3", "PCNA", "PSMA1"]
THR = -0.5


def dep_stats(v_p, v_o):
    return {"median_pancreas": float(np.median(v_p)), "median_other": float(np.median(v_o)),
            "frac_dep_pancreas": float((v_p < THR).mean()),
            "frac_dep_other": float((v_o < THR).mean()),
            "mw_p_two_sided": float(mannwhitneyu(v_p, v_o).pvalue), "n_pancreas": int(len(v_p))}


def dep_table(df, genes):
    panc = df.OncotreeLineage == "Pancreas"
    return {g: dep_stats(df.loc[panc, g].dropna().to_numpy(float),
                         df.loc[~panc, g].dropna().to_numpy(float)) for g in genes}


def expr_table(hpa, genes):
    out = {}
    for g in genes:
        a = hpa[hpa["Gene name"] == g]
        p = a[a.pancreatic]
        out[g] = {"median_ntpm_pancreatic_lines": float(p.nTPM.median()),
                  "median_ntpm_all_lines": float(a.nTPM.median()),
                  "frac_zero_pancreatic": float((p.nTPM == 0).mean()),
                  "frac_zero_all": float((a.nTPM == 0).mean()),
                  "n_pancreatic_lines": int(p["Cell line"].nunique())}
    return out


def main():
    df = pd.read_csv("results/depmap_panel_subset.csv")
    hpa = pd.read_csv("results/hpa_pancreatic_lines.tsv", sep="\t")
    out = {
        "design": __doc__.strip().splitlines()[0],
        "sources": {"depmap": "Public 24Q4 CRISPRGeneEffect.csv + Model.csv (figshare 27993248)",
                    "hpa": "rna_celline.tsv.zip (last-modified 2025-11-05)"},
        "threshold": THR,
        "n_models": int(len(df)), "n_pancreas_models": int((df.OncotreeLineage == "Pancreas").sum()),
        "dependency": dep_table(df, GENES), "dependency_controls": dep_table(df, CTRL),
        "expression": expr_table(hpa, GENES),
    }
    json.dump(out, open("results/depmap_hpa_panel_audit.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
