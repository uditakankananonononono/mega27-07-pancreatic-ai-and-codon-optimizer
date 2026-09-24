"""Extract committed subsets for the PDAC-driver feasibility audit.

Sources (fetched 2026-09-25, cached /tmp): DepMap Public 24Q4 CRISPRGeneEffect.csv
+ Model.csv (figshare 27993248, files 51064667/51065297); Human Protein Atlas
rna_celline.tsv.zip. Writes results/depmap_panel_subset.csv (4 PDAC drivers +
3 pan-essential controls x 1,178 models with lineage) and
results/hpa_pancreatic_lines.tsv (HPA nTPM for the 4 drivers in pancreatic
cell lines + all-line medians).
"""
import pandas as pd

GENES = ["KRAS", "TP53", "CDKN2A", "SMAD4"]
CTRL = ["RPS3", "PCNA", "PSMA1"]
PANC_KEY = ("PANC", "Capan", "CAPAN", "AsPC", "BxPC", "CFPAC", "HPAF", "HuP",
            "MIA PaCa", "PL45", "PaTu", "Panc ", "SW 1990", "SU.86")


def is_pancreatic(name):
    return any(k in name for k in PANC_KEY)


def main():
    hdr = pd.read_csv("/tmp/CRISPRGeneEffect.csv", nrows=0).columns
    sym = {c: c.split(" (")[0] for c in hdr[1:]}
    cols = [hdr[0]] + [c for c in hdr[1:] if sym[c] in set(GENES + CTRL)]
    df = pd.read_csv("/tmp/CRISPRGeneEffect.csv", usecols=cols).rename(columns={hdr[0]: "ModelID", **sym})
    model = pd.read_csv("/tmp/Model.csv", usecols=["ModelID", "OncotreeLineage"])
    df = model.merge(df, on="ModelID", how="inner")
    df.to_csv("results/depmap_panel_subset.csv", index=False)

    hpa = pd.read_csv("/tmp/rna_celline.tsv.zip", sep="\t", compression="zip",
                      usecols=["Gene name", "Cell line", "nTPM"])
    hpa = hpa[hpa["Gene name"].isin(GENES)]
    hpa["pancreatic"] = hpa["Cell line"].map(is_pancreatic)
    hpa.to_csv("results/hpa_pancreatic_lines.tsv", sep="\t", index=False)
    print("depmap", df.shape, "; hpa rows", len(hpa), "; pancreatic lines",
          hpa.loc[hpa.pancreatic, "Cell line"].nunique())


if __name__ == "__main__":
    main()
