#!/usr/bin/env python3
"""Figure for the Open Targets panel audit (reads committed JSON only)."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

a = json.load(open("results/opentargets_panel_audit.json"))
DIS = ["MONDO_0005192", "MONDO_0007254", "EFO_0004251"]
DLAB = {"MONDO_0005192": "exocrine pancreatic\ncarcinoma (n=11,589)",
        "MONDO_0007254": "breast cancer\n(n=17,252)",
        "EFO_0004251": "myeloproliferative\ndisorder (n=13,488)"}
GENES = ["KRAS", "TP53", "CDKN2A", "SMAD4", "BRCA1", "JAK2"]
COL = {"KRAS": "#c0392b", "TP53": "#c0392b", "CDKN2A": "#c0392b",
       "SMAD4": "#c0392b", "BRCA1": "#2980b9", "JAK2": "#2980b9"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.4))
w, x = 0.26, range(len(GENES))
for i, dis in enumerate(DIS):
    vals = [a["per_gene"][g]["diseases"][dis].get("rank",
            a["diseases"][dis]["n_targets"]) for g in GENES]
    ax1.bar([xi + (i - 1) * w for xi in x], vals, w,
            label=DLAB[dis], color=["#c0392b", "#f1b0b7", "#7f8c8d"][i])
ax1.set_yscale("log"); ax1.invert_yaxis()
ax1.set_xticks(list(x)); ax1.set_xticklabels(GENES, fontsize=8)
ax1.set_ylabel("association rank (lower = stronger)")
ax1.legend(fontsize=6.5, frameon=False)
ax1.set_title("A. Panel genes top the pancreatic table", fontsize=9)

import csv
bx, by = [], []
with open("results/opentargets_assoc_rows.csv") as f:
    for r in csv.DictReader(f):
        if r["disease"] != "MONDO_0005192":
            continue
        bx.append(float(r["europepmc"])); by.append(float(r["genetic_max"]))
ax2.scatter(bx, by, s=1, c="#bbbbbb", alpha=0.25, rasterized=True,
            label="all targets")
for g in GENES:
    d = a["per_gene"][g]["diseases"]["MONDO_0005192"]
    ax2.scatter([d["europepmc"]], [d["genetic_max"]], s=28, c=COL[g],
                marker="D" if COL[g] == "#2980b9" else "o", zorder=5)
    ax2.annotate(g, (d["europepmc"], d["genetic_max"]), fontsize=6.5,
                 xytext=(3, 3), textcoords="offset points")
ax2.set_xlabel("literature score (Europe PMC)")
ax2.set_ylabel("max genetic-datasource score")
ax2.set_title("B. Pancreatic association is genetically backed", fontsize=9)
fig.tight_layout()
fig.savefig("paper/figs/fig_opentargets.pdf")
print("wrote paper/figs/fig_opentargets.pdf")
