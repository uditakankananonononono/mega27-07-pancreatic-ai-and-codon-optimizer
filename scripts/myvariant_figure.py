#!/usr/bin/env python3
"""Figure for the MyVariant allele-level audit (reads committed results only)."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

a = json.load(open("results/myvariant_panel_audit.json"))
al = pd.read_csv("results/myvariant_missense_alleles.csv")
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.3))
ax = axs[0]
for grp, col in (("panel", "#c0392b"), ("passenger", "#7f8c8d")):
    ax.hist(al[al.group == grp].am.dropna(), bins=25, range=(0, 1), alpha=0.6, color=col,
            label=f"{grp} (n={int((al.group == grp).sum())})", density=True)
ax.axvline(0.564, ls="--", c="k", lw=0.8)
ax.set_xlabel("AlphaMissense score"); ax.set_ylabel("density")
ax.legend(fontsize=7, frameon=False)
ax.set_title("A. Driver vs passenger (AUROC %.2f)" % a["q1_positive_control"]["am"]["auroc_panel_vs_passenger"], fontsize=9)

ax = axs[1]
q2 = a["q2_recurrence"]; rep = q2["cohort_split_replication"]
sets = [("pooled", {s: q2["panel"][s]["auroc_recurrent_vs_singleton"] for s in ("am", "revel", "cadd")}),
        ("MSK-IMPACT", {s: rep["msk_impact"][s]["auroc_recurrent_vs_singleton"] for s in ("am", "revel", "cadd")}),
        ("WES/WGS", {s: rep["wes_cohorts"][s]["auroc_recurrent_vs_singleton"] for s in ("am", "revel", "cadd")})]
w = 0.26
for i, (sc, col) in enumerate((("am", "#8e44ad"), ("revel", "#2980b9"), ("cadd", "#27ae60"))):
    ax.bar([k + (i - 1) * w for k in range(3)], [v[sc] for _, v in sets], w, color=col,
           label={"am": "AlphaMissense", "revel": "REVEL", "cadd": "CADD"}[sc])
ax.axhline(0.5, c="k", lw=0.8)
ax.set_ylim(0.4, 0.75); ax.set_xticks(range(3)); ax.set_xticklabels([n for n, _ in sets], fontsize=8)
ax.set_ylabel("AUROC recurrent vs singleton"); ax.legend(fontsize=7, frameon=False)
ax.set_title("B. Within-panel recurrence", fontsize=9)

ax = axs[2]
pg = q2["per_gene_am"]; genes = ["KRAS", "TP53", "SMAD4", "CDKN2A"]
ax.bar(genes, [pg[g]["frac_patient_weighted_am_lp"] for g in genes], color="#c0392b")
ax.set_ylim(0, 1.05); ax.set_ylabel("patient-weighted frac. AM likely-path.")
ax.set_title("C. Missense carriers above AM threshold", fontsize=9)
plt.tight_layout()
plt.savefig("paper/figs/fig_myvariant.pdf")
print("wrote paper/figs/fig_myvariant.pdf")
