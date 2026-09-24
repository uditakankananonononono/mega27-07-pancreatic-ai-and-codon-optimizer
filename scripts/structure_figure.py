#!/usr/bin/env python3
"""Figure for the structural-context audit (committed results only)."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

a = json.load(open("results/structure_panel_audit.json"))
m = pd.read_csv("results/structure_panel_alleles.csv")
fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.2))
ax = axs[0]
for g, col in (("TP53", "#c0392b"), ("SMAD4", "#e67e22"), ("CDKN2A", "#8e44ad"), ("KRAS", "#2980b9")):
    s = m[m.gene == g]
    ax.scatter(s.rsa, s.am, s=8, alpha=0.6, c=col, label=f"{g} ({len(s)})")
ax.axhline(0.564, ls="--", c="k", lw=0.8)
ax.set_xlabel("relative solvent accessibility"); ax.set_ylabel("AlphaMissense")
ax.legend(fontsize=6.5, frameon=False, loc="lower left")
ax.set_title("A. Score vs exposure (rho %.2f)" % a["spearman_rsa_vs_am"], fontsize=9)
ax = axs[1]
t = a["low_am_x_contact"]["table_rows_highAM_lowAM_cols_no_yes"]
b = a["low_am_x_buried"]["table_rows_highAM_lowAM_cols_no_yes"]
vals = [[t[0][1] / sum(t[0]), t[1][1] / sum(t[1])], [b[0][1] / sum(b[0]), b[1][1] / sum(b[1])]]
x = [0, 1]
ax.bar([i - 0.18 for i in x], [v[0] for v in vals], 0.36, color="#c0392b", label="AM >= 0.564")
ax.bar([i + 0.18 for i in x], [v[1] for v in vals], 0.36, color="#8e44ad", label="AM < 0.564")
ax.set_xticks(x); ax.set_xticklabels(["partner/ligand contact", "buried (RSA<0.2)"], fontsize=8)
ax.set_ylabel("fraction of alleles"); ax.set_ylim(0, 1); ax.legend(fontsize=7, frameon=False)
ax.set_title("B. Low-AM alleles avoid contacts", fontsize=9)
plt.tight_layout(); plt.savefig("paper/figs/fig_structure.pdf"); print("ok")
