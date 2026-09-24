#!/usr/bin/env python3
"""Figure for the GTEx panel-background audit (from committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = json.load(open("results/gtex_panel_background.json"))
G = D["genes"]
ORDER = ["KRAS", "TP53", "CDKN2A", "SMAD4", "GAPDH", "ACTB", "PRSS1", "INS"]
COLS = {"panel": "#c0392b", "housekeeping": "#7f8c8d", "pancreas_specific": "#2980b9"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.6))
x = np.arange(len(ORDER)); w = 0.38
pan = [max(G[g]["pancreas_median_tpm"], 1e-3) for g in ORDER]
bld = [max(G[g]["whole_blood_median_tpm"], 1e-3) for g in ORDER]
ax1.bar(x - w/2, pan, w, color="#e67e22", label="Pancreas (target organ)")
ax1.bar(x + w/2, bld, w, color="#2c3e50", label="Whole blood (sampling compartment)")
ax1.set_yscale("log"); ax1.set_ylabel("median TPM (log scale)")
ax1.set_xticks(x); ax1.set_xticklabels(ORDER, rotation=45, ha="right", fontsize=8)
ax1.legend(fontsize=7, frameon=False); ax1.set_title("Normal background: target organ vs blood", fontsize=9)

taus = [G[g]["tau"] for g in ORDER]
ax2.bar(x, taus, color=[COLS[G[g]["group"]] for g in ORDER])
for i, g in enumerate(ORDER):
    ax2.text(i, taus[i] + 0.02, f"r{G[g]['pancreas_rank']}", ha="center", fontsize=7)
ax2.set_ylim(0, 1.12); ax2.set_ylabel("tau tissue specificity")
ax2.set_xticks(x); ax2.set_xticklabels(ORDER, rotation=45, ha="right", fontsize=8)
ax2.set_title("Specificity index (label: pancreas rank / 54)", fontsize=9)
import matplotlib.patches as mp
ax2.legend(handles=[mp.Patch(color=COLS[k], label=k.replace("_", " ")) for k in COLS], fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig("paper/figs/fig_gtex_background.pdf")
print("wrote paper/figs/fig_gtex_background.pdf")
