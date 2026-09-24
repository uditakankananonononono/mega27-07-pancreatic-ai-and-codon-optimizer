#!/usr/bin/env python3
"""Figure for the CDKN2A ARF-frame + cancerhotspots audit (committed results only)."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

a = json.load(open("results/cdkn2a_arf_frame.json"))
c = a["codon_position_control_shared_region"]
fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.2))
ax = axs[0]
st = c["strata_tables"]
labs, hi, lo = [], [], []
for k in ("1", "2"):
    t = st[k]   # rows high/low AM, cols ARF unaltered/altered
    labs.append(f"p16 codon pos {k}")
    hi.append(t[0][1] / max(1, sum(t[0]))); lo.append(t[1][1] / max(1, sum(t[1])))
x = range(len(labs))
ax.bar([i - 0.18 for i in x], hi, 0.36, color="#c0392b", label="AM >= 0.564")
ax.bar([i + 0.18 for i in x], lo, 0.36, color="#8e44ad", label="AM < 0.564")
for i, k in enumerate(("1", "2")):
    t = st[k]
    ax.text(i - 0.18, hi[i] + 0.02, f"n={sum(t[0])}", ha="center", fontsize=7)
    ax.text(i + 0.18, lo[i] + 0.02, f"n={sum(t[1])}", ha="center", fontsize=7)
ax.set_xticks(list(x)); ax.set_xticklabels(labs, fontsize=8); ax.set_ylim(0, 1.15)
ax.set_ylabel("frac. alleles altering p14ARF"); ax.legend(fontsize=7, frameon=False, loc="center right")
ax.set_title("A. ARF effect is set by codon position", fontsize=9)
ax = axs[1]
h = a["cancerhotspots"]
v = [h["auroc_hotspot_vs_not_am"], h["auroc_hotspot_vs_not_revel"], h["auroc_hotspot_vs_not_cadd"]]
ax.bar(["AlphaMissense", "REVEL", "CADD"], v, color=["#8e44ad", "#2980b9", "#27ae60"])
ax.axhline(0.5, c="k", lw=0.8); ax.set_ylim(0.4, 0.7)
ax.set_ylabel("AUROC hotspot vs non-hotspot residue")
ax.set_title("B. cancerhotspots.org residues (panel)", fontsize=9)
plt.tight_layout(); plt.savefig("paper/figs/fig_cdkn2a_arf.pdf"); print("ok")
