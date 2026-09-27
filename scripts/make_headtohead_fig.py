#!/usr/bin/env python3
"""Regenerate the head-to-head figure from archived summary, INCLUDING CodonOpt."""
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
S = json.load(open(ROOT / "results/icor_headtohead.json"))["summary"]

LABELS = {"original": "WT", "icor": "ICOR", "gensmart": "GenSmart", "HFC": "HFC",
          "BFC": "BFC", "URC": "URC", "ERC": "ERC", "ours_cnn": "CodonOpt"}

fig, ax = plt.subplots(figsize=(6.2, 4.6))
ax.axvspan(0.78, 0.95, color="#ffe9c9", alpha=0.6, zorder=0)
ax.text(0.795, -1.15, "CAI in ICOR band", color="#cc7a00", fontsize=9)
for k, v in S.items():
    if k == "ours_cnn":
        ax.scatter(v["cai_mean"], v["pred_expr_z_mean"], s=140, color="#d62728",
                   marker="*", zorder=5, label="CodonOpt (ours)")
        ax.annotate("CodonOpt", (v["cai_mean"], v["pred_expr_z_mean"]),
                    xytext=(v["cai_mean"] - 0.16, v["pred_expr_z_mean"] - 0.12), fontsize=11,
                    fontweight="bold", color="#d62728")
    else:
        ax.scatter(v["cai_mean"], v["pred_expr_z_mean"], s=55, color="gray", zorder=4)
        ax.annotate(LABELS.get(k, k), (v["cai_mean"], v["pred_expr_z_mean"]),
                    xytext=(v["cai_mean"] + 0.006, v["pred_expr_z_mean"] + 0.04), fontsize=9)
ax.set_xlabel("CAI (E. coli weights)")
ax.set_ylabel("Predicted expression (z, ExpressionCNN)")
ax.set_ylim(-1.35, 2.55)
fig.tight_layout()
fig.savefig(ROOT / "paper/figs/fig_icor_headtohead.pdf")
print("figure regenerated with CodonOpt at CAI", S["ours_cnn"]["cai_mean"], "z", S["ours_cnn"]["pred_expr_z_mean"])
