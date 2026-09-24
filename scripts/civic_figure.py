#!/usr/bin/env python3
"""Stacked disease-distribution bars for the CIViC evidence audit (accepted items)."""
import json, re
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PANCREAS = re.compile(r"pancrea", re.I)
BREAST_OVARIAN = re.compile(r"breast|ovarian", re.I)
HEMATOLOGIC = re.compile(r"myelo|leukemi|polycythemia|thrombocythemia|lymph", re.I)
GENES = [("KRAS", PANCREAS), ("TP53", PANCREAS), ("CDKN2A", PANCREAS), ("SMAD4", PANCREAS),
         ("BRCA1", BREAST_OVARIAN), ("JAK2", HEMATOLOGIC)]

fig, ax = plt.subplots(figsize=(6.4, 3.2))
names, tgt, other, nodis = [], [], [], []
for sym, pat in GENES:
    items = json.load(open(f"data/civic/{sym}.json"))["items"]
    acc = [e for e in items if e["status"] == "ACCEPTED"]
    k = sum(1 for e in acc if e["disease"] and pat.search(e["disease"]["name"] or ""))
    nd = sum(1 for e in acc if not (e["disease"] and e["disease"].get("name")))
    names.append(sym); tgt.append(k); nodis.append(nd); other.append(len(acc) - k - nd)
ax.bar(names, tgt, color="#b2182b", label="target-disease evidence (pancreatic / breast+ovarian / hematologic)")
ax.bar(names, other, bottom=tgt, color="#9ecae1", label="other-disease evidence")
ax.bar(names, nodis, bottom=[t + o for t, o in zip(tgt, other)], color="#d9d9d9", label="no disease recorded")
for i, (t, o, n) in enumerate(zip(tgt, other, nodis)):
    tot = t + o + n
    ax.text(i, tot + 3, f"{100*t/tot:.0f}%", ha="center", fontsize=8)
ax.set_ylabel("accepted evidence items")
ax.set_ylim(0, max(t + o + n for t, o, n in zip(tgt, other, nodis)) * 1.28)
ax.set_title("CIViC clinical-evidence concentration: panel vs tissue-specific controls", fontsize=9)
ax.legend(fontsize=6.5, loc="upper right")
fig.tight_layout()
fig.savefig("paper/figs/fig_civic.pdf")
print("wrote paper/figs/fig_civic.pdf")
