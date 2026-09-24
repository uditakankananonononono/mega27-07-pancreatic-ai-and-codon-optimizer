#!/usr/bin/env python3
"""Normal-tissue expression background of the PDAC detection panel (GTEx v8).

Question: does the 4-gene panel (KRAS/TP53/CDKN2A/SMAD4) carry any
expression-level tissue specificity that a detection assay could exploit, and
what is the normal background in pancreas (target organ) and whole blood
(liquid-biopsy sampling compartment)?

Metrics per gene (from data/gtex/<SYM>.json raw per-sample TPM vectors):
  median TPM per tissue (54 tissues), tau tissue-specificity index
  tau = sum_i(1 - x_i/x_max)/(n-1)  (0 = ubiquitous, 1 = single-tissue),
  pancreas rank among tissues, whole-blood median TPM, global median.

Positive controls: PRSS1/INS must be pancreas-top and high-tau;
housekeeping GAPDH/ACTB must be low-tau. Outputs:
  results/gtex_panel_background.json, results/gtex_panel_tissue_medians.csv"""
import json, os, statistics, csv

GENES = ["KRAS", "TP53", "CDKN2A", "SMAD4", "GAPDH", "ACTB", "PRSS1", "INS"]

genes = {}
rows = []
for sym in GENES:
    d = json.load(open(f"data/gtex/{sym}.json"))
    med = {t["tissueSiteDetailId"]: statistics.median(t["tpm"]) for t in d["tissues"]}
    nsamp = {t["tissueSiteDetailId"]: t["n_samples"] for t in d["tissues"]}
    xs = list(med.values())
    xmax = max(xs)
    tau = sum(1.0 - x / xmax for x in xs) / (len(xs) - 1) if xmax > 0 else None
    ranked = sorted(med.items(), key=lambda kv: -kv[1])
    prank = next(i + 1 for i, (t, _) in enumerate(ranked) if t == "Pancreas")
    genes[sym] = {
        "group": d["group"], "gencodeId": d["gencodeId"],
        "n_tissues": len(med), "n_samples": sum(nsamp.values()),
        "tau": round(tau, 4),
        "pancreas_median_tpm": round(med["Pancreas"], 4),
        "pancreas_rank": prank,
        "whole_blood_median_tpm": round(med.get("Whole_Blood", 0.0), 4),
        "cells_ebv_median_tpm": round(med.get("Cells_EBV-transformed_lymphocytes", 0.0), 4),
        "top_tissue": ranked[0][0], "top_tissue_median_tpm": round(ranked[0][1], 3),
        "global_median_tpm": round(statistics.median(xs), 4),
        "tissue_medians": {k: round(v, 4) for k, v in med.items()},
    }
    for t, v in med.items():
        rows.append({"gene": sym, "group": d["group"], "gencodeId": d["gencodeId"],
                     "tissueSiteDetailId": t, "median_tpm": round(v, 4),
                     "n_samples": nsamp[t]})

panel = [g for g in GENES if genes[g]["group"] == "panel"]
summary = {
    "source": "GTEx Portal API v2, dataset gtex_v8 (per-sample TPM, 17,382 samples x 54 tissues per gene)",
    "n_genes": len(GENES), "n_tissues_per_gene": 54,
    "n_tissue_gene_records": sum(g["n_tissues"] for g in genes.values()),
    "panel_max_tau": max(genes[g]["tau"] for g in panel),
    "panel_any_pancreas_top5": any(genes[g]["pancreas_rank"] <= 5 for g in panel),
    "controls_pancreas_top": {g: genes[g]["top_tissue"] for g in ("PRSS1", "INS")},
    "housekeeping_max_tau": max(genes[g]["tau"] for g in ("GAPDH", "ACTB")),
    "genes": genes,
}
json.dump(summary, open("results/gtex_panel_background.json", "w"), indent=1)
with open("results/gtex_panel_tissue_medians.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

for g in GENES:
    v = genes[g]
    print(f"{g:7s} tau={v['tau']:.3f} pancreas={v['pancreas_median_tpm']:>9.3f} (rank {v['pancreas_rank']:2d}) "
          f"blood={v['whole_blood_median_tpm']:>8.3f} top={v['top_tissue']} ({v['top_tissue_median_tpm']:.1f})")
print("records", summary["n_tissue_gene_records"])
