"""Pan-cancer CDKN2A landscape across fetched cBioPortal cohorts.

Per study: CDKN2A/KRAS/TP53/SMAD4 alteration frequencies with Wilson CIs;
KRAS-CDKN2A co-alteration odds ratio (Fisher exact) where both vary.
Meta summary: distribution of CDKN2A prevalence across cancer types;
PDAC cohorts flagged for comparison with the item's PDAC-core finding.
Output: results/pancan_cdkn2a.json + paper/pancan_table.tex
"""
import json, pathlib
import numpy as np
from scipy.stats import fisher_exact

ROOT = pathlib.Path("data/pancan")
GENES = ["CDKN2A", "KRAS", "TP53", "SMAD4"]

def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    m = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (round(float((c - m) / d), 4), round(float((c + m) / d), 4))

rows = {}
for d in sorted(ROOT.iterdir()):
    if not d.is_dir():
        continue
    muts = json.loads((d / "mutations.json").read_text())
    samples = json.loads((d / "samples.json").read_text())
    n = len(samples)
    if n < 25:
        continue
    alt = {g: {m["sampleId"] for m in muts if m["hugo"] == g} for g in GENES}
    row = {"n_samples": n}
    for g in GENES:
        k = len(alt[g])
        lo, hi = wilson(k, n)
        row[g] = {"k": k, "freq": round(k / n, 4), "ci95": [lo, hi]}
    a, b = alt["KRAS"], alt["CDKN2A"]
    if 0 < len(a) < n and 0 < len(b) < n:
        tab = [[len(a & b), len(a - b)], [len(b - a), n - len(a | b)]]
        orr, p = fisher_exact(tab)
        row["kras_cdkn2a_coalt"] = {"odds_ratio": round(float(orr), 3), "fisher_p": float(f"{p:.2e}")}
    rows[d.name] = row

freqs = np.array([r["CDKN2A"]["freq"] for r in rows.values()])
pdac = [k for k in rows if "paad" in k or "panc" in k or "paac" in k]
out = {
    "n_studies": len(rows),
    "cdkn2a_freq_median": round(float(np.median(freqs)), 4),
    "cdkn2a_freq_iqr": [round(float(np.percentile(freqs, 25)), 4), round(float(np.percentile(freqs, 75)), 4)],
    "cdkn2a_freq_max": {"study": max(rows, key=lambda k: rows[k]["CDKN2A"]["freq"]),
                        "freq": round(float(freqs.max()), 4)},
    "pdac_cohorts": {k: rows[k]["CDKN2A"] for k in pdac},
    "studies": rows,
}
json.dump(out, open("results/pancan_cdkn2a.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "studies"}, indent=1)[:800])
# paper table: top 12 by CDKN2A frequency + PDAC cohorts
top = sorted(rows.items(), key=lambda kv: -kv[1]["CDKN2A"]["freq"])[:12]
with open("paper/pancan_table.tex", "w") as fh:
    fh.write("\\begin{tabular}{lrrrr}\n\\hline\nStudy & $n$ & CDKN2A freq & 95\\% CI & KRAS-CDKN2A OR \\\\\n\\hline\n")
    for k, r in top:
        co = r.get("kras_cdkn2a_coalt", {})
        fh.write(f"{k.replace('_',' ')} & {r['n_samples']} & {r['CDKN2A']['freq']:.3f} & "
                 f"[{r['CDKN2A']['ci95'][0]:.3f},{r['CDKN2A']['ci95'][1]:.3f}] & {co.get('odds_ratio','-')} \\\\\n")
    fh.write("\\hline\n\\end{tabular}\n")
print("studies:", len(rows), "table rows:", len(top))
