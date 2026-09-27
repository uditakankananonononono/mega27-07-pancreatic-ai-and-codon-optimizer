"""Driver-discovery scan (verdict #7): score ALL genes in the PAAD mutation
matrix, with the 4 canonical drivers (KRAS, TP53, CDKN2A, SMAD4) as positive
controls, and report novel candidates above background.

Scoring (transparent, reproducible):
  recurrence  = fraction of cohort samples with >=1 nonsynonymous mutation
  damage_frac = fraction of a gene's mutations that are Nonsense/Frame/Splice
  score       = recurrence * (1 + damage_frac)  (frequency weighted by damage)
Known limits stated in the paper: no gene-length normalization (long genes
accumulate passenger missense), no replication-timing/covariate model - this is
a recurrence screen, not MutSigCV; candidates require orthogonal validation.
"""
import json, collections, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
muts = json.load(open(ROOT / "data" / "paad" / "mutations.json"))
samples = {s["sampleId"] if isinstance(s, dict) else s for s in json.load(open(ROOT / "data" / "paad" / "samples.json"))}
n_samples = len(samples)

NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
          "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation", "Translation_Start_Site"}
DAMAGING = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation"}
CONTROLS = {"KRAS", "TP53", "CDKN2A", "SMAD4"}

per_gene_samples = collections.defaultdict(set)
per_gene_types = collections.defaultdict(list)
for m in muts:
    if m.get("mutationType") in NONSYN:
        per_gene_samples[m["hugo"]].add(m["sampleId"])
        per_gene_types[m["hugo"]].append(m["mutationType"])

rows = []
for g, ss in per_gene_samples.items():
    types = per_gene_types[g]
    rec = len(ss) / n_samples
    dmg = sum(t in DAMAGING for t in types) / len(types)
    rows.append({"gene": g, "n_samples": len(ss), "recurrence": round(rec, 4),
                 "damage_frac": round(dmg, 3), "n_mut": len(types),
                 "score": round(rec * (1 + dmg), 4),
                 "known_driver": g in CONTROLS})
rows.sort(key=lambda r: -r["score"])

ctrl_ranks = {r["gene"]: i + 1 for i, r in enumerate(rows) if r["gene"] in CONTROLS}
novel = [r for r in rows if not r["known_driver"] and r["n_samples"] >= 5][:20]
out = {"n_samples": n_samples, "n_genes_scored": len(rows), "n_mutations": len(muts),
       "positive_control_ranks": ctrl_ranks,
       "positive_control_check": all(v <= 10 for v in ctrl_ranks.values()),
       "top10_overall": rows[:10],
       "novel_candidates_min5samples": novel}
json.dump(out, open(ROOT / "results" / "driver_discovery.json", "w"), indent=1)
print("controls:", ctrl_ranks)
print("top10:", [(r["gene"], r["n_samples"], r["damage_frac"]) for r in rows[:10]])
print("novel top12:", [(r["gene"], r["n_samples"], r["damage_frac"]) for r in novel[:12]])
