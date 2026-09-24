#!/usr/bin/env python3
"""Trinucleotide / CpG context of every missense SNV in results/myvariant_variants.csv
from the UCSC Genome Browser REST API (hg19 reference, one span query per gene;
raw spans cached in data/myvariant/ucsc_<gene>.json, untracked).
Purpose: CpG transitions are the most mutable somatic class, so allele
recurrence can reflect mutability instead of selection; the audit uses this
table as a covariate. Writes results/myvariant_allele_context.csv."""
import csv, json, os, time, urllib.request
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data", "myvariant")
RES = os.path.join(HERE, "..", "results")
API = "https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom={c};start={s};end={e}"
COMP = str.maketrans("ACGT", "TGCA")


def span(gene, chrom, lo, hi):
    dest = os.path.join(RAW, f"ucsc_{gene}.json")
    if not os.path.exists(dest):
        url = API.format(c=chrom, s=lo - 2, e=hi + 1)
        for a in range(4):
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=180) as r:
                    json.dump(json.load(r), open(dest, "w")); break
            except Exception:
                if a == 3: raise
                time.sleep(5)
    return json.load(open(dest))


def main():
    d = pd.read_csv(os.path.join(RES, "myvariant_variants.csv"))
    m = d[(d.mutationType == "Missense_Mutation") & d.hgvs.notna()].drop_duplicates("hgvs")
    rows = []
    for gene, g in m.groupby("gene"):
        chrom = "chr" + str(g.chr.iloc[0]).replace("chr", "")
        lo, hi = int(g.start.min()), int(g.start.max())
        sp = span(gene, chrom, lo, hi)
        seq, s0 = sp["dna"].upper(), sp["start"]
        for r in g.itertuples():
            i = int(r.start) - 1 - s0          # 0-based index of the variant base
            tri = seq[i - 1:i + 2]
            ok = tri[1:2] == r.ref
            cpg = (r.ref == "C" and tri[2:3] == "G") or (r.ref == "G" and tri[0:1] == "C")
            ts = (r.ref, r.alt) in {("C", "T"), ("G", "A"), ("T", "C"), ("A", "G")}
            rows.append({"hgvs": r.hgvs, "gene": gene, "trinuc": tri, "ref_matches": int(ok),
                         "cpg": int(cpg), "transition": int(ts), "cpg_transition": int(cpg and ts)})
    with open(os.path.join(RES, "myvariant_allele_context.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(len(rows), "ref match", sum(r["ref_matches"] for r in rows), "cpg_ts", sum(r["cpg_transition"] for r in rows))


if __name__ == "__main__":
    main()
