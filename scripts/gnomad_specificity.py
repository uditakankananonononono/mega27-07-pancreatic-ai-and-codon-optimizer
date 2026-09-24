#!/usr/bin/env python3
"""Two-sided specificity audit of the PDAC mutation-detection panel.

Somatic side: TCGA-PAAD mutation calls (data/paad/mutations.json, cBioPortal,
committed provenance via scripts/fetch_paad.py) - which KRAS/TP53/CDKN2A/SMAD4
alleles actually occur in pancreatic cancer and how often.
Population side: gnomAD r4 (results/gnomad_panel_variants.csv) - are those
exact alleles present in ~1.5M healthy-population chromosomes?
A blood-based panel allele is specific iff its population frequency is ~0.
Also: germline LoF carrier background per panel gene (specificity floor for
truncating-variant calls). Writes results/gnomad_panel_specificity.json and
paper/figs/fig_gnomad_specificity.pdf."""
import json, os, csv, re
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
LOF = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"}
AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G",
       "His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S",
       "Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Ter":"*"}
MUT_OK = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
          "In_Frame_Del", "In_Frame_Ins", "Splice_Site"}


def hgvsp_to_short(hgvsp):
    """p.Gly12Asp -> G12D ; p.Arg248Gln -> R248Q ; '' or odd -> ''."""
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)", hgvsp or "")
    if not m:
        return ""
    a, pos, b = m.groups()
    if a not in AA3 or (b not in AA3 and b != "*"):
        return ""
    return f"{AA3[a]}{pos}{AA3.get(b, '*')}"


def somatic_alleles(mut_path):
    muts = json.load(open(mut_path))
    genes = {"KRAS", "TP53", "CDKN2A", "SMAD4"}
    samples = set()
    counts = {}
    for m in muts:
        if m["hugo"] in genes:
            samples.add(m["sampleId"])
        if m["hugo"] in genes and m.get("proteinChange") and m.get("mutationType") in MUT_OK:
            key = (m["hugo"], m["proteinChange"])
            counts.setdefault(key, set()).add(m["sampleId"])
    return {k: len(v) for k, v in counts.items()}, len(samples)


def best_af(row):
    vals = [x for x in (row["af_exome"], row["af_genome"]) if x != ""]
    return max((float(x) for x in vals), default=None)


def main(make_fig=True):
    gdf = pd.read_csv(os.path.join(HERE, "..", "results", "gnomad_panel_variants.csv"),
                      keep_default_na=False)
    gdf["allele"] = gdf["hgvsp"].map(hgvsp_to_short)
    som, n_samples = somatic_alleles(os.path.join(HERE, "..", "data", "paad", "mutations.json"))

    rows = []
    for (gene, allele), cnt in sorted(som.items(), key=lambda kv: -kv[1]):
        if cnt < 2:
            continue
        sub = gdf[(gdf["gene"] == gene) & (gdf["allele"] == allele)]
        if len(sub) == 0:
            rows.append({"gene": gene, "allele": allele, "paad_samples": cnt,
                         "gnomad_status": "absent_from_record", "max_af": None, "ac": None, "an": None})
            continue
        afs = [best_af(r) for _, r in sub.iterrows()]
        afs = [a for a in afs if a is not None]
        acs = [int(r["ac_exome"] or 0) + int(r["ac_genome"] or 0) for _, r in sub.iterrows()]
        ans = [int(r["an_exome"] or 0) + int(r["an_genome"] or 0) for _, r in sub.iterrows()]
        rows.append({"gene": gene, "allele": allele, "paad_samples": cnt,
                     "gnomad_status": "observed" if max(acs) > 0 else "ac0",
                     "max_af": max(afs) if afs else 0.0, "ac": max(acs), "an": max(ans)})

    lof_bg = {}
    for gene, sub in gdf.groupby("gene"):
        lof = sub[sub["consequence"].isin(LOF)]
        afs = [best_af(r) or 0.0 for _, r in lof.iterrows()]
        lof_bg[gene] = {"n_lof_pass": int(sum(1 for _, r in lof.iterrows()
                                             if r["af_exome"] != "" or r["af_genome"] != "")),
                        "sum_af": float(np.sum(afs))}
    per_gene = gdf.groupby("gene").size().to_dict()

    tested = [r for r in rows]
    ac0 = [r for r in tested if r["gnomad_status"] in ("ac0", "absent_from_record")]
    flagged = [r for r in tested if (r["max_af"] or 0) > 1e-4]
    out = {
        "source": "gnomAD r4 GraphQL API (gnomad.broadinstitute.org); TCGA-PAAD via cBioPortal (paad_tcga_pan_can_atlas_2018)",
        "n_paad_samples_with_panel_gene_mutation": n_samples,
        "n_gnomad_variants": int(len(gdf)), "variants_per_gene": {k: int(v) for k, v in per_gene.items()},
        "n_somatic_alleles_tested": len(tested),
        "n_ac0_or_absent": len(ac0),
        "frac_ac0_or_absent": len(ac0) / max(1, len(tested)),
        "n_flagged_af_gt_1e-4": len(flagged), "flagged": flagged,
        "max_af_over_hotspots": max((r["max_af"] or 0.0) for r in tested),
        "lof_background": lof_bg,
        "allele_table": rows,
    }
    json.dump(out, open(os.path.join(HERE, "..", "results", "gnomad_panel_specificity.json"), "w"), indent=1)

    if make_fig:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        kras = [r for r in rows if r["gene"] == "KRAS"]
        fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.9))
        x = np.arange(len(kras))
        ax[0].bar(x, [r["paad_samples"] for r in kras], color="#c0392b")
        ax2 = ax[0].twinx()
        floor = 1e-7
        obs = [(i, max(r["max_af"] or floor, floor)) for i, r in enumerate(kras)
               if r["gnomad_status"] == "observed"]
        ab = [i for i, r in enumerate(kras) if r["gnomad_status"] != "observed"]
        if obs:
            ax2.plot([i for i, _ in obs], [a for _, a in obs], "ko", ms=4)
        if ab:
            ax2.plot(ab, [floor] * len(ab), "k^", ms=4, mfc="none")
        ax2.set_yscale("log"); ax2.set_ylim(5e-8, 1e-3)
        ax[0].set_xticks(x); ax[0].set_xticklabels([r["allele"] for r in kras], rotation=90, fontsize=7)
        ax[0].set_ylabel("PAAD samples (bars)"); ax2.set_ylabel("gnomAD r4 AF (dots)")
        genes = ["KRAS", "TP53", "CDKN2A", "SMAD4"]
        ax[1].bar(genes, [lof_bg[g]["sum_af"] for g in genes], color="#7f8c8d")
        ax[1].set_ylabel("sum AF of pass LoF variants", fontsize=8)
        ax[1].set_yscale("log")
        ax[1].tick_params(axis="x", labelsize=7, rotation=20)
        fig.tight_layout()
        os.makedirs(os.path.join(HERE, "..", "paper", "figs"), exist_ok=True)
        fig.savefig(os.path.join(HERE, "..", "paper", "figs", "fig_gnomad_specificity.pdf"))
    print(json.dumps({k: out[k] for k in ("n_gnomad_variants", "n_somatic_alleles_tested",
          "n_ac0_or_absent", "frac_ac0_or_absent", "n_flagged_af_gt_1e-4",
          "max_af_over_hotspots", "lof_background", "n_paad_samples_with_panel_gene_mutation")}, indent=1))


if __name__ == "__main__":
    main()
