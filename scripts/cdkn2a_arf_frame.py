#!/usr/bin/env python3
"""Dissect the CDKN2A AlphaMissense exception (myvariant_panel_audit.json).

For every CDKN2A somatic missense SNV (results/myvariant_missense_alleles.csv):
 (1) rebuild the p16INK4a (NM_000077) and p14ARF (NM_058195) coding sequences
     from UCSC hg19 refGene exon coordinates + hg19 sequence (UCSC REST API),
     verify both translate to full-length proteins (156 / 132 aa, Met start,
     single terminal stop), and compute each allele's consequence in BOTH
     reading frames;
 (2) look up every panel residue in cancerhotspots.org (Chang et al., single-
     residue hotspots, API /api/hotspots/single/byGene) to see which recurrent
     residues AlphaMissense scores below its likely-pathogenic threshold.
Tests whether low-AlphaMissense CDKN2A alleles preferentially alter p14ARF
(the untested explanation stated in the MyVariant section).
Raw API JSON -> data/myvariant/ (untracked). Writes results/cdkn2a_arf_frame.json
and results/cdkn2a_arf_alleles.csv."""
import csv, json, os, time, urllib.request
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data", "myvariant")
RES = os.path.join(HERE, "..", "results")
UCSC = "https://api.genome.ucsc.edu/getData"
HOT = "https://www.cancerhotspots.org/api/hotspots/single/byGene/{g}"
AM_LP = 0.564
CODON = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
COMP = str.maketrans("ACGT", "TGCA")


def get(url, dest):
    dest = os.path.join(RAW, dest)
    if not os.path.exists(dest):
        for a in range(4):
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=120) as r:
                    json.dump(json.load(r), open(dest, "w")); break
            except Exception:
                if a == 3: raise
                time.sleep(4)
    return json.load(open(dest))


def cds_map(tx, refgene, seqd):
    """Return list of genomic 1-based positions in coding order + the coding sequence."""
    r = next(x for x in refgene if x["name"] == tx)
    s0, seq = seqd["start"], seqd["dna"].upper()
    starts = [int(v) for v in r["exonStarts"].strip(",").split(",")]
    ends = [int(v) for v in r["exonEnds"].strip(",").split(",")]
    pos = []
    for s, e in zip(starts, ends):
        lo, hi = max(s, r["cdsStart"]), min(e, r["cdsEnd"])
        pos += list(range(lo + 1, hi + 1))          # 1-based
    assert r["strand"] == "-"
    pos = sorted(pos, reverse=True)
    cds = "".join(seq[p - 1 - s0] for p in pos).translate(COMP)
    return pos, cds


def translate(s):
    return "".join(CODON[s[i:i + 3]] for i in range(0, len(s) - 2, 3))


def consequence(pos, cds, gpos, alt_plus):
    if gpos not in pos:
        return "noncoding", None, None, None
    i = pos.index(gpos)
    alt = alt_plus.translate(COMP)
    c = i // 3
    ref_codon = cds[3 * c:3 * c + 3]
    mut = list(ref_codon); mut[i % 3] = alt; mut = "".join(mut)
    ra, aa = CODON[ref_codon], CODON[mut]
    kind = "synonymous" if ra == aa else ("nonsense" if aa == "*" else "missense")
    return kind, f"{ra}{c + 1}{aa}", ra, i % 3 + 1


def main():
    os.makedirs(RAW, exist_ok=True)
    lo, hi = 21967000, 21995000
    rg = get(f"{UCSC}/track?genome=hg19;track=refGene;chrom=chr9;start={lo};end={hi}", "ucsc_refgene_cdkn2a.json")["refGene"]
    sq = get(f"{UCSC}/sequence?genome=hg19;chrom=chr9;start={lo};end={hi}", "ucsc_seq_cdkn2a_locus.json")
    p16_pos, p16 = cds_map("NM_000077", rg, sq)
    arf_pos, arf = cds_map("NM_058195", rg, sq)
    p16_aa, arf_aa = translate(p16), translate(arf)
    checks = {"p16_len_aa": len(p16_aa) - 1, "arf_len_aa": len(arf_aa) - 1,
              "p16_ok": p16_aa[0] == "M" and p16_aa[-1] == "*" and "*" not in p16_aa[:-1],
              "arf_ok": arf_aa[0] == "M" and arf_aa[-1] == "*" and "*" not in arf_aa[:-1],
              "shared_coding_bp": len(set(p16_pos) & set(arf_pos))}
    al = pd.read_csv(os.path.join(RES, "myvariant_missense_alleles.csv"))
    al = al[al.gene == "CDKN2A"].copy()
    rows, mism = [], 0
    for r in al.itertuples():
        g = int(r.hgvs.split("g.")[1][:-3]); ref, alt = r.hgvs[-3], r.hgvs[-1]
        k16, c16, _, cpos16 = consequence(p16_pos, p16, g, alt)
        karf, carf, _, _ = consequence(arf_pos, arf, g, alt)
        cb = r.protein if isinstance(r.protein, str) else ""
        agree = int(c16 is not None and c16[1:] == cb[1:])   # position+alt must match cBioPortal p16 call
        mism += 1 - agree
        rows.append({"hgvs": r.hgvs, "protein_cbio": cb, "p16_calc": c16 or "", "p16_kind": k16,
                     "arf_change": carf or "", "arf_kind": karf, "n_pat": int(r.n_pat),
                     "am": r.am, "revel": r.revel, "cadd": r.cadd, "clinvar": r.clinvar if isinstance(r.clinvar, str) else "",
                     "p16_matches_cbio": agree, "p16_codon_pos": cpos16 or ""})
    df = pd.DataFrame(rows)
    df["arf_altered"] = df.arf_kind.isin(["missense", "nonsense"]).astype(int)
    df["low_am"] = (df.am < AM_LP).astype(int)
    df = df.dropna(subset=["am"])
    n_before = len(df)
    # keep only alleles that are missense in NM_000077 and agree with the cBioPortal p16 call
    df = df[(df.p16_kind == "missense") & (df.p16_matches_cbio == 1)].copy()
    excluded = n_before - len(df)
    t = pd.crosstab(df.low_am, df.arf_altered).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
    orr, p = stats.fisher_exact(t.values)
    # patient-weighted
    pw = {k: float((v.n_pat * v.arf_altered).sum() / v.n_pat.sum()) for k, v in df.groupby("low_am")}
    # within the ARF-overlapping region only (position confound: region differs in structure)
    inreg = df[df.arf_kind != "noncoding"]
    t2 = pd.crosstab(inreg.low_am, inreg.arf_altered).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
    or2, p2 = stats.fisher_exact(t2.values)
    # codon-position control: AM substitution severity and ARF-frame effect both depend
    # on which p16 codon position is hit; Mantel-Haenszel across codon positions
    from statsmodels.stats.contingency_tables import StratifiedTable
    tabs, strata = [], {}
    for cp, sub in inreg.groupby("p16_codon_pos"):
        tt = pd.crosstab(sub.low_am, sub.arf_altered).reindex(index=[0, 1], columns=[0, 1], fill_value=0).values
        strata[str(cp)] = tt.tolist()
        if (tt.sum(0) > 0).all() and (tt.sum(1) > 0).all():
            tabs.append(tt.astype(float) + 0.0)
    if tabs:
        st = StratifiedTable(tabs)
        mh = {"or_mh": float(st.oddsratio_pooled), "p_mh": float(st.test_null_odds(correction=True).pvalue),
              "n_informative_strata": len(tabs)}
    else:
        mh = {"or_mh": None, "p_mh": None, "n_informative_strata": 0}
    mh["strata_tables"] = strata
    # is low AlphaMissense itself a codon-position-1 phenomenon (whole gene)?
    c12 = df[df.p16_codon_pos.isin([1, 2])]
    t3 = pd.crosstab(c12.p16_codon_pos, c12.low_am).reindex(index=[1, 2], columns=[0, 1], fill_value=0)
    or3, p3 = stats.fisher_exact(t3.values)
    mh["low_am_by_codon_pos_table_rows_pos1_pos2_cols_high_low"] = t3.values.tolist()
    mh["low_am_pos1_vs_pos2_fisher_or"] = float(or3); mh["low_am_pos1_vs_pos2_fisher_p"] = float(p3)
    mh["frac_arf_altered_by_codon_pos_shared_region"] = {str(k): float(v.arf_altered.mean()) for k, v in inreg.groupby("p16_codon_pos")}
    # AM within shared region vs p16-only region
    am_shared, am_p16only = df[df.arf_kind != "noncoding"].am, df[df.arf_kind == "noncoding"].am
    mw = stats.mannwhitneyu(am_shared, am_p16only, alternative="two-sided")
    out = {"transcript_checks": checks, "n_alleles": int(len(df)), "n_p16_mismatch_vs_cbio": int(mism), "n_excluded_not_p16_missense": int(excluded),
           "codon_position_control_shared_region": mh,
           "n_in_arf_coding": int((df.arf_kind != "noncoding").sum()),
           "arf_kind_counts": df.arf_kind.value_counts().to_dict(),
           "table_low_am_x_arf_altered": t.values.tolist(), "fisher_or": float(orr), "fisher_p": float(p),
           "patient_weighted_frac_arf_altered": {"high_am": pw.get(0), "low_am": pw.get(1)},
           "shared_region_only": {"table": t2.values.tolist(), "fisher_or": float(or2), "fisher_p": float(p2)},
           "am_median_shared_region": float(am_shared.median()), "am_median_p16_only_region": float(am_p16only.median()),
           "am_shared_vs_p16only_mwu_p": float(mw.pvalue),
           "low_am_alleles": df[df.low_am == 1].sort_values("n_pat", ascending=False)[
               ["protein_cbio", "arf_change", "arf_kind", "n_pat", "am", "revel", "clinvar"]].head(12).to_dict("records")}
    # hotspots across all 4 panel genes
    alla = pd.read_csv(os.path.join(RES, "myvariant_missense_alleles.csv"))
    hot = {}
    for gname in ("KRAS", "TP53", "CDKN2A", "SMAD4"):
        hs = get(HOT.format(g=gname), f"cancerhotspots_{gname}.json")
        hot[gname] = {int(h["aminoAcidPosition"]["start"]) for h in hs if h.get("type") == "single residue"}
    pan = alla[alla.group == "panel"].dropna(subset=["am"]).copy()
    pan["pos"] = pan.protein.astype(str).str.extract(r"^[A-Z](\d+)")[0].astype(float)
    pan["hotspot"] = [int(p in hot[g]) if p == p else 0 for g, p in zip(pan.gene, pan.pos)]
    hs_low = pan[(pan.hotspot == 1) & (pan.am < AM_LP)]
    hstats = {"n_hotspot_residues": {g: len(v) for g, v in hot.items()},
              "n_panel_alleles": int(len(pan)), "n_hotspot_alleles": int(pan.hotspot.sum()),
              "hotspot_alleles_below_am_lp": int(len(hs_low)),
              "hotspot_patients_below_am_lp": int(hs_low.n_pat.sum()),
              "hotspot_patients_total": int(pan[pan.hotspot == 1].n_pat.sum()),
              "below_lp_hotspot_alleles": hs_low.sort_values("n_pat", ascending=False)[["gene", "protein", "n_pat", "am", "revel", "cadd", "clinvar"]].to_dict("records")}
    for sc in ("am", "revel", "cadd"):
        s = pan.dropna(subset=[sc])
        u = stats.mannwhitneyu(s[s.hotspot == 1][sc], s[s.hotspot == 0][sc], alternative="two-sided")
        hstats[f"auroc_hotspot_vs_not_{sc}"] = float(u.statistic / ((s.hotspot == 1).sum() * (s.hotspot == 0).sum()))
        hstats[f"p_{sc}"] = float(u.pvalue)
    out["cancerhotspots"] = hstats
    json.dump(out, open(os.path.join(RES, "cdkn2a_arf_frame.json"), "w"), indent=1, default=str)
    df.to_csv(os.path.join(RES, "cdkn2a_arf_alleles.csv"), index=False)
    print(json.dumps(out, indent=1, default=str)[:4500])


if __name__ == "__main__":
    main()
