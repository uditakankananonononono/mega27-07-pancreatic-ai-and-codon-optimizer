#!/usr/bin/env python3
"""Allele-level deleteriousness + ClinVar audit of PDAC panel somatic variants.

Input: results/myvariant_variants.csv (scripts/myvariant_fetch.py).
Questions:
 Q1 (positive control, must pass): somatic missense alleles in the 4 panel
    drivers score as more deleterious than missense in 8 long passenger genes
    (AlphaMissense, REVEL, CADD). If this fails, annotations are void.
 Q2: within panel genes, does predicted deleteriousness track recurrence
    (patients carrying the allele)? Negative control: same test in passengers.
 Q3: what fraction of panel-positive patients carry only alterations that are
    neither truncating, ClinVar P/LP, nor AlphaMissense likely-pathogenic?
    That is the fraction a mutation-based detection call would rest on
    unsupported alleles.
Writes results/myvariant_panel_audit.json."""
import json, os
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
AM_LP = 0.564   # AlphaMissense published likely-pathogenic threshold (Cheng et al. 2023)
TRUNC = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
         "Nonstop_Mutation", "Translation_Start_Site"}


def auc(pos, neg):
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    return float(u.statistic / (len(pos) * len(neg))), float(u.pvalue)


def main():
    d = pd.read_csv(os.path.join(RES, "myvariant_variants.csv"))
    for c in ("alphamissense", "revel", "cadd_phred"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    out = {"n_rows": int(len(d)), "n_patients_by_study": d.groupby("study").patientId.nunique().to_dict(),
           "n_snv_ids": int(d.hgvs.dropna().nunique()),
           "myvariant_found_fraction_snv": float(d[d.hgvs.notna()].mv_found.mean())}
    mis = d[d.mutationType == "Missense_Mutation"].copy()
    # unique allele table (one row per HGVS id), recurrence = distinct patients
    al = (mis.dropna(subset=["hgvs"]).groupby(["hgvs", "gene", "group"])
          .agg(n_pat=("patientId", lambda s: s.astype(str).nunique()), am=("alphamissense", "first"),
               revel=("revel", "first"), cadd=("cadd_phred", "first"), clinvar=("clinvar", "first"),
               protein=("proteinChange", "first")).reset_index())
    q1 = {}
    for sc in ("am", "revel", "cadd"):
        p = al[(al.group == "panel")][sc].dropna(); n = al[(al.group == "passenger")][sc].dropna()
        a, pv = auc(p.values, n.values)
        q1[sc] = {"panel_median": float(p.median()), "passenger_median": float(n.median()),
                  "n_panel_alleles": int(len(p)), "n_passenger_alleles": int(len(n)),
                  "auroc_panel_vs_passenger": a, "mwu_p": pv}
    q1["passed"] = all(q1[s]["auroc_panel_vs_passenger"] > 0.7 and q1[s]["mwu_p"] < 1e-6 for s in ("am", "revel", "cadd"))
    out["q1_positive_control"] = q1

    q2 = {}
    for grp, genes in (("panel", ["KRAS", "TP53", "CDKN2A", "SMAD4"]), ("passenger", None)):
        sub = al[al.group == grp] if genes is None else al[al.gene.isin(genes)]
        res = {}
        for sc in ("am", "revel", "cadd"):
            s = sub.dropna(subset=[sc])
            rho, p = stats.spearmanr(s.n_pat, s[sc])
            rec, sing = s[s.n_pat >= 3][sc], s[s.n_pat == 1][sc]
            a, pa = auc(rec.values, sing.values) if len(rec) >= 3 and len(sing) >= 3 else (None, None)
            res[sc] = {"spearman_rho_recurrence": float(rho), "p": float(p), "n_alleles": int(len(s)),
                       "n_recurrent_ge3": int(len(rec)), "n_singleton": int(len(sing)),
                       "auroc_recurrent_vs_singleton": a, "p_auroc": pa}
        q2[grp] = res
    per_gene = {}
    for g in ("KRAS", "TP53", "CDKN2A", "SMAD4"):
        s = al[(al.gene == g)].dropna(subset=["am"])
        rho, p = stats.spearmanr(s.n_pat, s.am) if len(s) > 3 else (np.nan, np.nan)
        rec, sing = s[s.n_pat >= 3].am, s[s.n_pat == 1].am
        a, pa = auc(rec.values, sing.values) if len(rec) >= 3 and len(sing) >= 3 else (None, None)
        per_gene[g] = {"n_alleles": int(len(s)), "n_patients_missense": int(s.n_pat.sum()),
                       "am_median": float(s.am.median()), "frac_alleles_am_lp": float((s.am >= AM_LP).mean()),
                       "frac_patient_weighted_am_lp": float((s.n_pat * (s.am >= AM_LP)).sum() / s.n_pat.sum()),
                       "spearman_rho": None if np.isnan(rho) else float(rho), "p": None if np.isnan(p) else float(p),
                       "n_recurrent_ge3": int(len(rec)), "auroc_recurrent_vs_singleton": a, "p_auroc": pa}
    q2["per_gene_am"] = per_gene
    top = al[al.group == "panel"].sort_values("n_pat", ascending=False).head(15)
    q2["top_panel_alleles"] = [{"gene": r.gene, "protein": r.protein, "n_pat": int(r.n_pat),
                                "am": None if pd.isna(r.am) else float(r.am), "clinvar": r.clinvar if isinstance(r.clinvar, str) else ""}
                               for r in top.itertuples()]
    # CpG-mutability control: recurrence can reflect mutability (CpG transitions)
    ctx = pd.read_csv(os.path.join(RES, "myvariant_allele_context.csv"))[["hgvs", "cpg_transition"]]
    alc = al.merge(ctx, on="hgvs", how="left")
    ctrl = {}
    for grp in ("panel", "passenger"):
        sub = alc[(alc.group == grp)].dropna(subset=["cpg_transition"])
        rho_c, p_c = stats.spearmanr(sub.n_pat, sub.cpg_transition)
        cres = {"n_alleles": int(len(sub)), "frac_cpg_transition": float(sub.cpg_transition.mean()),
                "spearman_recurrence_vs_cpg_ts": float(rho_c), "p": float(p_c)}
        for sc in ("am", "revel", "cadd"):
            s2 = sub.dropna(subset=[sc])
            r_score_cpg = stats.spearmanr(s2[sc], s2.cpg_transition)
            # partial Spearman: rank-residualize both on cpg_transition
            rn, rs, rc = (stats.rankdata(s2.n_pat), stats.rankdata(s2[sc]), s2.cpg_transition.values.astype(float))
            X = np.column_stack([np.ones(len(rc)), rc])
            res_n = rn - X @ np.linalg.lstsq(X, rn, rcond=None)[0]
            res_s = rs - X @ np.linalg.lstsq(X, rs, rcond=None)[0]
            pr, pp = stats.pearsonr(res_n, res_s)
            non = s2[s2.cpg_transition == 0]
            rn2, pn2 = stats.spearmanr(non.n_pat, non[sc])
            cres[sc] = {"spearman_score_vs_cpg_ts": float(r_score_cpg[0]), "p_score_vs_cpg": float(r_score_cpg[1]),
                        "partial_spearman_given_cpg": float(pr), "p_partial": float(pp),
                        "noncpg_only_rho": float(rn2), "noncpg_only_p": float(pn2), "n_noncpg": int(len(non))}
        ctrl[grp] = cres
    q2["cpg_control"] = ctrl
    # AlphaMissense ceiling: share of panel singletons already above the LP threshold
    ps = al[(al.group == "panel") & (al.n_pat == 1)].am.dropna()
    q2["am_ceiling"] = {"panel_singletons_frac_am_lp": float((ps >= AM_LP).mean()),
                        "panel_singletons_frac_am_ge_0_9": float((ps >= 0.9).mean()), "n": int(len(ps))}
    # cohort-split replication: MSK targeted (pdac_msk_2024) vs the six WES/WGS cohorts
    rep = {}
    pm = mis[(mis.group == "panel")].dropna(subset=["hgvs"])
    for name, part in (("msk_impact", pm[pm.study == "pdac_msk_2024"]), ("wes_cohorts", pm[pm.study != "pdac_msk_2024"])):
        a2 = (part.groupby("hgvs").agg(n_pat=("patientId", lambda s: s.astype(str).nunique()),
                                       am=("alphamissense", "first"), revel=("revel", "first"),
                                       cadd=("cadd_phred", "first")).reset_index())
        thr = 3 if name == "msk_impact" else 2
        r = {"n_alleles": int(len(a2)), "recurrent_threshold": thr}
        for sc in ("am", "revel", "cadd"):
            s2 = a2.dropna(subset=[sc])
            rec, sing = s2[s2.n_pat >= thr][sc], s2[s2.n_pat == 1][sc]
            a, pa = auc(rec.values, sing.values)
            r[sc] = {"auroc_recurrent_vs_singleton": a, "p": pa, "n_rec": int(len(rec)), "n_sing": int(len(sing))}
        rep[name] = r
    q2["cohort_split_replication"] = rep
    out["q2_recurrence"] = q2

    # Q3 patient-level support of panel calls
    pan = d[d.group == "panel"].copy()
    pan["trunc"] = pan.mutationType.isin(TRUNC)
    pan["clinvar_plp"] = pan.clinvar == "P/LP"
    pan["am_lp"] = (pan.mutationType == "Missense_Mutation") & (pan.alphamissense >= AM_LP)
    pan["supported"] = pan.trunc | pan.clinvar_plp | pan.am_lp
    pan["pid"] = pan.study + "|" + pan.patientId.astype(str)
    pt = pan.groupby("pid").agg(any_supported=("supported", "any"), any_clinvar=("clinvar_plp", "any"),
                                n_alt=("supported", "size"), study=("study", "first"))
    n_all = int(d.assign(pid=d.study + "|" + d.patientId.astype(str)).pid.nunique())
    uns = pan[~pan.supported]
    q3 = {"n_panel_positive_patients": int(len(pt)),
          "frac_panel_positive_with_supported_allele": float(pt.any_supported.mean()),
          "n_panel_positive_unsupported_only": int((~pt.any_supported).sum()),
          "frac_panel_positive_with_clinvar_plp": float(pt.any_clinvar.mean()),
          "unsupported_alteration_types": uns.mutationType.value_counts().to_dict(),
          "unsupported_by_gene": uns.gene.value_counts().to_dict(),
          "n_unsupported_alterations": int(len(uns)), "n_panel_alterations": int(len(pan)),
          "by_study_frac_supported": pt.groupby("study").any_supported.mean().to_dict(),
          "n_patients_any_fetched_gene_mutation": n_all,
          "note": "denominator = patients with >=1 panel alteration; panel-negative patients are not in the mutation table"}
    # clinvar annotation of the unsupported missense
    q3["unsupported_missense_clinvar"] = uns[uns.mutationType == "Missense_Mutation"].clinvar.fillna("none").value_counts().to_dict()
    out["q3_patient_support"] = q3
    json.dump(out, open(os.path.join(RES, "myvariant_panel_audit.json"), "w"), indent=1, default=str)
    al.to_csv(os.path.join(RES, "myvariant_missense_alleles.csv"), index=False)
    print(json.dumps(q2["cohort_split_replication"], indent=1)); print(json.dumps({"cpg": q2["cpg_control"], "ceiling": q2["am_ceiling"]}, indent=1))
    print(json.dumps({"q1": q1, "q2panel": q2["panel"], "q2pass": q2["passenger"], "per_gene": per_gene,
                      "q3": {k: v for k, v in q3.items()}}, indent=1, default=str))


if __name__ == "__main__":
    main()
