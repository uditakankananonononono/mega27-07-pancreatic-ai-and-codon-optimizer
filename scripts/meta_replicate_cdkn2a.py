"""Multi-cohort replication + Mantel-Haenszel meta-analysis of CDKN2A co-mutation.

Extends the 2-cohort replication to every cBioPortal PDAC cohort with mutation
data (9 independent cohorts + TCGA discovery). Pooled OR via Mantel-Haenszel,
heterogeneity via Cochran's Q and I^2. Verdicts per covariate: REPLICATED
(pooled p<0.05, same direction in >=2/3 cohorts), FALSIFIED (pooled ns or
direction flips), or INCONCLUSIVE.
"""
import json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from scipy.stats import fisher_exact, chi2
from pancreatic_ai import features

STUDIES = ["paad_qcmg_uq_2016", "pdac_msk_2024", "paad_icgc", "paad_utsw_2015",
           "paad_cptac_2021", "paac_jhu_2014", "pancreas_msk_2024", "paad_msk_2025", "paad_tcga"]

def per_study(study):
    d = pathlib.Path("data") / study
    muts = json.loads((d / "mutations.json").read_text())
    samples = json.loads((d / "samples.json").read_text())
    per = {s: {m["hugo"] for m in muts if m["sampleId"] == s
               and m["mutationType"] in features.NONSYNONYMOUS} for s in samples}
    y = np.array([1 if "CDKN2A" in per[s] else 0 for s in samples])
    out = {"study": study, "n": len(samples), "cdkn2a_mut_n": int(y.sum())}
    for label, gs in [("KRAS", ["KRAS"]), ("TP53_pathway", features.PATHWAYS["TP53"]),
                      ("SWI_SNF", features.PATHWAYS["SWI_SNF"]),
                      ("RTK_RAS_nonKRAS", [g for g in features.PATHWAYS["RTK_RAS"] if g != "KRAS"])]:
        g = np.array([1 if per[s] & set(gs) else 0 for s in samples])
        a = int(((y==1)&(g==1)).sum()); b = int(((y==1)&(g==0)).sum())
        c = int(((y==0)&(g==1)).sum()); dd = int(((y==0)&(g==0)).sum())
        odds, p = fisher_exact([[a,b],[c,dd]])
        out[label] = {"table": [[a,b],[c,dd]],
                      "odds_ratio": round(float(odds),3) if np.isfinite(odds) else None,
                      "p": float(f"{p:.2e}")}
    return out

def mantel_haenszel(tables):
    # tables: list of [[a,b],[c,d]] with 0.5 correction where needed
    num = den = 0.0; var = 0.0; ors = []
    for t in tables:
        a,b,c,d = [float(x) for x in (t[0][0],t[0][1],t[1][0],t[1][1])]
        if min(a,b,c,d) == 0: a,b,c,d = a+0.5,b+0.5,c+0.5,d+0.5
        n = a+b+c+d
        num += a*d/n; den += b*c/n
        ors.append((a*d)/(b*c))
    or_mh = num/den
    # Cochran's Q on log-OR with inverse-variance weights
    logs=[]; w=[]
    for t in tables:
        a,b,c,d = [float(x) for x in (t[0][0],t[0][1],t[1][0],t[1][1])]
        if min(a,b,c,d) == 0: a,b,c,d = a+0.5,b+0.5,c+0.5,d+0.5
        lo = np.log((a*d)/(b*c)); v = 1/a+1/b+1/c+1/d
        logs.append(lo); w.append(1/v)
    logs=np.array(logs); w=np.array(w)
    mu = (w*logs).sum()/w.sum()
    Q = float((w*(logs-mu)**2).sum()); k=len(logs)
    I2 = max(0.0,(Q-(k-1))/Q) if Q>0 else 0.0
    # MH z-test
    p = None
    # variance of log MH-OR (Robins-Greenland)
    R=S=P1=P2=0.0
    for t in tables:
        a,b,c,d = [float(x) for x in (t[0][0],t[0][1],t[1][0],t[1][1])]
        if min(a,b,c,d) == 0: a,b,c,d = a+0.5,b+0.5,c+0.5,d+0.5
        n=a+b+c+d
        R += a*d/n; S += b*c/n
        P1 += (a+d)/n * (a*d/n); P2 += (b+c)/n * (b*c/n)
    var = P1/(2*R*R) + P2/(2*S*S)
    from scipy.stats import norm
    z = np.log(or_mh)/np.sqrt(var)
    p = float(2*(1-norm.cdf(abs(z))))
    return or_mh, p, Q, I2

res = [per_study(s) for s in STUDIES]
meta = {}
for label in ["KRAS","TP53_pathway","SWI_SNF","RTK_RAS_nonKRAS"]:
    tables=[r[label]["table"] for r in res if sum(map(sum, r[label]["table"]))>0]
    ors=[r[label]["odds_ratio"] for r in res if r[label]["odds_ratio"]]
    same_dir = sum(1 for o in ors if (o>1)==(np.median(ors)>1))
    or_mh,p,Q,I2 = mantel_haenszel(tables)
    verdict = ("REPLICATED" if p<0.05 and same_dir>=2*len(ors)/3 else
               "FALSIFIED" if p>=0.05 or same_dir<len(ors)/2 else "INCONCLUSIVE")
    meta[label]={"pooled_or_mh":round(float(or_mh),3),"pooled_p":float(f"{p:.2e}"),
                 "cochrans_Q":round(Q,2),"I2":round(I2,3),"k_cohorts":len(tables),
                 "same_direction":f"{same_dir}/{len(ors)}","verdict":verdict}
out={"discovery":"paad_tcga_pan_can_atlas_2018","per_study":res,"meta":meta}
pathlib.Path("results/cdkn2a_meta_replication.json").write_text(json.dumps(out,indent=1))
print(json.dumps(meta,indent=1))
