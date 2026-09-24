"""Feature construction for TCGA-PAAD driver-mutation prediction.

Research angle: predict driver mutation status (KRAS, TP53, CDKN2A, SMAD4)
from the *rest* of the tumor's mutational landscape + clinical covariates -
the target gene itself is excluded to prevent trivial leakage. Mirrors the
precision-medicine aim of PANCREAS.AI (mutation-aware therapy selection)
using public TCGA-PAAD data instead of a private image set.
"""
import json
import pathlib
from collections import Counter, defaultdict

import numpy as np

DATA = pathlib.Path(__file__).resolve().parents[2] / "data" / "paad"

NONSYNONYMOUS = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Splice_Region",
    "Translation_Start_Site", "Nonstop_Mutation",
}
TRUNCATING = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
              "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation"}
MISSENSE = {"Missense_Mutation"}
DRIVERS = ["KRAS", "TP53", "CDKN2A", "SMAD4"]
N_GENES = 300
TYPE_INDEX = {"missense": 0, "truncating": 1, "other": 2}


def load_raw(data_dir=None):
    d = data_dir or DATA
    muts = json.loads((d / "mutations.json").read_text())
    samples = json.loads((d / "samples.json").read_text())
    patient = json.loads((d / "clinical_patient.json").read_text())
    return muts, samples, patient


def sample_patient_map(data_dir=None):
    """Map sampleId -> patientId (TCGA sample barcode prefix)."""
    muts, samples, _ = load_raw(data_dir)
    return {s: "-".join(s.split("-")[:3]) for s in samples}


def build_targets(muts, samples):
    """Binary nonsynonymous-mutation status per driver gene."""
    status = {s: {g: 0 for g in DRIVERS} for s in samples}
    for m in muts:
        if m["hugo"] in DRIVERS and m["mutationType"] in NONSYNONYMOUS:
            status[m["sampleId"]][m["hugo"]] = 1
    return status


def clinical_features(data_dir=None):
    """patientId -> (age_years, is_male). Missing -> np.nan."""
    _, _, patient = load_raw(data_dir)
    per_patient = defaultdict(dict)
    for rec in patient:
        per_patient[rec["patientId"]][rec["clinicalAttributeId"]] = rec["value"]
    out = {}
    for pid, attrs in per_patient.items():
        age = np.nan
        for key in ("AGE", "AGE_AT_INDEX", "AGE_AT_DIAGNOSIS"):
            if key in attrs:
                try:
                    age = float(attrs[key]); break
                except ValueError:
                    pass
        sex = attrs.get("SEX", "")
        out[pid] = (age, 1.0 if sex.lower().startswith("m") else 0.0)
    return out


def build_matrix(muts, samples, n_genes=N_GENES, exclude=None):
    """samples x (n_genes x 3) mutation-type count matrix, plus gene list."""
    exclude = set(exclude or [])
    counts = Counter(m["hugo"] for m in muts
                     if m["mutationType"] in NONSYNONYMOUS and m["hugo"] not in exclude)
    top_genes = [g for g, _ in counts.most_common(n_genes)]
    gidx = {g: i for i, g in enumerate(top_genes)}
    X = np.zeros((len(samples), n_genes, 3), dtype=np.float32)
    sidx = {s: i for i, s in enumerate(samples)}
    for m in muts:
        g, t = m["hugo"], m["mutationType"]
        if g not in gidx or t not in NONSYNONYMOUS:
            continue
        ti = TYPE_INDEX["missense" if t in MISSENSE else "truncating" if t in TRUNCATING else "other"]
        X[sidx[m["sampleId"]], gidx[g], ti] += 1.0
    X = np.log1p(X)
    return X, top_genes


def build_dataset(target_gene, data_dir=None):
    """Full dataset for one driver: X (n, 300, 3), clinical (n, 2), y (n,)."""
    muts, samples, _ = load_raw(data_dir)
    X, genes = build_matrix(muts, samples, exclude=[target_gene])
    status = build_targets(muts, samples)
    clin = clinical_features(data_dir)
    sp = sample_patient_map(data_dir)
    C = np.array([clin.get(sp[s], (np.nan, 0.0)) for s in samples], dtype=np.float32)
    age = C[:, 0]
    C[:, 0] = np.where(np.isnan(age), np.nanmean(age), age) / 100.0  # scale
    C = np.hstack([C, extra_features(muts, samples, exclude=[target_gene])])
    y = np.array([status[s][target_gene] for s in samples], dtype=np.int64)
    return X, C, y, genes


# Curated pathway membership per TCGA PanCancer Atlas PAAD marker paper
# (Raphael et al. 2017 / Cancer Cell 2017 PAAD integrated analysis).
PATHWAYS = {
    "RTK_RAS": ["KRAS", "BRAF", "MAP2K4", "MAPK1", "NF1", "RASA1", "EGFR", "ERBB2"],
    "TP53": ["TP53", "MDM2", "MDM4", "ATM", "CHEK2"],
    "CELL_CYCLE": ["CDKN2A", "CDKN2B", "CCND1", "CDK4", "CDK6", "RB1"],
    "TGF_BETA": ["SMAD4", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "ACVR1B"],
    "SWI_SNF": ["ARID1A", "ARID1B", "ARID2", "SMARCA4", "SMARCB1", "PBRM1"],
    "DNA_REPAIR": ["BRCA1", "BRCA2", "PALB2", "MLH1", "MSH2", "MSH6", "POLD1", "POLE"],
}


def extra_features(muts, samples, exclude=None):
    """Per-sample global + pathway features (target gene excluded):
    [log1p(burden), frac_truncating, log1p(n_genes), pathway mutated (6x)]."""
    exclude = set(exclude or [])
    per_sample = defaultdict(list)
    for m in muts:
        if m["mutationType"] in NONSYNONYMOUS and m["hugo"] not in exclude:
            per_sample[m["sampleId"]].append(m)
    out = np.zeros((len(samples), 9), dtype=np.float32)
    for i, s in enumerate(samples):
        ms = per_sample.get(s, [])
        burden = len(ms)
        trunc = sum(1 for m in ms if m["mutationType"] in TRUNCATING)
        genes = {m["hugo"] for m in ms}
        row = [np.log1p(burden), trunc / burden if burden else 0.0, np.log1p(len(genes))]
        for pname, pgenes in PATHWAYS.items():
            row.append(1.0 if genes & set(pgenes) else 0.0)
        out[i] = row
    return out
