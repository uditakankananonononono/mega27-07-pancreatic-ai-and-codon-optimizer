"""Tests for the MyVariant allele-level audit (committed results only, hermetic)."""
import csv, json, os

A = "results/myvariant_panel_audit.json"
V = "results/myvariant_variants.csv"
C = "results/myvariant_allele_context.csv"
L = "results/myvariant_missense_alleles.csv"


def test_files_consistent():
    for p in (A, V, C, L):
        assert os.path.exists(p)
    a = json.load(open(A))
    assert sum(1 for _ in open(V)) - 1 == a["n_rows"]
    ids = {r["hgvs"] for r in csv.DictReader(open(V)) if r["hgvs"]}
    assert len(ids) == a["n_snv_ids"]
    assert a["myvariant_found_fraction_snv"] > 0.95


def test_reference_context_matches_every_allele():
    rows = list(csv.DictReader(open(C)))
    assert rows and all(r["ref_matches"] == "1" for r in rows)


def test_positive_control_passes():
    q1 = json.load(open(A))["q1_positive_control"]
    assert q1["passed"] is True
    for s in ("am", "revel", "cadd"):
        assert q1[s]["auroc_panel_vs_passenger"] > 0.8


def test_alphamissense_ceiling_negative_is_preserved():
    q2 = json.load(open(A))["q2_recurrence"]
    assert q2["panel"]["am"]["p_auroc"] > 0.05
    assert q2["am_ceiling"]["panel_singletons_frac_am_lp"] > 0.8
    for split in ("msk_impact", "wes_cohorts"):
        assert q2["cohort_split_replication"][split]["am"]["p"] > 0.05


def test_cadd_signal_not_claimed_as_replicated():
    rep = json.load(open(A))["q2_recurrence"]["cohort_split_replication"]
    assert rep["msk_impact"]["cadd"]["p"] < 0.05
    assert rep["wes_cohorts"]["cadd"]["p"] > 0.05   # replication failed; paper must say so
    tex = open("paper/myvariant_sec.tex").read()
    assert "does not replicate" in tex


def test_cdkn2a_exception_and_support():
    a = json.load(open(A))
    pg = a["q2_recurrence"]["per_gene_am"]
    assert pg["CDKN2A"]["frac_patient_weighted_am_lp"] < min(pg[g]["frac_patient_weighted_am_lp"] for g in ("KRAS", "TP53", "SMAD4"))
    q3 = a["q3_patient_support"]
    assert 0.9 < q3["frac_panel_positive_with_supported_allele"] <= 1.0
