"""Tests for the Open Targets panel audit (reads committed results only)."""
import csv, json, os

A = "results/opentargets_panel_audit.json"
C = "results/opentargets_assoc_rows.csv"

def test_result_files_exist_and_consistent():
    assert os.path.exists(A) and os.path.exists(C)
    a = json.load(open(A))
    n_csv = sum(1 for _ in open(C)) - 1
    n_json = sum(d["n_targets"] for d in a["diseases"].values())
    assert n_csv == n_json == 42329

def test_controls_calibrate_the_metric():
    a = json.load(open(A))
    d = a["per_gene"]
    # JAK2 must top MPN and sit far down the pancreatic table
    assert d["JAK2"]["diseases"]["EFO_0004251"]["rank"] == 1
    assert d["JAK2"]["diseases"]["MONDO_0005192"]["rank"] > 100
    # BRCA1 must rank top-5 for breast cancer
    assert d["BRCA1"]["diseases"]["MONDO_0007254"]["rank"] <= 5

def test_panel_concentrates_on_pancreatic():
    a = json.load(open(A))
    for g in ["KRAS", "TP53", "CDKN2A", "SMAD4"]:
        assert a["per_gene"][g]["diseases"]["MONDO_0005192"]["rank"] <= 10
    assert a["panel_all4_top_rank_probability"] < 1e-9

def test_genetic_backing_and_drugs():
    a = json.load(open(A))
    bg = a["pancreatic_background"]
    assert bg["lit_only_n"] + 1 > 0 and bg["lit_only_share"] > 0.15
    assert bg["fisher_top100_vs_rest_litonly_p"] < 0.001
    kr = a["per_gene"]["KRAS"]["drug_candidates"]
    assert "SOTORASIB" in kr["pancreatic_indication_drugs"]
    assert a["per_gene"]["CDKN2A"]["drug_candidates"]["count"] == 0
