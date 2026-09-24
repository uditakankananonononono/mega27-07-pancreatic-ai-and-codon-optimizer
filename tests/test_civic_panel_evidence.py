"""Hermetic tests for the CIViC evidence audit; read only committed result files."""
import csv, json

R = json.load(open("results/civic_panel_evidence.json"))
ROWS = list(csv.DictReader(open("results/civic_panel_evidence_rows.csv")))

def test_controls_concentrate():
    assert R["genes"]["BRCA1"]["target_share"] > 0.5
    assert R["genes"]["JAK2"]["target_share"] > 0.8

def test_panel_is_diffuse():
    pooled = R["panel_pooled_accepted"]
    assert pooled["share"] < 0.10
    assert pooled["pancreatic_hits"] == 18 and pooled["n"] == 433
    assert R["fisher_panel_vs_brca1_breast_ovarian"]["p"] < 1e-10
    assert R["fisher_panel_vs_jak2_hematologic"]["p"] < 1e-10

def test_row_integrity():
    assert len(ROWS) == R["gene_eid_rows"] == 1153
    assert len({r["eid"] for r in ROWS}) == R["unique_eids"] == 1146
    acc = [r for r in ROWS if r["status"] == "ACCEPTED"]
    assert sum(int(r["pancreatic"]) for r in acc if r["group"] == "panel") == 18

def test_per_gene_counts():
    pg = R["per_gene_pancreatic_accepted"]
    assert pg["KRAS"] == {"k": 13, "n": 200}
    assert pg["TP53"] == {"k": 0, "n": 196}
    assert pg["CDKN2A"] == {"k": 2, "n": 30}
    assert pg["SMAD4"] == {"k": 3, "n": 7}
    assert R["genes"]["TP53"]["no_disease_accepted"] == 102
