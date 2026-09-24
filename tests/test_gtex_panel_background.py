"""Hermetic checks on the committed GTEx panel-background audit.
Reads results/gtex_panel_background.json only; no network."""
import json, os

R = os.path.join(os.path.dirname(__file__), "..", "results", "gtex_panel_background.json")
D = json.load(open(R))
G = D["genes"]
PANEL = ["KRAS", "TP53", "CDKN2A", "SMAD4"]


def test_coverage_complete():
    assert D["n_genes"] == 8
    assert D["n_tissue_gene_records"] == 432
    for g, v in G.items():
        assert v["n_tissues"] == 54 and v["n_samples"] == 17382
        assert 0.0 <= v["tau"] <= 1.0


def test_positive_controls_validate_metric():
    # PRSS1/INS are known pancreas-specific: the tau/rank machinery must recover that
    for g in ("PRSS1", "INS"):
        assert G[g]["top_tissue"] == "Pancreas"
        assert G[g]["pancreas_rank"] == 1
        assert G[g]["tau"] >= 0.95
        blood = max(G[g]["whole_blood_median_tpm"], 1e-3)
        assert G[g]["pancreas_median_tpm"] / blood >= 100.0


def test_panel_has_no_pancreatic_expression_specificity():
    # core negative: no panel gene is enriched in normal pancreas
    for g in PANEL:
        assert G[g]["pancreas_rank"] > 30
        assert G[g]["top_tissue"] != "Pancreas"


def test_leukocyte_background_exists():
    # blood (the liquid-biopsy compartment) expresses the panel genes
    for g in ("KRAS", "TP53", "SMAD4"):
        assert G[g]["whole_blood_median_tpm"] > 1.0
    # KRAS is higher in blood than in the target organ
    assert G["KRAS"]["whole_blood_median_tpm"] > G["KRAS"]["pancreas_median_tpm"]


def test_cdkn2a_silent_in_normal_pancreas_and_blood():
    assert G["CDKN2A"]["pancreas_median_tpm"] < 1.0
    assert G["CDKN2A"]["whole_blood_median_tpm"] < 1.0


def test_housekeeping_blood_background_high():
    for g in ("GAPDH", "ACTB"):
        assert G[g]["whole_blood_median_tpm"] > 1000.0
