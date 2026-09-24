import sys, os, json
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from gnomad_specificity import hgvsp_to_short, best_af

RES = os.path.join(os.path.dirname(__file__), "..", "results")


def test_hgvsp_parse():
    assert hgvsp_to_short("p.Gly12Asp") == "G12D"
    assert hgvsp_to_short("p.Arg248Gln") == "R248Q"
    assert hgvsp_to_short("p.Arg80Ter") == "R80*" or hgvsp_to_short("p.Arg80Ter") == ""
    assert hgvsp_to_short("") == "" and hgvsp_to_short(None) == ""
    assert hgvsp_to_short("p.Gly12_Gly13dup") == ""


def test_best_af():
    assert best_af({"af_exome": 1e-6, "af_genome": ""}) == 1e-6
    assert best_af({"af_exome": "", "af_genome": ""}) is None
    assert best_af({"af_exome": 1e-6, "af_genome": 2e-6}) == 2e-6


def test_committed_json_consistent_with_csv():
    out = json.load(open(os.path.join(RES, "gnomad_panel_specificity.json")))
    df = pd.read_csv(os.path.join(RES, "gnomad_panel_variants.csv"), keep_default_na=False)
    assert len(df) == out["n_gnomad_variants"] == 9410
    assert set(df["gene"].unique()) == {"KRAS", "TP53", "CDKN2A", "SMAD4"}
    # headline bounds from committed data
    assert out["n_flagged_af_gt_1e-4"] == 0
    assert out["max_af_over_hotspots"] < 3e-5
    assert out["n_ac0_or_absent"] == 7
    # KRAS G12D row present and specific
    g12d = [r for r in out["allele_table"] if r["gene"] == "KRAS" and r["allele"] == "G12D"][0]
    assert g12d["paad_samples"] == 49 and g12d["max_af"] < 2e-5
    # LoF background ordering preserved
    lb = out["lof_background"]
    assert lb["SMAD4"]["sum_af"] > lb["TP53"]["sum_af"] > 0


def test_allele_table_covers_recurrent_kras():
    out = json.load(open(os.path.join(RES, "gnomad_panel_specificity.json")))
    kras = {r["allele"] for r in out["allele_table"] if r["gene"] == "KRAS"}
    assert {"G12D", "G12V", "G12R", "Q61H"} <= kras
