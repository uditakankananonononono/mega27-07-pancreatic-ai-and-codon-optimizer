"""Tests for the CDKN2A p14ARF frame test + cancerhotspots check (committed results only)."""
import csv, json, os

A = "results/cdkn2a_arf_frame.json"
R = "results/cdkn2a_arf_alleles.csv"


def test_transcripts_rebuilt_full_length():
    tc = json.load(open(A))["transcript_checks"]
    assert tc["p16_len_aa"] == 156 and tc["arf_len_aa"] == 132
    assert tc["p16_ok"] and tc["arf_ok"] and tc["shared_coding_bp"] > 150


def test_allele_table_matches_json():
    a = json.load(open(A))
    rows = [r for r in csv.DictReader(open(R)) if r["p16_kind"] == "missense" and r["p16_matches_cbio"] == "1"]
    assert len(rows) == a["n_alleles"]


def test_frame_geometry_holds():
    f = json.load(open(A))["codon_position_control_shared_region"]["frac_arf_altered_by_codon_pos_shared_region"]
    assert f["1"] > 0.9 and f["2"] == 0.0


def test_arf_explanation_not_claimed():
    a = json.load(open(A))
    assert a["fisher_p"] < 0.01                      # raw association exists
    c = a["codon_position_control_shared_region"]
    assert c["p_mh"] > 0.05                           # vanishes after codon-position stratification
    assert c["low_am_pos1_vs_pos2_fisher_p"] < 0.001
    tex = open("paper/cdkn2a_arf_sec.tex").read()
    assert "not supported" in tex


def test_hotspot_ceiling_reproduced():
    h = json.load(open(A))["cancerhotspots"]
    assert h["p_am"] > 0.05 and h["p_revel"] < 0.01 and h["p_cadd"] < 0.01
