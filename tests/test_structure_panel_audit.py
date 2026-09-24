"""Tests for the structural-context audit (committed results only)."""
import json, os
import pandas as pd

A = "results/structure_panel_audit.json"


def test_numbering_verified():
    a = json.load(open(A))
    for g, s in a["structures"].items():
        assert s["uniprot_identity"] > 0.95, g
    m = pd.read_csv("results/structure_panel_alleles.csv")
    assert len(m) == a["n_alleles_mapped"] and (m.refaa == m.aa).all()


def test_burial_control_direction():
    a = json.load(open(A))
    assert a["spearman_rsa_vs_am"] < 0 and a["p_rsa_am"] < 0.05


def test_interface_hypothesis_falsified():
    a = json.load(open(A))
    assert a["low_am_x_contact"]["or"] < 1 and a["low_am_x_contact"]["p"] < 0.05
    assert a["logit"]["contact"]["coef"] < 0 and a["logit"]["rsa"]["coef"] > 0
    assert "falsified" in open("paper/structure_sec.tex").read()


def test_residue_table_contacts():
    r = pd.read_csv("results/structure_panel_residues.csv")
    assert set(r.gene) == {"KRAS", "TP53", "CDKN2A", "SMAD4"}
    assert r.rsa.between(0, 1).all() and r.partner_contact.sum() > 50
