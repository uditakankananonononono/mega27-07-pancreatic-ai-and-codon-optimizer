import sys, os
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from depmap_hpa_audit import dep_stats, dep_table, expr_table, THR


def test_dep_stats_detects_selective_dependency():
    rng = np.random.default_rng(1)
    vp = rng.normal(-1.5, 0.2, 40)
    vo = rng.normal(-0.2, 0.2, 400)
    s = dep_stats(vp, vo)
    assert s["median_pancreas"] < -1.2 and s["frac_dep_pancreas"] > 0.9
    assert s["mw_p_two_sided"] < 1e-10 and s["n_pancreas"] == 40


def test_dep_stats_tolerated_tsg_null():
    rng = np.random.default_rng(2)
    s = dep_stats(rng.normal(0.1, 0.1, 50), rng.normal(0.1, 0.1, 500))
    assert s["frac_dep_pancreas"] == 0.0 and s["mw_p_two_sided"] > 0.05


def test_dep_table_uses_lineage_flag():
    df = pd.DataFrame({"OncotreeLineage": ["Pancreas"] * 3 + ["Lung"] * 3,
                       "G": [-2, -2, -2, 0, 0, 0]})
    t = dep_table(df, ["G"])
    assert t["G"]["median_pancreas"] == -2 and t["G"]["frac_dep_other"] == 0.0


def test_expr_table_zero_fraction():
    hpa = pd.DataFrame({"Gene name": ["G"] * 6,
                        "Cell line": ["PANC-1", "BxPC-3", "HeLa", "A549", "K562", "U2OS"],
                        "nTPM": [0.0, 10.0, 5.0, 5.0, 5.0, 5.0],
                        "pancreatic": [True, True, False, False, False, False]})
    t = expr_table(hpa, ["G"])["G"]
    assert t["frac_zero_pancreatic"] == 0.5 and t["n_pancreatic_lines"] == 2
    assert t["median_ntpm_all_lines"] == 5.0


def test_threshold_strict():
    s = dep_stats(np.array([THR, THR - 0.01]), np.zeros(10))
    assert s["frac_dep_pancreas"] == 0.5
