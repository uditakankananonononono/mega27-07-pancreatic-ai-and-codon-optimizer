"""Hermetic tests for the codonopt CLI (bundled data only, tiny iteration budget)."""
import json, pathlib, subprocess, sys, tempfile

REPO = pathlib.Path(__file__).resolve().parents[1]
ENV = {"PYTHONPATH": str(REPO / "src")}
import os
ENV = {**os.environ, **ENV}

PROT = ">t1\nMKTAYIAKQRQISFVKSHFSRQDILDLWQ\n>t2\nMAEGEITTFTALTEKFNLPPGNYKKPK\n"


def run_cli(*args):
    return subprocess.run([sys.executable, "-m", "codon_optimizer.cli", *args],
                          capture_output=True, text=True, env=ENV, timeout=600)


def test_optimize_end_to_end(tmp_path):
    f = tmp_path / "in.fasta"; f.write_text(PROT)
    out = tmp_path / "out"
    r = run_cli("optimize", str(f), "-o", str(out), "--iters", "20", "--batch", "12")
    assert r.returncode == 0, r.stderr
    m = json.loads((out / "codonopt_metrics.json").read_text())
    assert set(m["genes"]) == {"t1", "t2"}
    for g, row in m["genes"].items():
        dna = (out / f"{g}_mo.fasta").read_text().splitlines()[1]
        assert len(dna) == row["len_nt"] == len(PROT.split("\n")[1 if g == "t1" else 3]) * 3
        assert row["tai"] >= row["tai_floor"] - 1e-9
        assert 0.35 <= row["gc"] <= 0.70


def test_metrics_roundtrip(tmp_path):
    f = tmp_path / "in.fasta"; f.write_text(PROT)
    out = tmp_path / "out"
    assert run_cli("optimize", str(f), "-o", str(out), "--iters", "5", "--batch", "8").returncode == 0
    r = run_cli("metrics", str(out / "t1_mo.fasta"))
    assert r.returncode == 0, r.stderr
    rows = json.loads(r.stdout)
    assert "t1_mo" in rows and rows["t1_mo"]["tai"] > 0
