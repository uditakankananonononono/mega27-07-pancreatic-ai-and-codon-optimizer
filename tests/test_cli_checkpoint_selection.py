import json,sys,subprocess,os
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from codon_optimizer.cli import load_model,main
def test_no_silent_default():
 with pytest.raises(SystemExit):main(['metrics','unused.fasta'])
def test_explicit_provenance_match():
 m=load_model(ROOT/'results_expr_unique_locus_model.pt',ROOT/'results/cli_unique_locus_audit.json');assert m.selection_metadata['scope'].startswith('unique-locus');assert len(m.selection_metadata['checkpoint_sha256'])==64
def test_mismatched_ledger_refused():
 with pytest.raises(ValueError,match='does not match'):load_model(ROOT/'results_expr_unique_locus_model.pt',ROOT/'results/cli_checkpoint_provenance.json')
