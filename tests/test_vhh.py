import re
from pathlib import Path

import pytest

from drugagent.modules import vhh as vh


def test_library_generation():
    seqs = vh.generate_vhh_library(500, seed=1)
    assert len(seqs) == 500
    assert len(set(seqs)) == 500
    for s in seqs:
        assert 90 < len(s) < 160
        assert set(s) <= set("ACDEFGHIKLMNPQRSTVWY")
        # CDR3 region (tail) should be enriched in aromatic/polar
        tail = s[-17:]
        assert len(tail) >= 7


def test_library_diversity():
    seqs = vh.generate_vhh_library(1000, seed=2)
    tails = {s[-15:] for s in seqs}
    assert len(tails) > 300  # diverse CDR3s


def test_fasta_roundtrip(tmp_path):
    seqs = vh.generate_vhh_library(50, seed=3)
    p = vh.save_library(seqs, tmp_path / "lib.fasta")
    loaded = vh.load_fasta(p)
    assert loaded == seqs


@pytest.mark.slow
def test_model_one_vhh():
    seqs = vh.generate_vhh_library(1, seed=4)
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        r = vh.model_vhh_one((0, seqs[0], d))
        assert r["ok"], r.get("error")
        assert r["plddt"] > 0
        assert (Path(d) / "vhh_0.pdb").exists()
