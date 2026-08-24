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


# ---------------------------------------------------------------- R10: G7

def test_dock_vhh_candidates_parallel(tmp_path, monkeypatch):
    """R10/G7: candidates are all docked, results carry docking fields,
    PDBQT conversion is idempotent (second call reuses the file)."""
    from drugagent.modules import vhh as vh
    ok = [{"idx": i, "plddt": 90 - i} for i in range(8)]
    for i in range(8):
        (tmp_path / f"vhh_{i}.pdb").write_text(
            "ATOM      1  N   UNK A   1    0.0  0.0  0.0  1.0  0.0\n"
            "END\n")
    calls = {"pdbqt": 0, "dock": 0}

    def fake_to_pdbqt(pdb, out):
        calls["pdbqt"] += 1
        out.write_text("REMARK x\n" + "ATOM      1  C   UNK A   1 "
                       "0.0 0.0 0.0 0.0 0.0\n" * 10 + "\n")
        return out

    def fake_dock_one(args):
        calls["dock"] += 1
        return {"ok": True, "score": -4.0 - (args[0] and 0)}

    monkeypatch.setattr(vh, "to_pdbqt", fake_to_pdbqt, raising=False)
    monkeypatch.setattr(vh, "dock_one", fake_dock_one)
    out = vh.dock_vhh_candidates(ok, tmp_path, "rec.pdbqt",
                           {"center": [0, 0, 0], "xsize": 12,
                            "ysize": 12, "zsize": 12}, n_jobs=1)
    assert len(out) == 8
    assert calls["dock"] == 8
    assert all(r.get("ok") for r in out)
    # idempotency: pdbqt files exist -> no reconversion on second run
    before = calls["pdbqt"]
    vh.dock_vhh_candidates(ok, tmp_path, "rec.pdbqt",
                           {"center": [0, 0, 0], "xsize": 12,
                            "ysize": 12, "zsize": 12}, n_jobs=1)
    assert calls["pdbqt"] == before
    # empty input
    assert vh.dock_vhh_candidates([], tmp_path, "r", {}) == []
