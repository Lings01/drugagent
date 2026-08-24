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

    def fake_to_pdbqt(pdb, out, **kw):
        calls["pdbqt"] += 1
        calls.setdefault("flex_args", []).append(kw.get("flex"))
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


# --------------------------------------------------------------- R11: G9/G10

def test_fast_vhh_defaults_widened():
    """R11/G9: fast mode widens the hit surface (measured pLDDT of modeled
    VHHs concentrates at 30-35; the old 45 gate left 1 candidate)."""
    from drugagent.config import resolve_defaults
    fast = resolve_defaults({"fast": True})
    assert fast.vhh_plddt_min == 35.0
    assert fast.vhh_screen_n == 80
    full = resolve_defaults({})
    assert full.vhh_plddt_min == 50.0
    assert full.vhh_screen_n == 1000
    # options override still wins
    ovr = resolve_defaults({"fast": True, "vhh_plddt_min": 40.0})
    assert ovr.vhh_plddt_min == 40.0
    assert ovr.vhh_dock_flex is False


def test_dock_vhh_candidates_rigid_default(tmp_path, monkeypatch):
    """R11/G10: default docking is RIGID — make_rigid_pdbqt zeroes the
    torsions; flex=True leaves them alone."""
    from drugagent.modules import vhh as vh
    ok = [{"idx": 0, "plddt": 60.0}]
    (tmp_path / "vhh_0.pdb").write_text(
        "ATOM      1  N   UNK A   1    0.0  0.0  0.0  1.0  0.0\nEND\n")
    rigid_calls = {"n": 0}

    def fake_to_pdbqt(pdb, out, **kw):
        out.write_text("REMARK x\n" + "ATOM      1  C   UNK A   1 "
                       "0.0 0.0 0.0 0.0 0.0\n" * 10
                       + "ROOT\nENDROOT\nTORSDOF 231\n")
        return out

    def fake_rigid(pdbqt):
        rigid_calls["n"] += 1
        text = pdbqt.read_text().replace("TORSDOF 231", "TORSDOF 0")
        pdbqt.write_text(text)
        return pdbqt

    def fake_dock_one(args):
        return {"ok": True, "score": -5.0}

    monkeypatch.setattr(vh, "to_pdbqt", fake_to_pdbqt)
    monkeypatch.setattr(vh, "make_rigid_pdbqt", fake_rigid)
    monkeypatch.setattr(vh, "dock_one", fake_dock_one)
    out = vh.dock_vhh_candidates(ok, tmp_path, "rec.pdbqt",
                                 {"center": [0, 0, 0], "xsize": 12,
                                  "ysize": 12, "zsize": 12}, n_jobs=1)
    assert len(out) == 1 and out[0]["ok"]
    assert rigid_calls["n"] == 1                      # rigid by default
    assert vh._pdbqt_is_flex(tmp_path / "vhh_0.pdbqt") is False
    # explicit flex=True leaves the torsions active
    (tmp_path / "vhh_0.pdbqt").unlink()
    rigid_calls["n"] = 0
    vh.dock_vhh_candidates(ok, tmp_path, "rec.pdbqt",
                           {"center": [0, 0, 0], "xsize": 12,
                            "ysize": 12, "zsize": 12}, n_jobs=1, flex=True)
    assert rigid_calls["n"] == 0
    assert vh._pdbqt_is_flex(tmp_path / "vhh_0.pdbqt") is True


def test_dock_vhh_candidates_stale_pdbqt_reconverted(tmp_path, monkeypatch):
    """R11/G10: a stale PDBQT whose flex/rigid mode mismatches the request
    is reconverted instead of silently reused."""
    from drugagent.modules import vhh as vh
    ok = [{"idx": 0, "plddt": 60.0}]
    (tmp_path / "vhh_0.pdb").write_text(
        "ATOM      1  N   UNK A   1    0.0  0.0  0.0  1.0  0.0\nEND\n")
    # stale FLEXIBLE pdbqt (active torsions)
    (tmp_path / "vhh_0.pdbqt").write_text(
        "REMARK stale\n" + "ATOM      1  C   UNK A   1 "
        "0.0 0.0 0.0 0.0 0.0\n" * 5
        + "ROOT\nENDROOT\nTORSDOF 3\n")
    assert vh._pdbqt_is_flex(tmp_path / "vhh_0.pdbqt") is True
    recon = {"n": 0}

    def fake_to_pdbqt(pdb, out, **kw):
        recon["n"] += 1
        out.write_text("REMARK fresh\n" + "ATOM      1  C   UNK A   1 "
                       "0.0 0.0 0.0 0.0 0.0\n" * 5
                       + "ROOT\nENDROOT\nTORSDOF 7\n")
        return out

    def fake_rigid(pdbqt):
        pdbqt.write_text(pdbqt.read_text().replace("TORSDOF 7",
                                                   "TORSDOF 0"))
        return pdbqt

    def fake_dock_one(args):
        return {"ok": True, "score": -5.0}

    monkeypatch.setattr(vh, "to_pdbqt", fake_to_pdbqt)
    monkeypatch.setattr(vh, "make_rigid_pdbqt", fake_rigid)
    monkeypatch.setattr(vh, "dock_one", fake_dock_one)
    pocket = {"center": [0, 0, 0], "xsize": 12, "ysize": 12, "zsize": 12}
    # rigid request vs stale flexible file -> reconvert once
    vh.dock_vhh_candidates(ok, tmp_path, "rec.pdbqt", pocket, n_jobs=1)
    assert recon["n"] == 1
    assert vh._pdbqt_is_flex(tmp_path / "vhh_0.pdbqt") is False
    # matching mode now -> no reconversion
    vh.dock_vhh_candidates(ok, tmp_path, "rec.pdbqt", pocket, n_jobs=1)
    assert recon["n"] == 1
    # flex request vs rigid file -> reconvert again
    vh.dock_vhh_candidates(ok, tmp_path, "rec.pdbqt", pocket, n_jobs=1,
                           flex=True)
    assert recon["n"] == 2
    assert vh._pdbqt_is_flex(tmp_path / "vhh_0.pdbqt") is True


def test_make_rigid_pdbqt(tmp_path):
    """R11/G10: make_rigid_pdbqt zeroes the torsion count, idempotent."""
    from drugagent.modules.target_prep import make_rigid_pdbqt
    p = tmp_path / "x.pdbqt"
    p.write_text("REMARK x\nATOM      1  C   UNK A   1 0 0 0 0 0 0 C\n"
                 "ROOT\nATOM      1  C   UNK A   1 0 0 0 0 0 0 C\n"
                 "ENDROOT\nTORSDOF 231\n")
    make_rigid_pdbqt(p)
    assert "TORSDOF 0" in p.read_text()
    make_rigid_pdbqt(p)  # idempotent
    assert p.read_text().count("TORSDOF") == 1


def test_pdbqt_is_flex_missing_file(tmp_path):
    from drugagent.modules import vhh as vh
    assert vh._pdbqt_is_flex(tmp_path / "nope.pdbqt") is False
