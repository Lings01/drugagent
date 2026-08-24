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


def test_screen_vhh_end_to_end_fake(tmp_path, monkeypatch):
    """R11/G9: screen_vhh runs the full fake path (filter -> log -> dock),
    catching name bugs in the log line, and honours the config gate."""
    from drugagent.modules import vhh as vh
    state = {
        "options": {"fast": True},
        "target_prep": {
            "receptor_pdbqt": str(tmp_path / "rec.pdbqt"),
            "pocket": {"center": [0, 0, 0], "xsize": 12,
                       "ysize": 12, "zsize": 12},
        },
    }
    (tmp_path / "rec.pdbqt").write_text("REMARK rec\n")
    monkeypatch.setattr(vh, "load_fasta",
                        lambda p: ["ACGT" * 30] * 5)
    monkeypatch.setattr(vh, "model_vhh_one",
                        lambda a: {"ok": True, "plddt": 40.0,
                                   "idx": a[0], "model": "x.pdb"})
    seen = {}

    def fake_dock(ok, model_dir, rec, pocket, **kw):
        seen.update(kw)
        return [dict(r, ok=True, score=-6.0, rmsd_lb=0.0, rmsd_ub=0.0,
                     top_pose_pdbqt=None) for r in ok]

    monkeypatch.setattr(vh, "dock_vhh_candidates", fake_dock)
    out = vh.screen_vhh(state, tmp_path, n=5, n_jobs=1)
    assert out["n_docked"] == 5
    assert out["results"][0]["plddt"] == 40.0
    assert seen.get("flex") is False
    # R13: full pLDDT distribution for the report histogram
    assert out["plddt_all"] == [40.0] * 5
    # fast gate (35) keeps pLDDT 40; full gate (50) would drop it
    state["options"] = {}
    monkeypatch.setattr(vh, "model_vhh_one",
                        lambda a: {"ok": True, "plddt": 40.0,
                                   "idx": a[0], "model": "x.pdb"})
    out2 = vh.screen_vhh(state, tmp_path, n=5, n_jobs=1)
    assert out2["n_docked"] == 0


def _fake_pdb(tmp_path, n_res=40, low_range=(18, 24), low_bf=40.0,
              high_bf=80.0, name="v.pdb"):
    """40-residue VHH-like PDB: resSeq 1..n, 1 CA + 1 CB atom each,
    B-factor = pLDDT (high everywhere, low in low_range)."""
    lines = []
    serial = 0
    for res in range(1, n_res + 1):
        bf = low_bf if low_range[0] <= res <= low_range[1] else high_bf
        for atom in ("CA", "CB"):
            serial += 1
            lines.append(
                f"ATOM  {serial:>5d} {atom:>4s} UNK A{res:>4d}    "
                f"{float(res) % 10:8.3f}{float(res) % 7:8.3f}"
                f"{float(res) % 13:8.3f}{1.00:6.2f}{bf:6.2f}          C")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\nEND\n")
    return p


def test_vhh_cdr_fragments_extraction(tmp_path):
    """R11/G10-v2: low-pLDDT runs become fragments, padded by 2, and
    the whole-structure run is rejected."""
    from drugagent.modules.vhh import vhh_cdr_fragments
    pdb = _fake_pdb(tmp_path)
    frags = vhh_cdr_fragments(pdb)
    assert len(frags) == 1
    frag = frags[0]
    # run 18..24 padded by 2 -> 16..26
    assert [n for _, n in frag] == list(range(16, 27))
    assert all(ch == "A" for ch, _ in frag)
    # everything low -> single run covers the structure -> no fragments
    pdb2 = _fake_pdb(tmp_path, name="v2.pdb",
                     low_range=(1, 40), low_bf=30.0)
    assert vhh_cdr_fragments(pdb2) == []
    # nothing low -> no fragments
    pdb3 = _fake_pdb(tmp_path, name="v3.pdb", low_range=(18, 18),
                     low_bf=40.0)
    # run of 1 < min_res=4 -> rejected
    assert vhh_cdr_fragments(pdb3) == []


def test_write_fragment_pdb(tmp_path):
    """R11/G10-v2: fragment PDB contains only the selected residues."""
    from drugagent.modules.vhh import _write_fragment_pdb, vhh_cdr_fragments
    pdb = _fake_pdb(tmp_path)
    frags = vhh_cdr_fragments(pdb)
    out = tmp_path / "frag.pdb"
    _write_fragment_pdb(pdb, frags[0], out)
    got = {int(l[22:26]) for l in out.read_text().splitlines()
           if l.startswith(("ATOM", "HETATM"))}
    assert got == set(range(16, 27))
    assert len(out.read_text().splitlines()) == 23  # 11 res x 2 atoms + END


def test_dock_vhh_candidates_cdr_only(tmp_path, monkeypatch):
    """R11/G10-v2: cdr_only docks fragments (composite = best score);
    falls back to full-VHH docking when no fragment qualifies."""
    from drugagent.modules import vhh as vh
    pdb = _fake_pdb(tmp_path, name="vhh_0.pdb")
    ok = [{"idx": 0, "plddt": 60.0}]
    frag_pdb = tmp_path / "vhh_0_frag0.pdb"
    calls = []

    def fake_frag_dock(model_dir, src, idx, rec, pocket, cpu, flex):
        calls.append((idx, str(src)))
        return {"ok": True, "score": -9.0 - len(calls), "rmsd_lb": 0.0,
                "rmsd_ub": 0.0, "top_pose_pdbqt": None}

    full_calls = []

    def fake_full_dock(model_dir, src, idx, rec, pocket, cpu, flex):
        full_calls.append(idx)
        return {"ok": True, "score": -5.0, "rmsd_lb": 0.0, "rmsd_ub": 0.0,
                "top_pose_pdbqt": None}

    monkeypatch.setattr(vh, "_dock_ligand",
                        lambda md, src, idx, rec, pocket, cpu, flex:
                        fake_frag_dock(md, src, idx, rec, pocket, cpu, flex)
                        if idx.endswith("frag0") else
                        fake_full_dock(md, src, idx, rec, pocket, cpu, flex))
    out = vh.dock_vhh_candidates(ok, tmp_path, "rec.pdbqt",
                                 {"center": [0, 0, 0], "xsize": 12,
                                  "ysize": 12, "zsize": 12},
                                 n_jobs=1, cdr_only=True)
    assert len(out) == 1 and out[0]["ok"]
    assert out[0]["n_fragments"] >= 1
    assert len(out[0]["fragment_scores"]) == out[0]["n_fragments"]
    assert out[0]["score"] == min(out[0]["fragment_scores"])
    # all fragment tags docked, no full dock
    assert all(c[0].endswith("frag0") for c in calls)
    assert not full_calls
    # fallback: pdb with no low-pLDDT run -> full-VHH path
    _fake_pdb(tmp_path, name="vhh_1.pdb", low_range=(18, 18), low_bf=40.0)
    ok2 = [{"idx": 1, "plddt": 60.0}]
    out2 = vh.dock_vhh_candidates(ok2, tmp_path, "rec.pdbqt",
                                  {"center": [0, 0, 0], "xsize": 12,
                                   "ysize": 12, "zsize": 12},
                                  n_jobs=1, cdr_only=True)
    assert out2[0]["score"] == -5.0
    assert "1" in full_calls


def test_fast_cdr_only_default():
    """R11/G10-v2: fast mode docks fragments by default; full mode does
    not; options override wins."""
    from drugagent.config import resolve_defaults
    assert resolve_defaults({"fast": True}).vhh_dock_cdr_only is True
    assert resolve_defaults({}).vhh_dock_cdr_only is False
    assert resolve_defaults(
        {"fast": True, "vhh_dock_cdr_only": False}).vhh_dock_cdr_only is False


def _big_low_pdb(tmp_path, name="big.pdb", n_res=80, low_range=(15, 50)):
    """80-residue PDB with a long (36-res) low-pLDDT run 15..50 (45% of
    the structure — below the 0.6 whole-structure rejection)."""
    lines = []
    serial = 0
    for res in range(1, n_res + 1):
        bf = 40.0 if low_range[0] <= res <= low_range[1] else 80.0
        for atom in ("CA", "CB"):
            serial += 1
            lines.append(
                f"ATOM  {serial:>5d} {atom:>4s} UNK A{res:>4d}    "
                f"{float(res) % 10:8.3f}{float(res) % 7:8.3f}"
                f"{float(res) % 13:8.3f}{1.00:6.2f}{bf:6.2f}          C")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\nEND\n")
    return p


def test_cdr_fragments_chunks_long_runs(tmp_path):
    """R12/G10-v3: a 36-res low run (padded 38) is split into <=20-res
    chunks so no single fragment dominates the dock runtime."""
    from drugagent.modules.vhh import vhh_cdr_fragments
    pdb = _big_low_pdb(tmp_path)
    frags = vhh_cdr_fragments(pdb)
    assert len(frags) >= 2
    for f in frags:
        assert len(f) <= 20
    # chunks are contiguous and cover the padded run
    nums = [n for _, n in frags[0]]
    assert nums == list(range(nums[0], nums[-1] + 1))
    # nothing dropped: union of fragments == padded run 13..52
    alln = sorted({n for f in frags for _, n in f})
    assert alln == list(range(13, 53))


def test_frag_box_adaptive(tmp_path):
    """R12/G10-v3: fragment box is pocket-centered and sized to the
    fragment diameter (+ margin), clipped to [12, 30] A."""
    from drugagent.modules import vhh as vh
    pdb = _big_low_pdb(tmp_path)
    frag = [("A", k) for k in range(20, 30)]  # 10-res fragment
    pocket = {"center": [5.0, -1.5, 14.1], "xsize": 25.84,
              "ysize": 25.84, "zsize": 25.84}
    box = vh._frag_box(pdb, frag, pocket)
    assert box["center"] == pocket["center"]
    # coords in the fake pdb are small (<= 12 A spread) -> box ~ diameter+8
    assert 12.0 <= box["xsize"] <= 30.0
    assert box["xsize"] == box["ysize"] == box["zsize"]
    # a tiny fragment -> clipped to the 12 A floor
    frag_small = [("A", 20)]
    box2 = vh._frag_box(pdb, frag_small, pocket)
    assert box2["xsize"] == 12.0
    # unknown coords (residue not in pdb) -> fallback to pocket box
    box3 = vh._frag_box(pdb, [("B", 1)], pocket)
    assert box3["xsize"] == pocket["xsize"]


def test_dock_vhh_candidates_no_frag_skip(tmp_path, monkeypatch):
    """R12/G10-v3: with full_fallback=False (fast default), a candidate
    without CDR fragments is skipped (ok=False, no_frag, NaN score)
    instead of paying a ~100-min full-VHH dock."""
    import numpy as np
    from drugagent.modules import vhh as vh
    # no low-pLDDT run -> no fragments
    _fake_pdb(tmp_path, name="vhh_2.pdb", low_range=(18, 18), low_bf=40.0)
    ok = [{"idx": 2, "plddt": 60.0}]
    full_calls = []

    def fake_full_dock(model_dir, src, idx, rec, pocket, cpu, flex):
        full_calls.append(idx)
        return {"ok": True, "score": -5.0, "rmsd_lb": 0.0, "rmsd_ub": 0.0,
                "top_pose_pdbqt": None}

    monkeypatch.setattr(vh, "_dock_ligand", fake_full_dock)
    out = vh.dock_vhh_candidates(ok, tmp_path, "rec.pdbqt",
                                 {"center": [0, 0, 0], "xsize": 12,
                                  "ysize": 12, "zsize": 12},
                                 n_jobs=1, cdr_only=True,
                                 full_fallback=False)
    assert len(out) == 1
    assert out[0]["ok"] is False
    assert out[0]["no_frag"] == "full_fallback_off"
    assert np.isnan(out[0]["score"])
    assert out[0]["n_fragments"] == 0
    assert not full_calls  # no full dock paid
    # with full_fallback=True the full dock IS paid (default behavior)
    out2 = vh.dock_vhh_candidates(ok, tmp_path, "rec.pdbqt",
                                  {"center": [0, 0, 0], "xsize": 12,
                                   "ysize": 12, "zsize": 12},
                                  n_jobs=1, cdr_only=True,
                                  full_fallback=True)
    assert out2[0]["ok"] and out2[0]["score"] == -5.0
    assert "2" in full_calls


# --------------------------------------------------------------------------- #
# R13: design validation cache (G8 at design granularity)
# --------------------------------------------------------------------------- #
def _design_pdb(tmp_path, name, res="GLY"):
    lines = []
    for i in range(1, 6):
        lines.append(
            f"ATOM  {i:>5d} {res:>4s} {res[:3]:>4s} A{i:>4d}    "
            f"{float(i):8.3f}{float(i):8.3f}{float(i):8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\nEND\n")
    return p


def test_score_designs_cache(tmp_path, monkeypatch):
    """R13: design validation is cached in scored.json keyed by design
    mtime — re-runs skip already-validated designs (no ESMFold re-pay),
    and a changed design (mtime) is re-validated."""
    from drugagent.modules import vhh as vh
    d = tmp_path / "vhh_designs"
    d.mkdir()
    p0 = _design_pdb(d, "vhh_design_0.pdb")
    p1 = _design_pdb(d, "vhh_design_1.pdb")
    calls = []

    def fake_predict(seqs, **kw):
        calls.append(tuple(seqs))
        return {"mean_plddt": 70.0, "plddt": [80.0] * 20, "pdb": "X"}

    def fake_im(comp_pdb, plddt):
        return {"interface_plddt_mean": 60.0, "interface_plddt_min": 50.0,
                "n_interface_residues": 10}

    # stub out the chain/complex helpers (single-chain fake designs);
    # patch the binder module — score_designs imports them locally
    from drugagent.modules import binder as binder_mod
    monkeypatch.setattr(binder_mod, "_chain_ids", lambda p: ["A"])
    monkeypatch.setattr(binder_mod, "_ca_sequence", lambda p, c: "GGLGG")
    monkeypatch.setattr(binder_mod, "_make_complex",
                        lambda a, b, o: open(o, "w").write("MODEL") or o)
    out1 = vh.score_designs([p0, p1], d, tmp_path / "clean.pdb",
                            predict_fn=fake_predict, im_fn=fake_im)
    assert len(out1) == 2 and calls and len(calls) == 2
    assert (d / "scored.json").is_file()
    assert out1[0]["complex_plddt"] == 70.0

    # second pass: fully cached -> zero new predictions
    calls.clear()
    out2 = vh.score_designs([p0, p1], d, tmp_path / "clean.pdb",
                            predict_fn=fake_predict, im_fn=fake_im)
    assert not calls
    assert out2[0]["complex_plddt"] == 70.0
    assert out2[0]["interface_plddt_mean"] == 60.0

    # changed design (new mtime) -> only that one re-validated
    import time as _t
    _t.sleep(0.01)
    p1.write_text(p1.read_text() + "REMARK changed\n")
    st1 = p1.stat().st_mtime
    p1_utime = st1
    import os
    os.utime(p1, (st1, st1 + 1.0))
    calls.clear()
    out3 = vh.score_designs([p0, p1], d, tmp_path / "clean.pdb",
                            predict_fn=fake_predict, im_fn=fake_im)
    assert len(calls) == 1  # only design_1 re-scored


def test_design_vhh_glob_skips_byproducts(tmp_path, monkeypatch):
    """R13: track B validates only top-level vhh_design_N.pdb, not the
    scoring by-products (vhh_design_N_binder/_complex.pdb) which the
    old glob picked up and re-validated (spawning more by-products)."""
    from drugagent.modules import vhh as vh
    from drugagent.modules import binder as binder_mod
    from drugagent.modules import esmfold_run as esm
    wd = tmp_path / "04_vhh"
    wd.mkdir()
    outdir = wd / "vhh_designs"
    outdir.mkdir()
    (tmp_path / "clean.pdb").write_text("ATOM      1  CA GLY A  1\nEND\n")
    state = {"options": {"fast": True},
             "target_prep": {"clean_pdb": str(tmp_path / "clean.pdb"),
                             "pocket": {"center": [0, 0, 0], "xsize": 10,
                                        "ysize": 10, "zsize": 10}}}
    for name in ("vhh_design_0.pdb", "vhh_design_1.pdb",
                 "vhh_design_0_binder.pdb", "vhh_design_0_complex.pdb"):
        (outdir / name).write_text(_design_pdb_text())
    calls = []

    def fake_predict(seqs, **kw):
        calls.append(tuple(seqs))
        return {"mean_plddt": 65.0, "plddt": [80.0] * 10, "pdb": "X"}

    monkeypatch.setattr(vh, "run_cmd", lambda *a, **k: 0)
    monkeypatch.setattr(binder_mod, "rf_repo", lambda: tmp_path)
    monkeypatch.setattr(binder_mod, "rf_python", lambda: "/bin/true")
    monkeypatch.setattr(binder_mod, "pocket_hotspots",
                        lambda *a: ["A1", "A2"])
    monkeypatch.setattr(binder_mod, "_chain_ids", lambda p: ["A"])
    monkeypatch.setattr(binder_mod, "_ca_sequence", lambda p, c: "GGLGG")
    monkeypatch.setattr(binder_mod, "_make_complex",
                        lambda a, b, o: open(o, "w").write("MODEL") or o)
    monkeypatch.setattr(esm, "predict", fake_predict)
    monkeypatch.setattr(esm, "interface_metrics",
                        lambda comp, plddt: {"interface_plddt_mean": 55.0,
                                             "interface_plddt_min": 40.0,
                                             "n_interface_residues": 5})
    out = vh.design_vhh(state, wd, n_designs=2)
    assert out["n_designs"] == 2, "by-products counted as designs"
    assert len(calls) == 2, "by-products were validated"
    names = [Path(d["design"]).name for d in out["designs"]]
    assert names == ["vhh_design_0.pdb", "vhh_design_1.pdb"]


def _design_pdb_text():
    lines = []
    for i in range(1, 6):
        lines.append(
            f"ATOM  {i:>5d} GLY   GLY A{i:>4d}    "
            f"{float(i):8.3f}{float(i):8.3f}{float(i):8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}")
    return "\n".join(lines) + "\nEND\n"
