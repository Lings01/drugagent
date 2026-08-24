import numpy as np
import pytest

from drugagent.modules import target_prep as tp
from drugagent.utils import centroid_from_pdb, pdb_ligands


def test_pdb_ligands(hivp_pdb):
    ligs = pdb_ligands(hivp_pdb)
    assert "A77" in ligs  # saquinavir


def test_analyze_completeness(hivp_pdb):
    r = tp.analyze_completeness(hivp_pdb)
    assert r["has_ligand"]
    assert "A77" in r["ligands"]
    assert r["n_chains"] == 2  # HIVp dimer
    assert r["multimer"]
    # liganded structure: no apo_target issue
    assert all(i["type"] != "apo_target" for i in r["issues"])


def test_analyze_completeness_apo_target_issue(tiny_pdb):
    """R11/R5: an apo structure carries an info-level apo_target issue so
    the agent (and the report) knows the flexible-target criterion may
    apply after MD."""
    r = tp.analyze_completeness(tiny_pdb)
    assert not r["has_ligand"]
    apo = [i for i in r["issues"] if i["type"] == "apo_target"]
    assert len(apo) == 1
    assert apo[0]["severity"] == "info"
    assert "柔性靶点" in apo[0]["suggestion"]


def test_clean_pdb_removes_waters_ions(hivp_pdb, tmp_path):
    out = tmp_path / "clean.pdb"
    # 1HVI: 3 waters, distances to ligand centroid ~2.1/2.4/2.7 A, so a
    # 2.2 A cutoff keeps 1 and drops 2
    stats = tp.clean_pdb(hivp_pdb, out, keep_resnames=["A77"],
                         keep_waters=True, water_keep_dist=2.2)
    txt = out.read_text()
    het = [l for l in txt.splitlines() if l.startswith("HETATM")]
    resnames = {l[17:20].strip() for l in het}
    assert "A77" in resnames
    assert "K" not in resnames
    assert stats["kept_atoms"] > 0
    assert stats["dropped_atoms"] > 0  # at least the far waters dropped
    assert "HOH" in resnames  # the close water is kept


def test_pocket_from_ligand(hivp_pdb):
    p = tp.pocket_from_ligand(hivp_pdb, "A77")
    c = centroid_from_pdb(hivp_pdb, resnames=["A77"])
    assert all(abs(a - b) < 1e-6 for a, b in zip(p["center"], c))
    assert 16 <= p["xsize"] <= 40
    assert p["method"].startswith("ligand_centroid")


def test_find_pocket_prefers_ligand(hivp_pdb):
    p = tp.find_pocket(hivp_pdb, "A77", brain=None)
    assert p["ligand"] == "A77"
    assert p["site_id"] == "S1"


def test_find_pocket_without_ligand(tiny_pdb):
    p = tp.find_pocket(tiny_pdb, None, brain=None)
    assert "center" in p and p["xsize"] > 0
    assert p["method"] in ("grid_cavity", "protein_centroid_fallback")


def test_grid_pockets_runs(tiny_pdb):
    pockets = tp.grid_pockets(tiny_pdb, grid=0.75)
    assert isinstance(pockets, list)
    for p in pockets:
        assert p["xsize"] > 0
        assert p["volume_A3"] > 0


def test_remove_extract_res(hivp_pdb, tmp_path):
    rec = tmp_path / "rec.pdb"
    tp._remove_res(hivp_pdb, rec, ["A77"])
    rec_lines = [l for l in rec.read_text().splitlines()
                 if l.startswith(("ATOM", "HETATM"))]
    assert not any(l[17:20].strip() == "A77" for l in rec_lines)
    lig = tmp_path / "lig.pdb"
    tp._extract_res(hivp_pdb, lig, ["A77"])
    txt = lig.read_text()
    het = [l for l in txt.splitlines() if l.startswith(("ATOM", "HETATM"))]
    assert len(het) > 20
    assert all(l[17:20].strip() == "A77" for l in het)


@pytest.mark.slow
def test_prepare_target_end_to_end(hivp_pdb, tmp_path, monkeypatch):
    def _stub_pdbqt(pdb, out, **kw):
        out.write_text("stub")
        return out
    monkeypatch.setattr(tp, "to_pdbqt", _stub_pdbqt)
    state = {"project_dir": str(tmp_path / "proj"),
             "target": {"kind": "pdb_file", "value": str(hivp_pdb)},
             "options": {"no_llm": True}}
    out = tp.prepare_target(state)
    prep = out["target_prep"]
    assert Path_exists(prep["clean_pdb"])
    assert prep["pocket"]["ligand"] == "A77"
    assert prep["ligand_pdbqt"].endswith(".pdbqt")
    assert prep["receptor_pdbqt"].endswith(".pdbqt")


def Path_exists(p) -> bool:
    from pathlib import Path
    return Path(p).exists()


# --------------------------------------------------------------------------- #
# structural pitfall detection + repair
# --------------------------------------------------------------------------- #
def _atom(serial, name, resname, chain, resseq, x=1.0, y=2.0, z=3.0,
          alt=" ", element="C"):
    # exact PDB columns: name 13-16, altloc 17, resname 18-20, chain 22,
    # resseq 23-26, xyz 31-54, element 77-78 (79 cols total)
    return (f"ATOM  {serial:5d} {name:<4s}{alt}{resname:3s} {chain}{resseq:4d}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  0.00  0.00    +0.000 {element}  ")


def _het(serial, name, resname, chain, resseq, element="Z"):
    return (f"HETATM{serial:5d} {name:<4s} {resname:3s} {chain}{resseq:4d}   "
            f"{1.0:8.3f}{2.0:8.3f}{3.0:8.3f}  0.00  0.00    +0.000 {element}  ")


def test_detect_multi_model_and_altloc(tmp_path):
    lines = ["MODEL      1"]
    for i in range(1, 21):
        lines.append(_atom(i, "CA ", "ALA", "A", i))
    lines.append("ENDMDL")
    lines.append("MODEL      2")
    for i in range(1, 21):
        lines.append(_atom(i, "CA ", "ALA", "A", i, x=5.0))
    lines.append("ENDMDL")
    lines.append("END")
    p = tmp_path / "nmr.pdb"
    p.write_text("\n".join(lines) + "\n")
    issues = tp.detect_structure_issues(p)
    types = {i["type"] for i in issues}
    assert "multiple_models" in types
    mm = [i for i in issues if i["type"] == "multiple_models"][0]
    assert mm["severity"] == "error"
    assert mm["fix"] == "dedupe_models"


def test_detect_metals_nucleic_and_gaps(tmp_path):
    lines = []
    for i in (1, 2, 3, 4, 5):
        lines.append(_atom(i, "CA ", "ALA", "A", i))
    # gap: residues 6-8 missing
    for i in (9, 10):
        lines.append(_atom(5 + i, "CA ", "ALA", "A", i))
    lines.append(_het(1, "ZN  ", "ZN ", "A", 1, element="ZN"))
    # a nucleotide (2-letter code)
    lines.append(_het(2, "P1  ", "DA ", "B", 1, element="P"))
    lines.append("END")
    p = tmp_path / "mix.pdb"
    p.write_text("\n".join(lines) + "\n")
    issues = tp.detect_structure_issues(p)
    types = {i["type"] for i in issues}
    assert "metals" in types
    assert "nucleic_acids" in types
    assert "missing_residues" in types
    met = [i for i in issues if i["type"] == "metals"][0]
    assert "ZN" in met["detail"]


def test_detect_disordered_termini(tmp_path):
    # 12 residues; first 4 backbone-only (disordered N-terminus)
    lines = []
    for i in range(1, 13):
        n_atoms = 2 if i <= 4 else 5  # backbone-only vs full
        for j in range(n_atoms):
            lines.append(_atom((i - 1) * 6 + j + 1, "CA ", "ALA", "A", i))
    lines.append("END")
    p = tmp_path / "dis.pdb"
    p.write_text("\n".join(lines) + "\n")
    issues = tp.detect_structure_issues(p)
    assert any(i["type"] == "disordered_termini" for i in issues)


def test_repair_dedupe_models_and_altloc(tmp_path):
    lines = ["MODEL      1"]
    for i in range(1, 11):
        lines.append(_atom(i, "CA ", "ALA", "A", i))
    # altloc B atom (should be dropped)
    lines.append(_atom(11, "CA ", "ALA", "A", 1, alt="B", x=9.0))
    lines.append("ENDMDL")
    lines.append("MODEL      2")
    for i in range(1, 11):
        lines.append(_atom(i, "CA ", "ALA", "A", i, x=5.0))
    lines.append("ENDMDL")
    lines.append("END")
    src = tmp_path / "in.pdb"
    src.write_text("\n".join(lines) + "\n")
    out = tmp_path / "out.pdb"
    res = tp.repair_structure(src, out,
                              actions=["dedupe_models", "keep_altloc_a"])
    txt = out.read_text().splitlines()
    atoms = [l for l in txt if l.startswith("ATOM")]
    # only model 1 (10 atoms) + altloc A (altloc B dropped) -> 10 atoms
    assert len(atoms) == 10
    assert "MODEL      2" not in txt
    assert res["applied"] == ["dedupe_models", "keep_altloc_a"]
    assert sum(res["removed_atoms"].values()) > 0


def test_repair_drop_metals(tmp_path):
    lines = [_atom(1, "CA ", "ALA", "A", 1),
             _atom(2, "CA ", "ALA", "A", 2),
             _het(1, "ZN  ", "ZN ", "A", 1, element="ZN"),
             "END"]
    src = tmp_path / "in.pdb"
    src.write_text("\n".join(lines) + "\n")
    out = tmp_path / "out.pdb"
    tp.repair_structure(src, out, actions=["drop_metals"])
    txt = out.read_text()
    assert "ZN " not in [l[17:20].strip() for l in txt.splitlines()
                         if l.startswith("HETATM")]


def test_analyze_includes_issues(hivp_pdb):
    r = tp.analyze_completeness(hivp_pdb)
    assert "issues" in r  # 1HVI is clean-ish; just ensure key present
