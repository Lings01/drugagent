import re
from pathlib import Path

import pytest

from drugagent.config import TOOLS
from drugagent.modules import binder as bd
from drugagent.modules import target_prep as tp


def test_pocket_hotspots(hivp_pdb):
    pocket = tp.pocket_from_ligand(hivp_pdb, "A77")
    hs = bd.pocket_hotspots(hivp_pdb, pocket, max_res=10)
    assert len(hs) <= 10
    assert all(re.fullmatch(r"A\d+", h) for h in hs)
    assert len(set(hs)) == len(hs)


def test_ca_sequence(hivp_pdb):
    seq = bd._ca_sequence(hivp_pdb, "A")
    assert 95 < len(seq) < 200  # HIVp chain ~99-110
    assert set(seq) <= set("ACDEFGHIKLMNPQRSTVWY")


def test_make_complex(hivp_pdb, tmp_path):
    b = tmp_path / "binder.pdb"
    b.write_text(hivp_pdb.read_text())  # reuse as pseudo-binder
    out = tmp_path / "complex.pdb"
    bd._make_complex(hivp_pdb, b, out)
    chains = set()
    for line in out.read_text().splitlines():
        if line.startswith("ATOM"):
            chains.add(line[21])
    assert chains == {"A", "B"}


def test_rename_chain(hivp_pdb, tmp_path):
    out = tmp_path / "renamed.pdb"
    bd._rename_chain(hivp_pdb, out, "Z")
    for line in out.read_text().splitlines():
        if line.startswith("ATOM"):
            assert line[21] == "Z"
            break


@pytest.mark.slow
def test_rfdesign_one(hivp_pdb, tmp_path):
    if not (TOOLS / "RFdiffusion" / "models" / "Base_ckpt.pt").exists():
        pytest.skip("RFdiffusion weights not downloaded")
    pocket = tp.pocket_from_ligand(hivp_pdb, "A77")
    designs = bd.rfdesign(hivp_pdb, pocket, tmp_path / "rf",
                          n_designs=1, length=(60, 70), hotspots=None)
    assert len(designs) >= 1
    txt = designs[0].read_text()
    assert "ATOM" in txt


def test_design_chain_extraction(tmp_path):
    """RF inpaint outputs: chain A = all-GLY design region, chain B = named target."""
    lines = []
    for ch, res, n in (("A", "GLY", 5), ("B", "ALA", 4)):
        for i in range(1, n + 1):
            lines.append(f"ATOM  {10*i:5d}  N   {res} {ch}{i:4d}    "
                         f"  1.000   2.000   3.000  1.00  0.00          N  ")
            lines.append(f"ATOM  {10*i+1:5d}  CA  {res} {ch}{i:4d}    "
                         f"  2.000   2.000   3.000  1.00  0.00           C  ")
    pdb = tmp_path / "design_multi.pdb"
    pdb.write_text("\n".join(lines) + "END\n")
    assert bd._chain_ids(pdb) == ["A", "B"]
    assert bd._design_chain(pdb) == "A"
    out = tmp_path / "binder_only.pdb"
    bd._extract_chain(pdb, "A", out)
    got = [l[21] for l in out.read_text().splitlines() if l.startswith("ATOM")]
    assert all(c == "A" for c in got) and len(got) == 10
    # single-chain PDB: the only chain is the design chain
    pdb1 = tmp_path / "design_one.pdb"
    pdb1.write_text("\n".join(l for l in lines if l[21] == "A") + "END\n")
    assert bd._design_chain(pdb1) == "A"


def test_pdbqt_element_and_graph_sanitizer(tmp_path):
    """obabel quirks: bogus element 'A' on ALA CA; missing ROOT/TORSDOF that
    this vina build requires for flex ligands."""
    from drugagent.modules.target_prep import (
        _add_molecule_graph, _sanitize_pdbqt_elements)
    pdbqt = tmp_path / "vhh.pdbqt"
    lines = [
        "REMARK  Name = vhh.pdb",
        "REMARK  2 active torsions:",
        "REMARK  status: ('A' for Active; 'I' for Inactive)",
        "REMARK    1  A    between atoms: CA_2  and  C_3",
        "REMARK    2  A    between atoms: C_3  and  CB_4",
    ]
    for i, name in enumerate(["CA", "C", "O", "CA"], start=1):
        # exact obabel line shape (79 cols, element at 77-78)
        lines.append(
            f"ATOM  {i:5d}  {name:<2s}  ALA A{i % 2:4d}    "
            f"       1.000   2.000   3.000  0.00  0.00    +0.000 A ")
    pdbqt.write_text("\n".join(lines) + "\n")
    _sanitize_pdbqt_elements(pdbqt, add_graph=True)
    out = pdbqt.read_text().splitlines()
    atoms = [l for l in out if l.startswith("ATOM")]
    # bogus 'A' replaced; valid elements (C/N/O/S) untouched
    assert all(l[76:78].strip() != "A" for l in atoms)
    assert all(l[76:78].strip() in {"C", "N", "O", "S"} for l in atoms)
    assert any(l == "ROOT" for l in out)
    assert any(l == "ENDROOT" for l in out)
    assert any(l == "TORSDOF 2" for l in out)
    # rigid receptor: no graph keywords added
    rec = tmp_path / "rec.pdbqt"
    rec.write_text("\n".join(lines) + "\n")
    _sanitize_pdbqt_elements(rec, add_graph=False)
    rec_out = rec.read_text().splitlines()
    assert not any(l.startswith("ROOT") for l in rec_out)


def test_ca_sequence_full_names(tmp_path):
    """_ca_sequence must use the full three-letter code, not name[0]
    (the old bug mapped ASP->A->ALA, GLU->G->GLY, LYS->L->LEU, ...)."""
    from drugagent.modules.binder import _ca_sequence
    # one residue per interesting code, chain A
    lines = []
    serial = 1
    for i, name in enumerate(["ALA", "ASP", "GLU", "LYS", "ARG", "GLN",
                              "TRP", "TYR", "ASN", "GLY"], start=1):
        x = i * 3.8
        lines.append(
            f"ATOM  {serial:5d}  CA  {name} A{i:4d}   "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          C\n")
        serial += 1
    pdb = tmp_path / "seq.pdb"
    pdb.write_text("".join(lines) + "END\n")
    # ALA ASP GLU LYS ARG GLN TRP TYR ASN GLY -> A D E K R Q W Y N G
    assert _ca_sequence(pdb) == "ADEKRQWYNG"


def test_binder_score_designs_cache(tmp_path, monkeypatch):
    """R14: binder design validation is cached (scored.json) keyed by
    design+target mtime + alt sequence — re-runs skip ESMFold, changed
    inputs re-score."""
    from drugagent.modules import esmfold_run as esm
    wd = tmp_path / "03_binder"
    wd.mkdir()

    def _pdb(residue_chain="A", n=3):
        lines = []
        for i in range(1, n + 1):
            lines.append(
                f"ATOM  {i:>5d}  CA  GLY {residue_chain}{i:>4d}    "
                f"{float(i):8.3f}{float(i):8.3f}{float(i):8.3f}"
                f"{1.00:6.2f}{0.00:6.2f}")
        return "\n".join(lines) + "\nEND\n"

    (tmp_path / "target.pdb").write_text(_pdb("A", 4))
    d0 = wd / "binder_design_0.pdb"
    d0.write_text(_pdb("B", 3))
    calls = []

    def fake_predict(seq, **kw):
        calls.append(seq)
        return {"mean_plddt": 60.0, "plddt": [70.0] * 10,
                "min_plddt": 40.0, "pdb": "X"}

    monkeypatch.setattr(esm, "predict", fake_predict)
    monkeypatch.setattr(esm, "interface_metrics",
                        lambda comp, p: {"interface_plddt_mean": 50.0,
                                         "interface_plddt_min": 30.0,
                                         "n_interface_residues": 4})
    out1 = bd.score_designs([d0], tmp_path / "target.pdb", wd, {})
    assert out1[0]["interface_plddt_mean"] == 50.0
    assert len(calls) == 2  # mono + complex

    calls.clear()
    out2 = bd.score_designs([d0], tmp_path / "target.pdb", wd, {})
    assert not calls, "cached design was re-validated"
    assert out2[0]["interface_plddt_mean"] == 50.0
    assert out2[0]["mono_plddt"] == 60.0

    import os
    st = d0.stat()
    os.utime(d0, (st.st_atime, st.st_mtime + 5.0))
    calls.clear()
    bd.score_designs([d0], tmp_path / "target.pdb", wd, {})
    assert len(calls) == 2, "mtime-changed design must be re-validated"
