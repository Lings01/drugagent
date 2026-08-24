import numpy as np
import pytest
from pathlib import Path

from drugagent.config import TOOLS, resolve_defaults
from drugagent.modules import screening as sc


def test_standardize_sdf(small_sdf, tmp_path):
    df = sc.standardize_sdf(small_sdf, tmp_path / "ok.sdf", n_jobs=2)
    assert len(df) == 5
    assert (df.status == "ok").sum() >= 4
    assert df[df.status == "ok"]["smiles"].notna().all()


def test_compute_features(small_sdf, tmp_path):
    df = sc.standardize_sdf(small_sdf, tmp_path / "ok.sdf", n_jobs=2)
    ok = df[df.status == "ok"]
    fdf = sc.compute_features(ok)
    assert len(fdf) == len(ok)
    for col in ("MW", "LogP", "QED", "SA"):
        assert col in fdf.columns
    assert fdf.MW.between(50, 700).all()


def test_physchem_filter(small_sdf, tmp_path):
    df = sc.standardize_sdf(small_sdf, tmp_path / "ok.sdf", n_jobs=2)
    ok = df[df.status == "ok"]
    fdf = sc.compute_features(ok)
    kept = sc.physchem_filter(fdf)
    assert len(kept) >= 3  # all our small mols should pass
    assert len(kept) <= len(fdf)


def test_ml_prefilter_fallback(small_sdf, tmp_path):
    df = sc.standardize_sdf(small_sdf, tmp_path / "ok.sdf", n_jobs=2)
    fdf = sc.compute_features(df[df.status == "ok"])
    out = sc.ml_prefilter(fdf, model_path=None, n_keep=3)
    assert len(out) <= 3
    assert "ml_score" in out.columns
    assert out.ml_score.is_monotonic_decreasing


def test_decide_hits_with_reference():
    scores = list(np.random.normal(-8, 1.5, 500))
    h = sc.decide_hits(scores, ref_score=-9.0, n_wanted=20, brain=None)
    assert h["threshold"] == pytest.approx(-7.0, abs=1e-6)
    assert h["n_hits"] <= 20
    assert h["n_hits"] >= 5
    assert h["rationale"].startswith("known ligand")


def test_decide_hits_without_reference():
    rng = np.random.default_rng(0)
    scores = list(rng.normal(-7, 1.2, 300))
    h = sc.decide_hits(scores, ref_score=None, n_wanted=20, brain=None)
    p5 = float(np.percentile(scores, 5))
    assert h["threshold"] == pytest.approx(p5, abs=1e-3)
    assert h["n_hits"] >= 5


def test_write_ligand_pdbqt(tmp_path):
    out = tmp_path / "mol.pdbqt"
    ok = sc.write_ligand_pdbqt("c1ccccc1", out, seed=1)
    assert ok
    txt = out.read_text()
    assert "ROOT" in txt  # ligands keep the ROOT wrapper (vina requires it)
    assert txt.count("ATOM") > 0


@pytest.mark.slow
def test_dock_one_with_vina(hivp_pdb, tmp_path):
    import shutil
    if shutil.which("vina") is None and not (
            (TOOLS / "vina" / "bin" / "vina").exists()):
        pytest.skip("vina not installed")
    from drugagent.modules import target_prep as tp
    rec_pdb = tmp_path / "rec.pdb"
    tp._remove_res(hivp_pdb, rec_pdb, ["A77"])
    # rigid receptor: no AD4 graph keywords (this vina build rejects
    # ROOT/TORSDOF in --receptor files)
    rec = tp.to_pdbqt(rec_pdb, tmp_path / "rec.pdbqt", flex=False)
    lig = tmp_path / "aspirin.pdbqt"
    assert sc.write_ligand_pdbqt("CC(=O)Oc1ccccc1C(=O)O", lig, seed=3)
    pocket = tp.pocket_from_ligand(hivp_pdb, "A77")
    r = sc.dock_one((str(rec), str(lig), str(tmp_path / "aspirin"),
                     pocket, 4, 4))
    assert r["ok"]
    assert -15 < r["score"] < 0


# --------------------------------------------------------------------------- #
# R2: receptor conformation selection + flexible docking
# --------------------------------------------------------------------------- #
import numpy as np  # noqa: E402
from drugagent.modules import screening as sc  # noqa: E402


def _rec_pdb_lines(n_res=4):
    """Minimal receptor PDB: residues 1..n on chain A along the x axis
    (3.8 A spacing). Exact 79-col PDB format (R0-verified template)."""
    lines = []
    serial = 1
    for res in range(1, n_res + 1):
        for name in ("N", "CA", "C", "O", "CB"):
            x = res * 3.8
            y = 0.0 if name in ("N", "CA", "C", "O") else 1.5
            lines.append(
                f"ATOM  {serial:5d} {name:<4s} ALA A{res:4d}   "
                f"{x:8.3f}{y:8.3f}{0.0:8.3f}  1.00  0.00          C\n")
            serial += 1
    lines.append("END\n")
    return lines


def _lig_pdb_lines(center, n_atoms=4):
    """Tiny ligand around a given center (79-col format)."""
    lines = []
    for i in range(n_atoms):
        x = center[0] + 0.5 * i
        lines.append(
            f"HETATM{i + 1:5d} {('C' + str(i + 1)):<4s} LIG A 999   "
            f"{x:8.3f}{center[1]:8.3f}{center[2]:8.3f}  1.00  0.00          C\n")
    lines.append("END\n")
    return lines


def test_flex_sidechain_selection(tmp_path):
    # ligand near residue 2 (x ~ 7.6): residues 1..3 within 5 A, 4 not
    rec = tmp_path / "rec.pdb"
    rec.write_text("".join(_rec_pdb_lines(4)))
    lig = tmp_path / "lig.pdb"
    lig.write_text("".join(_lig_pdb_lines((7.6, 0.5, 0.5))))
    out = tmp_path / "flex.pdbqt"
    r = sc.flex_sidechain_pdbqt(rec, lig, out, cutoff=5.0)
    assert r["n_residues"] == 3, r
    txt = out.read_text()
    # AD4 flex-residue format: one BEGIN_RES/END_RES block per residue with
    # a torsion graph
    assert txt.count("BEGIN_RES") == 3
    assert txt.count("END_RES") == 3
    assert "ROOT" in txt and "ENDROOT" in txt
    assert "BRANCH" in txt and "ENDBRANCH" in txt
    # side chains only: CB present; backbone N/C/O absent; CA allowed only
    # as the immobile graph root (exactly one CA line per residue)
    atom_lines = [l for l in txt.splitlines() if l.startswith(("ATOM", "HETATM"))]
    names = [l[12:16].strip() for l in atom_lines]
    assert "CB" in names
    assert "N" not in names and "C" not in names and "O" not in names
    assert names.count("CA") == 3  # one root per selected residue


def test_flex_sidechain_no_ligands_outside(tmp_path):
    rec = tmp_path / "rec.pdb"
    rec.write_text("".join(_rec_pdb_lines(4)))
    lig = tmp_path / "lig.pdb"
    lig.write_text("".join(_lig_pdb_lines((500.0, 0.0, 0.0))))
    out = tmp_path / "flex.pdbqt"
    r = sc.flex_sidechain_pdbqt(rec, lig, out, cutoff=5.0)
    assert r["n_residues"] == 0


def test_consensus_stats():
    r = sc.consensus_stats([-7.0, -7.4, -8.1])
    assert abs(r["mean"] - (-7.5)) < 0.01
    assert r["min"] == -8.1
    assert r["n"] == 3


def test_cluster_representatives_smoke():
    """R8: adaptive — the expected rep count is derived from the
    trajectory's actual cluster populations (the smoke run is rebuilt
    across rounds; the 1HVI 2 ns run is one dominant cluster)."""
    base = Path(__file__).resolve().parent.parent / \
        "projects/agent_smoke_0821_0404/05_md"
    if not (base / "analysis" / "clusters_r1.pdb").is_file():
        pytest.skip("smoke MD cluster reps not present")
    from drugagent.modules.mdsim import _parse_xvg
    t, y = _parse_xvg(base / "analysis" / "cluster_idx_r1.xvg")
    assert y, "empty cluster index"
    yr = np.round(np.asarray(y, dtype=float)).astype(int)
    min_pop = 0.004
    expected = sum(1 for c in set(yr.tolist())
                   if (np.sum(yr == c) / len(yr)) >= min_pop)
    expected = min(expected, 3)  # max_n
    reps = sc.cluster_representatives(base, rep=1, max_n=3, min_pop=min_pop)
    assert len(reps) == expected
    assert len(reps) >= 1
    for rep in reps:
        assert rep["population"] > 0
        pdb = rep["pdb"]
        assert pdb.is_file()
        txt = pdb.read_text()
        # protein only: no water residue names
        assert "HOH" not in txt and "SOL" not in txt
        assert "ATOM" in txt


def test_kabsch_identity_and_rotation():
    """Kabsch: identity for same frame; recovers a known rotation."""
    rng = np.random.default_rng(7)
    pts = rng.normal(0, 5, size=(50, 3))
    R, t = sc.kabsch(pts, pts)
    out = pts @ R.T + t
    assert np.max(np.abs(out - pts)) < 1e-6
    # rotate by 30 deg about z
    th = np.deg2rad(30)
    Rz = np.array([[np.cos(th), -np.sin(th), 0],
                   [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    moved = pts @ Rz.T + np.array([1.0, -2.0, 3.0])
    R2, t2 = sc.kabsch(pts, moved)
    out2 = pts @ R2.T + t2
    assert np.max(np.abs(out2 - moved)) < 1e-6


def test_align_pdb_to_reference(tmp_path):
    """align_pdb_to_reference moves a translated copy back onto the ref."""
    rec = tmp_path / "ref.pdb"
    rec.write_text("".join(_rec_pdb_lines(4)))
    # shifted copy: +10 in x, different serials fine (same coords+shift)
    shifted = []
    for l in _rec_pdb_lines(4):
        if l.startswith("ATOM"):
            x = float(l[30:38]) + 10.0
            shifted.append(l[:30] + f"{x:8.3f}" + l[38:])
        else:
            shifted.append(l)
    mov = tmp_path / "mov.pdb"
    mov.write_text("".join(shifted))
    out = tmp_path / "aligned.pdb"
    n = sc.align_pdb_to_reference(mov, rec, out)
    assert n == 4  # one CA per residue, all matched
    # aligned file: CA of residue 1 should match ref within 0.01 A
    ref_atoms = sc._atoms(rec)
    mov_atoms = sc._atoms(out)
    for a in ref_atoms:
        if a["name"] == "CA":
            b = next(x for x in mov_atoms
                     if x["name"] == "CA"
                     and x["resseq"] == a["resseq"])
            d = float(np.sqrt(sum((a["xyz"][k] - b["xyz"][k]) ** 2
                                  for k in range(3))))
            assert d < 0.05, (a["resseq"], d)


def test_align_renumbered_chain(tmp_path):
    """Chain B of the model is the same sequence as ref chain B but renumbered
    (+600 offset); alignment must still recover the frame."""
    lines = []
    serial = 1
    # chain A: 4 residues 1-4; chain B: 4 residues 601-604 (offset 600)
    for chain, start in (("A", 1), ("B", 601)):
        for i in range(4):
            for name in ("N", "CA", "C"):
                x = (10.0 if chain == "B" else 0.0) + i * 3.8
                lines.append(
                    f"ATOM  {serial:5d} {name:<4s} ALA {chain}{start + i:4d}   "
                    f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          C\n")
                serial += 1
    ref = tmp_path / "ref.pdb"
    ref_lines = []
    serial = 1
    for chain, start in (("A", 1), ("B", 1)):
        for i in range(4):
            for name in ("N", "CA", "C"):
                x = (10.0 if chain == "B" else 0.0) + i * 3.8
                ref_lines.append(
                    f"ATOM  {serial:5d} {name:<4s} ALA {chain}{start + i:4d}   "
                    f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          C\n")
                serial += 1
    ref.write_text("".join(ref_lines) + "END\n")
    # shifted model
    mov_lines = []
    for l in lines:
        if l.startswith("ATOM"):
            x = float(l[30:38]) + 7.0
            mov_lines.append(l[:30] + f"{x:8.3f}" + l[38:])
        else:
            mov_lines.append(l)
    mov = tmp_path / "mov.pdb"
    mov.write_text("".join(mov_lines) + "END\n")
    out = tmp_path / "out.pdb"
    n = sc.align_pdb_to_reference(mov, ref, out)
    assert n >= 8  # both chains matched
    mov2 = sc._atoms(out)
    ref2 = sc._atoms(ref)
    for a in ref2:
        b = next(x for x in mov2
                 if x["name"] == a["name"] and x["chain"] == a["chain"]
                 and abs(x["serial"] - a["serial"]) <= 1
                 and x["resseq"] in (a["resseq"], a["resseq"] + 600))
        d = float(np.sqrt(sum((a["xyz"][k] - b["xyz"][k]) ** 2
                              for k in range(3))))
        assert d < 0.05, (a["chain"], a["resseq"], d)

# ---------------------------------------------------------------- R7: pool


def _ca_pdb(path, coords, box=None):
    lines = []
    if box is not None:
        lines.append("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1"
                     % (box[0], box[1], box[2]))
    for i, (x, y, z) in enumerate(coords, start=1):
        lines.append(f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return path


def test_dock_one_cpu_cap(monkeypatch, tmp_path):
    """R10: vina --cpu is capped at the physical core count (with a
    warning) — oversubscription slows the run and RSS grows ~linearly,
    especially with --flex."""
    rec = tmp_path / "r.pdbqt"
    rec.write_text("REMARK test\nATOM      1  C   LIG A  1   0.0 0.0 0.0"
                   "\nENDMODEL\n")
    lig = tmp_path / "l.pdbqt"
    lig.write_text("REMARK test\nATOM      1  C   LIG A  1   1.0 0.0 0.0"
                   "\nENDMODEL\n")
    pocket = {"center": (0.0, 0.0, 0.0), "xsize": 10.0, "ysize": 10.0,
              "zsize": 10.0}
    captured = {}

    class _P:
        returncode = 1
        stdout = b"no table"

    def fake_run_cmd(cmd, **kw):
        captured["cmd"] = cmd
        return _P()

    monkeypatch.setattr(sc, "run_cmd", fake_run_cmd)
    import os
    monkeypatch.setattr(os, "cpu_count", lambda: 8)
    prefix = str(tmp_path / "pose")
    r = sc.dock_one((str(rec), str(lig), prefix, pocket, 4, 64))
    assert r["ok"] is False  # fake vina -> no score table
    i = captured["cmd"].index("--cpu")
    assert captured["cmd"][i + 1] == "8"


def test_ca_rmsd_pbc_frames(tmp_path):
    """R9: two C-alpha sets differing by PER-ATOM integer box vectors
    (simulation-box frame vs another wrap) must give the TRUE
    conformational RMSD, not the muddled plain-Kabsch number."""
    box = np.array([10.0, 10.0, 10.0])  # A
    n = 24
    A = np.zeros((n, 3))
    for i in range(n):
        A[i] = (i * 0.38, (i % 5) * 0.2, 0.0)
    # physically coherent wrap: shift the whole cloud by a NON-integer
    # box vector and wrap per atom (atoms crossing a boundary jump a full
    # box — mixed per-atom integer vectors, spatially coherent like a
    # real pdb2gmx wrap)
    shift = np.array([6.0, 4.5, 0.0])
    B = (A + shift) % box
    pa = _ca_pdb(tmp_path / "a.pdb", A, box=box)
    pb = _ca_pdb(tmp_path / "b.pdb", B, box=box)
    d = sc._ca_rmsd(pa, pb)
    assert d is not None and d < 0.05, d
    # sanity: plain Kabsch on the same clouds is badly muddled
    R, t = sc.kabsch(A, B)
    raw = float(np.sqrt(np.mean(((A @ R.T + t - B) ** 2).sum(1))))
    assert raw > 1.0, raw
    # no CRYST1: plain behavior (same-frame clouds still fine)
    pc = _ca_pdb(tmp_path / "c.pdb", A)
    pd_ = _ca_pdb(tmp_path / "d.pdb", B)
    d2 = sc._ca_rmsd(pc, pd_)
    assert d2 is not None


def _cluster_artifacts(anadir, rep, frac_c1, n=10):
    """cluster_idx xvg (frac_c1 fraction in cluster 1) + 2-model clusters pdb."""
    idx = []
    for i in range(n):
        c = 1 if i < round(n * frac_c1) else 2
        idx.append(f"  {i * 10.0:8.3f}    {c}")
    (anadir / f"cluster_idx_r{rep}.xvg").write_text(
        "# cluster index\n" + "\n".join(idx) + "\n")
    p0 = [(i * 0.5, 0.0, 0.0) for i in range(10)]           # straight
    p1 = ([(i * 0.5, 0.0, 0.0) for i in range(5)]
          + [(4 * 0.5, j * 0.5, 0.0) for j in range(5)])    # L-bend
    models = []
    for coords in (p0, p1):
        models.append("MODEL 1\n" if not models else "MODEL 2\n")
        for i, (x, y, z) in enumerate(coords, start=1):
            models.append(
                f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
        models.append("ENDMDL\n")
    (anadir / f"clusters_r{rep}.pdb").write_text("".join(models))


def test_rep_files_prefer_merged_xtc(tmp_path):
    repdir = tmp_path / "md_rep1"
    repdir.mkdir()
    (repdir / "md.tpr").write_bytes(b"x")
    (repdir / "md.xtc").write_bytes(b"x")
    tpr, xtc = sc._rep_files(repdir)
    assert xtc.name == "md.xtc"
    (repdir / "md_all.xtc").write_bytes(b"x")
    tpr, xtc = sc._rep_files(repdir)
    assert xtc.name == "md_all.xtc"


def test_dump_group_system_aware(monkeypatch, tmp_path):
    import drugagent.modules.mdsim as md

    def fake_index(tpr, gmx, out):
        return {"System", "DNA", "Water"}
    monkeypatch.setattr(md, "_index_groups", fake_index)
    assert sc._dump_group(tmp_path / "x.tpr", "gmx", tmp_path) == "DNA"

    def fake_index2(tpr, gmx, out):
        return {"Backbone", "Protein", "System"}
    monkeypatch.setattr(md, "_index_groups", fake_index2)
    assert sc._dump_group(tmp_path / "x.tpr", "gmx", tmp_path) == "Backbone"

    def fake_index3(tpr, gmx, out):
        raise RuntimeError("no gmx")
    monkeypatch.setattr(md, "_index_groups", fake_index3)
    assert sc._dump_group(tmp_path / "x.tpr", "gmx", tmp_path) == "Protein"


def test_pool_representatives_dedup(tmp_path):
    anadir = tmp_path / "analysis"
    anadir.mkdir()
    # replica 1: 70/30; replica 2: 60/40 — r2c1's model is a 0.05 A shift
    # of r1c1 (duplicate), r2c2 is an L-bend (new shape)
    _cluster_artifacts(anadir, 1, 0.7)
    _cluster_artifacts(anadir, 2, 0.6)
    # overwrite r2 models: c1 = near-duplicate of r1c1, c2 = L-bend
    p0 = [(i * 0.5, 0.0, 0.0) for i in range(10)]
    p0s = [(x + 0.05, y, z) for (x, y, z) in p0]
    pl = ([(i * 0.5, 0.0, 0.0) for i in range(5)]
          + [(4 * 0.5, j * 0.5, 0.0) for j in range(5)])
    models = []
    for k, coords in enumerate((p0s, pl)):
        models.append(f"MODEL {k + 1}\n")
        for i, (x, y, z) in enumerate(coords, start=1):
            models.append(
                f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
        models.append("ENDMDL\n")
    (anadir / "clusters_r2.pdb").write_text("".join(models))

    # dedup_nm=0.3: the 0.05 A-shift duplicate is far inside the threshold,
    # the L-bend (RMSD ~0.77 A against the straight model) is far outside
    pool = sc.pool_representatives(tmp_path, max_n=3, min_pop=0.05,
                                   dedup_nm=0.3)
    assert len(pool) == 2, pool
    # population order is preserved; the duplicate (r2c1, 0.6) is dropped
    assert pool[0]["rep"] == 1 and pool[0]["cluster"] == 1
    assert pool[0]["population"] == 0.7
    assert not any(r["rep"] == 2 and r["cluster"] == 1 for r in pool)
    # the L-bend (r2c2, 0.4) fills the remaining slot; r1c2 is the same
    # shape as r2c2 (identical model) so it dedups against it too
    assert pool[1]["rep"] == 2 and pool[1]["cluster"] == 2
    assert pool[1]["population"] == 0.4
    for r in pool:
        assert Path(r["pdb"]).is_file()

@pytest.mark.slow
def test_pool_representatives_dna(dna_pdb, tmp_path):
    """@slow R7 e2e: 1BNA 2-replica 0.5 ns MD -> analysis -> pooled
    representatives. The DNA-only tpr has no 'Protein' group, so the
    full-atom dump must go through the system-aware 'DNA' group and
    yield full_atoms=True (base/sugar phosphate atoms present)."""
    from drugagent.modules import mdsim as md
    env = md.gromacs()
    w = tmp_path / "md"
    w.mkdir()
    sysinfo = md.build_system(dna_pdb, w, env, is_ligand=False)
    reps = md.run_replicas(sysinfo, w, env, ns=0.5, reps=2)
    md.analyze_replicas(reps, w, env, is_ligand=False)
    pool = sc.pool_representatives(w, max_n=3, min_pop=0.05, gmx=env["gmx"])
    assert pool, "no cluster representatives produced"
    assert pool[0]["population"] >= pool[-1]["population"]
    for r in pool:
        assert Path(r["pdb"]).is_file()
        txt = Path(r["pdb"]).read_text()
        assert "ATOM" in txt
        # full-atom DNA dump: phosphate/base atoms beyond {N, CA, C, O}
        assert r["full_atoms"], txt[:200]
        assert "OP1" in txt or "P" in txt.split()
    # dedup invariant: any two pooled reps differ by >= 1.0 A Cα RMSD
    import itertools
    for a, b in itertools.combinations(pool, 2):
        d = sc._ca_rmsd(a["pdb"], b["pdb"])
        assert d is None or d >= 1.0, (a["pdb"], b["pdb"], d)

# ---------------------------------------------------------------- R8: cluster-mean reps


@pytest.mark.slow
def test_cluster_mean_pdb(tmp_path, dna_pdb):
    """R8: the cluster representative is the MEAN structure over ALL
    cluster frames (MDAnalysis), not an arbitrary first-frame snapshot."""
    import numpy as np
    import MDAnalysis as mda
    from drugagent.modules import mdsim as md
    env = md.gromacs()
    gmx = env["gmx"]
    w = tmp_path / "md"
    w.mkdir()
    info = md.build_system(Path(dna_pdb), w, env, is_ligand=False)
    rdir = w / "md_rep1"
    rdir.mkdir()
    (rdir / "md.mdp").write_text(md._mdp_md(0.1))
    import os
    e = {**os.environ, "OMP_NUM_THREADS": "4"}
    # R5: equilibrate first — production straight from EM can start on
    # unsettled waters (step-0 NaN blow-up -> mdrun segfault)
    eq = md.run_equilibration(info, w, env, nvt_ps=2.0, npt_ps=2.0)
    md.run_cmd([gmx, "grompp", "-f", str(rdir / "md.mdp"),
                "-c", eq["eq_gro"], "-t", eq["eq_cpt"],
                "-p", info["top"],
                "-o", str(rdir / "md.tpr"), "-maxwarn", "1"],
               log_file=rdir / "grompp.log", env=e)
    md.run_cmd([gmx, "mdrun", "-deffnm", str(rdir / "md"), "-ntmpi", "1",
                "-ntomp", "4"], log_file=rdir / "md.log", env=e)
    anadir = w / "analysis"
    anadir.mkdir()
    md.run_cmd([gmx, "cluster", "-s", str(rdir / "md.tpr"),
                "-f", str(rdir / "md.xtc"),
                "-method", "gromos", "-cutoff", "1.2",
                "-cl", str(anadir / "clusters_r1"),
                "-clid", str(anadir / "cluster_idx_r1.xvg")],
               log_file=anadir / "cluster.log", stdin="DNA\nDNA\n")
    t, y = md._parse_xvg(anadir / "cluster_idx_r1.xvg")
    yr = np.round(np.asarray(y, dtype=float)).astype(int)
    c = int(np.bincount(yr).argmax())
    idx = np.where(yr == c)[0]
    out = anadir / f"reps_r1_c{c}_mean.pdb"
    grp = sc._dump_group(rdir / "md.tpr", gmx, rdir)
    ok = sc._cluster_mean_pdb(rdir / "md.tpr", rdir / "md.xtc", idx,
                              grp, out)
    assert ok and out.is_file()
    # reference: per-frame mean of the CAs, straight from MDAnalysis
    u = mda.Universe(str(rdir / "md.tpr"), str(rdir / "md.xtc"))
    # compare on the same group _cluster_mean_pdb averaged over:
    # CAs for protein systems, the full nucleic group for NA systems
    ma_sel = {"Protein": "name CA", "DNA": "nucleic",
              "RNA": "nucleic", "System": "all"}.get(grp, "name CA")
    ca = u.atoms.select_atoms(ma_sel)
    # R8: _cluster_mean_pdb applies per-atom minimum image against the
    # cluster's first frame before accumulating (xtc coordinates are
    # unwrapped; an atom diffusing > half a box would otherwise smear
    # the mean across box images). Mirror that exactly.
    u.trajectory[idx[0]]
    ref0 = ca.positions.copy()
    box = np.asarray(u.dimensions[:3], dtype=float)
    acc = np.zeros_like(ca.positions)
    for k in idx:
        u.trajectory[k]
        p = ca.positions.copy()
        d = p - ref0
        d -= box * np.round(d / box)
        acc += ref0 + d
    ref = (acc / len(idx)).copy()
    # produced PDB: atom order == the MA selection order (the writer
    # iterates the selection). Compare in the MDA internal scale: this
    # GROMACS xtc build is read by MDA 10x (A as nm), and the PDB writer
    # emits the internal values verbatim as A, so both sides are in the
    # same units and the factor cancels.
    got = np.array([
        [float(l[30:38]), float(l[38:46]), float(l[46:54])]
        for l in out.read_text().splitlines()
        if l.startswith("ATOM")])
    assert len(got) == len(ca), \
        f"PDB has {len(got)} atoms, selection has {len(ca)}"
    d = float(np.abs(got - ref).max())
    assert d < 5e-4, f"max deviation from cluster mean: {d} nm"
    # and it must differ from the raw first frame of the cluster
    # (otherwise the mean added nothing) — unless the cluster is 1 frame
    if len(idx) > 1:
        first = ref * 0.0
        u.trajectory[idx[0]]
        first[:] = ca.positions
        diff = np.linalg.norm(ref - first)
        assert diff > 1e-4 or len(idx) == 1

def test_mi_kabsch():
    """R8: MD-rep PDBs live in the simulation-box frame while crystal
    PDBs live in the crystallographic frame — the two differ by
    per-atom integer box vectors, which a plain Kabsch averages into a
    muddled fit (190 A on the 1HVI dimer). _mi_kabsch must recover the
    true placement (model->target transform) to < 0.5 A per atom."""
    import numpy as np
    rng = np.random.default_rng(7)
    # a dimer-like point cloud, 200 points
    P = rng.normal(scale=3.0, size=(200, 3))
    P[:, 0] += np.arange(200) * 0.15  # elongate
    # model (rep) frame = crystal frame shifted and per-atom wrapped
    # into the box (the PBC mismatch) — wrap vectors are
    # position-coherent (adjacent atoms share one). The box is slightly
    # larger than the cloud (minimum image unambiguous) and shifted so
    # the core (first 10) shares one box vector while the far end
    # crosses a boundary.
    span = P.max(axis=0) - P.min(axis=0)
    box = span / 0.92
    Qraw = P
    Qw = None
    shift = 0.0
    for sh in np.linspace(0.0, 1.0, 101):
        cand = Qraw + sh * box[0]
        cw = cand % box
        ck = np.round((cand - cw) / box)
        if (ck[:10] == ck[0]).all() and len({tuple(x) for x in ck}) >= 4:
            Qraw, Qw, shift = cand, cw, sh
            break
    assert Qw is not None, "no suitable box placement found"
    # score a MODEL -> TARGET transform: aligned = Qw @ r.T + t,
    # residual modulo the box lattice
    def score(r, t):
        al = Qw @ r.T + t
        diff = al - P
        diff -= box * np.round(diff / box)
        return float(np.sqrt((diff ** 2).sum(1).mean()))
    # ground truth: identity rotation, translation -shift*box[0]
    assert score(np.eye(3), -shift * box[0]) < 0.01, "ground truth broken"
    Rp, tp = sc.kabsch(Qw, P)
    plain = score(Rp, tp)
    Rm, tm = sc._mi_kabsch(Qw, P, box)
    mi_rmsd = score(Rm, tm)
    assert mi_rmsd < 0.5, f"MI RMSD {mi_rmsd}"
    assert mi_rmsd < plain, f"MI ({mi_rmsd}) should beat plain ({plain})"
    # no box (None) -> identical to plain kabsch
    Rn, tn = sc._mi_kabsch(Qw, P, None)
    assert np.allclose(Rn, Rp) and np.allclose(tn, tp)




# ---------------------------------------------------------------- R10: G6 — library fallback

def test_resolve_library_direct_hit(tmp_path):
    from drugagent.modules import screening as sc
    (tmp_path / "dtp.sdf").write_bytes(b"x" * 2_000_000)
    lib, used = sc.resolve_library({"library": "dtp"}, resolve_defaults({"fast": True}),
                                   base=tmp_path)
    assert lib == tmp_path / "dtp.sdf"
    assert used == "dtp"


def test_resolve_library_fallback_when_missing(tmp_path):
    from drugagent.modules import screening as sc
    (tmp_path / "chembl35_small.sdf").write_bytes(b"x" * 2_000_000)
    # dtp missing -> fall back to the local ChEMBL subsample
    lib, used = sc.resolve_library({"library": "dtp"}, resolve_defaults({"fast": True}),
                                   base=tmp_path)
    assert lib.name == "chembl35_small.sdf"
    assert "fallback" in used and "dtp" in used


def test_resolve_library_fallback_when_corrupt(tmp_path):
    from drugagent.modules import screening as sc
    (tmp_path / "dtp.sdf").write_bytes(b"x" * 100)  # 0-byte-ish download
    (tmp_path / "chembl35_small.sdf").write_bytes(b"x" * 2_000_000)
    lib, used = sc.resolve_library({"library": "dtp"}, resolve_defaults({"fast": True}),
                                   base=tmp_path)
    assert lib.name == "chembl35_small.sdf"


def test_resolve_library_custom_path_no_fallback(tmp_path):
    from drugagent.modules import screening as sc
    with pytest.raises(FileNotFoundError):
        sc.resolve_library(
            {"library": "custom",
             "library_path": str(tmp_path / "missing.sdf")},
            resolve_defaults({"fast": True}), base=tmp_path)
