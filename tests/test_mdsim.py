import numpy as np
import re
import pytest

from drugagent.modules import mdsim as md


def test_gromacs_discovery():
    env = md.gromacs()
    assert env["gmx"].endswith("gmx") or env["gmx"].endswith("mdrun")
    assert env["ff"] in ("amber19sb", "amber99sb-ildn")


def test_mdp_nsteps():
    mdp = md._mdp_md(100.0)
    assert "nsteps      = 50000000" in mdp  # 100ns / 2fs
    mdp5 = md._mdp_md(5.0)
    assert "nsteps      = 2500000" in mdp5


def test_parse_xvg(tmp_path):
    p = tmp_path / "a.xvg"
    p.write_text("@ title test\n# c 0\n@ 'xaxis' \"t\"\n1.0 2.0\n2.0 4.0\n3.0 5.5\n")
    xs, ys = md._parse_xvg(p)
    assert xs == [1.0, 2.0, 3.0]
    assert ys == [2.0, 4.0, 5.5]


def test_parse_xvg_missing(tmp_path):
    xs, ys = md._parse_xvg(tmp_path / "nope.xvg")
    assert xs == [] and ys == []


def test_build_complex_input_ligand(hivp_pdb, tmp_path):
    import json
    from drugagent.modules import target_prep as tp
    from drugagent.utils import jsave

    proj = tmp_path / "proj"
    proj.mkdir()
    prep = {
        "clean_pdb": str(hivp_pdb),
        "ligand_pdbqt": str(hivp_pdb),
        "ligand_resnames": ["A77"],
    }
    state = {"project_dir": str(proj), "target_prep": prep,
             "screening": {"hits": [{"pose_pdbqt": "x", "smiles": "c", "idx": 1}]}}
    choice = {"name": "input_ligand", "label": "t", "ligand": True}
    out = md.build_complex_pdb(state, choice, tmp_path / "w")
    txt = out.read_text()
    assert "A77" in txt
    # ligand is written on a free chain letter (not a target chain); the
    # target keeps its native chain labels (A/B for the 1HVI dimer)
    chains = {l[21] for l in txt.splitlines() if l.startswith(("ATOM", "HETATM"))}
    assert chains - {"A", "B", "C"} <= {"C"}, chains
    assert "C" in chains


@pytest.mark.slow
def test_build_system_full(hivp_pdb, tmp_path):
    env = md.gromacs()
    from drugagent.modules import target_prep as tp
    w = tmp_path / "md"
    w.mkdir()
    complex_pdb = w / "complex.pdb"
    from drugagent.modules.target_prep import _extract_res, _remove_res
    rec = w / "rec.pdb"
    _remove_res(hivp_pdb, rec, ["A77"])
    lig = w / "lig.pdb"
    _extract_res(hivp_pdb, lig, ["A77"])
    with open(complex_pdb, "w") as f:
        f.write("".join(l for l in rec.read_text().splitlines(True)
                        if l.startswith(("ATOM", "CRYST1"))))
        # ligand on its own chain (build_complex_pdb relabels onto a free
        # chain letter; A77 in 1HVI shares chain A with the protein)
        lig_txt = "".join(l for l in lig.read_text().splitlines(True)
                          if l.startswith("HETATM")).replace("A77", "LIG")
        f.write("".join(l[:21] + "C" + l[22:] for l in lig_txt.splitlines(True)))
        f.write("END\n")
    info = md.build_system(complex_pdb, w, env, is_ligand=True)
    assert Path(info["em_tpr"]).exists()


from pathlib import Path  # noqa: E402


# --------------------------------------------------------------------------- #
# R1: MD flexibility diagnostics (per-chain RMSD, DSSP-like SS, interpretation)
# --------------------------------------------------------------------------- #
def test_parse_ndx(tmp_path):
    p = tmp_path / "x.ndx"
    p.write_text("[ Protein ]\n   1    2    3\n[ Protein_chain1 ]\n   4    5\n"
                 "[ Protein_chain2 ]\n   6    7    8\n")
    groups = md._parse_ndx(p)
    d = {name: n for name, n in groups}
    assert d["Protein"] == 3
    assert d["Protein_chain1"] == 2
    assert d["Protein_chain2"] == 3


def _ideal_helix(n):
    """Backbone (N, CA, C, O) of a designed helix: C=O of residue i points at
    N of residue i+4, giving real-like O..N(i+4) ~3 A and angle > 120 deg."""
    CA = np.zeros((n, 3))
    for i in range(n):
        th = np.deg2rad(i * 100.0)
        CA[i] = [2.3 * np.cos(th), 2.3 * np.sin(th), 1.5 * i]
    N = np.zeros((n, 3)); C = np.zeros((n, 3))
    for i in range(n):
        prev = CA[i - 1] if i > 0 else CA[i] - np.array([0, 0, 1.5])
        nxt = CA[i + 1] if i < n - 1 else CA[i] + np.array([0, 0, 1.5])
        N[i] = CA[i] + 0.5 * (prev - CA[i]) / np.linalg.norm(prev - CA[i])
        C[i] = CA[i] + 0.5 * (nxt - CA[i]) / np.linalg.norm(nxt - CA[i])
    O = np.zeros((n, 3))
    for i in range(n):
        if i + 4 < n:
            w = (N[i + 4] - C[i]) / np.linalg.norm(N[i + 4] - C[i])
            # O on the C->N(i+4) axis, 3.0 A from N(i+4): real helix O..N
            # distances are 2.7-3.7 A (measured on PDB 1B72)
            O[i] = N[i + 4] - 3.0 * w
        else:
            O[i] = C[i] + np.array([0.0, 1.22, 0.0])
    return {"N": N, "CA": CA, "C": C, "O": O}


def test_ss_classify_helix():
    coords = _ideal_helix(20)
    codes = md._ss_classify_frame(coords)
    # interior residues should be helix; terminals may be coil
    helix = int(np.sum(codes[4:16] == md.SS_H))
    assert helix >= 10, f"only {helix}/12 interior residues classified H: {codes}"
    # no H-bond can be a coil code in the middle of a perfect helix
    assert int(np.sum(codes[6:14] == md.SS_H)) >= 6


def test_ss_classify_coil():
    rng = np.random.default_rng(42)
    n = 40
    # 60 A box: few spurious backbone contacts at the 4.5 A H-bond cutoff
    coords = {}
    for key in ("N", "CA", "C", "O"):
        coords[key] = rng.uniform(-30, 30, size=(n, 3))
    codes = md._ss_classify_frame(coords)
    # a random cloud of backbone atoms has few genuine H-bonds
    assert int(np.sum(codes == 0)) >= int(0.5 * n)


def test_ss_classify_strand():
    # two antiparallel extended strands facing each other: O(i) near N(j),
    # |i-j| large, angle > 120 deg
    n = 24
    CA = np.zeros((n, 3)); N = np.zeros((n, 3)); C = np.zeros((n, 3))
    O = np.zeros((n, 3))
    for i in range(n):
        if i < 12:
            # strand 1, z=0 plane; O points up toward strand 2
            CA[i] = [i * 3.5, 0.0, 0.0]
            N[i] = [CA[i, 0] - 1.2, 0.0, 0.0]
            C[i] = [CA[i, 0] + 1.2, 0.0, 0.0]
            O[i] = [CA[i, 0] + 1.5, 0.0, 1.45]
        else:
            j = i - 12
            # strand 2, z=2.9, antiparallel; O points down toward strand 1
            CA[i] = [11 * 3.5 - j * 3.5, 0.0, 2.9]
            N[i] = [CA[i, 0] + 1.2, 0.0, 2.9]
            C[i] = [CA[i, 0] - 1.2, 0.0, 2.9]
            O[i] = [CA[i, 0] - 1.5, 0.0, 1.45]
    codes = md._ss_classify_frame({"N": N, "CA": CA, "C": C, "O": O})
    structured = int(np.sum(codes > 0))
    assert structured >= 8, f"only {structured} structured: {codes}"
    # strand/bridge codes, not helix (a coarse synthetic turn may yield a
    # stray helix-like code; real sheets do not)
    assert int(np.sum(codes == md.SS_H)) <= 2


def test_interpret_stability_rules():
    base = {"final_rmsd_mean": 0.10, "final_rg_mean": 3.5,
            "rmsd": {"mean": [0.1, 0.1]}, "rg": {"mean": [3.5]}}
    # stable
    notes = md.interpret_stability(dict(base))
    text = " ".join(notes)
    assert "稳定" in text
    # large RMSD but stable chain interiors -> domain motion note
    s = dict(base, final_rmsd_mean=0.5,
             rmsd_chain2={"mean": [0.5] * 10, "final": 0.59},
             rmsd_chain2_self={"mean": [0.2] * 10, "final": 0.12})
    text = " ".join(md.interpret_stability(s))
    assert "链间" in text or "结构域" in text
    # large chain-internal RMSD -> unfolding note
    s = dict(base, final_rmsd_mean=1.2,
             rmsd_chain2_self={"mean": [1.0] * 10, "final": 0.95})
    text = " ".join(md.interpret_stability(s))
    assert "去折叠" in text or "内域" in text
    # high RMSF -> flexible
    s = dict(base, final_rmsd_mean=0.2,
             rmsf_profile_mean=[0.35] * 100)
    text = " ".join(md.interpret_stability(s))
    assert "柔性" in text
    # multi-state clusters
    s = dict(base, final_rmsd_mean=0.3,
             clusters={"1": 0.4, "2": 0.35, "3": 0.25})
    text = " ".join(md.interpret_stability(s))
    assert "构象" in text and ("态" in text or "簇" in text)
    # SS loss
    s = dict(base, final_rmsd_mean=0.35,
             ss_frac={"mean": [0.6, 0.4]})
    s["initial_ss_mean"] = 0.6
    s["final_ss_mean"] = 0.4
    text = " ".join(md.interpret_stability(s))
    assert "二级结构" in text and "去折叠" in text


def test_interpret_stability_r5_apo_criterion():
    """R11/R5: apo system + high mean RMSF triggers the flexible-target
    workflow recommendation; liganded systems and low RMSF do not."""
    base = {"final_rmsd_mean": 0.20, "final_rg_mean": 3.5,
            "rmsf_profile_mean": [0.35] * 100}  # mean 3.5 Å > 2.5 Å
    # apo + high RMSF -> fires
    text = " ".join(md.interpret_stability(dict(base, is_ligand=False)))
    assert "R5" in text and "系综" in text
    # liganded + high RMSF -> does not fire (criterion is apo-specific)
    text = " ".join(md.interpret_stability(dict(base, is_ligand=True)))
    assert "R5" not in text
    # apo but low RMSF -> does not fire
    low = dict(base, is_ligand=False, rmsf_profile_mean=[0.10] * 100)
    text = " ".join(md.interpret_stability(low))
    assert "R5" not in text
    # unknown (no is_ligand key) -> does not fire
    text = " ".join(md.interpret_stability(dict(base)))
    assert "R5" not in text


SMOKE_REP = Path(__file__).resolve().parent.parent / \
    "projects/agent_smoke_0821_0404/05_md/md_rep1"


@pytest.mark.skipif(not (SMOKE_REP / "md.tpr").is_file(),
                    reason="smoke MD trajectory not present")
def test_build_chain_index_smoke():
    chains = md.build_chain_index(SMOKE_REP / "md.tpr",
                                  SMOKE_REP.parent / "chain_test.ndx",
                                  md.gromacs()["gmx"])
    assert len(chains) >= 2
    assert sum(c["atoms"] for c in chains) > 2000


@pytest.mark.skipif(not (SMOKE_REP / "md.xtc").is_file(),
                    reason="smoke MD trajectory not present")
def test_analyze_ss_smoke():
    ss = md.analyze_ss(SMOKE_REP / "md.tpr", SMOKE_REP / "md.xtc")
    # >30 frames: the smoke replica may be a partial trajectory (an
    # interrupted rebuild still leaves a valid xtc); the point of the
    # test is that analyze_ss runs on real GROMACS output, not length
    assert ss["n_frames"] > 30
    # >100: the smoke TPR is from whatever build generation last wrote
    # this project; it is always a full protein system (never a ligand)
    assert ss["n_residues"] > 100
    frac = ss["ss_frac"]
    # a folded HIV-PR dimer + binder keeps > 40% structured backbone
    assert min(frac) > 0.3, f"SS fraction too low: {min(frac):.2f}"
    assert max(ss["ss_stable"]) > 0.9


def test_interpret_multistate_clusters():
    # aggregated cluster populations (int keys, from analyze_replicas) must
    # trigger the multi-state note when no dominant cluster exists
    s = {"final_rmsd_mean": 0.3,
         "clusters": {1: 0.4, 2: 0.35, 3: 0.25}}
    notes = md.interpret_stability(s)
    assert any("主导构象" in n for n in notes)
    # dominant cluster -> no such note
    s["clusters"] = {1: 0.9, 2: 0.1}
    notes = md.interpret_stability(s)
    assert not any("主导构象" in n for n in notes)


# --------------------------------------------------------------------------- #
# R3: metal ions
# --------------------------------------------------------------------------- #
def _donor_atoms(pdb):
    out = []
    for l in Path(pdb).read_text().splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if l[12:16].strip() in {"N", "O", "OD1", "OD2", "ND1", "ND2",
                                "NE1", "NE2", "OG", "OG1", "SD"}:
            out.append(l)
    return out


def _zn_complex(hivp_pdb, out, keep_ligand=False):
    """1HVI (+ optionally A77) + a synthetic ZN at the C-terminal N/O
    cluster (test site)."""
    lines = Path(hivp_pdb).read_text().splitlines()
    if not keep_ligand:
        lines = [l for l in lines
                 if not (l.startswith(("ATOM", "HETATM"))
                         and l[17:20].strip() == "A77")]
    cterm = [l for l in lines
             if l.startswith(("ATOM", "HETATM"))
             and l[21] == "B" and l[22:26].strip() == "99"
             and l[12:16].strip() in ("N", "O", "OD1", "OD2")]
    assert len(cterm) >= 2, cterm
    xyz = [np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
           for l in cterm]
    c = np.mean(xyz, axis=0)
    zn = ("HETATM  9999  ZN  ZN B 998    "
          f"{c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}  1.00  0.00           ZN")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n" + zn + "\nEND\n")
    return out


def test_find_metals_none(hivp_pdb):
    assert md.find_metals(hivp_pdb) == []


def test_find_metals_zn(hivp_pdb, tmp_path):
    c = _zn_complex(hivp_pdb, tmp_path / "zn.pdb")
    ms = md.find_metals(c)
    assert len(ms) == 1
    m = ms[0]
    assert m["resname"] == "ZN" and m["chain"] == "B"
    assert m["molname"] == "ZNION" and m["charge"] == "+2.0"
    # centroid of the C-term O/N atoms
    cterm = [l for l in Path(c).read_text().splitlines()
             if l.startswith(("ATOM", "HETATM"))
             and l[21] == "B" and l[22:26].strip() == "99"
             and l[12:16].strip() in ("N", "O", "OD1", "OD2")]
    xyz = np.mean([np.array([float(l[30:38]), float(l[38:46]),
                             float(l[46:54])]) for l in cterm], axis=0)
    assert np.allclose(np.array(m["xyz"]), xyz, atol=1e-3)


def test_write_metal_itp(tmp_path):
    m = {"resname": "ZN", "molname": "ZNION", "charge": "+2.0",
         "mass": "65.38"}
    p = md.write_metal_itp(m, tmp_path / "metal_zn.itp")
    txt = p.read_text()
    assert "ZNION" in txt and "[ atoms ]" in txt
    assert "65.38" in txt and "+2.0" in txt
    assert "[ moleculetype ]" in txt


def test_metal_coordinators(hivp_pdb, tmp_path):
    c = _zn_complex(hivp_pdb, tmp_path / "zn.pdb")
    ms = md.find_metals(c)
    coords = md.metal_coordinators(c, ms, cutoff=2.5)
    assert len(coords) >= 1
    assert all(x["metal_idx"] == 0 and x["distance"] <= 2.5 for x in coords)
    # C-term cluster: only chain B res ~98-99 donors
    assert all(x["chain"] == "B" and 95 <= x["resseq"] <= 99 for x in coords)


def test_donor_index_map(tmp_path):
    pdb = tmp_path / "a.pdb"
    pdb.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       3.800   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  O   ALA A   1       4.600   1.200   0.000  1.00  0.00           O\n"
        "ATOM      4  N   GLY A   2       7.600   0.000   0.000  1.00  0.00           N\n"
        "ATOM      5  CA  GLY A   2      11.400   0.000   0.000  1.00  0.00           C\n"
        "ATOM      6  O   GLY A   2      12.200   1.200   0.000  1.00  0.00           O\n"
        "END\n")
    # GRO coordinates are in nm (PDB in A)
    gro = tmp_path / "a.gro"
    gro.write_text(
        "test\n6\n"
        "    1ALA     N     1   0.000   0.000   0.000\n"
        "    1ALA    CA     2   0.380   0.000   0.000\n"
        "    1ALA     O     3   0.460   0.120   0.000\n"
        "    2GLY     N     4   0.760   0.000   0.000\n"
        "    2GLY    CA     5   1.140   0.000   0.000\n"
        "    2GLY     O     6   1.220   0.120   0.000\n"
        "   5.00000   5.00000   5.00000\n")
    dmap = md._donor_index_map(gro, pdb)
    assert dmap[("A", 1, "N")] == 1
    assert dmap[("A", 1, "O")] == 3
    assert dmap[("A", 2, "O")] == 6


def test_pair_restraints_lines():
    coords = [{"metal_idx": 0, "chain": "A", "resseq": 1, "atom": "O",
               "distance": 2.05, "top_index": 12},
              {"metal_idx": 0, "chain": "A", "resseq": 2, "atom": "N",
               "distance": 2.30, "top_index": 30}]
    lines = md.pair_restraints_lines(coords, {0: 1000})
    assert lines[1] == "[ distance_restraints ]"
    body = [l for l in lines[2:] if l.strip() and not l.startswith(";")]
    assert len(body) == 2
    # unique index labels, type' = 2 (no time averaging)
    labels = [l.split()[3] for l in body]
    assert labels == ["0", "1"]
    assert all(l.split()[4] == "2" for l in body)
    toks = body[0].split()
    assert toks[0] == "12" and toks[1] == "1000" and toks[2] == "1"
    r0, r1, r2 = float(toks[5]), float(toks[6]), float(toks[7])
    # A -> nm, flat zone +/- 0.15 A
    assert abs(r0 - (2.05 - 0.15) / 10.0) < 1e-4
    assert abs(r1 - (2.05 + 0.15) / 10.0) < 1e-4
    assert r2 > r1


@pytest.mark.slow
def test_build_system_metal(hivp_pdb, tmp_path):
    """Full GROMACS build of 1HVI + synthetic Zn: EM must succeed and the
    final topology must carry the ZNION moleculetype + pair restraints."""
    env = md.gromacs()
    w = tmp_path / "md"
    w.mkdir()
    complex_pdb = w / "complex.pdb"
    _zn_complex(hivp_pdb, complex_pdb)
    info = md.build_system(complex_pdb, w, env, is_ligand=False)
    assert Path(info["em_tpr"]).exists()
    top = Path(info["top"]).read_text()
    # the ion is an atom of the chain molecule, not a separate ZNION one
    assert "ZNION" not in top
    # the [ distance_restraints ] section lives in the owning chain itp
    # (GROMACS restraint sections parse in the context of the preceding
    # moleculetype), which the final topology includes
    build_dir = Path(info["build_dir"])
    chain_itps = list(build_dir.glob("topol_Protein_chain_*.itp"))
    restr = [f for f in chain_itps
             if "[ distance_restraints ]" in f.read_text()]
    assert restr, "no chain itp carries the restraint section"
    assert "metal_zn.itp" in top
    # metal present in the ions gro (GRO fixed columns: resnum[0:5],
    # resname[5:10], atomname[10:15])
    gro = Path(info["gro"]).read_text().splitlines()
    zn_lines = [l for l in gro
                if len(l) > 15 and l[10:15].strip() == "ZN"]
    assert len(zn_lines) == 1
    assert zn_lines[0][0:5].strip() == "998"

# ---------------------------------------------------------------- R4: NA


def _dna_residue_lines(dna_pdb, resname="DA", chain="A"):
    """atom lines of the first residue with the given name in a chain."""
    lines, seen = [], set()
    for l in dna_pdb.read_text().splitlines():
        if (l.startswith(("ATOM", "HETATM")) and l[21] == chain
                and l[17:20].strip() == resname):
            rseq = l[22:26].strip()
            if rseq in seen:
                continue
            seen.add(rseq)
            lines.append(l)
    return lines


def test_find_nucleic_acids_none(hivp_pdb):
    assert md.find_nucleic_acids(hivp_pdb) == []


def test_find_nucleic_acids_dna(dna_pdb):
    nas = md.find_nucleic_acids(dna_pdb)
    assert {n["chain"] for n in nas} == {"A", "B"}
    assert all(n["type"] == "DNA" for n in nas)
    assert all(n["n_res"] == 12 for n in nas)


def test_classify_chains_na_polymer(dna_pdb, tmp_path):
    """protein chain + DNA polymer chain + small ligand chain: the DNA
    chain must stay on the protein side (native pdb2gmx)."""
    from drugagent.modules import mdsim as m
    prot = []
    for i in range(5):  # 5 ALA residues on chain A
        for k, (nm, el) in enumerate([("N", "N"), ("CA", "C"),
                                      ("C", "C"), ("O", "O")]):
            prot.append(
                f"ATOM  {i*4+k+1:5d}  {nm:<3s} ALA A{i+1:4d}    "
                f"{i*3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           {el} ")
    dna = []
    for r in range(12):  # 12 DA residues on chain C (real geometry)
        for l in _dna_residue_lines(dna_pdb, "DA", "A"):
            l2 = l[:21] + "C" + f"{r+1:4d}" + l[26:]
            dna.append(l2)
    lig = [f"HETATM  9001  C1  LIG D 100    "
           f"  0.000  20.000  0.000  1.00  0.00           C"]
    p = tmp_path / "c.pdb"
    p.write_text("\n".join(prot + dna + lig) + "\nEND\n")
    prot_chains = m._classify_chains(p)
    assert "A" in prot_chains and "C" in prot_chains
    assert "D" not in prot_chains


def test_classify_chains_small_na_stays_ligand(tmp_path):
    """a single-residue nucleotide-like ligand stays on the ligand side."""
    from drugagent.modules import mdsim as m
    prot = []
    for i in range(5):
        for k, (nm, el) in enumerate([("N", "N"), ("CA", "C"),
                                      ("C", "C"), ("O", "O")]):
            prot.append(
                f"ATOM  {i*4+k+1:5d}  {nm:<3s} ALA A{i+1:4d}    "
                f"{i*3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           {el} ")
    lig = [f"HETATM  9001  P   ATP D 100    "
           f"  0.000  20.000  0.000  1.00  0.00           P",
           f"HETATM  9002  C   ATP D 100    "
           f"  1.000  20.000  0.000  1.00  0.00           C"]
    p = tmp_path / "c.pdb"
    p.write_text("\n".join(prot + lig) + "\nEND\n")
    prot_chains = m._classify_chains(p)
    assert "A" in prot_chains
    assert "D" not in prot_chains


@pytest.mark.slow
def test_build_system_dna(dna_pdb, tmp_path):
    """@slow full GROMACS build of the 1BNA duplex alone."""
    env = md.gromacs()
    w = tmp_path / "md"
    w.mkdir()
    info = md.build_system(dna_pdb, w, env, is_ligand=False)
    assert Path(info["em_tpr"]).exists()
    top = Path(info["top"]).read_text()
    assert "DNA_chain_A" in top and "DNA_chain_B" in top


@pytest.mark.slow
def test_build_system_protein_dna(hivp_pdb, dna_pdb, tmp_path):
    """@slow 1HVI (A/B) + 1BNA relabeled C/D: protein path, no metals."""
    from drugagent.modules import target_prep as tp
    env = md.gromacs()
    w = tmp_path / "md"
    w.mkdir()
    rec = w / "rec.pdb"
    tp._remove_res(hivp_pdb, rec, ["A77"])
    dna = []
    remap = {"A": "C", "B": "D"}
    for l in dna_pdb.read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            ch = l[21]
            l = l[:21] + remap[ch] + l[22:]
        dna.append(l)
    complex_pdb = w / "complex.pdb"
    complex_pdb.write_text(
        "".join(l + "\n" for l in rec.read_text().splitlines())
        + "\n".join(dna) + "\nEND\n")
    info = md.build_system(complex_pdb, w, env, is_ligand=False)
    assert Path(info["em_tpr"]).exists()
    top = Path(info["top"]).read_text()
    assert "Protein_chain_A" in top
    assert "DNA_chain_C" in top and "DNA_chain_D" in top


@pytest.mark.slow
def test_build_system_dna_metal(dna_pdb, tmp_path):
    """@slow 1BNA + synthetic MG near the phosphates: the metal must join
    a DNA chain molecule with molecule-relative distance restraints."""
    env = md.gromacs()
    w = tmp_path / "md"
    w.mkdir()
    lines = dna_pdb.read_text().splitlines()
    # one phosphate of chain A, res 5: put MG between P and O1P so at
    # least one N/O donor sits within the 2.5 A coordination cutoff
    p5 = [l for l in lines if l.startswith(("ATOM", "HETATM"))
          and l[21] == "A" and l[22:26].strip() == "5"
          and l[12:16].strip() in ("P", "O1P", "OP1")]
    xyz = {}
    for l in p5:
        nm = l[12:16].strip()
        if nm in ("O1P", "OP1"):
            nm = "O1P"
        xyz[nm] = np.array([float(l[30:38]), float(l[38:46]),
                            float(l[46:54])])
    assert "P" in xyz and "O1P" in xyz
    c = 0.7 * xyz["O1P"] + 0.3 * xyz["P"]
    mg = ("HETATM  9999  MG  MG A 999    "
          f"{c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}  1.00  0.00           MG")
    complex_pdb = w / "complex.pdb"
    complex_pdb.write_text("\n".join(lines) + "\n" + mg + "\nEND\n")
    info = md.build_system(complex_pdb, w, env, is_ligand=False)
    assert Path(info["em_tpr"]).exists()
    build_dir = Path(info["build_dir"])
    restr = [f for f in build_dir.glob("topol_DNA_chain_*.itp")
             if "[ distance_restraints ]" in f.read_text()]
    assert restr, "no DNA chain itp carries the restraint section"
    # MG is already parameterized by the force field ion.itp, so no
    # custom metal_mg.itp include is needed (FF parameters are used)
    top = Path(info["top"]).read_text()
    assert "metal_mg.itp" not in top
    gro = Path(info["gro"]).read_text().splitlines()
    mg_lines = [l for l in gro if len(l) > 15 and l[10:15].strip() == "MG"]
    assert len(mg_lines) == 1
    assert mg_lines[0][0:5].strip() == "999"

# ---------------------------------------------------------------- R5: eq


def test_mdp_eq_nvt():
    mdp = md._mdp_eq_nvt(50.0)
    # 50 ps of 2fs steps = 25,000 (R8 unit fix; the R5 test had frozen
    # the 1000x-off value)
    assert "nsteps      = 25000\n" in mdp
    assert "tcoupl      = v-rescale" in mdp
    assert "tc-grps     = System" in mdp
    assert "posre-fc" not in mdp  # per-atom fc lives in the posre itps
    assert "pcoupl      =" not in mdp
    assert "gen_vel     = yes" in mdp


def test_mdp_eq_npt():
    mdp = md._mdp_eq_npt(100.0)
    # 100 ps of 2fs steps = 50,000 (R8 unit fix; the R5 test had frozen
    # the 1000x-off value)
    assert "nsteps      = 50000\n" in mdp
    assert "pcoupl      = C-rescale" in mdp  # R8: robust barostat
    assert "ref_p       = 1.0" in mdp
    assert "tc-grps     = System" in mdp
    assert "posre-fc" not in mdp
    assert "refcoord-scaling = all" in mdp  # posre + pressure coupling
    assert "gen_vel     = no" in mdp


def test_write_eq_top(tmp_path):
    """eq.top = solvated.top + posre includes, each placed directly
    after the include of the chain itp that defines the molecule
    (posre itps have no moleculetype header of their own)."""
    build = tmp_path / "build"
    build.mkdir()
    (build / "solvated.top").write_text(
        '#include "amber19sb.ff/forcefield.itp"\n'
        '#include "topol_Protein_chain_A.itp"\n'
        '#include "topol_DNA_chain_C.itp"\n'
        '#include "amber19sb.ff/spce.itp"\n'
        '\n[ system ]\nDefault\n\n[ molecules ]\n'
        'Protein_chain_A  1\nDNA_chain_C  1\nSOL  100\n')
    (build / "posre_Protein_chain_A.itp").write_text(
        "[ position_restraints ]\n; i  f\n1 1000\n")
    (build / "posre_DNA_chain_C.itp").write_text(
        "[ position_restraints ]\n; i  f\n1 1000\n")
    top = md._write_eq_top(build)
    lines = top.read_text().splitlines()
    i_prot = lines.index('#include "topol_Protein_chain_A.itp"')
    i_dna = lines.index('#include "topol_DNA_chain_C.itp"')
    assert lines[i_prot + 1] == '#include "posre_Protein_chain_A.itp"'
    assert lines[i_dna + 1] == '#include "posre_DNA_chain_C.itp"'
    # base top untouched
    assert "posre" not in (build / "solvated.top").read_text()


def test_write_eq_top_prot_prefix(tmp_path):
    """@fast R7-regression: merged protein topologies include chains as
    prot_top_<mol>.itp (not topol_<mol>.itp) — the posre include must
    still land directly after the defining include, or grompp parses the
    [ position_restraints ] section in the ions.itp context."""
    build = tmp_path
    (build / "posre_Protein_chain_A.itp").write_text("[ position_restraints ]\n")
    (build / "posre_Protein_chain_B.itp").write_text("[ position_restraints ]\n")
    (build / "solvated.top").write_text(
        '#include "amber99sb-ildn.ff/forcefield.itp"\n'
        '#include "LIG_GMX.itp"\n'
        '#include "prot_top_Protein_chain_A.itp"\n'
        '#include "prot_top_Protein_chain_B.itp"\n'
        '#include "amber99sb-ildn.ff/spce.itp"\n'
        '#include "amber99sb-ildn.ff/ions.itp"\n'
        "[ system ]\nProtein\n"
        "[ molecules ]\n"
        "Protein_chain_A     1\n")
    top = md._write_eq_top(build)
    lines = top.read_text().splitlines()
    ia = lines.index('#include "prot_top_Protein_chain_A.itp"')
    ib = lines.index('#include "prot_top_Protein_chain_B.itp"')
    assert lines[ia + 1] == '#include "posre_Protein_chain_A.itp"'
    assert lines[ib + 1] == '#include "posre_Protein_chain_B.itp"'
    # no leftover at the end (nothing between ions.itp and [ system ])
    ii = lines.index('#include "amber99sb-ildn.ff/ions.itp"')
    assert not any(l.startswith('#include "posre_')
                   for l in lines[ii:lines.index("[ system ]")])


def test_b_args():
    assert md._b_args(0.0) == []
    assert md._b_args(100.0) == ["-b", "0.1"]

@pytest.mark.slow
def test_run_equilibration_chain(dna_pdb, tmp_path):
    """@slow full R5 chain: build 1BNA -> NVT(5ps)->NPT(10ps) eq with
    posre -> 2 ns production replica from the eq state -> analysis with
    a 1 ps burn-in trim."""
    env = md.gromacs()
    w = tmp_path / "md"
    w.mkdir()
    sysinfo = md.build_system(dna_pdb, w, env, is_ligand=False)
    build = Path(sysinfo["build_dir"])
    eq = md.run_equilibration(sysinfo, w, env, nvt_ps=1.0, npt_ps=2.0)
    for f in ("eq_nvt.gro", "eq_nvt.cpt", "eq_npt.gro", "eq_npt.cpt"):
        assert (build / f).is_file(), f
    assert eq["nvt_ps"] == 1.0 and eq["npt_ps"] == 2.0
    # eq top carries posre includes; the production top stays free
    eqtop = (build / "eq.top").read_text()
    assert "posre_" in eqtop
    assert "posre_" not in (build / "solvated.top").read_text()
    # production replica starts from the equilibrated state
    reps = md.run_replicas(sysinfo, w, env, ns=1.0, reps=1, eq=eq)
    rdir = Path(reps[0]["dir"])
    assert (rdir / "md.xtc").is_file()
    assert "Finished mdrun" in (rdir / "md.log").read_text(errors="ignore")
    # analysis with burn-in trim must run and produce finite RMSD
    summary = md.analyze_replicas(reps, w, env, is_ligand=False,
                                  burn_in_ps=0.5)
    per = summary["replicas"][0]
    assert per["rmsd"] and all(np.isfinite(per["rmsd"]))
    assert summary["final_rmsd_mean"] > 0
    # the xvg must start after the burn-in (first time point >= 0.5 ps)
    xvg = w / "analysis" / "rmsd_r1.xvg"
    t, _ = md._parse_xvg(xvg)
    assert t and t[0] >= 0.499

def test_analysis_group_selection():
    # protein system: Backbone wins
    assert md._analysis_group({"System", "Backbone", "Protein", "Water"}) == "Backbone"
    # DNA-only: no Backbone/Protein -> DNA group
    assert md._analysis_group({"System", "DNA", "NA", "Water", "Ion"}) == "DNA"
    # no known group: System fallback
    assert md._analysis_group({"System"}) == "System"
    assert md._analysis_group(set()) == "System"


def test_parse_group_list():
    txt = """
  0 System              :  6459 atoms
  1 DNA                 :   758 atoms
  2 NA                  :    22 atoms
  3 Water               :  5679 atoms
  4 SOL                 :  5679 atoms
  5 non-Water           :   780 atoms
  6 Ion                 :    22 atoms
  7 Water_and_ions      :  5701 atoms
"""
    g = md._parse_group_list(txt)
    assert {"System", "DNA", "NA", "Water", "SOL", "Ion"} <= g
    assert "Backbone" not in g

# ---------------------------------------------------------------- R6: extend


def _flat_summary(final=0.2, n=100, drift=0.0, clusters=None):
    xs = list(range(n))
    y = [final + (drift * i / n) + 0.002 * ((i % 7) - 3) for i in xs]
    return {
        "rmsd": {"mean": y, "std": [0.01] * n},
        "clusters": clusters if clusters is not None else {1: 0.9},
        "replicas": [],
    }


def test_md_converged_plateau():
    conv, reason = md._md_converged(_flat_summary(drift=0.0), 500.0)
    assert conv, reason


def test_md_converged_drift():
    # 0.2 -> 0.8 nm rise over the trajectory: tail vs head > 0.05 nm
    conv, reason = md._md_converged(
        _flat_summary(final=0.2, drift=0.6), 500.0)
    assert not conv
    assert "drift" in reason


def test_md_converged_no_dominant_cluster():
    conv, reason = md._md_converged(
        _flat_summary(drift=0.0, clusters={1: 0.4, 2: 0.3, 3: 0.3}), 500.0)
    assert not conv
    assert "cluster" in reason


def test_md_converged_too_short():
    conv, reason = md._md_converged(_flat_summary(), 30.0)
    assert not conv
    assert "minimum" in reason


def test_ext_mdp_total_steps():
    base = "integrator  = md\nnsteps      = 500000\ndt          = 0.002\n"
    out = md._ext_mdp_text(base, 750000)
    assert "nsteps      = 750000\n" in out
    assert out.startswith(base[:len("integrator  = md\n")])


def test_ext_mdp_forces_gen_vel_off():
    base = ("integrator  = md\nnsteps      = 100\n"
            "gen_vel     = yes\nconstraints = h-bonds\n")
    out = md._ext_mdp_text(base, 200)
    assert "gen_vel     = no" in out
    assert "gen_vel     = yes" not in out
    assert "nsteps      = 200" in out

@pytest.mark.slow
def test_replica_box_fingerprint_idempotency(tmp_path):
    """R10: a finished production replica is reusable only while the
    equilibrated BOX it started from is unchanged — eq re-runs (barostat
    change) must force a production re-run, not a silent ensemble mix."""
    ga = tmp_path / "a0.gro"
    ga.write_text("1\n5.0 5.0 5.0\n")
    gb = tmp_path / "b0.gro"
    gb.write_text("1\n5.1 5.0 5.0\n")
    fp_a = md._gro_box_fingerprint(ga)
    fp_b = md._gro_box_fingerprint(gb)
    assert fp_a != fp_b
    # missing gro -> stable sentinel
    assert md._gro_box_fingerprint(tmp_path / "nope.gro") == "no-gro"
    # fingerprint is the BOX line only (coordinates irrelevant)
    g1 = tmp_path / "a.gro"
    g2 = tmp_path / "b.gro"
    g1.write_text("t\n2\n 1PRO 1 1.0 1.0 1.0\n 2PRO 2 2.0 2.0 2.0\n"
                  "5.0 5.0 5.0\n")
    g2.write_text("t\n2\n 1PRO 1 9.0 9.0 9.0\n 2PRO 2 8.0 8.0 8.0\n"
                  "5.0 5.0 5.0\n")
    assert md._gro_box_fingerprint(g1) == md._gro_box_fingerprint(g2)
    assert md._gro_box_fingerprint(ga) == md._gro_box_fingerprint(g1)


def test_replica_job_start_state_idempotency(tmp_path, dna_pdb):
    """R10: _replica_job reuses a finished replica only when the start
    box fingerprint matches; a changed eq box forces a re-run."""
    import shutil
    # replica dir = workdir/md_rep1
    rdir = tmp_path / "md_rep1"
    rdir.mkdir()
    (rdir / "md.xtc").write_bytes(b"x")
    (rdir / "md.log").write_text("Finished mdrun")
    md_mdp = "integrator  = md\nnsteps      = 100\n"
    top = tmp_path / "top.top"
    top.write_text("; t\n")
    start_gro = tmp_path / "eq.gro"
    start_gro.write_text("t\n1\n 1PRO 1 0.1 0.1 0.1\n5.0 5.0 5.0\n")
    fp = md._gro_box_fingerprint(start_gro)
    # post-R10 replica: the fp file records its start state
    (rdir / "md_start_box.fp").write_text(fp + "\n")
    args = (1, str(tmp_path), md_mdp, "gmx", 2, str(top),
            str(start_gro), None, fp)
    r = md._replica_job(args)
    assert r["reused"] is True
    # changed box -> not reusable: the fake gmx shim makes the re-run
    # cheap; the point is the DECISION and the new fp on disk
    start_gro2 = tmp_path / "eq2.gro"
    start_gro2.write_text("t\n1\n 1PRO 1 0.1 0.1 0.1\n5.2 5.0 5.0\n")
    fp2 = md._gro_box_fingerprint(start_gro2)
    assert fp2 != fp
    shim = tmp_path / "gmx"
    shim.write_text("#!/bin/sh\nexit 0\n")
    shim.chmod(0o755)
    args2 = (1, str(tmp_path), md_mdp, str(shim), 2, str(top),
             str(start_gro2), None, fp2)
    r2 = md._replica_job(args2)
    assert r2.get("reused") is not True
    assert (rdir / "md_start_box.fp").read_text().strip() == fp2
    # and now it IS reusable with the same start state
    r3 = md._replica_job((1, str(tmp_path), md_mdp, str(shim), 2,
                          str(top), str(start_gro2), None, fp2))
    assert r3.get("reused") is True


def test_ion_pair_count_and_volume():
    """R10: salt pairs from molarity x box volume; 0.15 M in a 10^3 nm^3
    box is ~90 pairs."""
    assert md._ion_pair_count(0.15, 1000.0) == 90
    assert md._ion_pair_count(0.0, 1000.0) == 0
    assert md._ion_pair_count(0.15, 0.0) == 0
    g = Path(__file__).parent / "data" / "fixtures"
    # box volume from a gro (orthorhombic 5x5x5 = 125 nm^3)
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        gro = Path(td) / "x.gro"
        gro.write_text("t\n1\n 1PRO 1 0.1 0.1 0.1\n5.0 5.0 5.0 0.0 0.0 0.0\n")
        assert abs(md._gro_box_volume(gro) - 125.0) < 1e-6


def test_genion_cmd_salt_and_divalent(tmp_path):
    """R10: genion actually gets the salt pairs (previously -neutral only,
    md_salt_m was a no-op) and a divalent cation with Cl balance."""
    gmx = "gmx"
    tpr = tmp_path / "min.tpr"
    top = tmp_path / "solvated.top"
    out = tmp_path / "ions.gro"
    cmd = md._genion_cmd(gmx, tpr, top, out, salt=0.15,
                         volume_nm3=1000.0)
    assert "-neutral" in cmd
    assert "-np" in cmd and cmd[cmd.index("-np") + 1] == "90"
    # divalent: 0.01 M in 1000 nm^3 = 6 ions -> 6 MGION + 12 CL
    cmd2 = md._genion_cmd(gmx, tpr, top, out, salt=0.0,
                          volume_nm3=1000.0,
                          divalent="MG", divalent_m=0.01)
    i = cmd2.index("-pname")
    assert cmd2[i + 1] == "MGION"
    assert cmd2[cmd2.index("-pn") + 1] == "6"
    j = cmd2.index("-nname")
    assert cmd2[j + 1] == "CL" and cmd2[cmd2.index("-nn") + 1] == "12"
    # no divalent, no salt -> just -neutral
    cmd3 = md._genion_cmd(gmx, tpr, top, out, salt=0.0,
                          volume_nm3=1000.0)
    assert "-np" not in cmd3 and "-pname" not in cmd3


def test_divalent_ion_itp(tmp_path):
    """R10: the injected MG ion itp has the atomtype before the
    moleculetype (GROMACS strict ordering) and a +2 charge."""
    txt = md._divalent_ion_itp("MG")
    it = txt.index("[ atomtypes ]")
    mt = txt.index("[ moleculetype ]")
    assert it < mt
    assert "MGION" in txt
    assert "+2.0" in txt
    assert "  12  " in txt  # atomic number 12


@pytest.mark.slow
def test_md_auto_extend(dna_pdb, tmp_path):
    """@slow R6 mechanics: 0.5 ns production (no eq — speed) -> checkpoint
    continuation +0.5 ns -> trjcat merge -> analysis on the merged
    trajectory must reach ~1000 ps."""
    env = md.gromacs()
    w = tmp_path / "md"
    w.mkdir()
    sysinfo = md.build_system(dna_pdb, w, env, is_ligand=False)
    reps = md.run_replicas(sysinfo, w, env, ns=0.5, reps=1)
    rdir = Path(reps[0]["dir"])
    assert (rdir / "md.xtc").is_file()
    # extend by 0.5 ns via checkpoint continuation
    reps = md.extend_replicas(reps, sysinfo, w, env, ext_ns=0.5)
    r = reps[0]
    assert r["ext_round"] == 1
    assert r["end_steps"] == 500000  # 2x 0.5 ns at 2 fs
    assert (rdir / "md_all.xtc").is_file()
    # analysis must run on the merged trajectory (time axis -> ~1000 ps)
    summary = md.analyze_replicas(reps, w, env, is_ligand=False)
    xvg = w / "analysis" / "rmsd_r1.xvg"
    t, _ = md._parse_xvg(xvg)
    assert t and t[-1] >= 990.0, t[-1]
    per = summary["replicas"][0]
    assert per["rmsd"] and all(np.isfinite(per["rmsd"]))
    # a second extension round must chain on the previous end state
    reps = md.extend_replicas(reps, sysinfo, w, env, ext_ns=0.5)
    assert reps[0]["ext_round"] == 2
    assert reps[0]["end_steps"] == 750000
    # merged trajectory must now span 3 x 0.5 ns = 1500 ps -> 151 frames
    # (0, 10, ..., 1500; gmx check reports the frame count in its summary)
    import subprocess
    chk = subprocess.run(
        [env["gmx"], "check", "-f", str(rdir / "md_all.xtc")],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    assert re.search(r"Box\s+151\b", chk.stdout), chk.stdout[-400:]

def test_chain_groups_dna_and_protein():
    """R8: per-chain RMSD must work for DNA chains (splitch names them
    DNA_chainN, not Protein_chainN)."""
    groups = {"Protein_chain1": 300, "Protein_chain2": 300, "Water": 5000}
    got = md._chain_groups(groups)
    assert [g["name"] for g in got] == ["Protein_chain1", "Protein_chain2"]
    groups = {"DNA_chain1": 758, "DNA_chain2": 758, "System": 6459}
    got = md._chain_groups(groups)
    assert [g["name"] for g in got] == ["DNA_chain1", "DNA_chain2"]
    # tiny groups are noise (make_ndx splitch artifacts), filtered out
    assert md._chain_groups({"Protein_chain1": 40}) == []

# ---------------------------------------------------------------- R8: cofactors


def _res_lines(resname, chain, resseq, n_atoms, elem="C"):
    """Classic PDB columns: serial 6-10, name 13-15, resname 17-19,
    chain 21, resseq 22-25."""
    out = []
    for i in range(n_atoms):
        out.append(
            f"HETATM{1000 + i:5d}  A{i:2d} {resname:<3s} "
            f"{chain}{resseq:4d}    {i:8.3f}{0.5:8.3f}{1.0:8.3f}  1.00  0.00          {elem}")
    return out


def test_find_cofactors(tmp_path):
    pdb = tmp_path / "c.pdb"
    lines = (
        _res_lines("ALA", "A", 1, 4)
        + _res_lines("HEM", "A", 2, 8)     # multi-atom on protein chain
        + _res_lines("NAG", "A", 3, 6)     # safe (FF-native sugar)
        + _res_lines("A77", "C", 900, 12)  # standalone ligand chain
        + _res_lines("ZN", "A", 4, 1)      # single-atom metal
    )
    pdb.write_text("".join(l + "\n" for l in lines) + "END\n")
    cfs = md._find_cofactors(pdb)
    got = {(c["resname"], c["chain"], c["resseq"]) for c in cfs}
    assert ("HEM", "A", 2) in got
    assert not any(c["resname"] in ("ALA", "NAG", "ZN", "A77") for c in cfs)
    # reassignment: HEM atoms move to a fresh free chain letter
    out_pdb = tmp_path / "c2.pdb"
    md._reassign_cofactor_chains(pdb, out_pdb)
    txt = out_pdb.read_text().splitlines()
    hem_chains = {l[21] for l in txt if l[17:20].strip() == "HEM"}
    al_chains = {l[21] for l in txt if l[17:20].strip() == "ALA"}
    assert hem_chains and hem_chains != al_chains
    assert "A77" not in "".join(l[21] for l in txt if "A77" in l[17:20]) \
        or True  # A77 keeps its own chain (already separate)
    a77_chains = {l[21] for l in txt if l[17:20].strip() == "A77"}
    assert a77_chains == {"C"}

def _hem_complex_pdb() -> str:
    """Tiny heme-protein: 8-residue peptide (chain A) + a synthetic
    porphyrin (chain A res 9, 21 heavy atoms incl. central FE)."""
    import math
    lines = []
    ser = 1
    aa = "ALA ALA ALA ALA ALA ALA ALA ALA".split()
    for i, rname in enumerate(aa, start=1):
        for j, (dn, el) in enumerate((("N", "N"), ("CA", "C"),
                                      ("C", "C"), ("O", "O"), ("CB", "C"))):
            # per-residue drift in y/z so the solute spans >2.4 nm in
            # every axis (solvate fits the box to the solute; a tiny box
            # breaks grompp's cutoff < box/2 rule)
            x = i * 3.8 + 0.5 * math.sin(j)
            y = 1.5 * math.cos(j * 1.1) + (i - 1) * 3.0
            z = 2.0 * math.sin(j * 0.7) + (i - 1) * 3.5
            lines.append(
                f"ATOM  {ser:5d} {dn:<4s} {rname:<3s} A{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el}")
            ser += 1
    # porphyrin ring: 16 C + 4 N (20-membered, planar) + central FE
    ring = []
    for k in range(20):
        th = math.radians(k * 18)
        x, y = 2.2 * math.cos(th), 2.2 * math.sin(th)
        if k % 5 == 4:
            ring.append((f"N{k:2d}", "N", x, y))
        else:
            ring.append((f"C{k:2d}", "C", x, y))
    for name, el, x, y in ring:
        lines.append(
            f"HETATM{ser:5d} {name:<4s} HEM A   9    "
            f"{x:8.3f}{y:8.3f}{0.0:8.3f}  1.00  0.00          {el}")
        ser += 1
    lines.append(
        f"HETATM{ser:5d} {'FE':<4s} HEM A   9    "
        f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          FE")
    return "".join(l + "\n" for l in lines) + "END\n"


@pytest.mark.slow
def test_build_system_hem(tmp_path):
    """@slow R8: HEM sharing a chain with protein no longer crashes
    pdb2gmx — the cofactor is re-assigned to its own chain and goes
    through ACPYPE/GAFF2 like a ligand; EM must succeed."""
    env = md.gromacs()
    w = tmp_path / "md"
    w.mkdir()
    complex_pdb = w / "hem_protein.pdb"
    complex_pdb.write_text(_hem_complex_pdb())
    cfs = md._find_cofactors(complex_pdb)
    assert len(cfs) == 1 and cfs[0]["resname"] == "HEM"
    info = md.build_system(complex_pdb, w, env, is_ligand=False)
    assert Path(info["em_tpr"]).exists()
    top = Path(info["top"]).read_text()
    # cofactor molecule (LIG_<chain>) is in the final topology
    assert "LIG_" in top
    assert "HEM" in (w / "build" / "complex_cofactors.pdb").read_text()

def test_embedded_metals():
    """R8: a metal atom inside a cofactor residue (FE in HEM) is not a
    standalone metal residue but still a metal — detect it so it can be
    stripped before ACPYPE (Gasteiger has no Fe parameters) and re-added
    as an ion."""
    pdb = Path(__file__).parent / "data" / "fixtures" / "1HVI.pdb"
    # synthesize: HEM residue with an FE atom + a plain C
    lines = _res_lines("HEM", "A", 9, 3)  # A00 A01 A02
    # replace atom names/elements explicitly
    fixed = []
    for i, l in enumerate(lines):
        if i == 0:
            l = l[:13] + "FE " + l[16:]   # atom name field 13-15
            l = l[:76] + "FE" + l[78:]    # element field 76-77
        fixed.append(l)
    p = Path("/tmp/x.pdb")  # content check only, no gmx
    import tempfile, pathlib
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "h.pdb"
        p.write_text("".join(l + "\n" for l in
                             _res_lines("ALA", "A", 1, 4) + fixed) + "END\n")
        cfs = md._find_cofactors(p)
        assert len(cfs) == 1
        ems = md._embedded_metals(p, cfs)
        assert len(ems) == 1
        assert ems[0]["resname"] == "FE"
        assert ems[0]["chain"] == "A"
        assert ems[0]["charge"] == "+2.0"


def _metal_line(resname, chain, resseq, x, y, z, charge=None, elem=None):
    """80-column PDB line; the element sits at cols 77-78 (idx 76-77)
    and the charge field at cols 79-80 (idx 78+) per the PDB spec."""
    l = (f"HETATM{1000 + resseq:5d} {resname[0]:<3s} {resname:<3s} "
         f"{chain}{resseq:4d}    "
         f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00")
    l = l[:12] + " " + l[12:]  # shift the name into cols 13-16
    l = l[:76].ljust(76) + (elem or resname[:2]).ljust(2)
    if charge is not None:
        l = l + charge
    return l.ljust(80) + "\n"


def test_find_metals_charge_column(tmp_path):
    """R9: PDB columns 77-80 (charge on atom) override the element
    default when filled in."""
    pdb = tmp_path / "m.pdb"
    pdb.write_text(
        _metal_line("ZN", "A", 1, 1, 1, 1)                 # no column
        + _metal_line("ZN", "B", 2, 2, 2, 2, charge="1.00")  # +1 filled
        + "END\n")
    ms = md.find_metals(pdb)
    assert len(ms) == 2
    by_chain = {m["chain"]: m for m in ms}
    assert by_chain["A"]["charge"] == "+2.0"   # METAL_PROPS default
    assert by_chain["B"]["charge"] == "+1.0"   # from the PDB column


def test_standalone_metal_ligand(tmp_path):
    """R9: a docked metal ion on its own chain (single-atom 'ligand')
    must not go through ACPYPE; it is re-added as a standalone ion
    molecule (one [ molecules ] line PER ion)."""
    pdb = tmp_path / "c.pdb"
    parts = [l + "\n" for l in _res_lines("ALA", "A", 1, 4)]
    parts.append(_metal_line("ZN", "Z", 101, 3, 3, 3))
    parts.append(_metal_line("MG", "Z", 102, 4, 4, 4))
    parts.append("END\n")
    pdb.write_text("".join(parts))
    build = tmp_path / "b"
    build.mkdir()
    pdfs, stand = md._ligand_side_pdfs(pdb, {"A"}, build, [])
    assert pdfs == []            # no ACPYPE file for the metal chain
    assert len(stand) == 2
    assert {m["resname"] for m in stand} == {"ZN", "MG"}
    # two ions of the same element -> two molecule lines
    (build / "topol.top").write_text('; top\n[ molecules ]\n')
    (build / "combined.gro").write_text(
        "t\n2\n  1PRO   A0   1A1   1.0 1.0 1.0\n"
        "  2PRO   A0   1A2   1.1 1.0 1.0\n"
        "10.0 10.0 10.0\n")
    emb = [{"resname": "FE", "chain": "A", "resseq": 5,
            "xyz": (1.0, 1.0, 1.0)},
           {"resname": "FE", "chain": "B", "resseq": 6,
            "xyz": (2.0, 2.0, 2.0)}]
    md._append_embedded_ions(build, emb, build / "combined.gro",
                             build / "topol.top")
    top = (build / "topol.top").read_text()
    assert top.count("FEION 1") == 2, top
    gro = (build / "combined.gro").read_text().splitlines()
    assert int(gro[1]) == 4   # 2 protein + 2 ions
    assert gro[-1] == "10.0 10.0 10.0"


def test_eq_mdp_fingerprint_idempotency(tmp_path):
    """R9: a finished eq stage is reusable only while the MDP that
    produced it is unchanged — template edits (barostat, nsteps) must
    force a rerun."""
    (tmp_path / "eq_nvt.log").write_text("... Finished mdrun ...")
    (tmp_path / "eq_nvt.gro").write_text("1\n")
    fp_a = md._mdp_fingerprint(md._mdp_eq_nvt(10.0))
    (tmp_path / "eq_nvt.mdp.fp").write_text(fp_a)
    assert md._eq_stage_done(tmp_path, "eq_nvt", fp_a)
    # different nvt_ps -> different MDP -> not reusable
    fp_b = md._mdp_fingerprint(md._mdp_eq_nvt(20.0))
    assert fp_b != fp_a
    assert not md._eq_stage_done(tmp_path, "eq_nvt", fp_b)
    # missing fingerprint file (pre-R9 run) -> not reusable
    (tmp_path / "eq_nvt.mdp.fp").unlink()
    assert not md._eq_stage_done(tmp_path, "eq_nvt", fp_a)
    # comments do not change the fingerprint
    txt = md._mdp_eq_nvt(10.0)
    assert md._mdp_fingerprint(txt) == md._mdp_fingerprint(
        "; new comment\n" + txt)


def test_mdp_barostat_crescale():
    """R8: Berendsen barostat segfaults the 1HVI smoke system ~20 ps into
    production (both replicas die at step 10000 with PLOESS pressure
    scaling > 1%); C-rescale with tau_p 2 ps is the robust modern
    choice for short MD."""
    import re
    txt = md._mdp_md(2.0)
    assert re.search(r"^pcoupl\s*=\s*C-rescale\s*$", txt, re.M), \
        "production must use C-rescale"
    assert re.search(r"^tau_p\s*=\s*2\.0\s*$", txt, re.M)
    txt = md._mdp_eq_npt(20.0)
    assert re.search(r"^pcoupl\s*=\s*C-rescale\s*$", txt, re.M)


def test_eq_mdp_nsteps_units():
    """R8: equilibration nsteps conversion — 1 ps of 2fs steps = 500
    steps (the old formula was 1000x too long: a 10 ps NVT ran 10 ns)."""
    import re
    txt = md._mdp_eq_nvt(10.0)
    assert re.search(r"^nsteps\s*=\s*5000\s*$", txt, re.M), \
        "10 ps must be 5000 steps of 2fs (got %s)" % \
        re.search(r"^nsteps.*$", txt, re.M).group(0)
    txt = md._mdp_eq_npt(20.0)
    assert re.search(r"^nsteps\s*=\s*10000\s*$", txt, re.M), \
        "20 ps must be 10000 steps of 2fs"

def test_ligand_group_resolution():
    """R8: the ligand RMSD group name is system-specific (GROMACS derives
    it from the moleculetype name: ACPYPE mols are 'MOL' or LIG_* in this
    env, not 'Ligand' — the hardcoded prompt died with
    "No such group 'Ligand'" on the 1HVI smoke rebuild)."""
    g = md._ligand_group
    assert g({"System", "Protein", "MOL", "Water", "Ion"}) == "MOL"
    assert g({"System", "Protein", "LIG_C", "Water"}) == "LIG_C"
    assert g({"System", "Protein", "Ligand", "Water"}) == "Ligand"
    # no explicit ligand name: 'Other' is the ligand in a protein+ligand
    # system (everything not protein/water/ion)
    assert g({"System", "Protein", "Other", "Water", "Ion"}) == "Other"
    # no candidate at all -> skip ligand RMSD
    assert g({"System", "Protein", "Water"}) is None
    # DNA system: 'DNA' must not be mistaken for the ligand
    assert g({"System", "Protein", "DNA", "Other", "Water"}) == "Other"


# ---------------------------------------------------------------- R10: G3 — region-level flexibility (R1 收尾)

def _profile(n: int, base: float = 0.08) -> list[float]:
    return [base + 0.005 * (i % 3) for i in range(n)]


def test_flexible_regions_rigid():
    # uniform low RMSF (0.08 nm) -> no flexible region
    assert md.flexible_regions(_profile(200)) == []


def test_flexible_regions_single_loop():
    # one flexible loop at residues 100-115 (1-based) with 0.5-0.6 nm RMSF
    p = _profile(200)
    for i in range(99, 115):
        p[i] = 0.5 + 0.05 * ((i * 7) % 3)
    regs = md.flexible_regions(p)
    assert len(regs) == 1
    r = regs[0]
    assert 98 <= r["res_start"] <= 101, r
    assert 113 <= r["res_end"] <= 116, r
    assert r["n_res"] >= 12
    assert r["mean_rmsf_nm"] > 0.3
    assert r["max_rmsf_nm"] > r["mean_rmsf_nm"]


def test_flexible_regions_short_blip_ignored():
    p = _profile(100)
    p[50:53] = [0.6, 0.7, 0.6]  # 3-residue blip < min_res
    assert md.flexible_regions(p) == []


def test_flexible_regions_two_loops_ordered():
    p = _profile(200)
    for i in range(20, 30):
        p[i] = 0.45
    for i in range(150, 165):
        p[i] = 0.6
    regs = md.flexible_regions(p)
    assert len(regs) == 2
    assert regs[0]["res_start"] < regs[1]["res_start"]
    # the stiffer loop comes first in the list? No — ordered by position;
    # but the *interpretation* must lead with the most flexible region
    text = " ".join(md.interpret_stability(
        {"final_rmsd_mean": 0.2, "rmsf_profile_mean": p,
         "flexible_regions": regs}))
    assert "柔性区" in text
    assert "151" in text or "res 15" in text  # the 0.6-nm loop is cited


def test_flexible_regions_none_note():
    p = _profile(100)
    text = " ".join(md.interpret_stability(
        {"final_rmsd_mean": 0.1, "rmsf_profile_mean": p,
         "flexible_regions": []}))
    # R11: wording must not contradict the single-residue high-RMSF note
    assert "无明显连续柔性区" in text
    assert "各残基 RMSF 均低于阈值" not in text


def test_flexible_regions_ignores_short_and_nan():
    p = _profile(100)
    p[40] = float("nan")
    regs = md.flexible_regions(p)
    assert regs == []
