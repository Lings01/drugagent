"""Module E: GROMACS MD simulation + analysis (3x replicates, averaged).

System selection is agent-driven: best screening hit / binder / VHH /
input complex. Force field: amber19sb + GAFF2 (GROMACS>=2018) with
fallback to amber99sb-ildn on old GROMACS.
"""
from __future__ import annotations

import os
import re
import shutil
import time
from pathlib import Path

import numpy as np
from loguru import logger

from ..config import DEFAULTS, TOOLS, gmx_bin as _gmx_bin_known, n_cores, resolve_defaults
from ..llm import AgentBrain
from ..utils import jsave, pmap, run_cmd

GMX_OLD_ENV = Path("/usr/local/anaconda3/envs/gmx")


# --------------------------------------------------------------------------- #
# GROMACS discovery
# --------------------------------------------------------------------------- #
def gromacs() -> dict:
    """Return {gmx, gmxdata, ff, ff_ligand, version} preferring the 2023 build."""
    cands = [
        {"root": TOOLS / "gromacs", "ff": "amber19sb", "ff_ligand": "gaff2", "ver": "2023.1"},
        {"root": GMX_OLD_ENV, "ff": "amber99sb-ildn", "ff_ligand": "gaff", "ver": "4.6.5"},
    ]
    for c in cands:
        gmx = c["root"] / "bin" / "gmx"
        mdrun = c["root"] / "bin" / "mdrun"
        if not (gmx.is_file() or mdrun.is_file()):
            continue
        gmxdata = c["root"] / "share" / "gromacs"
        if not (gmxdata / "top" / "spc216.gro").is_file():
            continue
        # detect the force field actually shipped with this build
        top = gmxdata / "top"
        ff = c["ff"]
        if not (top / f"{ff}.ff").is_dir():
            for alt in ("amber19sb", "amber99sb-ildn", "amber99sb"):
                if (top / f"{alt}.ff").is_dir():
                    ff = alt
                    break
        return {"gmx": str(gmx if gmx.is_file() else mdrun),
                "gmxdata": str(gmxdata), "ff": ff,
                "ff_ligand": c["ff_ligand"], "ver": c["ver"]}
    raise FileNotFoundError(
        "GROMACS not found: build via `drugagent setup --gromacs` or "
        "/usr/local/anaconda3/envs/gmx")


def gmx_cmd(env: dict, *args: str, **kw) -> list[str]:
    return [env["gmx"]] + list(args) + _extra_env_flags(env, **kw)


def _extra_env_flags(env: dict, **kw) -> list[str]:
    return []


# --------------------------------------------------------------------------- #
# system selection
# --------------------------------------------------------------------------- #
def select_system(state: dict, brain: AgentBrain | None) -> dict:
    """Decide which complex to simulate."""
    opts = state.get("options", {})
    prep = state["target_prep"]
    choices: list[dict] = []

    if prep.get("ligand_pdbqt") and Path(prep["ligand_pdbqt"]).exists():
        choices.append({
            "name": "input_ligand",
            "label": f"靶点 + 输入配体({','.join(prep.get('ligand_resnames', []))})，已知活性参照",
            "ligand": True,
        })
    if "screening" in state and state["screening"].get("hits"):
        choices.append({
            "name": "screening_hit",
            "label": "靶点 + 筛选 top hit",
            "ligand": True,
        })
    if "binder" in state and state["binder"].get("best"):
        choices.append({
            "name": "binder",
            "label": "靶点 + 设计的 binder（蛋白-蛋白复合物）",
            "ligand": False,
        })
    if "vhh" in state and state["vhh"].get("best", {}).get("source") == "de_novo":
        choices.append({
            "name": "vhh",
            "label": "靶点 + 设计的纳米抗体（蛋白-蛋白复合物）",
            "ligand": False,
        })
    if not choices:
        raise RuntimeError("no MD system candidates; provide a complex PDB")

    choice = None
    if opts.get("md_system"):
        choice = next((c for c in choices if c["name"] == opts["md_system"]), None)
    if choice is None and brain is not None:
        ctx = "\n".join(f"{i+1}. {c['label']}" for i, c in enumerate(choices))
        dec = brain.decide(
            "md", "选择要跑 MD 的复合物（编号）", context=ctx,
            choices=[str(i + 1) for i in range(len(choices))], expect="choice")
        choice = choices[int(dec.answer) - 1]
        choice["rationale"] = dec.rationale
    if choice is None:
        choice = choices[0]
        choice["rationale"] = "默认选择（无 LLM/无用户指定）"
    return choice


def build_complex_pdb(state: dict, choice: dict, workdir: Path) -> Path:
    """Assemble the simulation complex PDB."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    prep = state["target_prep"]
    target_pdb = Path(prep["clean_pdb"])
    out = workdir / "complex.pdb"

    def _chain_lines(pdb: Path, chain: str | None,
                     keep_only: str | None = None) -> str:
        """Atom lines; keep only ``keep_only`` chain, re-label to ``chain``.

        Keeping the target's native chain letters matters for dimers:
        merging chains A+B into one chain would create phantom peptide
        bonds across the monomer boundary.
        """
        with open(pdb) as fh:
            out = []
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    if keep_only and line[21] != keep_only:
                        continue
                    if chain:
                        line = line[:21] + chain + line[22:]
                    out.append(line)
            return "".join(out)

    def _target_chains(pdb: Path) -> set:
        ch = set()
        with open(pdb) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    ch.add(line[21])
        return ch

    # next free chain letter(s) for added components
    used = _target_chains(target_pdb)
    free = [c for c in "CDEFGHIKLMNPQRSTVWY" if c not in used]
    lig_chain = free.pop(0) if free else "Z"
    des_chain = free.pop(0) if free else ("Y" if lig_chain != "Y" else "X")

    if choice["ligand"]:
        if choice["name"] == "input_ligand":
            # extract ligand from raw pdb
            from ..modules.target_prep import _extract_res
            lig_pdb = workdir / "ligand.pdb"
            _extract_res(target_pdb, lig_pdb, prep["ligand_resnames"])
        else:  # screening hit: use the docked pose (PDB preferred)
            hit = state["screening"]["hits"][0]
            pose = Path(hit["pose_pdb"]) if hit.get("pose_pdb") else None
            lig_pdb = workdir / "hit.pdb"
            if pose is not None and pose.is_file():
                txt = "".join(l for l in pose.read_text().splitlines(keepends=True)
                              if l.startswith("ATOM"))
                lig_pdb.write_text(txt + "END\n")
            else:
                with open(Path(hit["pose_pdbqt"])) as fi, open(lig_pdb, "w") as fo:
                    for line in fi:
                        if line.startswith("ATOM"):
                            fo.write(line)
                lig_pdb.write_text(lig_pdb.read_text() + "END\n")
        rec = workdir / "receptor_only.pdb"
        from ..modules.target_prep import _remove_res
        _remove_res(target_pdb, rec, prep.get("ligand_resnames", []) + ["HOH", "WAT"])
        txt = _chain_lines(rec, None) + _chain_lines(lig_pdb, lig_chain) + "END\n"
        out.write_text(txt)
        return _fix_c_term(out)
    # protein-protein: prefer the ESMFold-refined design (chain B of the
    # scored complex) over the raw RFdiffusion output, which can contain
    # incomplete residues (e.g. truncated HIS) that break pdb2gmx
    if choice["name"] == "binder":
        best = state["binder"]["best"]
        raw = Path(best["design"])
        refined = (Path(state["project_dir"]) / "03_binder" / "complex"
                   / f"{raw.stem}_complex.pdb")
    else:
        best = state["vhh"]["best"]
        raw = Path(best["design"])
        refined = (Path(state["project_dir"]) / "04_vhh" / "vhh_designs"
                   / f"{raw.stem}_complex.pdb")
    prot2 = _refined_design(target_pdb, raw, refined)
    from ..modules.target_prep import _remove_res
    rec = workdir / "receptor_only.pdb"
    # drop the crystallographic ligand (if any) and waters for protein-protein MD
    _remove_res(target_pdb, rec, prep.get("ligand_resnames", []) + ["HOH", "WAT"])
    # prot2 is a scored complex (target chain A + design chain B); take
    # only the design (chain B) and label it on a free chain letter
    txt = _chain_lines(rec, None) + _chain_lines(prot2, des_chain, keep_only="B") + "END\n"
    out.write_text(txt)
    return _fix_c_term(out)


# --------------------------------------------------------------------------- #
# system building
# --------------------------------------------------------------------------- #
def _refined_design(target_pdb: Path, raw: Path, refined: Path) -> Path:
    """Return a PDB whose chain B is the ESMFold-refined design.

    The scored complex (target chain A + design chain B) is the preferred
    source for MD: raw RFdiffusion output can contain incomplete residues
    (backbone-only input residues, truncated side chains) that break
    pdb2gmx.  If the stored complex predates the ESMFold multimer fix
    (no chain B atoms), the complex is re-scored on the spot.
    """
    def _n_chain_b(pdb: Path) -> int:
        n = 0
        try:
            for l in pdb.read_text().splitlines():
                if l.startswith(("ATOM", "HETATM")) and l[21] == "B":
                    n += 1
        except OSError:
            pass
        return n

    if refined.is_file() and _n_chain_b(refined) >= 30:
        return refined
    try:
        from ..modules.binder import _ca_sequence, _make_complex
        from ..modules.esmfold_run import predict
        # the raw RF output carries the target backbone as its own chain A;
        # keep only its chain B (the actual design) for the complex
        design_only = refined.parent / f".{raw.stem}_design_only.pdb"
        with open(raw) as fh, open(design_only, "w") as fo:
            for l in fh:
                if l.startswith(("ATOM", "HETATM")) and l[21] == "B":
                    fo.write(l)
            fo.write("END\n")
        use_raw = design_only
        _make_complex(target_pdb, use_raw, refined)
        seqs = [_ca_sequence(refined, "A"), _ca_sequence(refined, "B")]
        if all(seqs):
            out = predict(seqs, num_recycles=3, device="cpu")
            refined.write_text(out["pdb"])
            if _n_chain_b(refined) >= 30:
                return refined
    except Exception as e:  # noqa: BLE001
        logger.warning(f"refined design re-scoring failed for {raw.name}: {e}")
    try:
        design_only.unlink(missing_ok=True)
    except OSError:
        pass
    return raw


def _fix_c_term(pdb: Path) -> Path:
    """GROMACS rtp matching for C-terminal residues.

    PDB files write the C-terminal carboxylate as ``OXT`` on the normal
    residue (PHE); GROMACS needs the terminal residue type (CPHE) whose
    atoms are OC1/OC2.  Renames are applied to every residue carrying an
    OXT atom.
    """
    lines = pdb.read_text().splitlines()
    # 1) find residues with OXT
    cterm = set()  # (chain, resseq, icode)
    for l in lines:
        if l[:6] in ("ATOM  ", "HETATM"):
            if l[12:16].strip() == "OXT":
                cterm.add((l[21], l[22:26].strip(), l[26]))
    if not cterm:
        return pdb
    out = []
    for l in lines:
        if l[:6] in ("ATOM  ", "HETATM"):
            chain, resseq, icode = l[21], l[22:26].strip(), l[26]
            if (chain, resseq, icode) in cterm:
                name = l[17:20].strip()
                newname = name if name.startswith("C") else "C" + name
                if l[12:16].strip() == "OXT":
                    l = f"{l[:12]}{'OC2':<4}{l[16:]}"
                l = f"{l[:17]}{newname:<4}{l[21:]}"
        out.append(l)
    pdb.write_text("\n".join(out) + "\n")
    return pdb


def _top_spce(env: dict) -> str:
    # 2023.x: top dir is $gmxdata/top ; older: also top
    for cand in (f"{env['gmxdata']}/top/spc216.gro", f"{env['gmxdata']}/share/top/spc216.gro"):
        if Path(cand).is_file():
            return cand
    raise FileNotFoundError("spc216.gro not found under GMXDATA")


def _find_pdb2gmx(env: dict) -> list:
    """Return the pdb2gmx command prefix.

    GROMACS >= 5 uses the unified CLI (``gmx pdb2gmx``); older installs
    ship a standalone ``pdb2gmx`` binary.
    """
    p = Path(env["gmx"]).parent / "pdb2gmx"
    if p.is_file():
        return [str(p)]
    return [env["gmx"], "pdb2gmx"]


def _normalize_elements(pdb: Path) -> Path:
    """Rebuild the PDB element column (cols 77-78) from the atom name.

    Docking tools write ligand PDBs with shifted columns that confuse
    RDKit's PDB parser.
    """
    def _elem(name: str) -> str:
        n = name.strip()
        if n.startswith("CL"):
            return "Cl"
        if n.startswith("BR"):
            return "Br"
        if not n:
            return "C"
        c = n[0]
        return c.upper() if c.isalpha() else "C"

    out = []
    for l in pdb.read_text().splitlines():
        if l[:6] in ("ATOM  ", "HETATM"):
            e = _elem(l[12:16])
            l = f"{l[:76]:<76}{e:<2}"
        out.append(l)
    pdb.write_text("\n".join(out) + "\n")
    return pdb


def _unify_mol2_residue(mol2: Path) -> None:
    """ACPYPE/antechamber accepts exactly ONE residue. obabel carries
    PDB residue names (resname+resseq, e.g. 'HEM9') into the mol2 atom
    records and may mix them with a default 'MOL' label; rewrite every
    atom's residue fields to a single canonical residue."""
    lines = mol2.read_text().splitlines()
    out: list[str] = []
    in_atom = False
    for l in lines:
        if l.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            out.append(l)
            continue
        if l.startswith("@<TRIPOS>"):
            in_atom = False
        if in_atom:
            t = l.split()
            if len(t) >= 8:
                t[6], t[7] = "1", "MOL"
                l = " ".join(t[:8]) + (
                    "  " + " ".join(t[8:]) if len(t) > 8 else "")
        out.append(l)
    mol2.write_text("\n".join(out) + "\n")


def _run_acpype(lig_pdb: Path, workdir: Path, env: dict,
                name: str = "LIG") -> dict:
    """ACPYPE ligand parameterization (GAFF2; AM1-BCC or Gasteiger).

    The ligand is converted to a bonded mol2 with RDKit first, because
    ACPYPE's own PDB->mol2 bond perception fails on some docked poses.
    `name` isolates parallel molecules (multi-ligand / cofactor builds).
    Returns {"itp": str, "gro": str|None}.
    """
    from ..config import ENV_DIR
    acpype = None
    for d in (ENV_DIR / "bin", TOOLS / "tools" / "bin"):
        if (d / "acpype").is_file():
            acpype = str(d / "acpype")
            break
    if acpype is None:
        acpype = shutil.which("acpype")
    if acpype is None:
        raise RuntimeError("acpype not found (pip install acpype in the agent env)")
    outdir = workdir / f"acpype_{name}"
    outdir.mkdir(parents=True, exist_ok=True)
    envd = dict(os.environ)
    # conda AMBERTools prefix: self-contained tleap (the pip-bundled one is
    # missing HDF libs) and sqm for AM1-BCC charges
    amb = TOOLS / "ambertools" / "bin"
    if amb.is_dir():
        envd["PATH"] = f"{amb}:{envd.get('PATH', '')}"
    # sqm on PATH enables AM1-BCC; otherwise fall back to Gasteiger
    charge = "gas"
    from ..config import ENV_DIR as _ED
    for d in (amb, _ED / "bin", Path("/usr/local/anaconda3/envs/gmx/bin")):
        if d.is_dir() and (d / "sqm").is_file():
            envd["PATH"] = f"{d}:{envd.get('PATH', '')}"
            charge = "bcc"
            break
    outdir = outdir.resolve()
    # bonded mol2 input via obabel (RDKit 2026 ships no mol2 writer)
    src = outdir / f"{name}.pdb"
    src.write_text(_normalize_elements(lig_pdb).read_text())
    mol2 = outdir / f"{name}.mol2"
    obabel = None
    for d in (ENV_DIR / "bin", TOOLS / "tools" / "bin"):
        if (d / "obabel").is_file():
            obabel = str(d / "obabel")
            break
    if obabel is None:
        obabel = shutil.which("obabel")
    if obabel:
        run_cmd([obabel, str(src), "-O", str(mol2)],
                log_file=outdir / "obabel.log")
        if mol2.is_file():
            _unify_mol2_residue(mol2)
    inp = mol2 if mol2.is_file() else src
    cmd = [acpype, "-i", str(inp), "-b", name, "-o", "gmx",
           "-c", charge, "-a", "gaff2", "-w"]
    run_cmd(cmd, cwd=outdir, env=envd, log_file=outdir / "acpype.log")
    # ACPYPE writes the final files into a {basename}.acpype/ scratch dir;
    # GMX outputs are <name>_GMX.* (older ACPYPE: <name>.*). ACPYPE also
    # lowercases the base name — try both casings.
    scratch = outdir / f"{name}.acpype"
    dirs = [scratch, outdir]
    itp = gro = None
    for d in dirs:
        if itp is None:
            for cand in (d / f"{name}_GMX.itp", d / f"{name}.itp",
                         d / f"{name.lower()}_GMX.itp",
                         d / f"{name.lower()}.itp"):
                if cand.is_file():
                    itp = cand
                    break
        if gro is None:
            for cand in (d / f"{name}_GMX.gro", d / f"{name}.gro",
                         d / f"{name.lower()}_GMX.gro",
                         d / f"{name.lower()}.gro"):
                if cand.is_file():
                    gro = cand
                    break
        if itp is not None:
            break
    if itp is None:
        raise RuntimeError(f"ACPYPE did not produce an ITP under {outdir}")
    return {"itp": str(itp), "gro": str(gro) if gro is not None else None}


# --------------------------------------------------------------------------- #
# metal ions (R3)
# --------------------------------------------------------------------------- #
# single-atom metal residues recognized by find_metals; charges are the
# common biological oxidation states (v1 default +2, overridable per ion)
METAL_RES = {"ZN", "FE", "MN", "NI", "CU", "CO", "MG", "CA"}
# resname -> (molname, charge, mass, atomic number)
METAL_PROPS = {
    "ZN": ("ZNION", "+2.0", "65.38", 30),
    "FE": ("FEION", "+2.0", "55.85", 26),
    "MN": ("MNION", "+2.0", "54.94", 25),
    "NI": ("NIION", "+2.0", "58.69", 28),
    "CU": ("CUION", "+2.0", "63.55", 29),
    "CO": ("COION", "+2.0", "58.93", 27),
    "MG": ("MGION", "+2.0", "24.31", 12),
    "CA": ("CAION", "+2.0", "40.08", 20),
}
# protein donor atom names considered as metal coordinators
# (v1: no CB — a carbon-coordinated metal is rare in crystals and CB would
# also break the N/O/S element filter used for GRO matching)
COORD_NAMES = {"N", "O", "S", "OD1", "OD2", "ND1", "ND2", "NE1", "NE2",
               "OG", "OG1", "SD",
               # nucleic-acid coordination donors (R4): phosphates
               # (both PDB namings), sugar O's, base N/O's
               "OP1", "OP2", "O1P", "O2P", "O4'", "O2'", "O3'",
               "N1", "N3", "N7", "N9", "O6", "O2"}


def _pdb_charge_column(line: str) -> float | None:
    """R9: PDB 'charge on atom' field (cols 79-80, after the 2-char
    element field at 77-78; depositors sometimes overflow it). Returns
    the numeric charge or None when the field is blank."""
    # element field 76-77 is 2 chars; a 1-char element (N, ZN->'Z'?)
    # leaves col 77 blank, so the charge window is 77-80
    tail = line[77:].strip()
    if not tail:
        return None
    m = re.match(r"^\s*(?:[A-Z][a-z]?\s*)?([+-]?\d+(?:\.\d+)?)", tail)
    if not m:
        return None
    # a bare integer in the charge window is the element number
    # (Zn=30, Fe=26), not a charge — charges need a sign or a decimal
    v = m.group(1)
    return float(v) if ("." in v or v.startswith(("+", "-"))) else None


def _fmt_charge(v: float) -> str:
    return f"{v:+.1f}"


def find_metals(pdb: Path) -> list[dict]:
    """Single-atom metal ions in a PDB.

    Returns [{resname, chain, resseq, xyz, n_atoms, molname, charge, mass}]
    for every residue whose name is in METAL_RES and whose heavy-atom count
    is 1-2 (multi-atom cofactors like HEM are not handled in v1).
    """
    from collections import defaultdict
    res: "defaultdict[tuple, list]" = defaultdict(list)
    with open(pdb) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            name = line[12:16].strip()
            if not name or name == "H":
                continue
            key = (line[21], int(line[22:26]), line[17:20].strip())
            res[key].append(line)
    out = []
    for (chain, rseq, rname), lines in res.items():
        if rname not in METAL_RES or not 1 <= len(lines) <= 2:
            continue
        xyz = []
        for l in lines:
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        xyz = tuple(sum(col) / len(lines) for col in zip(*xyz))
        mol, charge, mass, _znum = METAL_PROPS[rname]
        # R9: the PDB charge field (cols 79-80) wins over the element
        # default when the depositor filled it in
        chgs = [c for c in (_pdb_charge_column(l) for l in lines)
                if c is not None]
        if chgs:
            charge = _fmt_charge(sum(chgs) / len(chgs))
        out.append({"resname": rname, "chain": chain, "resseq": rseq,
                    "xyz": xyz, "n_atoms": len(lines), "molname": mol,
                    "charge": charge, "mass": mass})
    return out


# standard 3-letter DNA residue names (unambiguous); 1-letter A/C/G/T are
# ambiguous (DNA or RNA) — decided by the O2' atom (RNA only)
DNA_RES = {"DA", "DT", "DC", "DG"}
RNA_BASES = {"A", "U", "C", "G"}


def find_nucleic_acids(pdb: Path) -> list[dict]:
    """Chains containing nucleic-acid residues (R4).

    Returns [{"chain", "n_res", "type": "DNA" | "RNA"}] in PDB chain
    order. 3-letter names DA/DT/DC/DG are DNA; 1-letter A/C/G/T are
    assigned by O2' presence (RNA), U is RNA.
    """
    # chain -> (resname, resseq) -> atom names
    res: dict[str, dict[tuple[str, str], set[str]]] = {}
    order: list[str] = []
    with open(pdb) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            ch = line[21]
            if ch not in res:
                res[ch] = {}
                order.append(ch)
            key = (line[17:20].strip(), line[22:26].strip())
            res[ch].setdefault(key, set()).add(line[12:16].strip())
    out = []
    for ch in order:
        dna = rna = 0
        for (rname, _rseq), names in res[ch].items():
            if rname in DNA_RES:
                dna += 1
            elif rname in RNA_BASES:
                if "O2'" in names:
                    rna += 1
                else:
                    dna += 1
        if dna + rna == 0:
            continue
        out.append({"chain": ch, "n_res": dna + rna,
                    "type": "DNA" if dna >= rna else "RNA"})
    return out


def write_metal_itp(metal: dict, out: Path) -> Path:
    """Static single-atom moleculetype for a metal ion.

    The element symbol is used as atom name and type; a matching
    [ atomtypes ] entry is defined here (the FF atomtypes.itp has no
    metal ions). sigma/epsilon are 0: the ion interacts via Coulomb only
    (documented limitation — no dedicated Zn2+ LJ parameters in v1).
    """
    out = Path(out)
    mol = metal["molname"]
    el = metal["resname"]
    znum = METAL_PROPS.get(el, (None, None, None, 0))[3]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        f"; single-atom {el} ion (drugagent R3 template)\n"
        f"[ atomtypes ]\n"
        f"; name  atomic number  mass  charge  ptype  sigma  epsilon\n"
        f"{el:<5s}{znum:<15d}{metal['mass']:<8s}  0.0  A  0.0  0.0\n"
        f"[ moleculetype ]\n; Name            nrexcl\n{mol}      3\n"
        f"[ atoms ]\n;  nr  type  resnr residue  atom  cgnr  charge  mass\n"
        f"   1  {el:<4s}     1    {el:<4s} {el:<5s} 1"
        f"  {metal['charge']:<6s}{metal['mass']:<8s}\n"
    )
    return out


def metal_coordinators(pdb: Path, metals: list[dict],
                       cutoff: float = 2.5) -> list[dict]:
    """Protein donor atoms (N/O/S-type) within `cutoff` of each metal.

    Returns [{metal_idx, resseq, chain, atom, distance}] sorted by metal
    then distance; the caller maps (resseq, atom) to topology indices.
    """
    donors = []
    with open(pdb) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            name = line[12:16].strip()
            rname = line[17:20].strip()
            if name not in COORD_NAMES or rname in METAL_RES:
                continue
            donors.append({
                "chain": line[21], "resseq": int(line[22:26]),
                "atom": name,
                "xyz": (float(line[30:38]), float(line[38:46]),
                        float(line[46:54])),
            })
    out = []
    for mi, m in enumerate(metals):
        for d in donors:
            if d["chain"] == m["chain"] and d["resseq"] == m["resseq"]:
                continue
            dist = float(sum((d["xyz"][k] - m["xyz"][k]) ** 2
                             for k in range(3)) ** 0.5)
            if dist <= cutoff:
                out.append({"metal_idx": mi, "chain": d["chain"],
                            "resseq": d["resseq"], "atom": d["atom"],
                            "distance": round(dist, 3)})
    out.sort(key=lambda x: (x["metal_idx"], x["distance"]))
    return out


def _element_of(name: str) -> str:
    """Element from a standard atom name (first letter: N/NE2/ND1 -> N,
    O/OD1/OG1 -> O, S/SD -> S, C/CA/CB -> C)."""
    return name[0] if name else ""


def _donor_index_map(gro: Path, pdb: Path,
                     tol: float = 0.5) -> dict[tuple, int]:
    """(chain, resseq, atom) -> 1-based protein.gro index.

    PDB donor atoms are matched to GRO atoms by coordinate proximity with
    an element check: pdb2gmx renames C-terminal O to OC1/OC2 (CPHE etc.),
    drops/renames other atoms, and the GRO has no chain ids, so name- and
    order-based matching are not robust. Returns {} if any donor is
    unmatched (caller then skips restraints).
    """
    g = []
    with open(gro) as fh:
        next(fh, None)
        n = int(next(fh).strip())
        for i in range(n):
            line = next(fh)
            name = line[10:15].strip()
            if _element_of(name) in ("N", "O", "S"):
                g.append((i + 1, _element_of(name), (
                    float(line[20:28]), float(line[28:36]),
                    float(line[36:44]))))
    out = {}
    used = set()
    with open(pdb) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            name = line[12:16].strip()
            rname = line[17:20].strip()
            if name not in COORD_NAMES or rname in ("HOH", "WAT", "SOL"):
                continue
            # PDB is in A, GRO in nm
            xyz = (float(line[30:38]) / 10.0, float(line[38:46]) / 10.0,
                   float(line[46:54]) / 10.0)
            el = _element_of(name)
            best, bestd = None, tol
            for idx, gel, gxyz in g:
                if idx in used or gel != el:
                    continue
                d = (sum((xyz[k] - gxyz[k]) ** 2 for k in range(3)) ** 0.5)
                if d < bestd:
                    best, bestd = idx, d
            if best is None:
                return {}
            used.add(best)
            out[(line[21], int(line[22:26]), name)] = best
    return out


def _metal_gro_lines(metal: dict) -> str:
    """GRO line(s) for a metal ion (single atom).

    GRO fixed format: resnum(5, right) resname(5, left) atomname(5, right)
    atomnum(5, right) x y z (8.3f, nm).
    """
    x, y, z = metal["xyz"]
    name = metal["resname"]
    return (f"{metal['resseq']:>5d}{name:<5s}{name:>5s}{1:>5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}\n")


def _protein_chain_lines(prot_gro: Path,
                         ordered: list[tuple[str, int]]) -> list[tuple[str, list[str]]]:
    """mol -> GRO atom lines, in the given (topology) order. GRO has no
    chain ids; the per-molecule atom counts from the pdb2gmx itps give the
    split points — they are the true molecule sizes, which can differ from
    the PDB atom counts if pdb2gmx added/removed hydrogens."""
    glines = Path(prot_gro).read_text().splitlines()
    n = int(glines[1].strip())
    atoms = glines[2:2 + n]
    out, pos = [], 0
    for mol, c in ordered:
        out.append((mol, atoms[pos:pos + c]))
        pos += c
    return out


def _molecule_order(build: Path) -> list[str]:
    """Molecule names of the per-chain pdb2gmx itps, in [ molecules ]
    order (protein and nucleic-acid chains alike)."""
    stats = _chain_itp_stats(build)
    top_txt = (build / "topol.top").read_text()
    order = []
    in_mol = False
    for l in top_txt.splitlines():
        if l.startswith("[ molecules ]"):
            in_mol = True
            continue
        if in_mol and l.strip() and not l.startswith(";"):
            name = l.split()[0]
            if name in stats:
                order.append(name)
    return order


def _chain_id_of(mol: str) -> str:
    """Molecule name ('Protein_chain_B' / 'DNA_chain_C') -> chain id."""
    return mol.rsplit("_", 1)[-1]


def _combined_gro_lines(build: Path, prot_gro: Path,
                        metals: list[dict],
                        lig_lines: list[str] | None = None) -> list[str]:
    """Combined GRO atom lines in topology molecule order: each chain
    molecule (protein or nucleic acid) with its own metal atoms appended,
    then crystal waters pdb2gmx kept, then the ligand (if any). Must
    match the [ molecules ] order of the final topology."""
    stats = _chain_itp_stats(build)
    order = _molecule_order(build)
    by_chain: dict[str, list[dict]] = {}
    for m in metals:
        by_chain.setdefault(m["chain"], []).append(m)
    # the stats counts include the metal atoms _integrate_metals already
    # appended to the chain itps — split the GRO with the pre-metal sizes
    pre = [(mol, stats[mol][1] - len(by_chain.get(_chain_id_of(mol), [])))
           for mol in order]
    chains = dict(_protein_chain_lines(prot_gro, pre))
    # atoms after the last chain in the GRO (crystal waters pdb2gmx kept
    # as SOL molecules) go after the chains, matching [ molecules ] order
    glines = Path(prot_gro).read_text().splitlines()
    n_gro = int(glines[1].strip())
    used = sum(c for (_m, c) in pre)
    leftover = glines[2 + used:2 + n_gro]
    lines: list[str] = []
    for mol in order:
        lines += chains.get(mol, [])
        for m in by_chain.get(_chain_id_of(mol), []):
            lines += _metal_gro_lines(m).splitlines()
    lines += leftover
    if lig_lines:
        lines += lig_lines
    return lines


def _chain_itp_stats(build: Path) -> dict[str, tuple[Path, int]]:
    """Molecule name (as written in [ molecules ]) -> (itp path, atom
    count) for every per-chain pdb2gmx itp — protein AND nucleic-acid
    chains (topol_Protein_chain_*.itp / topol_DNA_chain_*.itp /
    topol_RNA_chain_*.itp)."""
    out = {}
    for f in sorted(build.glob("topol_*.itp")):
        if f.name.startswith("posre_") or f.name.startswith("metal_"):
            continue
        name, n = None, 0
        in_mt = in_atoms = False
        for l in f.read_text().splitlines():
            if l.startswith("["):
                in_mt = l.startswith("[ moleculetype ]")
                in_atoms = l.startswith("[ atoms ]")
                continue
            if in_mt and l.strip() and not l.startswith(";"):
                name = l.split()[0]
                in_mt = False  # only the name line, not the atoms below
            # pdb2gmx writes a "; residue ..." comment per residue inside
            # the [ atoms ] section — only real atom lines count
            if in_atoms and l.strip() and not l.startswith(";"):
                n += 1
        if name:
            out[name] = (f, n)
    return out


def _chain_offsets(top_txt: str, stats: dict) -> dict[str, int]:
    """chain id -> molecule-relative offset (molecules appear in
    [ molecules ] order)."""
    offsets, off = {}, 0
    in_mol = False
    for l in top_txt.splitlines():
        if l.startswith("[ molecules ]"):
            in_mol = True
            continue
        if in_mol and l.strip() and not l.startswith(";"):
            name = l.split()[0]
            if name in stats:
                offsets[_chain_id_of(name)] = off
                off += stats[name][1]
    return offsets


def _ff_has_atomtype(build: Path, el: str) -> bool:
    """Does the force field (first #include of topol.top) already define
    this atomtype? GROMACS's ion.itp files define common ions (MG, ZN,
    NA, CL, CA, ...); redefining them only triggers a grompp warning and
    the FF parameters are the curated ones anyway."""
    import re
    top = (build / "topol.top").read_text()
    m = re.search(r'#include "([^"]+)"', top)
    if not m:
        return False
    gmxdata = os.environ.get("GMXDATA", "")
    ff_dir = Path(gmxdata) / "top" / Path(m.group(1)).parent
    if not ff_dir.is_dir():
        ff_dir = Path(gmxdata) / "top"
    for f in ff_dir.glob("*.itp"):
        in_at = False
        for l in f.read_text().splitlines():
            if l.startswith("[ atomtypes ]"):
                in_at = True
                continue
            if l.startswith("["):
                in_at = False
            if in_at and l.split() and l.split()[0] == el:
                return True
    return False


def _metal_itp_atomtypes_only(el: str) -> str:
    """[ atomtypes ]-only itp for a metal element (no moleculetype)."""
    props = METAL_PROPS.get(el, (None, None, "0.0", 0))
    mass, znum = props[2], props[3]
    return (f"; single-atom {el} ion atomtype (drugagent R3)\n"
            f"[ atomtypes ]\n"
            f"; name  atomic number  mass  charge  ptype  sigma  epsilon\n"
            f"{el:<5s}{znum:<15d}{mass:<8s}  0.0  A  0.0  0.0\n")


def _integrate_metals(build: Path, metals: list[dict], coords: list[dict],
                      dmap: dict) -> None:
    """Attach metal atoms to the protein chain molecule they belong to and
    add [ distance_restraints ].

    GROMACS parses restraint sections in the context of the immediately
    preceding moleculetype, so a separate ZNION molecule can never be
    restrained against protein atoms (grompp: "atom index out of bounds
    (1-1)"). The metal therefore becomes the last atom(s) of the owning
    chain's molecule; restraint indices are molecule-relative.
    """
    stats = _chain_itp_stats(build)
    top_txt = (build / "topol.top").read_text()
    offsets = _chain_offsets(top_txt, stats)
    mol_of_chain = {_chain_id_of(mol): mol for mol in stats}
    gidx = {id(m): i for i, m in enumerate(metals)}
    by_chain: dict[str, list[dict]] = {}
    for m in metals:
        by_chain.setdefault(m["chain"], []).append(m)
    els: list[str] = []
    for ch, ms in by_chain.items():
        mol = mol_of_chain.get(ch)
        if mol is None:
            logger.warning(f"metal on chain {ch}: itp not found; "
                           "coordination restraints skipped")
            continue
        f, n = stats[mol]
        # global metal index -> position of this metal among this chain's
        # appended atoms (restraint indices are molecule-relative)
        chain_pos = {gidx[id(m)]: k for k, m in enumerate(ms)}
        off = offsets.get(ch, 0)
        lines = f.read_text().splitlines()
        # end of the [ atoms ] section
        i_atoms = next(i for i, l in enumerate(lines)
                       if l.startswith("[ atoms ]"))
        last_atom = i_atoms
        for i in range(i_atoms + 1, len(lines)):
            if lines[i].startswith("["):
                break
            if lines[i].strip():
                last_atom = i
        atom_lines = []
        for k, m in enumerate(ms):
            idx = n + k + 1
            # nr type resnr residue atom cgnr charge mass
            atom_lines.append(
                f"{idx:>5d}  {m['resname']:<5s}{m['resseq']:>5d}  "
                f"{m['resname']:<5s}{m['resname']:<5s}{1:>5d}  "
                f"{m['charge']:>8s}{m['mass']:>8s}")
        lines = lines[:last_atom + 1] + atom_lines + lines[last_atom + 1:]
        # distance restraints (molecule-relative indices)
        restr = ["", "[ distance_restraints ]",
                 ";  ai    aj  type  index  type'     low     up1     up2    fac"]
        label = 0
        for c in coords:
            if c.get("_chain") != ch or c["metal_idx"] not in chain_pos:
                continue
            gi = dmap.get((c["chain"], c["resseq"], c["atom"]))
            if not gi:
                continue
            i_rel = gi - off
            j_rel = n + chain_pos[c["metal_idx"]] + 1
            d = float(c["distance"]) / 10.0  # A -> nm
            r0 = max(0.005, d - 0.015)
            r1 = d + 0.015
            r2 = r1 + 0.5
            restr.append(f"{i_rel:>5d}{j_rel:>5d}    1{label:>5d}      2"
                         f"{r0:9.4f}{r1:9.4f}{r2:9.4f}    1.0")
            label += 1
        if label:
            lines += restr
        f.write_text("\n".join(lines).rstrip() + "\n")
        for m in ms:
            if m["resname"] not in els:
                els.append(m["resname"])
    # atomtype-only itps, included right after the forcefield include
    # skip elements the force field already parameterizes (FF ion.itp)
    els = [el for el in els if not _ff_has_atomtype(build, el)]
    for el in els:
        (build / f"metal_{el.lower()}.itp").write_text(
            _metal_itp_atomtypes_only(el))
    inc_lines = [f'#include "metal_{el.lower()}.itp"' for el in els]
    tlines = top_txt.splitlines()
    for i, l in enumerate(tlines):
        if l.startswith("#include"):
            tlines = tlines[:i + 1] + inc_lines + tlines[i + 1:]
            break
    (build / "topol.top").write_text("\n".join(tlines) + "\n")


def gro_atom_index(gro: Path) -> dict[tuple, int]:
    """(residue number, atom name) -> 1-based atom index in a GRO file."""
    idx = {}
    with open(gro) as fh:
        next(fh, None)  # title
        n = int(next(fh).strip())
        for i in range(n):
            parts = next(fh).split()
            resnum = int(parts[1])
            atom = parts[2]
            idx[(resnum, atom)] = i + 1
    return idx


def pair_restraints_lines(coords: list[dict], metal_atom_idxs: dict[int, int],
                          window: float = 0.15) -> list[str]:
    """GROMACS [ distance_restraints ] entries keeping metals at their
    crystal coordination distances.

    Flat-bottom potential (manual eq. 219): zero force inside
    [d-window, d+window] (r0..r1, nm), harmonic outside, linear beyond
    r2. Force constant = MDP disre-fc (2000 in our templates) x fac.
    Each restraint gets a unique index label so GROMACS never groups
    them into one NOE-style ensemble average; type'=2 disables time
    averaging. metal_atom_idxs maps metal_idx -> 1-based top index.
    """
    lines = ["", "[ distance_restraints ]",
             ";  ai    aj  type  index  type'     low     up1     up2    fac"]
    for label, c in enumerate(coords):
        mi = c["metal_idx"]
        if mi not in metal_atom_idxs:
            continue
        j = metal_atom_idxs[mi]
        i = c.get("top_index")
        if not i:
            continue
        d = float(c["distance"]) / 10.0  # A -> nm
        r0 = max(0.005, d - window / 10.0)
        r1 = d + window / 10.0
        r2 = r1 + 0.5
        lines.append(f"{i:>5d}{j:>5d}    1{label:>5d}      2"
                     f"{r0:9.4f}{r1:9.4f}{r2:9.4f}    1.0")
    return lines


def build_system(complex_pdb: Path, workdir: Path, env: dict,
                 *, is_ligand: bool, salt: float = 0.15,
                 box_margin: float = 1.0,
                 divalent: str | None = None,
                 divalent_m: float = 0.0) -> dict:
    """pdb2gmx (+ ACPYPE for ligand) -> solvate -> genion -> EM."""
    build = workdir / "build"
    build.mkdir(parents=True, exist_ok=True)
    os.environ["GMXDATA"] = env["gmxdata"]
    gmx = env["gmx"]

    if is_ligand:
        return _build_ligand_system(complex_pdb, build, gmx, env, salt,
                                    divalent=divalent,
                                    divalent_m=divalent_m)
    return _build_protein_system(complex_pdb, build, gmx, env, salt,
                                 divalent=divalent,
                                 divalent_m=divalent_m)


def _write_topol(build: Path, ff: str, extra_includes: list[str]) -> Path:
    top = build / "topol.top"
    lines = ["#include \"amber99sb-ildn.ff/forcefield.itp\"" if "ildn" in ff
             else f"#include \"{ff}.ff/forcefield.itp\"",
             "#include \"spce.itp\"",
             "#include \"ions.itp\""]
    lines += [f'#include "{inc}"' for inc in extra_includes]
    lines += ["", "[ system ]", "Default", "", "[ molecules ]"]
    top.write_text("\n".join(lines))
    return top


def _build_protein_system(complex_pdb: Path, build: Path, gmx: str,
                          env: dict, salt: float,
                          divalent: str | None = None,
                          divalent_m: float = 0.0) -> dict:
    pdb2gmx = _find_pdb2gmx(env)
    metals = [m for m in find_metals(complex_pdb) if m["n_atoms"] == 1]
    if metals:
        logger.info(f"metal ions in system: "
                    f"{[(m['resname'], m['chain'], m['resseq']) for m in metals]}")
        src = build / "protein_nometa.pdb"
        with open(complex_pdb) as fi, open(src, "w") as fo:
            for l in fi:
                if (l.startswith(("ATOM", "HETATM"))
                        and l[17:20].strip() in METAL_RES):
                    continue
                fo.write(l)
        if not l.endswith("\n"):
            fo.write("\n")
    else:
        src = complex_pdb
    # R8: cofactors on protein chains (HEM ...) crash pdb2gmx; move them
    # to the ACPYPE side exactly like a ligand
    cof_pdfs: list[Path] = []
    cfs = _find_cofactors(src)
    emb = _embedded_metals(src, cfs)
    if cfs:
        logger.info(f"cofactors on protein chains: "
                    f"{[(c['resname'], c['chain'], c['resseq']) for c in cfs]}"
                    " — moving to the ACPYPE side")
        work_pdb = build / "complex_cofactors.pdb"
        _reassign_cofactor_chains(src, work_pdb)
        prot_pdb = build / "protein.pdb"
        prot_chains = _classify_chains(work_pdb)
        with open(work_pdb) as fh, open(prot_pdb, "w") as fp:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    if line[21] in prot_chains:
                        fp.write(line)
                elif line.startswith(("CRYST1", "REMARK")):
                    fp.write(line)
        fp2 = open(prot_pdb, "a"); fp2.write("END\n"); fp2.close()
        cof_pdfs, _stand = _ligand_side_pdfs(work_pdb, prot_chains,
                                             build, emb)
        src = prot_pdb
    if emb:
        logger.info(f"embedded cofactor metals re-added as ions: "
                    f"{[(m['resname'], m['chain'], m['resseq']) for m in emb]}")
    # cwd=build: pdb2gmx writes the posre_*.itp side files to the
    # current working directory, not next to -p
    run_cmd(pdb2gmx + ["-f", str(src),
             "-o", str(build / "protein.gro"),
             "-p", str(build / "topol.top"),
             "-ff", env["ff"], "-water", "spce", "-ignh"],
            log_file=build / "pdb2gmx.log", cwd=build)
    cof_acps = [_run_acpype(p, build, env, name=f"LIG_{p.stem[-1]}")
                for p in cof_pdfs]
    cof_gro_lines: list[str] = []
    if cof_acps:
        top_txt = _merge_tops((build / "topol.top").read_text(),
                              [a["itp"] for a in cof_acps], env["ff"])
        (build / "topol.top").write_text(top_txt)
        for a in cof_acps:
            if a["gro"]:
                lg = Path(a["gro"]).read_text().splitlines()
                n = int(lg[1])
                cof_gro_lines += lg[2:2 + n]
            top_txt = (build / "topol.top").read_text()
            (build / "topol.top").write_text(
                top_txt.rstrip() + f"\n{_itp_molname(a['itp'])} 1\n")
    if metals:
        dmap = _donor_index_map(build / "protein.gro", src)
        if not dmap:
            logger.warning("donor atom counts mismatch; metal "
                           "coordination restraints skipped")
        coords = metal_coordinators(complex_pdb, metals)
        for c in coords:
            c["_chain"] = metals[c["metal_idx"]]["chain"]
        _integrate_metals(build, metals, coords, dmap)
        combined = build / "combined.gro"
        atom_lines = _combined_gro_lines(build, build / "protein.gro",
                                         metals, cof_gro_lines)
        combined.write_text("\n".join(
            ["protein+metals+cofactors" if cof_gro_lines else "protein+metals",
             str(len(atom_lines)),
             *atom_lines, "100000.0"]) + "\n")
        _append_embedded_ions(build, emb, combined, build / "topol.top")
        return _finish_system(build, gmx, env, salt,
                              combined_gro=combined,
                              divalent=divalent,
                              divalent_m=divalent_m)
    if cof_gro_lines:
        # combined gro: protein + cofactor atoms
        prot_gro_lines = (build / "protein.gro").read_text().splitlines()
        n_prot = int(prot_gro_lines[1])
        combined = build / "combined.gro"
        combined.write_text("\n".join(
            [f"protein+cofactors\n{n_prot + len(cof_gro_lines)}",
             *prot_gro_lines[2:2 + n_prot], *cof_gro_lines,
             "100000.0"]) + "\n")
        _append_embedded_ions(build, emb, combined, build / "topol.top")
        return _finish_system(build, gmx, env, salt,
                              combined_gro=combined,
                              base_top=build / "topol.top",
                              divalent=divalent,
                              divalent_m=divalent_m)
    return _finish_system(build, gmx, env, salt,
                          divalent=divalent,
                          divalent_m=divalent_m)


# a chain is protein if the majority of its residue names are standard
# amino acids (4-letter C-terminal names like CPHE are truncated by the
# 3-column read, so compare against both forms; HOH/waters stay with the
# protein chains)
PROT_RES = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS",
            "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP",
            "TYR","VAL","MSE","CPH","CGL","CGLU","CAL","HOH"}
# nucleic-acid chains with at least this many residues are polymers and
# are parameterized natively by pdb2gmx (amber99sb-ildn carries b-DNA /
# RNA parameters); shorter NA pieces (single-residue nucleotide ligands
# like NAD/ATP) go to the ligand side for ACPYPE/GAFF2
NA_POLYMER_MIN_RES = 10


# non-protein/NA residues that pdb2gmx handles natively (amber99sb-ildn
# carries these); anything else multi-atom is a cofactor (HEM, NAD, ATP,
# sugars not in the FF, ...) that must go to ACPYPE/GAFF2
SAFE_COFACTOR_RES = {"HOH", "WAT", "SOL", "TIP3", "SPC", "GOL", "NAG",
                     "BGC", "NGL", "GLC", "MAN", "GAL", "FUC", "XYL",
                     "EDO", "SO4", "PO4", "PSE", "ACT", "PEG", "DOD",
                     "CL", "NA", "K", "MG", "ZN", "CA", "FE", "MN", "NI",
                     "CU", "CO", "ZNK", "CLB", "NAK", "MGM", "CAL"}


def _find_cofactors(pdb: Path) -> list[dict]:
    """R8: multi-atom non-standard residues that crash pdb2gmx when they
    share a chain with protein ('inconsistent type' fatal: protein +
    'Other' mixed in one chain). Single-atom metals are handled by
    find_metals; FF-native residues (sugars, common ions) are safe."""
    from collections import defaultdict
    res: "defaultdict[tuple, list]" = defaultdict(list)
    prot_chains: set[str] = set()
    with open(pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                name = line[12:16].strip()
                if not name or name == "H":
                    continue
                rname = line[17:20].strip()
                res[(line[21], int(line[22:26]), rname)].append(line)
                if rname in PROT_RES:
                    prot_chains.add(line[21])
    out = []
    for (chain, rseq, rname), lines in sorted(
            res.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        # only cofactors sharing a chain with protein crash pdb2gmx;
        # non-standard residues on their own chains are already routed
        # to the ligand side by the chain split
        if rname in PROT_RES or rname in DNA_RES or rname in METAL_RES \
                or rname in SAFE_COFACTOR_RES or len(lines) <= 2 \
                or chain not in prot_chains:
            continue
        out.append({"resname": rname, "chain": chain, "resseq": rseq,
                    "n_atoms": len(lines)})
    return out


def _embedded_metals(pdb: Path, cfs: list[dict]) -> list[dict]:
    """R8: metal atoms inside a cofactor residue (the FE of HEM). GROMACS
    metals are single-atom residues, but in a cofactor the metal is just
    an atom of a multi-atom residue — and Gasteiger (our only charge
    method without sqm) has no parameters for Fe. Returns find_metals-
    style dicts so the ions can be stripped before ACPYPE and re-added
    as standalone ion molecules."""
    if not cfs:
        return []
    cof_keys = {(c["chain"], c["resseq"]) for c in cfs}
    out = []
    with open(pdb) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            key = (line[21], int(line[22:26]))
            if key not in cof_keys:
                continue
            name = line[12:16].strip()
            if name not in METAL_PROPS:
                continue
            mol, charge, mass, _z = METAL_PROPS[name]
            out.append({
                "resname": name, "chain": key[0], "resseq": key[1],
                "xyz": (float(line[30:38]), float(line[38:46]),
                        float(line[46:54])),
                "n_atoms": 1, "molname": mol, "charge": charge,
                "mass": mass, "embedded": True})
    return out


def _reassign_cofactor_chains(pdb: Path, out: Path) -> list[dict]:
    """Rewrite a PDB copy giving every cofactor residue its own fresh
    chain letter (first free letter not used by the file). Downstream
    chain classification then routes cofactors to the ligand side
    (one ACPYPE molecule per residue) instead of crashing pdb2gmx."""
    cfs = _find_cofactors(pdb)
    if not cfs:
        out.write_text(pdb.read_text())
        return cfs
    used = set()
    with open(pdb) as fh:
        for l in fh:
            if l.startswith(("ATOM", "HETATM")):
                used.add(l[21])
    free = [c for c in (chr(i) for i in range(65, 90)) if c not in used]
    mapping: dict[tuple, str] = {}
    for c in cfs:
        mapping[(c["chain"], c["resseq"], c["resname"])] = free[len(mapping)]
    with open(pdb) as fi, open(out, "w") as fo:
        for l in fi:
            if l.startswith(("ATOM", "HETATM")):
                key = (l[21], int(l[22:26]), l[17:20].strip())
                if key in mapping:
                    l = l[:21] + mapping[key] + l[22:]
            fo.write(l)
    return cfs


def _classify_chains(complex_pdb: Path) -> set[str]:
    """chain ids that belong on the protein side of the split."""
    res_count: dict[str, dict[str, int]] = {}
    with open(complex_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                rc = res_count.setdefault(line[21], {})
                name = line[17:20].strip()
                rc[name] = rc.get(name, 0) + 1
    prot_chains = set()
    na_poly = {n["chain"] for n in find_nucleic_acids(complex_pdb)
               if n["n_res"] >= NA_POLYMER_MIN_RES}
    for ch, rc in res_count.items():
        total = sum(rc.values())
        if not total:
            continue
        prot = sum(n for nm, n in rc.items() if nm in PROT_RES)
        if prot * 2 >= total:
            prot_chains.add(ch)
            continue
        # nucleic-acid polymer chain (R4): native pdb2gmx parameterization
        if ch in na_poly:
            prot_chains.add(ch)
    return prot_chains


def _append_embedded_ions(build: Path, emb: list[dict], combined: Path,
                          top: Path) -> None:
    """R8: re-add cofactor-embedded metals (e.g. HEM's FE) as standalone
    ion molecules after the cofactor was stripped for ACPYPE (Gasteiger
    has no Fe parameters). Writes <el>ION itps (atomtypes only when the
    FF lacks the element, R3 pattern), appends molecule lines and gro
    atom lines. v1: no coordination restraints (standalone molecule —
    GROMACS restraints must live inside one moleculetype)."""
    if not emb:
        return
    els = sorted({m["resname"] for m in emb})
    # GROMACS 2023 is strict: once ANY moleculetype has been parsed, no
    # [ atomtypes ] may follow. So the atomtype-only part is included
    # right after the forcefield (before all itps with moleculetypes),
    # the moleculetype part after the last include.
    types_inc = []
    mol_inc = []
    for el in els:
        props = METAL_PROPS[el]
        charge, mass = props[2], props[2]  # charge, mass
        charge = props[1]
        fn = f"embion_{el.lower()}.itp"
        lines = [f"; standalone {el} ion from a cofactor (drugagent R8)",
                 "[ moleculetype ]",
                 "; molname       nrexcl",
                 f"{el}ION           1",
                 "",
                 "[ atoms ]",
                 ";   id  type   res  resid  atom   cgnr  charge  mass",
                 f"    1  {el:<4s}  1  {el}ION  {el:<4s}     1  "
                 f"{charge:>7s}  {mass:>7s}"]
        (build / fn).write_text("\n".join(lines) + "\n")
        mol_inc.append(f'#include "{fn}"')
        if not _ff_has_atomtype(build, el):
            tfn = f"embion_{el.lower()}_types.itp"
            znum = props[3]
            (build / tfn).write_text(
                f"; {el} atomtype only (drugagent R8)\n"
                f"[ atomtypes ]\n"
                f";  name  atomic number  mass  charge  ptype  "
                f"sigma  epsilon\n"
                f"{el:<5s}{znum:<15d}{mass:<8s}  0.0  A  0.0  0.0\n")
            types_inc.append(f'#include "{tfn}"')
    tlines = top.read_text().splitlines()
    first = last = None
    for i, l in enumerate(tlines):
        if l.startswith("#include"):
            if first is None:
                first = i
            last = i
    if last is not None:
        tlines = (tlines[:first + 1] + types_inc
                  + tlines[first + 1:last + 1] + mol_inc
                  + tlines[last + 1:])
    # one [ molecules ] line PER ION (two HEMs -> two FE ions; a
    # per-element line would undercount the topology)
    for m in emb:
        tlines.append(f"{m['resname']}ION 1")
    top.write_text("\n".join(tlines) + "\n")
    # gro atom lines (after the existing atoms; count line updated)
    g = combined.read_text().splitlines()
    n = int(g[1])
    box = g[-1]  # the trailing box-size line must stay last
    add = []
    for k, m in enumerate(emb):
        n += 1
        x, y, z = (v / 10.0 for v in m["xyz"])  # PDB A -> gro nm
        add.append(f"{n:5d}{m['resname'] + 'ION':>5s}{m['resname']:>5s}"
                   f"{m['resseq']:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    g[1] = str(n)
    g = g[:-1] + add + [box]
    combined.write_text("\n".join(g) + "\n")


def _ligand_side_pdfs(split_pdb: Path, prot_chains: set[str],
                      build: Path, emb: list[dict] | None = None) -> list[Path]:
    """R8: one PDB per non-protein chain (each file = one ACPYPE molecule:
    the docked ligand, cofactors re-assigned by _reassign_cofactor_chains,
    single-residue nucleotide ligands). Embedded cofactor metals (`emb`)
    are excluded — they return as standalone ion molecules."""
    emb = emb or []
    from collections import defaultdict
    by_chain: dict[str, list[str]] = defaultdict(list)
    header: list[str] = []
    # match on (resseq, atom name): the chain letter was rewritten by
    # _reassign_cofactor_chains before this file is produced
    skip = {(m["resseq"], m["resname"]) for m in emb}
    with open(split_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                key = (int(line[22:26]), line[12:16].strip())
                if key in skip:
                    continue
                by_chain[line[21]].append(line)
            elif line.startswith(("CRYST1", "REMARK")):
                header.append(line)
    out = []
    standalone: list[dict] = []
    for ch in sorted(by_chain):
        if ch in prot_chains:
            continue
        # R9: a chain whose residues are ONLY single-atom metal ions
        # (a docked Zn2+/Mg2+ "ligand") must not go through ACPYPE
        # (single-atom mol2s are fragile); it returns for
        # standalone-ion treatment
        resnames = {l[17:20].strip() for l in by_chain[ch]}
        if resnames and resnames <= set(METAL_RES):
            seen = set()
            for l in by_chain[ch]:
                key = (l[17:20].strip(), l[22:26].strip())
                if key in seen:
                    continue
                seen.add(key)
                standalone.append({
                    "resname": l[17:20].strip(),
                    "chain": l[21], "resseq": int(l[22:26]),
                    "xyz": (float(l[30:38]), float(l[38:46]),
                            float(l[46:54])),
                    "n_atoms": 1, "standalone_ligand": True})
            continue
        f = build / f"ligand_{ch}.pdb"
        f.write_text("".join(header) + "".join(by_chain[ch]) + "END\n")
        out.append(f)
    return out, standalone


def _build_ligand_system(complex_pdb: Path, build: Path, gmx: str,
                         env: dict, salt: float,
                         divalent: str | None = None,
                         divalent_m: float = 0.0) -> dict:
    """protein (amber99sb-ildn) + ligands/cofactors (ACPYPE/GAFF2)."""
    pdb2gmx = _find_pdb2gmx(env)
    # R8: cofactors (HEM, ...) sharing a chain with protein crash
    # pdb2gmx ('inconsistent type'); re-assign them to fresh chain
    # letters so the split routes them to the ligand side (ACPYPE)
    work_pdb = complex_pdb
    cfs = _find_cofactors(complex_pdb)
    emb = _embedded_metals(complex_pdb, cfs)
    if cfs:
        logger.info(f"cofactors on protein chains: "
                    f"{[(c['resname'], c['chain'], c['resseq']) for c in cfs]}"
                    " — moving to the ACPYPE side")
        work_pdb = build / "complex_cofactors.pdb"
        _reassign_cofactor_chains(complex_pdb, work_pdb)
    if emb:
        logger.info(f"embedded cofactor metals re-added as ions: "
                    f"{[(m['resname'], m['chain'], m['resseq']) for m in emb]}")
    # split
    # build_complex_pdb writes the ligand on the first free chain letter
    # (not "L" when the target is a dimer); the ligand is the only
    # component whose chain id does not occur in the target protein, so
    # identify it as "every chain that is not a protein chain"
    prot_pdb = build / "protein.pdb"
    prot_chains = _classify_chains(work_pdb)
    with open(work_pdb) as fh, open(prot_pdb, "w") as fp:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                # single-atom metal ions are their own molecules (R3)
                if line[17:20].strip() in METAL_RES:
                    continue
                if line[21] in prot_chains:
                    fp.write(line)
            elif line.startswith(("CRYST1", "REMARK")):
                fp.write(line)
    fp2 = open(prot_pdb, "a"); fp2.write("END\n"); fp2.close()
    # one PDB per non-protein chain (ligand + re-assigned cofactors);
    # each becomes its own ACPYPE molecule
    lig_pdfs, stand_metals = _ligand_side_pdfs(work_pdb, prot_chains,
                                                build, emb)
    if stand_metals:
        sm = [(m["resname"], m["chain"], m["resseq"])
              for m in stand_metals]
        logger.info(f"standalone metal ligands re-added as ions: {sm}")
    if not lig_pdfs and not stand_metals:
        raise RuntimeError("ligand chain not found in complex PDB")

    # cwd=build: pdb2gmx writes the posre_*.itp side files to the
    # current working directory, not next to -p
    run_cmd(pdb2gmx + ["-f", str(prot_pdb),
             "-o", str(build / "protein.gro"),
             "-p", str(build / "prot_top.top"),
             "-ff", env["ff"], "-water", "spce", "-ignh"],
            log_file=build / "pdb2gmx.log", cwd=build)

    acps = [_run_acpype(p, build, env, name=f"LIG_{p.stem[-1]}")
            for p in lig_pdfs]
    # combined topology: protein + all ligand/cofactor itps + ions
    prot_top = (build / "prot_top.top").read_text()
    top_txt = _merge_tops(prot_top, [a["itp"] for a in acps], env["ff"])
    (build / "topol.top").write_text(top_txt)

    # combined gro: protein atoms + all ligand atoms
    prot_gro_lines = (build / "protein.gro").read_text().splitlines()
    n_prot = int(prot_gro_lines[1])
    lig_lines: list[str] = []
    for a in acps:
        if a["gro"]:
            lg = Path(a["gro"]).read_text().splitlines()
            n = int(lg[1])
            lig_lines += lg[2:2 + n]
    combined = build / "combined.gro"
    with open(combined, "w") as fh:
        fh.write(f"Combined protein+ligands\n{n_prot + len(lig_lines)}\n")
        fh.write("\n".join(prot_gro_lines[2:2 + n_prot]) + "\n")
        fh.write("\n".join(lig_lines) + "\n")
        fh.write("100000.0\n")
    # append ligand molecules to [ molecules ] (names from the itps)
    top_txt2 = (build / "topol.top").read_text()
    for a in acps:
        top_txt2 = top_txt2.rstrip() + f"\n{_itp_molname(a['itp'])} 1\n"
    (build / "topol.top").write_text(top_txt2)

    # metal ions (R3): appended to the owning protein chain molecule +
# [ distance_restraints ] (GROMACS restraint sections parse in the
# context of the preceding moleculetype)
    metals = [m for m in find_metals(complex_pdb) if m["n_atoms"] == 1]
    if metals:
        logger.info(f"metal ions in system: "
                    f"{[(m['resname'], m['chain'], m['resseq']) for m in metals]}")
        dmap = _donor_index_map(build / "protein.gro", prot_pdb)
        if not dmap:
            logger.warning("donor atom counts mismatch; metal "
                           "coordination restraints skipped")
        coords = metal_coordinators(complex_pdb, metals)
        for c in coords:
            c["_chain"] = metals[c["metal_idx"]]["chain"]
        _integrate_metals(build, metals, coords, dmap)
        # rebuild combined.gro in topology order (chain atoms + that
        # chain's metals, then ligand)
        atom_lines = _combined_gro_lines(build, build / "protein.gro",
                                         metals, lig_lines)
        combined.write_text("\n".join(
            ["Combined protein+metals+ligand", str(len(atom_lines)),
             *atom_lines, "100000.0"]) + "\n")

    _append_embedded_ions(build, emb, combined, build / "topol.top")
    # solvate + genion + EM (topology merge done in _finish_system)
    return _finish_system(build, gmx, env, salt,
                          start_gro=None, combined_gro=combined,
                          base_top=build / "topol.top",
                          divalent=divalent,
                          divalent_m=divalent_m)


def _itp_molname(itp: str | Path) -> str:
    """moleculetype name declared in an ACPYPE itp (first token of the
    first data line after [ moleculetype ])."""
    lines = Path(itp).read_text().splitlines()
    for i, l in enumerate(lines):
        if l.strip().startswith("[ moleculetype ]"):
            for l2 in lines[i + 1:]:
                if l2.strip() and not l2.lstrip().startswith(";"):
                    return l2.split()[0]
    raise RuntimeError(f"no moleculetype in {itp}")


def _merge_tops(prot_top: str, lig_itps: list[str], ff: str = "amber99sb-ildn") -> str:
    """Merge a pdb2gmx protein topology with ligand ITPs (ACPYPE, self-contained)."""
    lines = prot_top.splitlines()
    out: list[str] = []
    in_mol = False
    molecules: list[str] = []
    for line in lines:
        if line.startswith("[ molecules ]"):
            in_mol = True
            continue
        if in_mol:
            if line.strip():
                molecules.append(line)
            continue
        # cut the forcefield include (we rebuild)
        if line.startswith("#include") and ("forcefield" in line or line.strip() == '#include "spce.itp"'):
            continue
        out.append(line)
    head = "\n".join(out).rstrip()
    header = f'#include "{ff}.ff/forcefield.itp"'
    lig_incs = "".join(f'#include "{p}"\n' for p in lig_itps)
    return (header + "\n" + lig_incs + head
            + "\n\n[ molecules ]\n" + "\n".join(molecules) + "\n")


def _gro_box_volume(gro: str | Path) -> float:
    """Box volume (nm^3) from a GRO's last line (9 box vectors;
    GROMACS writes them transposed: v1x v2x v3x v1y v2y v3y ...)."""
    last = Path(gro).read_text().splitlines()[-1].split()
    v = [float(x) for x in last]
    if len(v) >= 9:
        # GROMACS writes the 3x3 box matrix ROW-major:
        # v1x v1y v1z v2x v2y v2z v3x v3y v3z (v1/v2/v3 = box vectors)
        a = (v[0], v[1], v[2])
        b = (v[3], v[4], v[5])
        c = (v[6], v[7], v[8])
        det = (a[0] * (b[1] * c[2] - b[2] * c[1])
               - a[1] * (b[0] * c[2] - b[2] * c[0])
               + a[2] * (b[0] * c[1] - b[1] * c[0]))
        return abs(det)
    # minimal gro: only the 3 diagonal lengths
    if len(v) >= 3:
        return abs(v[0] * v[1] * v[2])
    return 0.0


def _ion_pair_count(molarity: float, volume_nm3: float) -> int:
    """R10: number of ions for a molarity in a box volume
    (N = c [mol/L] * V [L]; 1 nm^3 = 1e-27 L; N_A = 6.022e23)."""
    if molarity <= 0 or volume_nm3 <= 0:
        return 0
    # 1 L = 1e27 nm^3, so N = c * V[nm^3] * 1e-27 L/nm^3 * N_A
    # = c * V * 6.022e-4 * 1e3 (the extra 1e3 is L per m^3:
    # 1 nm^3 = 1e-21 m^3 = 1e-24 L... worked out: 0.1 M in 1000 nm^3
    # is 9.0e-3 mol? no — 0.1 mol/L * 1e-24 L = 1e-25 mol * N_A =
    # 60 ions per 1000 nm^3 at 0.1 M; the factor 6.022e-4 * 1e3 =
    # 0.6022 per (M*nm^3) is the correct empirical constant here)
    return int(round(molarity * volume_nm3 * 6.022e-4 * 1000))


def _divalent_ion_itp(el: str) -> str:
    """R10: standalone divalent-cation ion for genion (MG/CA ...):
    atomtype (GROMACS ion.itp only has NA/CL/K) + 1-atom molecule.
    The cation's charge comes from METAL_PROPS (divalent metals are
    +2 in our usage)."""
    charge, mass = METAL_PROPS[el][1], METAL_PROPS[el][2]
    znum = METAL_PROPS[el][3]
    return ("; divalent counterion (drugagent R10)\n"
            "[ atomtypes ]\n"
            ";  name  atomic number  mass  charge  ptype  sigma  epsilon\n"
            f"{el:<5s}{znum:<15d}{mass:<8s}  0.0  A  0.0  0.0\n"
            "[ moleculetype ]\n; molname       nrexcl\n"
            f"{el}ION           1\n"
            "[ atoms ]\n"
            ";   id  type   res  resid  atom   cgnr  charge  mass\n"
            f"    1  {el:<4s}  1  {el}ION  {el:<4s}     1  "
            f"{charge:>7s}  {mass:>7s}\n")


def _genion_cmd(gmx: str, min_tpr: Path, solv_top: Path, ions_gro: Path,
                salt: float, volume_nm3: float,
                divalent: str | None = None,
                divalent_m: float = 0.0) -> list[str]:
    """R10: genion invocation. Previously genion ran with ONLY -neutral
    and the `salt` argument was silently ignored (md_salt_m was a no-op).
    Now: -neutral (counterions for the net charge) + salt pairs + an
    optional divalent cation with its Cl- balance (2 Cl- per divalent
    cation keeps the added ions charge-neutral)."""
    cmd = [gmx, "genion", "-s", str(min_tpr), "-p", str(solv_top),
           "-o", str(ions_gro), "-neutral"]
    n_pairs = _ion_pair_count(salt, volume_nm3)
    if n_pairs:
        cmd += ["-np", str(n_pairs)]
    if divalent:
        n_div = _ion_pair_count(divalent_m, volume_nm3)
        if n_div:
            cmd += ["-pname", f"{divalent}ION", "-pn", str(n_div),
                    "-nname", "CL", "-nn", str(2 * n_div)]
    return cmd


def _finish_system(build: Path, gmx: str, env: dict, salt: float,
                   start_gro: Path | None = None,
                   combined_gro: Path | None = None,
                   base_top: Path | None = None,
                   divalent: str | None = None,
                   divalent_m: float = 0.0) -> dict:
    # idempotent: if EM already finished, reuse; otherwise drop partial
    # state (genion rewrites solvated.top, so a stale solvated.gro from a
    # previous attempt no longer matches it)
    solv_top = build / "solvated.top"
    base_top = base_top or (build / "topol.top")
    solv_gro = start_gro or (build / "solvated.gro")
    if (build / "em.cpt").is_file() and (build / "em.tpr").is_file():
        return {"em_tpr": str(build / "em.tpr"),
                "top": str(solv_top),
                "gro": str(build / "ions.gro"),
                "build_dir": str(build)}
    for stale in ("solvated.gro", "solvated.top", "min.tpr", "ions.gro",
                  "em.tpr", "em.cpt", "em.log", "em.edr", "em.gro", "em.xtc"):
        (build / stale).unlink(missing_ok=True)
    base = combined_gro or (build / "protein.gro")
    # GROMACS >=5 command-line parser requires -p to exist first
    solv_top.write_text("")
    run_cmd([gmx, "solvate", "-cp", str(base),
             "-cs", _top_spce(env), "-p", str(solv_top),
             "-o", str(solv_gro)], log_file=build / "solvate.log")
    _merge_solvate_top(base_top, solv_top)
    # R10: divalent counterion (MG/CA): the ion must exist in the
    # topology BEFORE genion — atomtype (GROMACS ion.itp only ships
    # NA/CL/K) + standalone molecule, same include pattern as the R8
    # embedded ions
    if divalent:
        el = divalent.upper()
        itptxt = _divalent_ion_itp(el)
        (build / f"div_{el.lower()}.itp").write_text(itptxt)
        tlines = (build / "solvated.top").read_text().splitlines()
        for i, l in enumerate(tlines):
            if l.startswith("#include"):
                tlines = (tlines[:i + 1]
                          + [f'#include "div_{el.lower()}.itp"']
                          + tlines[i + 1:])
                break
        tlines.append(f"{el}ION 1")
        (build / "solvated.top").write_text("\n".join(tlines) + "\n")
    # tpr for genion
    min_tpr = build / "min.tpr"
    (build / "em0.mdp").write_text(_mdp_em())
    run_cmd([gmx, "grompp", "-f", str(build / "em0.mdp"),
             "-c", str(solv_gro), "-p", str(build / "solvated.top"),
             "-o", str(min_tpr)], log_file=build / "grompp_em0.log")
    ions_gro = build / "ions.gro"
    # genion (GROMACS >=5) interactively asks for the solvent group;
    # answer "SOL" on stdin. R10: salt pairs + divalent cation are
    # actually applied now (previously -neutral only, salt ignored)
    volume = _gro_box_volume(solv_gro)
    gcmd = _genion_cmd(gmx, min_tpr, build / "solvated.top", ions_gro,
                       salt, volume, divalent, divalent_m)
    run_cmd(gcmd, log_file=build / "genion.log", stdin="SOL\n")
    # EM
    em_tpr = build / "em.tpr"
    (build / "em.mdp").write_text(_mdp_em())
    run_cmd([gmx, "grompp", "-f", str(build / "em.mdp"),
             "-c", str(ions_gro), "-p", str(build / "solvated.top"),
             "-o", str(em_tpr)], log_file=build / "grompp_em.log")
    # GROMACS 2023: -ntomp and OMP_NUM_THREADS must agree; override the
    # env var for this process to match the requested thread count
    omp = str(max(4, n_cores() // 4))
    # thread_mpi build: distance-restraint init segfaults with multiple
    # thread ranks (same reason the production MD run uses -ntmpi 1)
    run_cmd([gmx, "mdrun", "-deffnm", str(build / "em"),
             "-ntmpi", "1", "-ntomp", omp],
            log_file=build / "em_mdrun.log",
            env=dict(os.environ, OMP_NUM_THREADS=omp))
    return {"em_tpr": str(em_tpr), "top": str(build / "solvated.top"),
            "gro": str(ions_gro), "build_dir": str(build)}


def _merge_solvate_top(base_top: Path, solv_top: Path) -> None:
    """Combine the pdb2gmx/merged base topology with the solvent line.

    ``gmx solvate -p`` (GROMACS >=5) only ever writes the SOL molecule
    line, so the full system topology has to be assembled by hand.
    """
    sol_lines = solv_top.read_text().splitlines()
    sol_mol = next((l for l in sol_lines if l.strip().startswith("SOL ")), None)
    keep: list[str] = []
    molecules: list[str] = []
    in_mol = False
    for l in base_top.read_text().splitlines():
        if l.startswith("[ molecules ]"):
            in_mol = True
            continue
        if in_mol:
            if l.strip():
                molecules.append(l)
            continue
        keep.append(l)
    header = "\n".join(keep).rstrip()
    out = header + "\n\n[ molecules ]\n" + "\n".join(molecules)
    if sol_mol:
        out += "\n" + sol_mol
    solv_top.write_text(out + "\n")


def _mdp_em() -> str:
    return """integrator  = steep
nsteps      = 10000
nstenergy   = 100
nstlog      = 1000
cutoff-scheme = Verlet
nstlist     = 10
rlist       = 1.0
disre       = Simple
disre-fc    = 2000
nstdisreout = 0
"""


def _mdp_md(ns: float) -> str:
    nsteps = int(ns * 1000 / 0.002)  # ps -> steps of 2fs
    return f"""integrator  = md
nsteps      = {nsteps}
dt          = 0.002
nstcomm     = 100
nstenergy   = 1000
nstlog      = 1000
nstxout     = 5000
nstxout-compressed = 5000
cutoff-scheme = Verlet
nstlist     = 10
rlist       = 1.0
coulombtype = PME
rcoulomb    = 1.0
fourierspacing = 0.16
disre       = Simple
disre-fc    = 2000
nstdisreout = 0
pcoupltype  = isotropic
tcoupl      = v-rescale
; "System" (not "Protein"): the default Protein group is derived from
; moleculetype names, which misses DNA_chain_* / ligand-only systems
tc-grps     = System
tau_t       = 0.02
ref_t       = 300
gen_vel     = yes
; R8: C-rescale (not Berendsen): Berendsen segfaulted the 1HVI smoke
; system ~20 ps into production on both replicas (PLOESS pressure
; scaling > 1% at the crash); C-rescale + tau_p 2 ps is robust
pcoupl      = C-rescale
ref_p       = 1.0
tau_p       = 2.0
compressibility = 4.5e-5
constraints = h-bonds
; dispersion correction
DispCorr    = EnerPres
"""


def _mdp_eq_nvt(ps: float) -> str:
    """NVT equilibration: heat to 300 K, position restraints on the
    protein (posre-fc per tc-group: Protein 1000, Non-Protein 0)."""
    nsteps = int(ps / 0.002)  # R8: ps of 2fs steps (1 ps = 500 steps;
    # the old ps*1000/0.002 was 1000x too long: 10 ps NVT ran 10 ns)
    return f"""integrator  = md
nsteps      = {nsteps}
dt          = 0.002
nstcomm     = 100
nstenergy   = 100
nstlog      = 1000
nstxout     = 25000
nstxout-compressed = 25000
cutoff-scheme = Verlet
nstlist     = 10
rlist       = 1.0
coulombtype = PME
rcoulomb    = 1.0
fourierspacing = 0.16
disre       = Simple
disre-fc    = 2000
nstdisreout = 0
pcoupltype  = isotropic
tcoupl      = v-rescale
; "System" (not "Protein"): the default Protein group is derived from
; moleculetype names, which misses DNA_chain_* / ligand-only systems
tc-grps     = System
tau_t       = 0.02
ref_t       = 300
gen_vel     = yes
; per-atom restraint force constants live in the posre_*.itp includes
constraints = h-bonds
"""


def _mdp_eq_npt(ps: float) -> str:
    """NPT equilibration: 300 K / 1 bar, protein restraints, box relax."""
    nsteps = int(ps / 0.002)  # R8: ps of 2fs steps (1 ps = 500 steps)
    return f"""integrator  = md
nsteps      = {nsteps}
dt          = 0.002
nstcomm     = 100
nstenergy   = 100
nstlog      = 1000
nstxout     = 25000
nstxout-compressed = 25000
cutoff-scheme = Verlet
nstlist     = 10
rlist       = 1.0
coulombtype = PME
rcoulomb    = 1.0
fourierspacing = 0.16
disre       = Simple
disre-fc    = 2000
nstdisreout = 0
pcoupltype  = isotropic
tcoupl      = v-rescale
; "System" (not "Protein"): the default Protein group is derived from
; moleculetype names, which misses DNA_chain_* / ligand-only systems
tc-grps     = System
tau_t       = 0.02
ref_t       = 300
gen_vel     = no
; per-atom restraint force constants live in the posre_*.itp includes
; R8: C-rescale (not Berendsen) — see production template
pcoupl      = C-rescale
; position restraints + pressure coupling: scale restraint reference
; coordinates with the box (otherwise grompp warns about artifacts)
refcoord-scaling = all
ref_p       = 1.0
tau_p       = 2.0
compressibility = 4.5e-5
constraints = h-bonds
"""


def _write_eq_top(build: Path) -> Path:
    """eq.top = solvated.top + posre_*.itp includes (position restraints
    during equilibration only; the production top stays restraint-free).

    pdb2gmx posre itps have no [ moleculetype ] header of their own —
    their [ position_restraints ] section is parsed in the context of the
    immediately preceding moleculetype — so each posre include must go
    directly after the include of the chain itp that defines the
    molecule (posre_DNA_chain_A.itp after topol_DNA_chain_A.itp)."""
    base = (build / "solvated.top").read_text()
    posres = sorted(f.name for f in build.glob("posre_*.itp"))
    if not posres:
        logger.warning("no posre_*.itp found in build dir; "
                       "equilibration runs without position restraints")
    lines = base.splitlines()
    out_lines = []
    for l in lines:
        out_lines.append(l)
        for p in posres:
            # posre_<mol>.itp follows the include that defines <mol>.
            # The defining include's prefix varies by build path
            # (topol_<mol>.itp for plain pdb2gmx, prot_top_<mol>.itp for
            # merged protein topologies), so match on the suffix.
            if not l.startswith('#include "'):
                continue
            inc = l[len('#include "'):-len('"')]
            if inc.endswith(f"_{p[len('posre_'):-len('.itp')]}.itp") \
                    and not inc.startswith("posre_"):
                out_lines.append(f'#include "{p}"')
    # any posre whose defining include was not found: after the last
    # include before [ system ] (fallback, with a warning)
    placed = {l[len('#include "'):-len('"')] for l in out_lines
              if l.startswith('#include "posre_')}
    leftover = [p for p in posres if p not in placed]
    if leftover:
        logger.warning(f"posre includes placed at end of topology "
                       f"(defining include not found): {leftover}")
        idx = 0
        for i, l in enumerate(out_lines):
            if l.startswith("#include"):
                idx = i
        out_lines = out_lines[:idx + 1] +             [f'#include "{p}"' for p in leftover] + out_lines[idx + 1:]
    out = build / "eq.top"
    out.write_text("\n".join(out_lines) + "\n")
    return out


def _mdp_fingerprint(mdp_text: str) -> str:
    """R9: hash of the *meaningful* MDP content (comments and blank
    lines stripped) — the MDP is part of a run's identity, not just its
    output artifacts."""
    import hashlib
    lines = [l.strip() for l in mdp_text.splitlines()
             if l.strip() and not l.strip().startswith(";")]
    return hashlib.sha1("|".join(lines).encode()).hexdigest()[:12]


def _gro_box_fingerprint(gro: str | Path) -> str:
    """R10: fingerprint of a GRO's box vector line (last line) — the
    equilibrated BOX is part of the production start state: reusing a
    production replica started from a different box (e.g. after the
    eq MDP/barostat changed and the eq stage re-ran) silently mixes
    ensembles."""
    import hashlib
    p = Path(gro)
    if not p.is_file():
        return "no-gro"
    try:
        last = p.read_text().splitlines()[-1]
    except Exception:
        return "no-gro"
    return hashlib.sha1(last.strip().encode()).hexdigest()[:12]


def _eq_stage_done(build: Path, name: str, mdp_fp: str) -> bool:
    """R9: an equilibration stage is reusable only when its artifacts
    exist AND the MDP fingerprint matches the one that produced them —
    otherwise a template change (e.g. barostat Berendsen -> C-rescale)
    is silently ignored and the old MDP keeps running."""
    log = Path(build) / f"{name}.log"
    fpf = Path(build) / f"{name}.mdp.fp"
    return ((Path(build) / f"{name}.gro").is_file()
            and log.is_file()
            and "Finished mdrun" in log.read_text(errors="ignore")
            and fpf.is_file()
            and fpf.read_text().strip() == mdp_fp)


def _b_args(burn_in_ps: float) -> list[str]:
    """gmx -b flag for trimming the production burn-in (ps -> ns).
    :g format — .1f would round 0.5 ps (0.0005 ns) to 0.0."""
    if burn_in_ps > 0:
        return ["-b", f"{burn_in_ps / 1000.0:g}"]
    return []


# --------------------------------------------------------------------------- #
# MD replicas
# --------------------------------------------------------------------------- #
def run_equilibration(sysinfo: dict, workdir: Path, env: dict,
                      *, nvt_ps: float = 50.0, npt_ps: float = 100.0) -> dict:
    """R5: NVT -> NPT equilibration (once per build; all replicas fork
    from the equilibrated state). Artifacts live in the build dir.

    Returns {"eq_gro", "eq_cpt", "eq_top"} (paths)."""
    build = Path(sysinfo["build_dir"])
    gmx = env["gmx"]
    omp = str(max(4, n_cores() // 4))
    eqtop = _write_eq_top(build)
    # idempotent: reuse a finished NPT equilibration (MDP must match)
    npt_fp = _mdp_fingerprint(_mdp_eq_npt(npt_ps))
    if _eq_stage_done(build, "eq_npt", npt_fp):
        logger.info("equilibration already complete; reusing")
        return {"eq_gro": str(build / "eq_npt.gro"),
                "eq_cpt": str(build / "eq_npt.cpt"),
                "eq_top": str(eqtop),
                "nvt_ps": None, "npt_ps": None}
    # NVT (heating, protein restraints) — idempotent: a finished NVT
    # (eq_nvt.gro + "Finished mdrun" + matching MDP fingerprint) is
    # reused; its checkpoint feeds the NPT stage below
    nvt_fp = _mdp_fingerprint(_mdp_eq_nvt(nvt_ps))
    nvt_done = _eq_stage_done(build, "eq_nvt", nvt_fp)
    if nvt_done:
        logger.info("NVT equilibration already complete; reusing")
    else:
        (build / "eq_nvt.mdp").write_text(_mdp_eq_nvt(nvt_ps))
        (build / "eq_nvt.mdp.fp").write_text(nvt_fp)
        # -r: position-restraint reference coordinates (same as -c)
        run_cmd([gmx, "grompp", "-f", str(build / "eq_nvt.mdp"),
                 "-c", str(build / "em.gro"), "-r", str(build / "em.gro"),
                 "-p", str(eqtop), "-o", str(build / "eq_nvt.tpr"),
                 "-maxwarn", "1"],
                log_file=build / "grompp_eq_nvt.log", env=dict(os.environ,
                OMP_NUM_THREADS=omp))
        run_cmd([gmx, "mdrun", "-deffnm", str(build / "eq_nvt"),
                 "-ntmpi", "1", "-ntomp", omp],
                log_file=build / "eq_nvt.log",
                env=dict(os.environ, OMP_NUM_THREADS=omp))
    # NPT (box relaxation, still restraints) — continues from the NVT state
    (build / "eq_npt.mdp").write_text(_mdp_eq_npt(npt_ps))
    (build / "eq_npt.mdp.fp").write_text(npt_fp)
    run_cmd([gmx, "grompp", "-f", str(build / "eq_npt.mdp"),
             "-c", str(build / "eq_nvt.gro"),
             "-r", str(build / "em.gro"),
             "-t", str(build / "eq_nvt.cpt"),
             "-p", str(eqtop), "-o", str(build / "eq_npt.tpr"),
             "-maxwarn", "1"],
            log_file=build / "grompp_eq_npt.log", env=dict(os.environ,
            OMP_NUM_THREADS=omp))
    run_cmd([gmx, "mdrun", "-deffnm", str(build / "eq_npt"),
             "-ntmpi", "1", "-ntomp", omp],
            log_file=build / "eq_npt.log",
            env=dict(os.environ, OMP_NUM_THREADS=omp))
    return {"eq_gro": str(build / "eq_npt.gro"),
            "eq_cpt": str(build / "eq_npt.cpt"),
            "eq_top": str(eqtop),
            "nvt_ps": nvt_ps, "npt_ps": npt_ps}


def _md_converged(summary: dict, length_ps: float, *,
                  min_len_ps: float = 50.0, drift_nm: float = 0.05,
                  min_cluster: float = 0.5) -> tuple[bool, str]:
    """R6 heuristic convergence gate for the auto-extension loop:

    (1) the trajectory must be at least `min_len_ps` long;
    (2) the RMSD (mean across replicas) must have plateaued — the mean of
        the last 20% of the series must not deviate from the mean of the
        preceding 80% by more than `drift_nm`;
    (3) clustering must show a dominant state (largest population >=
        `min_cluster`) — a system still hopping between states is not
        converged even if the RMSD mean is flat.

    Returns (converged, reason) — reason doubles as the log line."""
    if length_ps < min_len_ps:
        return False, (f"trajectory {length_ps:.0f} ps below minimum "
                       f"{min_len_ps:.0f} ps")
    mean = (summary.get("rmsd") or {}).get("mean") or []
    if len(mean) < 10:
        return False, f"only {len(mean)} RMSD samples"
    arr = np.asarray(mean, dtype=float)
    k = max(1, int(len(arr) * 0.2))
    tail, head = arr[-k:], arr[:-k]
    drift = abs(float(tail.mean()) - float(head.mean()))
    if drift >= drift_nm:
        return False, (f"RMSD still drifting (tail-head {drift * 10:.1f} A "
                       f"> {drift_nm * 10:.1f} A)")
    cl = summary.get("clusters") or {}
    largest = max((float(v) for v in cl.values()), default=0.0)
    if largest < min_cluster:
        return False, (f"no dominant cluster (largest {largest * 100:.0f}% "
                       f"< {min_cluster * 100:.0f}%)")
    return True, (f"RMSD plateaued (drift {drift * 10:.1f} A) with a "
                 f"dominant cluster ({largest * 100:.0f}%)")


def _ext_mdp_text(base_mdp: str, total_steps: int) -> str:
    """Production MDP text for a checkpoint continuation: nsteps is the
    TOTAL step count (original + extension) because mdrun -cpi reads the
    tpr nsteps as the run's end step; gen_vel is forced off so the
    continuation inherits the checkpoint velocities (no re-thermalization
    kick at the join)."""
    import re as _re
    out = _re.sub(r"^nsteps\s*=\s*\d+",
                  f"nsteps      = {int(total_steps)}",
                  base_mdp, flags=_re.M)
    if _re.search(r"^gen_vel\s*=\s*yes", out, flags=_re.M):
        out = _re.sub(r"^gen_vel\s*=\s*yes", "gen_vel     = no",
                      out, flags=_re.M)
    return out


def _latest_part_file(pattern: str, rdir: Path, kind: str) -> Path | None:
    """Highest-part-numbered file matching pattern (mdrun -cpi writes
    <deffnm>.part000N.<ext> continuation files)."""
    import re as _re
    cands = []
    for f in rdir.glob(pattern):
        m = _re.search(r"part(\d+)", f.name)
        n = int(m.group(1)) if m else 0
        cands.append((n, f.stat().st_mtime, f))
    if not cands:
        return None
    cands.sort()
    return cands[-1][2]


def _replica_job(args: tuple) -> dict:
    """One production replica (module-level on purpose: closures drag
    module globals through cloudpickle, and under pytest's stream capture
    the loguru logger holds a non-picklable captured stream — a plain
    data payload + importable function pickles cleanly)."""
    (i, workdir, md_mdp, gmx, cores_per, top, start_gro, start_cpt,
     start_box_fp) = args
    rdir = Path(workdir) / f"md_rep{i}"
    rdir.mkdir(parents=True, exist_ok=True)
    # idempotent: if this replica already ran to completion (trajectory
    # present + mdrun "Finished" in the log) AND its start state is
    # unchanged (R10: the equilibrated box fingerprint matches), skip
    # the expensive re-run.
    xtc = rdir / "md.xtc"
    mlog = rdir / "md.log"
    box_fp_file = rdir / "md_start_box.fp"
    done = (xtc.is_file() and mlog.is_file()
            and "Finished mdrun" in mlog.read_text(errors="ignore"))
    fp_ok = (box_fp_file.is_file()
             and box_fp_file.read_text().strip() == start_box_fp)
    logger.debug(f"replica {i} idempotency: done={done} fp_ok={fp_ok} "
                 f"xtc={xtc.is_file()} mlog={mlog.is_file()}")
    if done and fp_ok:
        logger.info(f"replica {i} already complete; skipping mdrun")
        return {"rep": i, "dir": str(rdir), "wall_h": 0.0,
                "reused": True}
    if done:
        logger.info(f"replica {i} start state changed (equilibrated box "
                    f"fingerprint mismatch); re-running")
    (rdir / "md.mdp").write_text(md_mdp)
    # Start MD from the equilibrated state (R5) when available,
    # otherwise from the EM-relaxed coordinates (em.gro). With a start
    # .cpt, grompp reads box + velocities from it (gen_vel in the mdp
    # is then ignored); without one, velocities are regenerated.
    # -maxwarn 1: the Berendsen barostat emits a "not a correct
    # ensemble" opinion warning that is harmless for short MD.
    grompp_cmd = [gmx, "grompp", "-f", str(rdir / "md.mdp"),
                  "-c", str(start_gro),
                  "-p", top, "-o", str(rdir / "md.tpr"),
                  "-maxwarn", "1"]
    if start_cpt is not None:
        grompp_cmd += ["-t", str(start_cpt)]
    (box_fp_file).write_text(start_box_fp)
    run_cmd(grompp_cmd, log_file=rdir / "grompp.log")
    t0 = time.time()
    run_cmd([gmx, "mdrun", "-deffnm", str(rdir / "md"),
             "-ntmpi", "1", "-ntomp", str(cores_per)],
            log_file=rdir / "mdrun.log",
            env=dict(os.environ, OMP_NUM_THREADS=str(cores_per)))
    return {"rep": i, "dir": str(rdir),
            "wall_h": round((time.time() - t0) / 3600, 2)}


def run_replicas(sysinfo: dict, workdir: Path, env: dict, *, ns: float,
                 reps: int, n_jobs: int | None = None,
                 save_ps: float = 10.0, eq: dict | None = None) -> list[dict]:
    """Production MD replicas. With `eq` (from run_equilibration) the
    replicas start from the equilibrated coordinates + velocities;
    otherwise from the EM output (pre-R5 behavior)."""
    build = Path(sysinfo["build_dir"])
    nstxout = max(1, int(save_ps * 1000 / 2))  # 2fs steps
    md_mdp = _mdp_md(ns)
    md_mdp = md_mdp.replace("nstxout     = 5000", f"nstxout     = {nstxout}")
    md_mdp = md_mdp.replace("nstxout-compressed = 5000",
                            f"nstxout-compressed = {nstxout}")
    gmx = env["gmx"]
    cores_total = n_jobs or n_cores()
    # cap at 32: for small systems (few k atoms) more OMP threads than
    # this oversubscribes and slows mdrun down (measured: 260 vs
    # 1300 steps/s for a 7 k-atom box at 64 vs 16 threads)
    cores_per = max(2, min(cores_total // reps, 32))
    start_gro = Path(eq["eq_gro"]) if eq else build / "em.gro"
    start_cpt = Path(eq["eq_cpt"]) if eq else None
    if start_cpt is not None:
        # with a start .cpt, grompp reads box + coordinates from it;
        # gen_vel=no makes it also use the equilibrated velocities
        # (with gen_vel=yes the -t velocities are ignored)
        md_mdp = md_mdp.replace("gen_vel     = yes", "gen_vel     = no")
    # R10: the equilibrated BOX is part of the start state identity
    box_src = start_cpt if start_cpt is not None else start_gro
    start_box_fp = _gro_box_fingerprint(box_src)
    payloads = [(i + 1, str(workdir), md_mdp, gmx, cores_per,
                 sysinfo["top"], str(start_gro),
                 str(start_cpt) if start_cpt is not None else None,
                 start_box_fp)
                for i in range(reps)]
    return pmap(_replica_job, payloads, n_jobs=reps)


def _extend_job(args: tuple) -> dict:
    """One replica extension round (module-level: picklable — see
    _replica_job for why closures are avoided here)."""
    (r, top, ext_ns, gmx, cores) = args
    import re as _re
    rdir = Path(r["dir"])
    round_i = int(r.get("ext_round", 0)) + 1
    deffnm = rdir / f"md_ext{round_i}"
    # current end state
    if round_i == 1:
        end_cpt = rdir / "md.cpt"
        base_steps = int(_re.search(r"^nsteps\s*=\s*(\d+)",
                                    (rdir / "md.mdp").read_text(),
                                    _re.M).group(1))
    else:
        end_cpt = Path(r["end_cpt"])
        base_steps = int(r["end_steps"])
    if not end_cpt.is_file():
        raise RuntimeError(f"end checkpoint missing: {end_cpt}")
    ext_steps = int(ext_ns * 1000 / 0.002)
    total_steps = base_steps + ext_steps
    (deffnm.with_suffix(".mdp")).write_text(
        _ext_mdp_text((rdir / "md.mdp").read_text(), total_steps))
    run_cmd([gmx, "grompp", "-f", str(deffnm.with_suffix(".mdp")),
             "-c", str(rdir / "md.gro"), "-t", str(end_cpt),
             "-p", top,
             "-o", str(deffnm.with_suffix(".tpr")), "-maxwarn", "1"],
            log_file=rdir / f"grompp_ext{round_i}.log")
    run_cmd([gmx, "mdrun", "-deffnm", str(deffnm),
             "-ntmpi", "1", "-ntomp", str(cores),
             "-cpi", str(end_cpt), "-noappend"],
            log_file=rdir / f"mdrun_ext{round_i}.log",
            env=dict(os.environ, OMP_NUM_THREADS=str(cores)))
    new_xtc = _latest_part_file(f"{deffnm.name}.part*.xtc", rdir, "xtc")
    if new_xtc is None:
        new_xtc = deffnm.with_suffix(".xtc")
    if not new_xtc.is_file():
        raise RuntimeError(f"extension xtc not found for {deffnm.name}")
    # merge with the current (possibly already extended) trajectory
    cur = Path(r.get("merged_xtc") or (rdir / "md.xtc"))
    merged = rdir / "md_all.xtc"
    run_cmd([gmx, "trjcat", "-f", str(cur), str(new_xtc),
             "-o", str(merged)], log_file=rdir / f"trjcat{round_i}.log")
    # new end state: the continuation run's final checkpoint
    new_cpt = _latest_part_file(f"{deffnm.name}.part*.cpt", rdir, "cpt")
    if new_cpt is None:
        new_cpt = deffnm.with_suffix(".cpt")
    r["ext_round"] = round_i
    r["end_cpt"] = str(new_cpt)
    r["end_steps"] = total_steps
    r["merged_xtc"] = str(merged)
    r["extended"] = True
    return r


def extend_replicas(reps: list[dict], sysinfo: dict, workdir: Path,
                    env: dict, *, ext_ns: float,
                    n_jobs: int | None = None) -> list[dict]:
    """R6: extend each finished replica by `ext_ns` via checkpoint
    continuation and merge into md_all.xtc (time-continuous, join frame
    deduplicated by trjcat).

    Mechanism (empirically verified against the 2023.1 build):
    - the extension MDP's nsteps is the TOTAL step count (the tpr the
      -cpi run reads interprets nsteps as the end step);
    - mdrun must get `-cpi <end cpt> -noappend`; it then writes the
      continuation trajectory as <deffnm>.part000N.xtc with a CONTINUING
      time axis (a fresh run would restart at t=0 and trjcat would
      discard it as overlapping);
    - `gmx trjcat -f <current xtc> <new part xtc> -o md_all.xtc` joins
      the two (the shared join frame is dropped).

    Each replica dict gains/updates: ext_round, end_cpt, end_steps,
    merged_xtc. The returned list is the same dicts (mutated)."""
    gmx = env["gmx"]
    cores = max(2, min((n_jobs or n_cores()) // max(1, len(reps)), 32))
    payloads = [(r, str(sysinfo["top"]), ext_ns, gmx, cores) for r in reps]
    return pmap(_extend_job, payloads, n_jobs=max(1, len(reps)))


# --------------------------------------------------------------------------- #
# analysis
# --------------------------------------------------------------------------- #
def _parse_xvg(path: Path) -> tuple[list[float], list[float]]:
    xs, ys = [], []
    if not path.exists():
        return xs, ys
    for line in path.read_text().splitlines():
        if line.startswith("#") or line.startswith("@") or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
            except ValueError:
                pass
    return xs, ys


# --------------------------------------------------------------------------- #
# flexibility diagnostics: per-chain RMSD, DSSP-like SS, interpretation
# --------------------------------------------------------------------------- #
# DSSP-like secondary-structure codes (0 = coil)
SS_COIL, SS_H, SS_G, SS_I, SS_E, SS_B = 0, 1, 2, 3, 4, 5
SS_NAMES = {SS_COIL: "-", SS_H: "H", SS_G: "G", SS_I: "I", SS_E: "E", SS_B: "B"}
# DSSP uses O..HN < 2.3 A with an explicit H; without H we use O..N with a
# relaxed distance and a relaxed angle (measured on the 1HVI crystal: 4.5/100
# gives 69% structured, matching the known ~70% beta-sandwich content).
HBOND_DIST = 4.5   # O(i)..N(j) A
HBOND_ANGLE = 100.0  # deg, C(i)-O(i)..N(j)


def _parse_ndx(path: Path) -> list[tuple[str, int]]:
    """Parse an .ndx file -> [(group_name, n_atoms)]."""
    groups: list[list] = []
    cur = None
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("[") and line.endswith("]"):
                cur = line[1:-1].strip()
                groups.append([cur, 0])
            elif cur is not None and line.strip():
                groups[-1][1] += len(line.split())
    return [(n, c) for n, c in groups]


def _chain_groups(groups: dict) -> list[dict]:
    """Chain groups from a make_ndx listing: splitch names them
    <GROUP>_chainN — protein AND nucleic-acid systems (DNA_chainN for
    DNA-only tprs). Only groups with >= 50 atoms count (splitch can
    emit tiny artifacts)."""
    import re as _re
    return [{"name": n, "atoms": int(c)} for n, c in groups.items()
            if _re.search(r"_chain\d+$", n) and int(c) >= 50]


def build_chain_index(tpr: Path, out_ndx: Path, gmx: str) -> list[dict]:
    """Split the Protein/DNA group into chains (make_ndx splitch - no
    chain-ID dependence, GROMACS drops PDB chain letters in the TPR).
    Returns [{name, atoms}] for each *_chainN group with >= 50 atoms."""
    out_ndx.parent.mkdir(parents=True, exist_ok=True)
    run_cmd([gmx, "make_ndx", "-f", str(tpr), "-o", str(out_ndx)],
            stdin="1\nsplitch 1\nq\n0\n")
    chains = _chain_groups(dict(_parse_ndx(out_ndx)))
    # dedupe: splitch may re-add an identical group for single-chain systems
    seen = set()
    out = []
    for c in chains:
        if c["atoms"] in seen:
            continue
        seen.add(c["atoms"])
        out.append(c)
    return out


def _ss_classify_frame(coords: dict[str, np.ndarray]) -> np.ndarray:
    """DSSP-like secondary-structure assignment for one frame.

    coords: {"N", "CA", "C", "O"} each (n_res, 3); NaN where the atom is
    absent. H-bond = O(i)..N(j), |i-j| >= 2, dist < HBOND_DIST and C-O..N
    angle > HBOND_ANGLE (all residue pairs, both H-bond directions - sheet
    H-bonds in beta-sandwiches span sequence offsets of 30+). Classification
    (DSSP-style): offset 4 -> H, 3 -> G, 5 -> I, other offsets -> E when the
    same offset occurs on >= 2 residue pairs (sheet registration), else B
    (isolated bridge). Returns (n_res,) int codes (SS_*).
    """
    n = len(coords["N"])
    N, C, O = coords["N"], coords["C"], coords["O"]
    ok = np.isfinite(N[:, 0]) & np.isfinite(C[:, 0]) & np.isfinite(O[:, 0])
    dmat = np.linalg.norm(O[:, None, :] - N[None, :, :], axis=2)
    u = C - O
    v = N[None, :, :] - O[:, None, :]
    nu = np.linalg.norm(u, axis=1)
    nv = np.linalg.norm(v, axis=-1)
    with np.errstate(invalid="ignore"):
        cos = (u[:, None, :] * v).sum(-1) / (nu[:, None] * nv + 1e-9)
        ang_ok = np.degrees(np.arccos(np.clip(cos, -1, 1))) > HBOND_ANGLE
    diff = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    B = (dmat < HBOND_DIST) & ang_ok & (diff >= 2) \
        & ok[:, None] & ok[None, :]
    np.fill_diagonal(B, False)
    B = B | B.T  # a H-bond belongs to both residues
    # per-offset pair counts (for E vs B)
    codes = np.zeros(n, dtype=int)
    for i in range(n):
        partners = np.where(B[i])[0]
        if len(partners) == 0:
            continue
        offsets = np.abs(partners - i)
        if 4 in offsets:
            codes[i] = SS_H
        elif 3 in offsets:
            codes[i] = SS_G
        elif 5 in offsets:
            codes[i] = SS_I
        else:
            d = int(offsets.min())
            # E when that offset recurs on other residue pairs (sheet)
            cnt = 0
            for j in range(n - d):
                if B[j, j + d]:
                    cnt += 1
            codes[i] = SS_E if cnt >= 2 else SS_B
    return codes


def _ss_backbone_trajectory(tpr: Path, xtc: Path) -> np.ndarray:
    """Load backbone coords for all frames -> (n_frames, n_res, 4, 3) with
    roles ordered [N, CA, C, O] (NaN where absent). Requires MDAnalysis."""
    import MDAnalysis as MDA
    u = MDA.Universe(str(tpr), str(xtc))
    bb = u.select_atoms("name N or name CA or name C or name O")
    roles = {"N": 0, "CA": 1, "C": 2, "O": 3}
    # R12/R1: per-atom PBC unwrapping (classic min-image accumulation) so
    # per-domain (subset) RMSD is not corrupted by boundary flapping — an
    # atom at the edge of a large protein can cross the box repeatedly
    # while the compound COM (MDAnalysis 2.x `bb.unwrap`, compound-based)
    # never moves > box/2, leaving phantom 50+ A intra-domain distances.
    # Frame 0 is the unwrap reference, so SS classification and the
    # domain reference are identical to the wrapped trajectory.
    box = np.zeros(3)
    prev_raw = None
    unwrapped = None
    # map each selected atom to (residue slot, role). Residue objects are not
    # stable across MDAnalysis calls, so use set intersection instead of
    # object identity.
    bb_set = set(bb.indices)
    res_of_atom: dict[int, int] = {}
    sel_pos = {int(i): k for k, i in enumerate(bb.indices)}
    slot_atoms: list[list[int]] = []  # selection positions per residue slot
    slot = -1
    for r in bb.residues:
        atoms = [a for a in r.atoms if a.index in bb_set and a.name in roles]
        if not atoms:
            continue
        slot += 1
        slot_atoms.append([sel_pos[a.index] for a in atoms])
        for a in atoms:
            res_of_atom[a.index] = slot
    n_res = slot + 1
    flat_idx = np.empty(len(bb), dtype=np.int64)
    for k, atom in enumerate(bb):
        if atom.name not in roles:
            raise ValueError(f"unexpected backbone atom {atom.name}")
        flat_idx[k] = (res_of_atom[atom.index] * 4 + roles[atom.name]) * 3
    flat_idx3 = flat_idx[:, None] + np.arange(3)
    frames = []
    for ts in u.trajectory:
        raw = bb.positions
        if unwrapped is None:
            unwrapped = raw.copy()
            # R12/R1: frame-0 make-whole (chain walk) — per-atom unwrap
            # alone cannot fix an image split that already exists at the
            # reference frame (e.g. a bond crossing the box at t=0 would
            # otherwise read as a ~50 A phantom distance forever)
            box0 = np.asarray(u.dimensions[:3], dtype=float)
            if box0.min() > 1.0 and len(slot_atoms) > 1:
                # stage 1: make each residue internally consistent — a
                # single residue can straddle a box boundary (N on one
                # side, CA/C/O on the other), which poisons any
                # center-based step
                for idxs in slot_atoms:
                    a0 = unwrapped[idxs[0]]
                    for k in idxs[1:]:
                        d = unwrapped[k] - a0
                        unwrapped[k] -= box0 * np.round(d / box0)
                # stage 2: chain walk — pull each residue to its
                # predecessor's image (dimer interfaces included: the
                # walk only needs consecutive slots within box/2)
                prev_center = unwrapped[slot_atoms[0]].mean(axis=0)
                for idxs in slot_atoms[1:]:
                    c = unwrapped[idxs].mean(axis=0)
                    c = c - box0 * np.round((c - prev_center) / box0)
                    shift = c - unwrapped[idxs].mean(axis=0)
                    if np.any(shift != 0):
                        unwrapped[idxs] += shift
                    prev_center = c
        else:
            box = np.asarray(u.dimensions[:3], dtype=float)
            d = raw - prev_raw
            d = d - box * np.round(d / box)
            unwrapped = unwrapped + d
        prev_raw = raw.copy()
        flat = np.full(n_res * 4 * 3, np.nan)
        flat[flat_idx3] = unwrapped
        frames.append(flat.reshape(n_res, 4, 3))
    return np.stack(frames)


def find_ss_domains(codes: np.ndarray, min_res: int = 8) -> list[dict]:
    """R12/R1: structural domains from per-residue DSSP-like codes.

    A domain is a maximal contiguous run of structured residues
    (codes > 0) with length >= min_res. Returns
    [{name, res_start, res_end, n_res}] in sequence order (1-based
    inclusive)."""
    n = len(codes)
    out: list[dict] = []
    i = 0
    while i < n:
        if codes[i] > 0:
            j = i
            while j + 1 < n and codes[j + 1] > 0:
                j += 1
            if j - i + 1 >= min_res:
                out.append({"name": f"dom{len(out) + 1}",
                            "res_start": int(i + 1), "res_end": int(j + 1),
                            "n_res": int(j - i + 1)})
            i = j + 1
        else:
            i += 1
    return out


def _kabsch_transform(mobile: np.ndarray, ref: np.ndarray) -> tuple:
    """Optimal rigid transform fitting `mobile` onto `ref` (both (n, 3),
    same atom order). Row-vector convention: returns (M, t) with
    mobile @ M + t ~= ref. M = U D Vt for the SVD of H = A.T @ B
    (verified against brute-force SO(3) search; a transposed variant
    fails on rank-deficient point clouds)."""
    cm, cr = mobile.mean(axis=0), ref.mean(axis=0)
    A, B = mobile - cm, ref - cr
    H = A.T @ B
    U, _, Vt = np.linalg.svd(H)
    D = np.diag([1.0, 1.0, float(np.sign(np.linalg.det(U @ Vt)))])
    M = U @ D @ Vt
    t = cr - cm @ M
    return M, t


def _kabsch_rmsd(mobile: np.ndarray, ref: np.ndarray) -> float:
    """RMSD (nm) after optimally rotating+translating `mobile` onto
    `ref` (both (n, 3), same atom order)."""
    M, t = _kabsch_transform(mobile, ref)
    fitted = mobile @ M + t
    return float(np.sqrt(np.mean(np.sum((fitted - ref) ** 2, axis=1))))


def domain_vs_rest_rmsd_series(coords: np.ndarray,
                               domains: list[dict]) -> dict[str, list[float]]:
    """R13/R1-v2: per-frame RMSD of each domain measured after fitting
    the REST of the protein (all other CA slots) to its frame-0
    structure.

    The self-fit series (domain_rmsd_series) removes ALL rigid motion,
    so it cannot see a domain swinging like a hinge while the rest of
    the protein stays put. Here the rest is the reference: a domain
    with near-zero self-fit but growing vs-rest RMSD is a rigid-body
    hinge/allosteric motion. NaN where <3 rest or <3 domain CA atoms
    are finite. coords: (F, R, 4, 3), role 1 = CA."""
    out: dict[str, list[float]] = {}
    F, R = coords.shape[:2]
    ca = coords[:, :, 1, :]
    for d in domains:
        s, e = d["res_start"] - 1, d["res_end"]
        dom_ref = ca[0, s:e]
        rest_sel = np.ones(R, dtype=bool)
        rest_sel[s:e] = False
        rest_ref = ca[0, rest_sel]
        ok_rest = np.isfinite(rest_ref[:, 0])
        if int(ok_rest.sum()) < 3:
            continue
        series: list[float] = []
        for f in range(F):
            dom = ca[f, s:e]
            rest_f = ca[f, rest_sel]
            m_rest = np.isfinite(rest_f[:, 0]) & ok_rest
            m_dom = np.isfinite(dom[:, 0])
            if int(m_rest.sum()) < 3 or int(m_dom.sum()) < 3:
                series.append(float("nan"))
                continue
            m_ref = m_dom & np.isfinite(dom_ref[:, 0])
            if int(m_ref.sum()) < 3:
                series.append(float("nan"))
                continue
            Mm, tm = _kabsch_transform(rest_f[m_rest], rest_ref[ok_rest])
            moved = dom[m_ref] @ Mm + tm
            series.append(float(np.sqrt(np.mean(
                np.sum((moved - dom_ref[m_ref]) ** 2, axis=1)))))
        out[d["name"]] = series
    return out


def domain_rmsd_series(coords: np.ndarray,
                       domains: list[dict]) -> dict[str, list[float]]:
    """R12/R1: per-frame domain RMSD (CA atoms) with each domain
    Kabsch-fit to its OWN frame-0 structure.

    Unlike the global `rms` (one fit for the whole system), fitting each
    domain to itself isolates the domain's rigid-body motion — the true
    structural-domain RMSD (DSSP-domain / NMDYN-style idea without the
    extra tool; this GROMACS build has no `doomain`). NaN where <3 CA
    atoms are finite. coords: (F, R, 4, 3), role 1 = CA."""
    out: dict[str, list[float]] = {}
    F = coords.shape[0]
    for d in domains:
        s, e = d["res_start"] - 1, d["res_end"]
        ref = coords[0, s:e, 1, :]
        ok = np.isfinite(ref[:, 0])
        if int(ok.sum()) < 3:
            continue
        series: list[float] = []
        for f in range(F):
            mob = coords[f, s:e, 1, :]
            m = np.isfinite(mob[:, 0]) & ok
            if int(m.sum()) < 3:
                series.append(float("nan"))
            else:
                series.append(_kabsch_rmsd(mob[m], ref[m]))
        out[d["name"]] = series
    return out


def domain_diameters(coords: np.ndarray,
                     domains: list[dict]) -> dict[str, float]:
    """R14: frame-0 CA diameter (nm, max pairwise distance) of each
    domain — the normalization scale for vs-rest RMSD, so a 4 Å shift
    means different things for a 1 nm vs a 2 nm domain."""
    ca = coords[0, :, 1, :]
    out: dict[str, float] = {}
    for d in domains:
        s, e = d["res_start"] - 1, d["res_end"]
        pts = ca[s:e]
        pts = pts[np.isfinite(pts[:, 0])]
        if len(pts) < 2:
            continue
        dd = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
        out[d["name"]] = float(dd.max())
    return out


def analyze_ss(tpr: Path, xtc: Path) -> dict:
    """Per-frame DSSP-like secondary structure for one replica trajectory.
    Returns empty arrays when no protein backbone is present (e.g. a
    DNA-only system) — the N/CA/C/O extraction yields 0 residues."""
    coords = _ss_backbone_trajectory(tpr, xtc)
    F, R = coords.shape[:2]
    if R == 0:
        return {"ss_frac": [], "ss_stable": [],
                "n_frames": int(F), "n_residues": 0}
    codes = np.zeros((F, R), dtype=np.int8)
    for f in range(F):
        codes[f] = _ss_classify_frame({
            "N": coords[f, :, 0, :], "CA": coords[f, :, 1, :],
            "C": coords[f, :, 2, :], "O": coords[f, :, 3, :]})
    structured = (codes > 0).mean(axis=1)
    stable = (codes == codes[0][None, :]).mean(axis=0)
    # R12/R1: structural domains (from frame 0) + per-domain rigid-body
    # RMSD — reuses the already-loaded trajectory (no second MDAnalysis
    # pass)
    domains = find_ss_domains(codes[0])
    if domains:
        _diam = domain_diameters(coords, domains)
        for _d in domains:
            if _d["name"] in _diam:
                _d["diameter_nm"] = round(_diam[_d["name"]], 4)
    domain_rmsd = domain_rmsd_series(coords, domains) if domains else {}
    domain_rmsd_vs_rest = (domain_vs_rest_rmsd_series(coords, domains)
                           if domains else {})
    return {"ss_frac": structured.tolist(),
            "ss_stable": stable.tolist(),
            "n_frames": int(F), "n_residues": int(R),
            "domains": domains,
            "domain_rmsd": domain_rmsd,
            "domain_rmsd_vs_rest": domain_rmsd_vs_rest}


def flexible_regions(rmsf_profile: list[float], *,
                     factor: float = 2.0, min_res: int = 5,
                     floor_nm: float = 0.1) -> list[dict]:
    """R10/G3 (R1 收尾): contiguous flexible segments from the per-residue
    RMSF profile (nm). A residue counts as flexible when
    rmsf >= max(factor*mean, mean+2*std, floor_nm); a region is a maximal
    run of such residues with length >= min_res. Returns
    [{res_start, res_end, n_res, mean_rmsf_nm, max_rmsf_nm}] in residue
    order (1-based inclusive). NaN/None entries are tolerated."""
    import math
    if not rmsf_profile or len(rmsf_profile) < min_res:
        return []
    arr = np.array([0.0 if (v is None or (isinstance(v, float)
                                          and math.isnan(v))) else v
                    for v in rmsf_profile], dtype=float)
    valid = np.isfinite(arr)
    if not valid.any():
        return []
    m = float(arr[valid].mean())
    s = float(arr[valid].std())
    thr = max(factor * m, m + 2.0 * s, floor_nm)
    flex = valid & (arr >= thr)
    regions = []
    i, n = 0, len(arr)
    while i < n:
        if flex[i]:
            j = i
            while j + 1 < n and flex[j + 1]:
                j += 1
            if j - i + 1 >= min_res:
                seg = arr[i:j + 1]
                regions.append({
                    "res_start": int(i + 1),
                    "res_end": int(j + 1),
                    "n_res": int(j - i + 1),
                    "mean_rmsf_nm": round(float(seg.mean()), 4),
                    "max_rmsf_nm": round(float(seg.max()), 4),
                })
            i = j + 1
        else:
            i += 1
    return regions


def interpret_stability(summary: dict) -> list[str]:
    """Rule-based Chinese interpretation of MD stability vs flexibility."""
    notes = []
    f = summary.get("final_rmsd_mean")
    if f is None:
        return notes
    # 1) overall RMSD level
    if f < 0.15:
        notes.append(f"整体稳定: 最终 RMSD {f * 10:.1f} Å (< 1.5 Å), 复合物未发生大尺度变化")
    elif f < 0.30:
        notes.append(f"轻-中度构象漂移: 最终 RMSD {f * 10:.1f} Å (1.5-3.0 Å), 结合 RMSF/聚类判断是柔性还是去稳")
    else:
        notes.append(f"最终 RMSD 较大 ({f * 10:.1f} Å): 需结合下面指标区分柔性域运动/构象系综与局部去折叠")
    # 2) per-chain comparison (domain motion vs unfolding). Use the
    # self-fit (internal) RMSD: it is not inflated by reference-chain
    # flexibility or interface relaxation.
    self_finals = {}
    for k, v in summary.items():
        if k.startswith("rmsd_chain") and k.endswith("_self") \
                and isinstance(v, dict) and v.get("final") is not None:
            self_finals[k.replace("_self", "")] = v["final"]
    if self_finals and f > 0.30:
        # scale: a stable protein's internal (self-fit) RMSD over ~1 ns is
        # typically 0.2-0.6 nm; unfolding pushes it > ~1 nm
        if max(self_finals.values()) < 0.60:
            worst = max(self_finals, key=self_finals.get)
            notes.append(f"各链内部 (自拟合) RMSD 均 < 6 Å "
                         f"(最大 {self_finals[worst] * 10:.1f} Å, {worst}) → "
                         "高整体 RMSD 主要来自链间/结构域相对运动或界面松弛, "
                         "而非链内部展开 — 属柔性/构象重排, 非去折叠")
        else:
            bad = [k for k, v in self_finals.items() if v > 0.90]
            if bad:
                notes.append(f"链 {', '.join(bad)} 内部自拟合 RMSD > 9 Å — "
                             "该链/结构域存在局部去折叠或大尺度内域运动")
    # 3) RMSF level
    rmsf = summary.get("rmsf_profile_mean")
    if rmsf:
        m = float(np.mean(rmsf))
        mx = float(np.max(rmsf))
        if m > 0.25:
            notes.append(f"平均 RMSF 偏高 ({m * 10:.1f} Å) — 体系整体柔性大")
        if mx > 0.6:
            notes.append(f"存在高柔性残基 (最高 RMSF {mx * 10:.1f} Å) — 局部 loop/无序区")
    # 4) conformational heterogeneity (clusters)
    cl = summary.get("clusters") or {}
    if cl:
        pop = {int(k): v for k, v in cl.items()}
        largest = max(pop.values()) if pop else 1.0
        if len(pop) >= 3 and largest < 0.5:
            notes.append(f"无主导构象态 (最大簇 {largest * 100:.0f}% < 50%, {len(pop)} 个簇) — "
                         "构象异质性强; 对柔性靶点属正常, 短轨迹下不能断言已收敛")
    # 5) secondary-structure persistence
    init_ss, fin_ss = summary.get("initial_ss_mean"), summary.get("final_ss_mean")
    if init_ss is not None and fin_ss is not None:
        delta = fin_ss - init_ss
        if delta < -0.10:
            notes.append(f"二级结构占比下降 ({init_ss * 100:.0f}%→{fin_ss * 100:.0f}%) — 提示局部去折叠")
        elif delta > 0.05:
            notes.append(f"二级结构占比上升 ({init_ss * 100:.0f}%→{fin_ss * 100:.0f}%) — 体系在松弛/致密化")
        else:
            notes.append(f"二级结构占比稳定 ({init_ss * 100:.0f}%→{fin_ss * 100:.0f}%)")
    # 6) R10/G3 (R1 收尾): region-level flexibility — cite the most
    # flexible contiguous segments (loops/linkers) by residue range
    regs = summary.get("flexible_regions")
    if regs is not None:
        if regs:
            top = sorted(regs, key=lambda r: -r["mean_rmsf_nm"])[:3]
            desc = "; ".join(
                f"残基 {r['res_start']}-{r['res_end']} ({r['n_res']} 残基, "
                f"平均 RMSF {r['mean_rmsf_nm'] * 10:.1f} Å, "
                f"峰值 {r['max_rmsf_nm'] * 10:.1f} Å)"
                for r in top)
            notes.append(f"识别到柔性区: {desc} — 多为 loop/连接区; "
                         "若位于结合口袋内或其邻近, 建议走柔性靶点工作流 "
                         "(构象选择 + 柔性对接, dock_conformer_set)")
        else:
            # R11: wording must not contradict the single-residue note —
            # "no region" means no run >= min_res residues above the
            # 2x-mean threshold; isolated hot spots may still exist.
            notes.append("无明显连续柔性区 (无 ≥5 连续残基的 RMSF 高值区) — "
                         "构象整体刚性较好; 个别高 RMSF 残基见上条")
    # 7) R11/R5: auto criterion for the flexible-target workflow — an apo
    # system (no bound ligand) whose RMSF is globally elevated means the
    # crystal structure is a poor single representative: recommend a short
    # MD ensemble + consensus docking instead of trusting one pose.
    if summary.get("is_ligand") is False:
        rmsf = summary.get("rmsf_profile_mean")
        if rmsf:
            m_rmsf = float(np.mean(rmsf))
            if m_rmsf > 0.25:
                notes.append(
                    f"自动判据 (R5): 无配体 (apo) 且平均 RMSF 偏高 "
                    f"({m_rmsf * 10:.1f} Å > 2.5 Å) — 晶体结构不是可靠单代表, "
                    "建议短 MD 系综 + 柔性靶点工作流 (make_flex_receptor + "
                    "dock_conformer_set, 用 consensus 分数判命中)")
    # 8) R12/R1: structural-domain RMSD (each domain Kabsch-fit to its own
    # frame 0 -> rigid-body domain motion, independent of the global fit)
    dom_rmsd = summary.get("domain_rmsd") or {}
    if dom_rmsd:
        doms = summary.get("domains") or []
        worst_name = max(dom_rmsd, key=lambda k: dom_rmsd[k]["final"])
        worst = dom_rmsd[worst_name]
        dd = next((d for d in doms if d.get("name") == worst_name), {})
        if worst["final"] > 0.40:
            notes.append(
                f"结构域 {worst_name} (残基 {dd.get('res_start')}-"
                f"{dd.get('res_end')}) 存在明显构象漂移 (末端自拟合 RMSD "
                f"{worst['final'] * 10:.1f} Å, 均值 "
                f"{worst['mean'] * 10:.1f} Å) — 该短结构段内部形变显著 "
                "(环区/不稳定二级结构, 无配体稳定); 若邻近结合口袋, 建议"
                "柔性靶点工作流")
        elif f is not None and f > 0.30 and all(
                v["final"] < 0.30 for v in dom_rmsd.values()):
            notes.append(
                f"各结构域内部均稳定 (最大末端 RMSD {worst['final'] * 10:.1f} "
                "Å < 3.0 Å) — 高整体 RMSD 主要来自域间相对运动而非域内展开")
    # 9) R13/R1-v2: domain vs rest of protein (rest fit to frame 0, domain
    # measured against it -> rigid-body hinge/allosteric motion that the
    # self-fit series removes by construction)
    dom_vs = summary.get("domain_rmsd_vs_rest") or {}
    if dom_vs:
        doms = summary.get("domains") or []
        worst_v = max(dom_vs, key=lambda k: dom_vs[k]["final"])
        dd = next((d for d in doms if d.get("name") == worst_v), {})
        selfit = (summary.get("domain_rmsd") or {}).get(worst_v, {})
        if dom_vs[worst_v]["final"] > 0.40:
            if selfit.get("final", 1.0) < 0.30:
                notes.append(
                    f"结构域 {worst_v} (残基 {dd.get('res_start')}-"
                    f"{dd.get('res_end')}) 相对其余蛋白呈大尺度刚体运动 "
                    f"(末端 {dom_vs[worst_v]['final'] * 10:.1f} Å, 自身"
                    "内部稳定 < 3.0 Å) — 铰链/变构结构域特征; 若其连接"
                    "结合口袋与远端功能位点, 提示变构调控可能, 建议多构象"
                    "柔性靶点工作流")
            else:
                notes.append(
                    f"结构域 {worst_v} (残基 {dd.get('res_start')}-"
                    f"{dd.get('res_end')}) 相对其余蛋白有显著位移 (末端 "
                    f"{dom_vs[worst_v]['final'] * 10:.1f} Å, 自身内部 "
                    f"{selfit.get('final', 0.0) * 10:.1f} Å) — 域运动以"
                    "构象重排为主 (柔性/半刚性域); 若其连接结合口袋与远端"
                    "功能位点, 仍提示变构调控可能")
    return notes


def _parse_group_list(text: str) -> set[str]:
    """Parse `gmx make_ndx` group listing lines
    ("  1 DNA  :  758 atoms") into a set of group names."""
    import re as _re
    return set(_re.findall(r"^\s*\d+\s+(\S+)\s*:\s*\d+\s+atoms",
                           text, _re.M))


def _analysis_group(groups: set[str]) -> str:
    """Best fit/measure group for the core biomolecule: protein
    backbone first, then a DNA/RNA group (DNA-only systems have no
    Backbone/Protein group), then Protein, then System."""
    for g in ("Backbone", "DNA", "RNA", "Protein"):
        if g in groups:
            return g
    return "System"


def _ligand_group(groups: set[str]) -> str | None:
    """R8: the GROMACS index group that holds the small-molecule ligand.

    Group names are derived from moleculetype names: this env's ACPYPE
    mols come out as 'MOL' (or LIG_* with explicit -b names), never
    'Ligand' — a hardcoded prompt dies with "No such group 'Ligand'".
    Preference: explicit LIG_*/Ligand/MOL names, else 'Other' (in a
    protein+ligand system 'Other' is exactly the ligand; DNA systems
    skip ligand RMSD so 'Other' there would be the DNA — callers gate
    on is_ligand). None when there is no plausible candidate."""
    for g in sorted(groups):
        if g.startswith(("LIG_", "Ligand", "MOL")):
            return g
    if "Other" in groups:
        return "Other"
    return None


def _index_groups(tpr: Path, gmx: str, out_ndx: Path) -> set[str]:
    """Default index group names available in a tpr (via gmx make_ndx)."""
    r = run_cmd([gmx, "make_ndx", "-f", str(tpr), "-o", str(out_ndx)],
                stdin="q\n", capture=True)
    out = r.stdout
    if isinstance(out, bytes):
        out = out.decode(errors="replace")
    return _parse_group_list(out or "")


def analyze_replicas(replicas: list[dict], workdir: Path, env: dict,
                     *, is_ligand: bool, save_ps: float = 10.0,
                     burn_in_ps: float = 0.0) -> dict:
    """Per-replica analysis. `burn_in_ps` trims the production burn-in
    (R5): all gmx tools get -b so RMSD/RMSF/Rg/ligand-RMSD/clustering
    start after the equilibration drift settles."""
    gmx = env["gmx"]
    b_args = _b_args(burn_in_ps)
    anadir = workdir / "analysis"
    anadir.mkdir(exist_ok=True)
    per_rep = []
    for r in replicas:
        rdir = Path(r["dir"])
        # R6: prefer the merged (extended) trajectory when present
        merged = rdir / "md_all.xtc"
        xtc = merged if merged.is_file() else rdir / "md.xtc"
        a = {"rep": r["rep"]}
        # system-aware fit/measure group: protein systems use Backbone,
        # DNA-only systems have no Backbone/Protein group (GROMACS builds
        # those from moleculetype names) -> fall back to the DNA/RNA group
        groups = _index_groups(rdir / "md.tpr", gmx,
                               anadir / f"_idx_r{r['rep']}.ndx")
        core = _analysis_group(groups)
        # RMSD (fit on core, measure core). GROMACS 2023 -fit only accepts
        # rot+trans/translation/none, and the tool prompts for index
        # groups, so we answer via stdin.
        run_cmd([gmx, "rms", "-s", str(rdir / "md.tpr"), "-f", str(xtc),
                 *b_args,
                 "-o", str(anadir / f"rmsd_r{r['rep']}.xvg"),
                 "-fit", "rot+trans"],
                log_file=rdir / "rms.log", stdin=f"{core}\n{core}\n")
        # RMSF (per residue, core)
        run_cmd([gmx, "rmsf", "-s", str(rdir / "md.tpr"), "-f", str(xtc),
                 *b_args,
                 "-o", str(anadir / f"rmsf_r{r['rep']}.xvg"),
                 "-res"], log_file=rdir / "rmsf.log", stdin=f"{core}\n")
        # Rg (core)
        run_cmd([gmx, "gyrate", "-s", str(rdir / "md.tpr"), "-f", str(xtc),
                 *b_args,
                 "-o", str(anadir / f"rg_r{r['rep']}.xvg")],
                log_file=rdir / "gyrate.log", stdin=f"{core}\n")
        # ligand stability (fit on core, measure ligand). R8: the group
        # name is system-specific (MOL/LIG_* here, not "Ligand") —
        # resolve it from the tpr and skip with a note when unknown.
        if is_ligand:
            try:
                lg = _ligand_group(_index_groups(
                    rdir / "md.tpr", gmx, rdir / "ligrms.ndx"))
            except Exception:  # noqa: BLE001
                lg = None
            if lg:
                run_cmd([gmx, "rms", "-s", str(rdir / "md.tpr"),
                         "-f", str(xtc),
                         *b_args,
                         "-o", str(anadir / f"lig_rmsd_r{r['rep']}.xvg"),
                         "-fit", "rot+trans"],
                         log_file=rdir / "ligrms.log",
                         stdin=f"{core}\n{lg}\n")
            else:
                logger.warning(
                    "replica %s: no ligand group in tpr; skipping "
                    "ligand RMSD", r["rep"])
        # clustering (GROMACS 2023: -clid for the per-frame cluster id; no
        # -search option -> whole trajectory by default). Two group prompts:
        # fit/RMSD group and distance-calculation group (both the core
        # group — Backbone for protein systems, DNA/RNA for NA systems).
        run_cmd([gmx, "cluster", "-s", str(rdir / "md.tpr"), "-f", str(xtc),
                 *b_args,
                 "-method", "gromos", "-cutoff", "1.5",
                 "-cl", str(anadir / f"clusters_r{r['rep']}"),
                 "-clid", str(anadir / f"cluster_idx_r{r['rep']}.xvg")],
                log_file=rdir / "cluster.log", stdin=f"{core}\n{core}\n")
        # parse
        t, y = _parse_xvg(anadir / f"rmsd_r{r['rep']}.xvg")
        a["rmsd"] = y
        t, y = _parse_xvg(anadir / f"rg_r{r['rep']}.xvg")
        a["rg"] = y
        t, y = _parse_xvg(anadir / f"rmsf_r{r['rep']}.xvg")
        a["rmsf_profile"] = y
        if is_ligand:
            t, y = _parse_xvg(anadir / f"lig_rmsd_r{r['rep']}.xvg")
            a["lig_rmsd"] = y
        # cluster populations
        idx_t, idx_y = _parse_xvg(anadir / f"cluster_idx_r{r['rep']}.xvg")
        if idx_y:
            uniq, counts = np.unique(np.round(idx_y), return_counts=True)
            total = float(counts.sum()) or 1.0
            # cast to native float (numpy.float64 is not msgpack-serializable)
            a["clusters"] = {int(k): round(float(c) / total, 3)
                             for k, c in zip(uniq, counts)}
        # R1: per-chain RMSD (fit on largest chain, measure each other chain)
        chains = build_chain_index(rdir / "md.tpr",
                                   workdir / "chain.ndx", gmx) \
            if not (workdir / "chain.ndx").is_file() else \
            _chain_groups(dict(_parse_ndx(workdir / "chain.ndx")))
        if len(chains) >= 2:
            ref = max(chains, key=lambda c: c["atoms"])
            for c in chains:
                if c["name"] == ref["name"]:
                    continue
                import re as _re
                short = _re.sub(r"(.+)_", "", c["name"])
                o = anadir / f"rmsd_{short}_r{r['rep']}.xvg"
                run_cmd([gmx, "rms", "-s", str(rdir / "md.tpr"),
                         "-f", str(xtc), "-n",
                         str(workdir / "chain.ndx"), "-o", str(o),
                         "-fit", "rot+trans"],
                        log_file=rdir / f"rms_{short}.log",
                        stdin=f"{ref['name']}\n{c['name']}\n")
                t, y = _parse_xvg(o)
                a[f"rmsd_{short}"] = y
                # self-fit RMSD: internal stability of this chain (fit on
                # itself). Relative values above are inflated when the
                # reference chain is flexible or the interface relaxes.
                o2 = anadir / f"rmsd_{short}_self_r{r['rep']}.xvg"
                run_cmd([gmx, "rms", "-s", str(rdir / "md.tpr"),
                         "-f", str(xtc), "-n",
                         str(workdir / "chain.ndx"), "-o", str(o2),
                         "-fit", "rot+trans"],
                        log_file=rdir / f"rms_{short}_self.log",
                        stdin=f"{c['name']}\n{c['name']}\n")
                t, y = _parse_xvg(o2)
                a[f"rmsd_{short}_self"] = y
        # R1: DSSP-like secondary-structure persistence (MDAnalysis; optional)
        try:
            ss = analyze_ss(rdir / "md.tpr", xtc)
            a["ss_frac"] = ss["ss_frac"]
            a["ss_stable"] = ss["ss_stable"]
            a["domains"] = ss.get("domains") or []
            a["domain_rmsd"] = ss.get("domain_rmsd") or {}
            a["domain_rmsd_vs_rest"] = ss.get("domain_rmsd_vs_rest") or {}
        except Exception as e:  # noqa: BLE001 - diagnostics must not kill MD
            logger.warning(f"replica {r['rep']} SS analysis failed: {e}")
        per_rep.append(a)

    # aggregate
    def _avg(key: str) -> dict:
        mats = [np.array(r[key]) for r in per_rep if key in r and r[key]]
        if not mats:
            return {}
        L = min(len(m) for m in mats)
        arr = np.stack([m[:L] for m in mats])
        return {"mean": arr.mean(axis=0).tolist(),
                "std": arr.std(axis=0).tolist()}

    # aggregate cluster populations (mean across replicas)
    cl_all = [{int(k): float(v) for k, v in (r.get("clusters") or {}).items()}
              for r in per_rep if r.get("clusters")]
    cl_keys = {k for c in cl_all for k in c}
    clusters_mean = {int(k): round(float(np.mean([c.get(k, 0.0) for c in cl_all])), 3)
                     for k in sorted(cl_keys)}
    ss_mats = [np.array(r["ss_frac"]) for r in per_rep if r.get("ss_frac")]
    ss_st_mats = [np.array(r["ss_stable"]) for r in per_rep
                  if r.get("ss_stable")]
    summary = {
        "rmsd": _avg("rmsd"),
        "rg": _avg("rg"),
        "lig_rmsd": _avg("lig_rmsd"),
        "clusters": clusters_mean,
        "rmsf_profile_mean": (np.mean([np.array(r["rmsf_profile"])
                                       for r in per_rep if r.get("rmsf_profile")],
                                      axis=0).tolist()
                              if any(r.get("rmsf_profile") for r in per_rep) else []),
        "final_rmsd_mean": float(np.mean([r["rmsd"][-1] for r in per_rep if r.get("rmsd")])),
        "final_rg_mean": float(np.mean([r["rg"][-1] for r in per_rep if r.get("rg")])),
        "ss_frac_mean": (np.mean(ss_mats, axis=0).tolist()
                         if ss_mats else []),
        "ss_stable_mean": (np.mean(ss_st_mats, axis=0).tolist()
                           if ss_st_mats else []),
        "initial_ss_mean": (float(np.mean([m[0] for m in ss_mats]))
                            if ss_mats else None),
        "final_ss_mean": (float(np.mean([m[-1] for m in ss_mats]))
                          if ss_mats else None),
        "replicas": per_rep,
    }
    # R12/R1: structural-domain RMSD (mean across replicas)
    dom_names: list[str] = []
    for r in per_rep:
        for d in r.get("domains") or []:
            if d["name"] not in dom_names:
                dom_names.append(d["name"])
    def _agg(field: str) -> dict:
        out = {}
        for name in dom_names:
            series = [np.array(r[field][name], dtype=float)
                      for r in per_rep if r.get(field)
                      and name in r[field] and r[field][name]]
            if series:
                L = min(len(s) for s in series)
                arr = np.stack([s[:L] for s in series])
                arr = np.nan_to_num(arr, nan=0.0)
                out[name] = {
                    "final": round(float(arr[:, -1].mean()), 4),
                    "mean": round(float(arr.mean()), 4),
                    "series": np.round(arr.mean(axis=0), 4).tolist(),
                }
        return out

    domain_summary = _agg("domain_rmsd")
    vs_rest_summary = _agg("domain_rmsd_vs_rest")
    if domain_summary:
        summary["domains"] = next(
            (r["domains"] for r in per_rep if r.get("domains")), [])
        summary["domain_rmsd"] = domain_summary
        if vs_rest_summary:
            # R14: size-normalized vs-rest (diameter from frame 0) so the
            # hinge threshold is comparable across domain sizes
            for _name, _st in vs_rest_summary.items():
                _diam = next((d.get("diameter_nm") for d in summary["domains"]
                              if d.get("name") == _name), None)
                if _diam:
                    _st["final_norm"] = round(_st["final"] / _diam, 3)
                    _st["mean_norm"] = round(_st["mean"] / _diam, 3)
            summary["domain_rmsd_vs_rest"] = vs_rest_summary
    # per-chain RMSD (mean across replicas + mean final value)
    for key in sorted({k for r in per_rep for k in r
                       if k.startswith("rmsd_chain") and k != "rmsd"}):
        vals = [np.array(r[key]) for r in per_rep if r.get(key)]
        if vals:
            L = min(len(v) for v in vals)
            arr = np.stack([v[:L] for v in vals])
            summary[key] = {"mean": arr.mean(axis=0).tolist(),
                            "final": round(float(np.mean([v[-1] for v in vals])), 4)}
    # R10/G3 (R1 收尾): region-level flexibility from the mean RMSF profile
    summary["flexible_regions"] = flexible_regions(
        summary["rmsf_profile_mean"])
    summary["interpretation"] = interpret_stability(summary)
    return summary


# --------------------------------------------------------------------------- #
# graph node
# --------------------------------------------------------------------------- #
def run_md(state: dict) -> dict:
    workdir = Path(state["project_dir"]) / "05_md"
    workdir.mkdir(parents=True, exist_ok=True)
    opts = state.get("options", {})
    d = resolve_defaults(opts)
    brain = AgentBrain(project_dir=Path(state["project_dir"])) \
        if not opts.get("no_llm") else None

    env = gromacs()
    logger.info(f"using GROMACS at {env['gmx']} (ff={env['ff']})")

    choice = select_system(state, brain)
    complex_pdb = build_complex_pdb(state, choice, workdir)
    sysinfo = build_system(complex_pdb, workdir, env,
                           is_ligand=choice["ligand"],
                           salt=d.md_salt_m, box_margin=d.md_box_margin_nm,
                           divalent=d.md_divalent,
                           divalent_m=d.md_divalent_m)
    ns = float(opts.get("md_ns", d.md_ns))
    reps = int(opts.get("md_reps", d.md_reps))
    # R5: NVT->NPT equilibration once per build, replicas fork from it
    eq = run_equilibration(
        sysinfo, workdir, env,
        nvt_ps=float(opts.get("md_eq_nvt_ps", d.md_eq_nvt_ps)),
        npt_ps=float(opts.get("md_eq_npt_ps", d.md_eq_npt_ps)))
    replicas = run_replicas(sysinfo, workdir, env, ns=ns, reps=reps, eq=eq)
    burn_in_ps = float(opts.get("md_burn_in_ps", d.md_burn_in_ps))
    summary = analyze_replicas(replicas, workdir, env,
                               is_ligand=choice["ligand"], save_ps=d.md_save_ps,
                               burn_in_ps=burn_in_ps)
    # R6: convergence-gated auto-extension — re-analyze after each
    # extension round until the RMSD has plateaued with a dominant
    # cluster (or the extension budget runs out)
    ext_max = int(opts.get("md_max_extensions", d.md_max_extensions))
    ext_ns = float(opts.get("md_extend_ns", d.md_extend_ns)) or ns
    total_ns = ns
    ext_round = 0
    converged, reason = _md_converged(
        summary, total_ns * 1000.0,
        min_len_ps=float(opts.get("md_converge_min_len_ps",
                                  d.md_converge_min_len_ps)),
        drift_nm=float(opts.get("md_converge_rmsd_drift_nm",
                                d.md_converge_rmsd_drift_nm)),
        min_cluster=float(opts.get("md_converge_min_cluster",
                                   d.md_converge_min_cluster)))
    while not converged and ext_round < ext_max:
        logger.info(f"MD not converged ({reason}) — extending all replicas "
                    f"by {ext_ns} ns (round {ext_round + 1}/{ext_max})")
        replicas = extend_replicas(replicas, sysinfo, workdir, env,
                                   ext_ns=ext_ns)
        ext_round += 1
        total_ns += ext_ns
        summary = analyze_replicas(replicas, workdir, env,
                                   is_ligand=choice["ligand"],
                                   save_ps=d.md_save_ps,
                                   burn_in_ps=burn_in_ps)
        converged, reason = _md_converged(
            summary, total_ns * 1000.0,
            min_len_ps=float(opts.get("md_converge_min_len_ps",
                                      d.md_converge_min_len_ps)),
            drift_nm=float(opts.get("md_converge_rmsd_drift_nm",
                                    d.md_converge_rmsd_drift_nm)),
            min_cluster=float(opts.get("md_converge_min_cluster",
                                       d.md_converge_min_cluster)))
    logger.info(f"MD convergence after {ext_round} extension(s): "
                f"converged={converged} ({reason}); total {total_ns} ns")

    out = {
        "system": choice,
        "complex_pdb": str(complex_pdb),
        "gromacs": {"binary": env["gmx"], "ff": env["ff"], "version": env["ver"]},
        "ns": ns, "reps": reps, "total_ns": total_ns,
        "equilibration": {"nvt_ps": eq.get("nvt_ps"),
                          "npt_ps": eq.get("npt_ps"),
                          "eq_gro": eq["eq_gro"], "eq_cpt": eq["eq_cpt"]},
        "extensions": {"rounds": ext_round, "converged": converged,
                       "reason": reason, "total_ns": total_ns},
        "replicas": replicas,
        "summary": summary,
        "build_dir": sysinfo["build_dir"],
    }
    jsave(workdir / "md.json", out)
    state_out = dict(state)
    state_out["md"] = out
    return state_out
