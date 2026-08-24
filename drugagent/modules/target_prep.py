"""Module A: target preparation.

Resolve PDB file / PDB ID / sequence -> completeness analysis -> agent
judgment -> cleaning -> pocket detection -> docking-ready PDBQT.
"""
from __future__ import annotations

import re
import shutil
from pathlib import Path

import numpy as np
from loguru import logger

from ..config import ENV_DIR, WEIGHTS
from ..llm import AgentBrain
from ..utils import (centroid_from_pdb, download_file, ensure_parent, fetch_pdb,
                     is_pdb_id,
                     jsave, pdb_extent, pdb_ligands, pdb_chains, run_cmd)

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "H": 1.20, "P": 1.80,
       "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}
SKIP_RES = {"HOH", "WAT", "H2O"}

METAL_RES = {"ZN", "MG", "FE", "NA", "K", "CA", "MN", "CO", "CU", "NI", "MO",
             "ZR", "CD", "LI", "CR", "SE", "RB"}
# standard 2-letter nucleotide codes (unambiguous); single-letter A/C/G/T/U
# are ambiguous with amino acids and only count when the residue carries P
NUC2_RES = {"DA", "DT", "DC", "DG", "RA", "RT", "RC", "RG"}


# --------------------------------------------------------------------------- #
# 1. resolve input
# --------------------------------------------------------------------------- #
def resolve_target(target: dict, workdir: Path, *, modeler: str = "esmfold") -> dict:
    """target: {"kind": pdb_file|pdb_id|fasta|sequence, "value": str}
    Returns {"pdb_path": Path, "kind": str, "note": str}"""
    kind = target["kind"]
    value = target["value"]
    workdir.mkdir(parents=True, exist_ok=True)
    if kind == "pdb_file":
        p = Path(value)
        if not p.exists():
            raise FileNotFoundError(f"PDB file not found: {p}")
        return {"pdb_path": p, "kind": kind, "note": f"using user file {p.name}"}
    if kind == "pdb_id":
        p = fetch_pdb(value.upper(), workdir / "raw")
        return {"pdb_path": p, "kind": kind, "note": f"fetched {value.upper()} from RCSB"}
    if kind == "fasta":
        p = Path(value)
        seq = _fasta_seq(p)
        pdb_path = _model_from_sequence(seq, workdir, modeler)
        return {"pdb_path": pdb_path, "kind": kind,
                "note": f"modeled from FASTA {p.name} via {modeler}"}
    if kind == "sequence":
        pdb_path = _model_from_sequence(value, workdir, modeler)
        return {"pdb_path": pdb_path, "kind": kind,
                "note": f"modeled from {len(value)}-aa sequence via {modeler}"}
    raise ValueError(f"unknown target kind: {kind}")


def _fasta_seq(path: Path) -> str:
    seqs = []
    cur = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line.strip().upper())
    if cur:
        seqs.append("".join(cur))
    if not seqs:
        raise ValueError(f"no sequence found in {path}")
    return max(seqs, key=len)


def _model_from_sequence(seq: str, workdir: Path, modeler: str) -> Path:
    from .esmfold_run import predict, write_pdb

    out = predict(seq, num_recycles=3)
    pdb_path = write_pdb(out, workdir / "modeled" / "model.pdb")
    jsave(workdir / "modeled" / "model_info.json",
          {"mean_plddt": out["mean_plddt"], "min_plddt": out["min_plddt"],
           "length": len(seq)})
    return pdb_path


# --------------------------------------------------------------------------- #
# 2. completeness analysis (rule based)
# --------------------------------------------------------------------------- #
def analyze_completeness(pdb_path: Path) -> dict:
    """Rule-based structural completeness report."""
    chains = pdb_chains(pdb_path)
    ligands = pdb_ligands(pdb_path)
    n_atoms = 0
    n_het = 0
    resnames: set[str] = set()
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                if line.startswith("HETATM"):
                    n_het += 1
                else:
                    n_atoms += 1
                    resnames.add(line[17:20].strip())
    n_res = len(resnames)
    report = {
        "file": str(pdb_path),
        "n_chains": len(chains),
        "chains": chains,
        "n_atoms": n_atoms,
        "n_hetatoms": n_het,
        "n_residue_types": n_res,
        "ligands": ligands,
        "has_ligand": len(ligands) > 0,
        "n_protein_residues_approx": n_atoms // 8,
    }
    # multimer check: identical-ish chains (crystallographic or biological)
    report["multimer"] = len(chains) > 1
    # structural pitfalls (multi-model, altloc, gaps, metals, nucleic acids...)
    try:
        report["issues"] = detect_structure_issues(pdb_path)
    except Exception as e:  # never let pre-check break the pipeline
        report["issues"] = [{"type": "precheck_error", "severity": "warn",
                             "detail": str(e), "suggestion": "manual review",
                             "fix": None}]
    # obvious problems
    problems = []
    if not chains:
        problems.append("no ATOM records")
    if n_res < 3:
        problems.append("very few residue types - possibly incomplete")
    return report


def agent_judge_target(report: dict, brain: AgentBrain | None,
                       context: str = "") -> dict:
    """LLM judgment on how to proceed; falls back to rules if no LLM."""
    ctx = (f"chains={report['chains']} ligands={report['ligands']} "
           f"n_res={report['n_residue_types']} multimer={report['multimer']}")
    issues = report.get("issues") or []
    if issues:
        ctx += "\n结构预检发现的问题:\n" + "\n".join(
            f"- [{i.get('severity')}] {i.get('type')}: {i.get('detail')}"
            + (f" (可自动修复: {i.get('fix')})" if i.get("fix") else "")
            for i in issues)
    if brain is None:
        decision = {
            "action": "clean_and_use",
            "rationale": "LLM unavailable; rule-based default: clean and use structure",
        }
        return decision
    dec = brain.decide(
        "target_prep",
        "该结构能否直接使用？选择 action: clean_and_use(清洗后直接用) / "
        "keep_ligand(清洗并保留配体作对接位点) / remodel(结构不完整需重新建模)。"
        "若结构含配体且残基完整，通常选 keep_ligand。",
        context=ctx + ("\n" + context if context else ""),
        choices=["clean_and_use", "keep_ligand", "remodel"],
        expect="choice",
    )
    return {"action": dec.answer, "rationale": dec.rationale}


# --------------------------------------------------------------------------- #
# 2b. structure issues: detect + repair
# --------------------------------------------------------------------------- #
def detect_structure_issues(pdb_path: Path) -> list[dict]:
    """Scan a PDB for structural pitfalls that break docking/MD.

    Returns a list of issues: {type, severity, detail, suggestion, fix}.
    severity: error (likely corrupts results) / warn / info.
    """
    issues: list[dict] = []
    n_models = 0
    altlocs: dict[str, int] = {}
    resseq: dict[str, dict[int, str]] = {}   # chain -> {resseq: resname}
    seqres: dict[str, list[tuple[int, str]]] = {}
    metals: dict[str, int] = {}
    nuc = 0
    p_residues: set[tuple[str, int]] = set()
    heavy: dict[tuple[str, int], int] = {}
    ssbond = 0
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("MODEL "):
                n_models += 1
            elif line.startswith("SSBOND"):
                ssbond += 1
            elif line.startswith("SEQRES"):
                ch = line[7]
                try:
                    first = int(line[8:13])
                except ValueError:
                    first = 1
                i = 14
                n = 0
                while i + 3 <= len(line.rstrip()):
                    rn = line[i:i + 3].strip()
                    if not rn:
                        break
                    n += 1
                    seqres.setdefault(ch, []).append((first + n - 1, rn))
                    i += 3
            elif line.startswith(("ATOM", "HETATM")):
                alt = line[16]
                if alt not in (" ", "A"):
                    altlocs[alt] = altlocs.get(alt, 0) + 1
                resname = line[17:20].strip()
                ch = line[21]
                try:
                    n = int(line[22:26])
                except ValueError:
                    continue
                # nucleotides may be ATOM (standard) or HETATM (some files)
                if resname in NUC2_RES or (
                        resname in {"A", "C", "G", "T", "U"}
                        and line[76:78].strip() == "P"):
                    nuc += 1
                if line.startswith("HETATM"):
                    if resname in METAL_RES:
                        metals[resname] = metals.get(resname, 0) + 1
                    continue
                if "P" in line[76:78]:
                    p_residues.add((ch, n))
                resseq.setdefault(ch, {})[n] = resname
                heavy[(ch, n)] = heavy.get((ch, n), 0) + 1

    if n_models > 1:
        issues.append({
            "type": "multiple_models", "severity": "error",
            "detail": f"{n_models} MODEL records (NMR ensemble / alternate conformations); "
                      "naive cleaning merges all models into one corrupt structure",
            "suggestion": "keep the best conformation (model 1 or ligand-bound)",
            "fix": "dedupe_models"})
    if altlocs:
        issues.append({
            "type": "alternate_locations", "severity": "warn",
            "detail": f"alternate locations present: "
                      + ", ".join(f"{k}:{v}" for k, v in sorted(altlocs.items())),
            "suggestion": "keep altloc A (drop B/C)", "fix": "keep_altloc_a"})

    # missing residues: gaps in ATOM numbering, and SEQRES minus ATOM
    for ch, resmap in resseq.items():
        if not resmap:
            continue
        nums = sorted(resmap)
        gaps = []
        for a, b in zip(nums, nums[1:]):
            if b - a > 1:
                gaps.append(f"{a}..{b - 1}")
        missing_from_seqres = []
        if ch in seqres:
            have = set(nums)
            missing_from_seqres = [f"{n}({r})" for n, r in seqres[ch]
                                   if n not in have][:10]
        if gaps or missing_from_seqres:
            detail = []
            if gaps:
                detail.append(f"numbering gaps: {', '.join(gaps[:8])}"
                              + (" ..." if len(gaps) > 8 else ""))
            if missing_from_seqres:
                detail.append(f"in SEQRES but not modeled: {', '.join(missing_from_seqres)}")
            issues.append({
                "type": "missing_residues", "severity": "warn",
                "detail": f"chain {ch}: " + "; ".join(detail),
                "suggestion": "unresolved loops may sit in the binding site; "
                              "flag or model before docking",
                "fix": None})

    if metals:
        issues.append({
            "type": "metals", "severity": "warn",
            "detail": "metal ions: "
                      + ", ".join(f"{k}x{v}" for k, v in sorted(metals.items())),
            "suggestion": "keep for MD if catalytically required (e.g. Zn "
                          "metalloprotease); drop for pure docking",
            "fix": "keep_metals / drop_metals"})
    if nuc > 0:
        issues.append({
            "type": "nucleic_acids", "severity": "warn",
            "detail": f"{nuc} nucleotide residues — MD handles NA chains "
                      "natively (amber99sb-ildn carries b-DNA/RNA parameters; "
                      "R4 keeps NA polymer chains on the protein side of the "
                      "split); short NA ligands go through ACPYPE/GAFF2",
            "suggestion": "keep NA chains for MD; for docking, the grid "
                          "box should cover the binding interface on the NA",
            "fix": None})
    if ssbond:
        issues.append({
            "type": "disulfide_bonds", "severity": "info",
            "detail": f"{ssbond} SSBOND records (pdb2gmx will apply them)",
            "suggestion": "none needed", "fix": None})

    # disordered termini: terminal runs of backbone-only residues (<4 heavy atoms)
    for ch, resmap in resseq.items():
        nums = sorted(resmap)
        if len(nums) < 6:
            continue
        def _bad(n: int) -> bool:
            return heavy.get((ch, n), 0) < 4
        n_start = 0
        for n in nums:
            if _bad(n):
                n_start += 1
            else:
                break
        n_end = 0
        for n in reversed(nums):
            if _bad(n):
                n_end += 1
            else:
                break
        if n_start >= 3 or n_end >= 3:
            issues.append({
                "type": "disordered_termini", "severity": "info",
                "detail": f"chain {ch}: ~{n_start} N-term / ~{n_end} C-term "
                          "terminal residues are backbone-only (disordered)",
                "suggestion": "trim for docking/MD if far from the pocket",
                "fix": "trim_disordered"})
    return issues


def repair_structure(pdb_path: Path, out_path: Path, *,
                     actions: list[str] | None = None) -> dict:
    """Apply structure repairs, write a fixed PDB, return what changed.

    actions: any of dedupe_models (keep model 1), keep_altloc_a,
    trim_disordered (drop backbone-only terminal runs), drop_metals,
    drop_hetatm (drop all HETATM).
    """
    actions = set(actions or [])
    in_model = 0
    removed: dict[str, int] = {}

    def _cut(kind: str) -> None:
        removed[kind] = removed.get(kind, 0) + 1

    # pass 1: find disordered terminal runs if trimming
    trim: set[tuple[str, int]] = set()
    if "trim_disordered" in actions:
        heavy: dict[tuple[str, int], int] = {}
        resmap: dict[str, list[int]] = {}
        with open(pdb_path) as fh:
            for line in fh:
                if line.startswith("ATOM"):
                    try:
                        key = (line[21], int(line[22:26]))
                    except ValueError:
                        continue
                    heavy[key] = heavy.get(key, 0) + 1
                    resmap.setdefault(key[0], []).append(key[1])
        for ch, nums in resmap.items():
            nums = sorted(set(nums))
            if len(nums) < 6:
                continue
            for n in nums:
                if heavy.get((ch, n), 0) < 4:
                    trim.add((ch, n))
                else:
                    break
            for n in reversed(nums):
                if heavy.get((ch, n), 0) < 4:
                    trim.add((ch, n))
                else:
                    break

    dedupe = "dedupe_models" in actions
    out_lines = []
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("MODEL ") or line.startswith("ENDMDL"):
                if line.startswith("MODEL "):
                    in_model += 1
                if dedupe:
                    continue  # output is a plain single-conformation PDB
            elif dedupe and in_model > 1 and not line.startswith(("ATOM", "HETATM")):
                continue  # skip non-atom lines of later models too
            elif line.startswith(("ATOM", "HETATM")):
                if "dedupe_models" in actions and in_model > 1:
                    _cut("other_models")
                    continue
                if "keep_altloc_a" in actions and line[16] not in (" ", "A"):
                    _cut("altloc")
                    continue
                try:
                    ch, n = line[21], int(line[22:26])
                except ValueError:
                    ch, n = line[21], None
                if n is not None and (ch, n) in trim:
                    _cut("disordered")
                    continue
                if line.startswith("HETATM"):
                    rn = line[17:20].strip()
                    if "drop_hetatm" in actions:
                        _cut("hetatm")
                        continue
                    if "drop_metals" in actions and rn in METAL_RES:
                        _cut("metals")
                        continue
            out_lines.append(line)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(out_lines))
    return {"out": str(out_path), "applied": sorted(actions),
            "removed_atoms": removed, "n_models_kept":
            1 if "dedupe_models" in actions else None}


# --------------------------------------------------------------------------- #
# 3. cleaning
# --------------------------------------------------------------------------- #
def clean_pdb(pdb_path: Path, out_path: Path, *,
              keep_resnames: list[str] | None = None,
              keep_waters: bool = False, water_keep_dist: float = 5.0,
              keep_chain: str | None = None) -> dict:
    """Write a cleaned PDB: ATOM (all chains or one) + selected HETATM + waters."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    keep = set(keep_resnames or [])
    keep_w = set(w.upper() for w in (keep_resnames or []) if w.upper() in SKIP_RES)
    center = None
    if keep_waters:
        try:
            center = centroid_from_pdb(pdb_path, resnames=list(keep) or None)
        except ValueError:
            center = None
    n_keep = n_drop = 0
    with open(pdb_path) as fh:
        lines_out = []
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                chain = line[21]
                if keep_chain and chain != keep_chain:
                    n_drop += 1
                    continue
                if line.startswith("ATOM"):
                    lines_out.append(line)
                    n_keep += 1
                elif resname in keep:
                    lines_out.append(line)
                    n_keep += 1
                elif resname.upper() in SKIP_RES and keep_waters and center is not None:
                    x, y, z = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    d = ((x - center[0]) ** 2 + (y - center[1]) ** 2
                         + (z - center[2]) ** 2) ** 0.5
                    if d <= water_keep_dist:
                        lines_out.append(line)
                        n_keep += 1
                    else:
                        n_drop += 1
                else:
                    n_drop += 1
            elif line.startswith(("MASTER", "END", "CONECT")):
                continue
            else:
                lines_out.append(line)
        lines_out.append("END\n")
    out_path.write_text("".join(lines_out))
    return {"kept_atoms": n_keep, "dropped_atoms": n_drop,
            "kept_resnames": sorted(keep), "out": str(out_path)}


# --------------------------------------------------------------------------- #
# 4. pocket detection
# --------------------------------------------------------------------------- #
def pocket_from_ligand(pdb_path: Path, lig_resname: str,
                       pad: float = 2.5, max_size: float = 40.0) -> dict:
    center = centroid_from_pdb(pdb_path, resnames=[lig_resname])
    extent = pdb_extent(pdb_path, center, resnames=[lig_resname])
    r = min(max(extent + pad, 8.0), max_size / 2)
    return {
        "center": list(center),
        "xsize": round(2 * r, 2), "ysize": round(2 * r, 2), "zsize": round(2 * r, 2),
        "method": f"ligand_centroid({lig_resname})",
        "ligand": lig_resname,
        "site_id": "S1",
    }


def grid_pockets(pdb_path: Path, grid: float = 0.5, pad: float = 5.0,
                 min_vol: float = 60.0, max_vol: float = 4000.0) -> list[dict]:
    """Simple grid-based cavity detection. Returns ranked pocket list."""
    xs, ys, zs = [], [], []
    atoms: list[tuple[float, float, float, float]] = []  # x,y,z,radius
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM"):
                el = line[12:16].strip()
                el1 = re.sub(r"[\d\-+]", "", el)
                r = VDW.get(el1[:2].upper() if el1[:2].upper() in VDW else el1[0].upper(), 1.8)
                atoms.append((float(line[30:38]), float(line[38:46]),
                              float(line[46:54]), r))
                xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
    if not atoms:
        return []
    xmin, ymin, zmin = min(xs) - pad, min(ys) - pad, min(zs) - pad
    xmax, ymax, zmax = max(xs) + pad, max(ys) + pad, max(zs) + pad
    nx, ny, nz = int((xmax - xmin) / grid) + 1, int((ymax - ymin) / grid) + 1, int((zmax - zmin) / grid) + 1
    # coarse occupancy: mark cells within vdW radius of any atom
    occ = np.zeros((nx, ny, nz), dtype=np.uint8)
    for ax, ay, az, r in atoms:
        rad_cells = int(r / grid) + 1
        ix0, ix1 = int((ax - xmin) / grid) - rad_cells, int((ax - xmin) / grid) + rad_cells + 1
        iy0, iy1 = int((ay - ymin) / grid) - rad_cells, int((ay - ymin) / grid) + rad_cells + 1
        iz0, iz1 = int((az - zmin) / grid) - rad_cells, int((az - zmin) / grid) + rad_cells + 1
        ix0, iy0, iz0 = max(0, ix0), max(0, iy0), max(0, iz0)
        ix1, iy1, iz1 = min(nx, ix1), min(ny, iy1), min(nz, iz1)
        if ix1 <= ix0 or iy1 <= iy0 or iz1 <= iz0:
            continue
        occ[ix0:ix1, iy0:iy1, iz0:iz1] = 1
    # flood fill from border through empty cells
    from collections import deque
    empty = occ == 0
    visited = np.zeros_like(occ, dtype=np.uint8)
    q = deque()
    for i in range(nx):
        for j in range(ny):
            for k in (0, nz - 1):
                if empty[i, j, k] and not visited[i, j, k]:
                    visited[i, j, k] = 1
                    q.append((i, j, k))
    for i in range(nx):
        for k in range(nz):
            for j in (0, ny - 1):
                if empty[i, j, k] and not visited[i, j, k]:
                    visited[i, j, k] = 1
                    q.append((i, j, k))
    for j in range(ny):
        for k in range(nz):
            for i in (0, nx - 1):
                if empty[i, j, k] and not visited[i, j, k]:
                    visited[i, j, k] = 1
                    q.append((i, j, k))
    while q:
        i, j, k = q.popleft()
        for di, dj, dk in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            ni, nj, nk = i + di, j + dj, k + dk
            if 0 <= ni < nx and 0 <= nj < ny and 0 <= nk < nz \
                    and empty[ni, nj, nk] and not visited[ni, nj, nk]:
                visited[ni, nj, nk] = 1
                q.append((ni, nj, nk))
    cav = empty & (~visited.astype(bool)) & ~occ.astype(bool)
    # label connected components (6-conn)
    labels = np.zeros_like(occ, dtype=np.int32)
    pockets = []
    lid = 0
    cidx = np.argwhere(cav)
    lab_map = {tuple(c): -1 for c in cidx}
    for c in cidx:
        i, j, k = int(c[0]), int(c[1]), int(c[2])
        if lab_map[(i, j, k)] != -1:
            continue
        lid += 1
        q = deque([(i, j, k)])
        lab_map[(i, j, k)] = lid
        cells = []
        while q:
            ci, cj, ck = q.popleft()
            cells.append((ci, cj, ck))
            for di, dj, dk in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                ni, nj, nk = ci + di, cj + dj, ck + dk
                key = (ni, nj, nk)
                if 0 <= ni < nx and 0 <= nj < ny and 0 <= nk < nz \
                        and lab_map.get(key, -1) == -1 and cav[ni, nj, nk]:
                    lab_map[key] = lid
                    q.append((ni, nj, nk))
        if len(cells) < 8:
            continue
        cells = np.array(cells, dtype=float)
        # center in Angstrom
        center = np.array([
            xmin + (cells[:, 0] + 0.5) * grid,
            ymin + (cells[:, 1] + 0.5) * grid,
            zmin + (cells[:, 2] + 0.5) * grid,
        ]).mean(axis=0)
        vol = len(cells) * grid ** 3
        # extent
        ext = (cells.max(axis=0) - cells.min(axis=0)).max() * grid
        if not (min_vol <= vol <= max_vol):
            continue
        r = min(max(ext / 2 + 2.0, 8.0), 20.0)
        pockets.append({
            "center": [round(float(v), 2) for v in center],
            "xsize": round(2 * r, 2), "ysize": round(2 * r, 2), "zsize": round(2 * r, 2),
            "method": "grid_cavity",
            "volume_A3": round(vol, 1),
            "n_cells": len(cells),
        })
    pockets.sort(key=lambda p: p["volume_A3"], reverse=True)
    return pockets


def find_pocket(target_pdb: Path, lig_resname: str | None,
                brain: AgentBrain | None = None) -> dict:
    """Choose docking site. Known ligand wins; otherwise top grid cavity."""
    if lig_resname:
        p = pocket_from_ligand(target_pdb, lig_resname)
        p["rationale"] = f"已知配体 {lig_resname} 所在位点，最可靠的对接位点"
        return p
    cands = grid_pockets(target_pdb)
    if not cands:
        # fallback: protein geometric center
        c = centroid_from_pdb(target_pdb)
        return {"center": list(c), "xsize": 24, "ysize": 24, "zsize": 24,
                "method": "protein_centroid_fallback", "rationale": "未检测到空腔，用蛋白中心兜底"}
    if brain is not None and len(cands) > 1:
        ctx = "; ".join(
            f"cavity{i+1}: center={p['center']} vol={p['volume_A3']}A3"
            for i, p in enumerate(cands[:5]))
        dec = brain.decide("target_prep", "选择最佳对接口袋（编号1-5）",
                           context=ctx, choices=[str(i + 1) for i in range(min(5, len(cands)))],
                           expect="choice")
        idx = int(dec.answer) - 1
        chosen = cands[idx]
        chosen["rationale"] = f"LLM 选择了空腔{idx+1}: {dec.rationale}"
        return chosen
    cands[0]["rationale"] = "按体积最大的网格空腔选择"
    return cands[0]


# --------------------------------------------------------------------------- #
# 5. pdbqt conversion
# --------------------------------------------------------------------------- #
def obabel_bin() -> str:
    p = ENV_DIR / "bin" / "obabel"
    if p.is_file():
        return str(p)
    import shutil as _s
    found = _s.which("obabel")
    if found:
        return found
    raise FileNotFoundError("obabel not found (run `drugagent setup`)")


def strip_root_sections(pdbqt: Path) -> Path:
    """Remove ROOT/ENDROOT/BRANCH/ENDBRAN lines (some vina builds reject them;
    a fully rigid molecule needs no explicit rigid-body tags)."""
    txt = pdbqt.read_text()
    if "ROOT" not in txt:
        return pdbqt
    lines = [l for l in txt.splitlines()
             if l.split()[:1] not in (["ROOT"], ["ENDROOT"], ["BRANCH"],
                                      ["ENDBRANCH"], ["TORSDOF"])]
    pdbqt.write_text("\n".join(lines) + "\n")
    return pdbqt


def to_pdbqt(pdb: Path, out: Path, *, keep_resnames: list[str] | None = None,
             partial_charges: str = "gasteiger", flex: bool = True) -> Path:
    """PDB -> PDBQT (obabel). `flex=False` for rigid receptors: this vina
    build rejects molecule-graph keywords (ROOT/TORSDOF) in rigid receptors
    but requires them in flex ligands."""
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.is_file() and out.stat().st_size > 0:
        # cache: deterministic conversion — but re-convert if the source PDB
        # is newer (stale PDBQTs from an earlier code version bite otherwise)
        if out.stat().st_mtime >= pdb.stat().st_mtime:
            return out
    cmd = [obabel_bin(), str(pdb), "-O", str(out), "--partial_charges", partial_charges]
    run_cmd(cmd, log_file=out.with_suffix(".log"))
    strip_root_sections(out)
    _sanitize_pdbqt_elements(out, add_graph=flex)
    return out


_KNOWN_ELEMENTS = {"C", "N", "O", "S", "P", "F", "CL", "BR", "I", "ZN", "FE",
                   "MG", "NA", "K", "CA", "MN", "MO", "NI", "CO", "CU", "LI",
                   "SE", "H"}


def _sanitize_pdbqt_elements(pdbqt: Path, *, add_graph: bool = True) -> None:
    """obabel sometimes drops the PDB element column and derives bogus
    element tags (e.g. 'A' for ALA CA atoms), which Vina rejects with
    'Unknown or inappropriate tag'. Re-derive from atom names when needed;
    add the AD4 molecule-graph keywords when this vina build needs them."""
    lines = pdbqt.read_text().splitlines()
    fixed = 0
    for i, l in enumerate(lines):
        if not l.startswith("ATOM") or len(l) < 78:
            continue
        el = l[76:78].strip()
        if el in _KNOWN_ELEMENTS:
            continue
        name = l[12:16].strip()
        if name.startswith("N"):
            el = "N"
        elif name.startswith("O"):
            el = "O"
        elif name.startswith("S"):
            el = "S"
        else:
            el = "C"
        lines[i] = l[:76] + f" {el} "
        fixed += 1
    if fixed:
        lines = [l for l in lines if l.strip()]
        pdbqt.write_text("\n".join(lines) + "\n")
    if add_graph:
        _add_molecule_graph(pdbqt)


def _add_molecule_graph(pdbqt: Path) -> None:
    """This vina build requires the AD4 molecule-graph keywords for flex
    ligands (ROOT/ENDROOT + TORSDOF); obabel omits them, producing
    'Unknown or inappropriate tag'. Wrap the ATOM lines when missing."""
    lines = pdbqt.read_text().splitlines()
    if any(l.startswith("ROOT") for l in lines):
        return
    torsion = 0
    for l in lines:
        if l.startswith("REMARK") and "active torsions" in l:
            try:
                torsion = int(l.split()[1])
            except (IndexError, ValueError):
                pass
    atoms = [l for l in lines if l.startswith("ATOM")]
    remarks = [l for l in lines if l.startswith("REMARK")]
    new = remarks + ["ROOT"] + atoms + ["ENDROOT", f"TORSDOF {torsion}"]
    pdbqt.write_text("\n".join(new) + "\n")


# --------------------------------------------------------------------------- #
# 6. graph node
# --------------------------------------------------------------------------- #
def prepare_target(state: dict) -> dict:
    """LangGraph node: full target preparation pipeline."""
    workdir = Path(state["project_dir"]) / "01_target"
    workdir.mkdir(parents=True, exist_ok=True)
    opts = state.get("options", {})
    brain = AgentBrain(project_dir=Path(state["project_dir"])) \
        if not opts.get("no_llm") else None

    # 1) resolve
    resolved = resolve_target(state["target"], workdir,
                              modeler=opts.get("modeler", "esmfold"))
    raw_pdb: Path = resolved["pdb_path"]

    # 2) analyze (+ structural pitfall pre-check)
    report = analyze_completeness(raw_pdb)
    issues = report.get("issues", [])
    ligands = report["ligands"]
    judge = agent_judge_target(report, brain)
    action = judge["action"]

    # auto-repair the safe fixes (multi-model merge / altloc duplication)
    # so the deterministic path never feeds a corrupt structure downstream;
    # the agent path can do more targeted repairs via repair_structure
    safe_fixes = {i.get("fix") for i in issues if i.get("fix")}
    if "dedupe_models" in safe_fixes or "keep_altloc_a" in safe_fixes:
        acts = [a for a in ("dedupe_models", "keep_altloc_a")
                if a in safe_fixes]
        repaired = workdir / "target_repaired.pdb"
        repair_stats = repair_structure(raw_pdb, repaired, actions=acts)
        report["repaired"] = repair_stats
        raw_pdb = repaired  # downstream uses the fixed file
        for i in issues:
            if i.get("fix") in acts:
                i["auto_fixed"] = True

    # 3) clean
    keep = ligands if action in ("keep_ligand", "clean_and_use") and ligands else None
    if action == "clean_and_use" and ligands:
        keep = ligands  # keep ligands but site from grid? keep simple: keep them
    clean_path = workdir / "target_clean.pdb"
    clean_stats = clean_pdb(raw_pdb, clean_path,
                            keep_resnames=ligands if ligands else None,
                            keep_waters=bool(ligands),
                            water_keep_dist=5.0)

    # 4) pocket
    lig = ligands[0] if (ligands and action != "remodel") else None
    pocket = find_pocket(clean_path, lig, brain)
    pocket["site_id"] = "S1"

    # 5) pdbqt (receptor without ligand atoms for docking; ligand kept separately)
    receptor_pdb = workdir / "receptor.pdb"
    if ligands:
        # receptor = clean pdb minus ligand residues
        _remove_res(clean_path, receptor_pdb, ligands)
    else:
        receptor_pdb = clean_path
    receptor_pdbqt = to_pdbqt(receptor_pdb, workdir / "receptor.pdbqt",
                              flex=False)  # rigid receptor: no graph keywords

    # ligand pdbqt (for positive control / MD)
    lig_pdbqt = None
    if ligands:
        lig_pdb = workdir / "ligand.pdb"
        _extract_res(clean_path, lig_pdb, ligands)
        try:
            lig_pdbqt = to_pdbqt(lig_pdb, workdir / "ligand.pdbqt")
        except RuntimeError as e:
            logger.warning(f"ligand pdbqt failed: {e}")

    out = {
        "raw_pdb": str(raw_pdb),
        "resolved": resolved,
        "completeness": report,
        "judgment": judge,
        "clean_pdb": str(clean_path),
        "clean_stats": clean_stats,
        "pocket": pocket,
        "receptor_pdb": str(receptor_pdb),
        "receptor_pdbqt": str(receptor_pdbqt),
        "ligand_pdb": str(lig_pdb) if lig_pdb is not None else None,
        "ligand_pdbqt": str(lig_pdbqt) if lig_pdbqt else None,
        "ligand_resnames": ligands,
    }
    jsave(workdir / "target_prep.json", out)
    state_out = dict(state)
    state_out["target_prep"] = out
    return state_out


def _remove_res(src: Path, dst: Path, resnames: list[str]) -> Path:
    drop = set(resnames)
    with open(src) as fh:
        lines = [l for l in fh
                 if not (l.startswith("HETATM") and l[17:20].strip() in drop)]
    ensure_parent(dst)
    dst.write_text("".join(lines))
    return dst


def _extract_res(src: Path, dst: Path, resnames: list[str]) -> Path:
    keep = set(resnames)
    with open(src) as fh:
        lines = [l for l in fh
                 if l.startswith("HETATM") and l[17:20].strip() in keep]
    lines.append("END\n")
    ensure_parent(dst)
    dst.write_text("".join(lines))
    return dst
