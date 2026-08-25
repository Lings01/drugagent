"""Module B: small-molecule virtual screening.

Pipeline: library -> standardize -> physchem/ML prefilter -> adaptive N ->
Vina parallel docking -> GNINA rescore -> agent hit criteria -> top hits.
"""
from __future__ import annotations
import re

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

from ..config import DEFAULTS, LIBRARIES, TOOLS, which_tool, resolve_defaults
from ..llm import AgentBrain
from ..utils import jsave, pmap, run_cmd

# --------------------------------------------------------------------------- #
# library sources
# --------------------------------------------------------------------------- #
LIBRARY_SOURCES = {
    # R18: default master library — local build, union of the two
    # provenance files below (scripts/merge_libraries.py). No download
    # needed when present.
    "nci_npatlas": {
        "url": None,
        "alt_urls": [],
        "size_hint": "master library: NCI-Open UNION NPAtlas, deduped by "
                     "InChIKey (~265k, ~1.4 GB; data/libraries/nci_npatlas.sdf)",
    },
    "nci_open": {
        "url": "https://dctd.cancer.gov/data-tools-biospecimens/data",
        "alt_urls": ["https://cellminer.cancer.gov/"],
        "size_hint": "NCI/DTP Open Chemical Repository (Dec-2010 release, "
                     "265,242 structures, 3D, NSC ids; 1.2 GB)",
    },
    "npatlas": {
        "url": "https://www.npatlas.com/",
        "alt_urls": [],
        "size_hint": "NPAtlas natural products (2024-09, 36,454 structures, "
                     "2D, InChI/InChIKey; 175 MB)",
    },
    "dtp": {
        # legacy: third-party mirror of the NCI/DTP open set — dead
        # (502 since R18); kept for --library dtp compat (falls back).
        "url": "http://www.dtpbase.org/download/All_Druglike_Compounds.sdf.gz",
        "alt_urls": [
            "http://www.dtpbase.org/download",
        ],
        "size_hint": "351k compounds (~1.5 GB) [mirror dead, use nci_npatlas]",
    },
    "chembl": {
        # date dir discovered at setup; fallback pattern
        "url": "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/",
        "size_hint": "ChEMBL molecules SDF (~2 GB)",
    },
    "pdbbind": {
        "url": "http://www.pdbbind.org.cn/pdbbind/DownloadData/PDBBind_v20201216.tar.gz",
        "size_hint": "PDBBind 2.015 (~30 MB)",
    },
    "custom": {"url": None, "size_hint": "user-provided SDF"},
}


def library_path(name: str) -> Path:
    return LIBRARIES / f"{name}.sdf"


def resolve_library(opts: dict, d, base: Path | None = None) -> tuple[Path, str]:
    """R10/G6: resolve the screening library with a fallback chain.

    The DTP/PDBBind mirrors are flaky (a 0-byte download is a known state
    of this box); a missing/corrupt < 1 MB library should not kill a run.
    Fall back to the local master library (nci_npatlas), then the ChEMBL 35
    subsample (50k mols), and report it in
    the returned label so state/report make the substitution explicit.
    Custom SDF paths (opts["library_path"]) do NOT fall back — a user
    pointing at a file means exactly that file.
    Returns (path, label)."""
    base = Path(base) if base else LIBRARIES
    lib_name = opts.get("library", d.default_library)
    if opts.get("library_path"):
        p = Path(opts["library_path"])
        if p.exists() and p.stat().st_size > 1_000_000:
            return p, lib_name
        raise FileNotFoundError(
            f"library not found (or < 1 MB): {p}")
    lib = base / f"{lib_name}.sdf"
    if lib.exists() and lib.stat().st_size > 1_000_000:
        return lib, lib_name
    for cand in ("nci_npatlas", "chembl35_small", "chembl35"):
        if cand == lib_name:
            continue
        c = base / f"{cand}.sdf"
        if c.exists() and c.stat().st_size > 1_000_000:
            logger.warning(
                f"library {lib_name} missing/corrupt ({lib}, "
                f"{lib.stat().st_size if lib.exists() else 0} bytes) — "
                f"falling back to {c.name} "
                f"({c.stat().st_size / 1e6:.0f} MB)")
            return c, f"{cand} (fallback for {lib_name})"
    raise FileNotFoundError(
        f"library not found: {lib} (run `drugagent setup --libraries`)")


# --------------------------------------------------------------------------- #
# standardization
# --------------------------------------------------------------------------- #
def _make_salt_remover():
    """Build a SaltRemover; falls back to identity when no salt pickle is
    available (RDKit >= 2024 changed the constructor signature)."""
    from rdkit.Chem.SaltRemover import SaltRemover
    import rdkit.Chem.SaltRemover as _sr_mod
    pkl = Path(_sr_mod.__file__).parent / "Salts.pkl"
    if pkl.exists():
        for kwargs in ({"defnFilename": str(pkl), "useMolBlockWts": True},
                       {"defnFilename": str(pkl)}, {}):
            try:
                return SaltRemover(**kwargs)
            except (TypeError, ValueError):
                continue
    try:
        return SaltRemover(useMolBlockWts=True)
    except TypeError:
        class _Noop:
            def StripMol(self, m):
                return m
        return _Noop()


def _parse_sdf_blocks(sdf_path: Path) -> list[tuple[str, dict, str]]:
    """Parse SDF blocks into (title, properties, molblock). Tolerant to
    vendor molfile dialects (e.g. ChEMBL) that RDKit\'s SDMolSupplier rejects."""
    out = []
    title, props, block_lines = "", {}, []
    i = 0
    with open(sdf_path) as fh:
        lines = fh.readlines()
    while i < len(lines):
        line = lines[i]
        if line.rstrip("\n") == "$$$$":
            if title or block_lines:
                out.append((title, props, "".join(block_lines)))
            title, props, block_lines = "", {}, []
            i += 1
            continue
        stripped = line.rstrip("\n")
        if not block_lines and not title and stripped:
            title = stripped
            # R18: the title line is part of the molblock — MolFromMolBlock
            # needs the canonical 3-line pre-header (title, comment, blank)
            # before the V2000/V3000 counts line, so keep it in the block
            # text (previously dropped; masked before because every
            # deployed library carried canonical_smiles and the molblock
            # path was never exercised).
            block_lines.append(line)
            i += 1
            continue
        if stripped.lstrip().startswith(">  <"):
            key = stripped.lstrip()[4:].rstrip(">").strip()
            i += 1
            val_lines = []
            while i < len(lines) and lines[i].rstrip("\n") != "$$$$":
                val_lines.append(lines[i].rstrip("\n"))
                i += 1
            props[key] = "\n".join(val_lines)
            continue
        if block_lines and not stripped.startswith("  ") and not title:
            pass
        block_lines.append(line)
        i += 1
    if title or block_lines:
        out.append((title, props, "".join(block_lines)))
    return out


def standardize_sdf(sdf_path: Path, out_path: Path, *,
                    max_heavy_atoms: int = 60, n_jobs: int = 1) -> pd.DataFrame:
    """Sanitize + 3D embed. Returns DataFrame with SMILES/status columns.

    Molecules are parsed from the block property ``canonical_smiles`` when
    present (robust to vendor molfile dialects), else from the molblock.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    remover = _make_salt_remover()
    blocks = _parse_sdf_blocks(sdf_path)

    def _one(item):
        i, (title, props, molblock) = item
        smi = props.get("canonical_smiles") or props.get("smiles")
        mol = Chem.MolFromSmiles(smi) if smi else None
        src = "smiles"
        if mol is None and molblock.strip():
            try:
                mol = Chem.MolFromMolBlock(molblock, removeHs=False,
                                           sanitize=False)
            except Exception:  # noqa: BLE001
                mol = None
            src = "molblock"
        if mol is None:
            return (i, None, "null")
        try:
            mol = remover.StripMol(mol)
            mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))  # sanitize round-trip
            if mol is None:
                return (i, None, "sanitize_fail")
            n_hv = mol.GetNumHeavyAtoms()
            if n_hv < 4 or n_hv > max_heavy_atoms:
                return (i, None, "size")
            molH = Chem.AddHs(mol)
            try:
                ok = AllChem.EmbedMolecule(molH, randomSeed=42 + i % 1000) == 0
                AllChem.MMFFOptimizeMolecule(molH)
            except Exception:  # noqa: BLE001
                ok = False
            if not ok:
                return (i, None, "embed_fail")
            return (i, Chem.MolToSmiles(mol), "ok")
        except Exception:  # noqa: BLE001
            return (i, None, "error")

    results = pmap(_one, enumerate(blocks), n_jobs=n_jobs)
    rows = [{"idx": i, "smiles": smi, "status": status}
            for i, smi, status in results]
    df = pd.DataFrame(rows)
    jsave(out_path.with_suffix(".stats.json"),
          {"n_total": len(df),
           "n_ok": int((df.status == "ok").sum()),
           "fails": df[~df.status.isin(["ok"])].status.value_counts().to_dict()})
    return df


def compute_features(df: pd.DataFrame, smiles_col: str = "smiles",
                     n_jobs: int = 1) -> pd.DataFrame:
    """Physchem + ML features (RDKit + Mordred)."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, QED
    try:
        from rdkit.Chem.SA_Score import sascore
    except ImportError:  # RDKit >= 2024 moved SA_Score out of core
        sascore = None

    RDLogger.DisableLog("rdApp.*")

    def _feats(row):
        if not row.get(smiles_col):
            return None
        m = Chem.MolFromSmiles(row[smiles_col])
        if m is None:
            return None
        f = {
            "MW": Descriptors.MolWt(m),
            "LogP": Descriptors.MolLogP(m),
            "HBA": Descriptors.NumHAcceptors(m),
            "HBD": Descriptors.NumHDonors(m),
            "TPSA": Descriptors.TPSA(m),
            "RotBonds": Descriptors.NumRotatableBonds(m),
            "AromRings": Descriptors.NumAromaticRings(m),
            "QED": QED.qed(m),
        }
        try:
            if sascore is not None:
                f["SA"] = sascore.calculateScore(m)
            else:
                f["SA"] = 3.0
        except Exception:  # noqa: BLE001
            f["SA"] = 5.0
        return f

    if len(df) == 0:
        return pd.DataFrame(columns=["idx", smiles_col, "MW", "LogP",
                                     "HBA", "HBD", "TPSA", "RotBonds",
                                     "AromRings", "QED", "SA"])
    idxs = list(df["idx"]) if "idx" in df.columns else list(range(len(df)))
    smiles = (list(df[smiles_col]) if smiles_col in df.columns
              else [None] * len(df))
    feats = df.apply(_feats, axis=1)
    rows = [(i, s, f) for i, s, f in zip(idxs, smiles, feats)
            if f is not None]
    fdf = pd.DataFrame([f for _, _, f in rows])
    if len(fdf):
        fdf.insert(0, "idx", [i for i, _, _ in rows])
        fdf.insert(1, smiles_col, [s for _, s, _ in rows])
    return fdf


def physchem_filter(fdf: pd.DataFrame) -> pd.DataFrame:
    """Soft Lipinski + synthetic-access filter (keeps a permissive set)."""
    m = (
        (fdf.MW >= 150) & (fdf.MW <= 650)
        & (fdf.LogP >= -2) & (fdf.LogP <= 6.5)
        & (fdf.HBA <= 14) & (fdf.HBD <= 7)
        & (fdf.TPSA <= 160)
        & (fdf.RotBonds <= 15)
        & (fdf.SA <= 5.5)
    )
    return fdf[m].copy()


def ml_prefilter(fdf: pd.DataFrame, model_path: Path | None = None,
                 n_keep: int | None = None) -> pd.DataFrame:
    """GBM/RF ranking on Mordred descriptors. Falls back to QED+SA blend."""
    if model_path is not None and model_path.exists():
        import joblib
        model = joblib.load(model_path)
        try:
            from mordred import Calculator, descriptors
            calc = Calculator(descriptors, ignore_3D=True)
            X = calc.on_df(fdf.set_index("idx"))
            X = X.loc[fdf["idx"]]
            X = X.replace([np.inf, -np.inf], np.nan).fillna(0.0)
            score = model.predict_proba(X)[:, 1] if getattr(model, "classes_", None) is not None \
                else model.predict(X)
            fdf = fdf.assign(ml_score=score)
            fdf = fdf.sort_values("ml_score", ascending=False)
            if n_keep:
                fdf = fdf.head(n_keep)
            return fdf
        except Exception as e:  # noqa: BLE001
            logger.warning(f"ML prefilter failed ({e}); using QED+SA fallback")
    fdf = fdf.assign(ml_score=fdf.QED * 0.7 + (5.0 - fdf.SA) / 5.0 * 0.3)
    fdf = fdf.sort_values("ml_score", ascending=False)
    if n_keep:
        fdf = fdf.head(n_keep)
    return fdf


# --------------------------------------------------------------------------- #
# docking
# --------------------------------------------------------------------------- #
def vina_bin() -> Path:
    p = which_tool("vina", [TOOLS / "vina" / "bin", ENV_BIN()])
    if p is None:
        raise FileNotFoundError("vina not found (run `drugagent setup`)")
    return Path(p)


def ENV_BIN() -> Path:
    from ..config import ENV_DIR
    return ENV_DIR / "bin"


def gnina_bin() -> Path | None:
    p = which_tool("gnina", [TOOLS / "gnina", ENV_BIN()])
    return Path(p) if p else None


def write_ligand_pdbqt(smi: str, out: Path, seed: int = 42) -> bool:
    """SMILES -> 3D PDBQT (RDKit PDBQTWriterLib, else Meeko).

    Returns False (no exception) when the ligand cannot be prepared.
    Ligands keep ROOT/BRANCH sections (this vina build requires them);
    receptors have them stripped in to_pdbqt.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return False
        # keep the largest fragment (salts/dimers from vendor libraries)
        frags = Chem.GetMolFrags(m, asMols=True)
        if len(frags) > 1:
            m = max(frags, key=lambda f: f.GetNumHeavyAtoms())
        mH = Chem.AddHs(m)
        if AllChem.EmbedMolecule(mH, randomSeed=seed) != 0:
            return False
        AllChem.MMFFOptimizeMolecule(mH)
        Chem.AssignStereochemistry(mH, cleanIt=True, force=True)
        out = Path(out)
        out.parent.mkdir(parents=True, exist_ok=True)
        try:
            from rdkit.Chem.PDBQTWriterLib import PDBQTWriter
            with open(out, "w") as fh:
                PDBQTWriter.WritePDBQT(mH, fh, {"hydrogens": True})
        except ImportError:
            # RDKit >= 2024: PDBQTWriterLib is not in the wheel; use Meeko
            from meeko import MoleculePreparation, PDBQTWriterLegacy
            setups = MoleculePreparation().prepare(mH)
            text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
            if not ok:
                raise RuntimeError(f"meeko pdbqt write failed: {err}")
            out.write_text(text)
        return True
    except Exception:  # noqa: BLE001
        return False
def dock_one(args: tuple) -> dict:
    """Vina one ligand. args = (receptor_pdbqt, lig_pdbqt, out_prefix, pocket,
    exh, n_jobs[, flex_pdbqt]) — 7th element enables Vina --flex."""
    rec, lig, prefix, pocket, exh, cpu = args[:6]
    flex = args[6] if len(args) > 6 else None
    out = {
        "smiles": Path(lig).stem,
        "score": np.nan, "rmsd_lb": np.nan, "rmsd_ub": np.nan,
        "ok": False, "top_pose_pdbqt": None,
    }
    # defense: this vina build rejects molecule-graph keywords (ROOT/
    # TORSDOF) in RIGID receptors; a stale/external receptor file may
    # still carry them (e.g. built with to_pdbqt flex=True)
    rec = Path(rec)
    if any(l.split()[:1] == ["ROOT"] for l in rec.read_text().splitlines()):
        clean = rec.with_name(rec.stem + "_nograph.pdbqt")
        kept = [l for l in rec.read_text().splitlines()
                if l.split()[:1] not in (["ROOT"], ["ENDROOT"],
                                         ["BRANCH"], ["ENDBRANCH"],
                                         ["TORSDOF"])]
        clean.write_text("\n".join(kept) + "\n")
        rec = clean
    # cache: reuse a previous successful dock (deterministic seed=42)
    logp = Path(str(prefix) + ".log")
    posep = Path(str(prefix) + ".pdbqt")
    if logp.is_file() and posep.is_file():
        try:
            txt = logp.read_text()
            for l in txt.splitlines():
                toks = l.split()
                if len(toks) >= 4 and toks[0] == "1":
                    out["score"] = float(toks[1])
                    out["rmsd_lb"] = float(toks[2])
                    out["rmsd_ub"] = float(toks[3])
                    out["ok"] = True
                    out["top_pose_pdbqt"] = str(posep)
                    return out
        except Exception:  # noqa: BLE001
            pass
    # R10: --cpu beyond the physical core count oversubscribes and
    # (worse with --flex) inflates RSS ~linearly per thread (measured:
    # 37 flex residues at 21 cpu -> ~23 GB). Cap at the core count and
    # warn when the caller asked for more.
    import os
    phys = os.cpu_count() or 1
    if cpu > phys:
        logger.warning(f"vina --cpu {cpu} > {phys} physical cores; "
                       f"capping (oversubscription slows the run and "
                       f"RSS grows ~linearly, especially with --flex)")
        cpu = phys
    cmd = [str(vina_bin()),
           "--receptor", str(rec), "--ligand", str(lig)]
    if flex:
        cmd += ["--flex", str(flex)]
    cmd += [
           "--center_x", str(pocket["center"][0]),
           "--center_y", str(pocket["center"][1]),
           "--center_z", str(pocket["center"][2]),
           "--size_x", str(pocket["xsize"]),
           "--size_y", str(pocket["ysize"]),
           "--size_z", str(pocket["zsize"]),
           "--exhaustiveness", str(exh),
           "--out", str(prefix + ".pdbqt"),
           "--cpu", str(cpu), "--seed", "42"]
    try:
        proc = run_cmd(cmd, check=False)
        txt = (proc.stdout.decode(errors="replace")
               if proc.stdout is not None else "")
        with open(prefix + ".log", "w") as fh:
            fh.write(txt)
    except Exception as e:  # noqa: BLE001
        out["error"] = str(e)
        return out
    try:
        # score table: first row whose first token is "1" and second is a float
        lines = []
        for l in txt.splitlines():
            toks = l.split()
            if len(toks) >= 4 and toks[0] == "1":
                try:
                    float(toks[1])
                    lines.append(l)
                    break
                except ValueError:
                    continue
        if lines:
            parts = lines[0].split()
            # 1  -XX.X  rmsd_lb  rmsd_ub
            out["score"] = float(parts[1])
            out["rmsd_lb"] = float(parts[2])
            out["rmsd_ub"] = float(parts[3])
            out["ok"] = True
            out["top_pose_pdbqt"] = prefix + ".pdbqt"
    except Exception as e:  # noqa: BLE001
        out["error"] = f"log parse: {e}"
    return out


def gnina_rescore(rec_pdbqt: str, pose_pdbqt: str, pocket: dict) -> dict:
    """GNINA scoring-only on a pose. Returns {gnina_score, cn, ok}."""
    g = gnina_bin()
    if g is None:
        return {"gnina_score": np.nan, "cn": np.nan, "ok": False}
    out = {"gnina_score": np.nan, "cn": np.nan, "ok": False}
    # GNINA ELF needs CUDA libs from the env's nvidia pip packages
    import os as _os
    from ..config import ENV_DIR
    envd = None
    nv = ENV_DIR / "lib" / "python3.12" / "site-packages" / "nvidia"
    if nv.is_dir():
        libs = sorted(str(d) for d in nv.glob("*/lib") if d.is_dir())
        if libs:
            envd = dict(_os.environ,
                        LD_LIBRARY_PATH=_os.pathsep.join(
                            libs + [_os.environ.get("LD_LIBRARY_PATH", "")]))
    cmd = [str(g), "score",
           "--receptor", rec_pdbqt, "--complex", pose_pdbqt,
           "--center_x", str(pocket["center"][0]),
           "--center_y", str(pocket["center"][1]),
           "--center_z", str(pocket["center"][2]),
           "--size_x", str(pocket["xsize"]),
           "--size_y", str(pocket["ysize"]),
           "--size_z", str(pocket["zsize"]),
           "--exhaustiveness", "16", "--cpu", "4"]
    try:
        proc = run_cmd(cmd, check=False, env=envd)
        txt = (proc.stdout or b"").decode(errors="replace")
        import re as _re
        m = _re.search(r"GNINA Score:?\s*(-?[\d.]+)", txt)
        c = _re.search(r"CN:?\s*(-?[\d.]+)", txt)
        if m:
            out["gnina_score"] = float(m.group(1))
            out["cn"] = float(c.group(1)) if c else np.nan
            out["ok"] = True
    except Exception as e:  # noqa: BLE001
        out["error"] = str(e)
    return out


# --------------------------------------------------------------------------- #
# hit criteria (agent judgment)
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# R2: receptor conformation selection + flexible docking
# --------------------------------------------------------------------------- #
BACKBONE_ATOMS = {"N", "CA", "C", "O"}
SOLVENT_RES = {"HOH", "WAT", "SOL", "TIP3", "SPC", "NA", "CL", "K", "MG",
               "ZN", "CA", "NA+", "CL-", "K+"}


def consensus_stats(scores: list[float]) -> dict:
    """Robust consensus over per-conformer docking scores (kcal/mol, lower
    = better). Mean is the consensus value: one bad conformer should not
    dominate, but the best-case (min) is reported for reference."""
    vals = [float(x) for x in scores]
    return {"mean": round(float(np.mean(vals)), 3),
            "min": round(min(vals), 3),
            "max": round(max(vals), 3),
            "n": len(vals)}


def _pdb_serial(line: str) -> int:
    """Atom serial from the first integer after the record name; tolerant of
    both PDB (5-col field) and AD4-PDBQT (6-col field) layouts."""
    m = re.match(r"\s*\d+", line[4:14])
    return int(m.group(0)) if m else 0


def _atoms(pdb: Path) -> list[dict]:
    """Parse ATOM/HETATM lines -> [{name, resname, chain, resseq, xyz, serial,
    line}]. Works for both PDB and AD4-PDBQT serial layouts (the serial is
    read from a tolerant window)."""
    out = []
    with open(pdb) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                out.append({
                    "name": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "chain": line[21],
                    "resseq": int(line[22:26]),
                    "xyz": (float(line[30:38]), float(line[38:46]),
                            float(line[46:54])),
                    "serial": _pdb_serial(line),
                    "line": line.rstrip("\n"),
                })
            except (ValueError, IndexError):
                continue
    return out


def _renumber_pdbqt_line(line: str, serial: int) -> str:
    """Replace the atom serial of a PDBQT line (AD4 layout: 7-char field at
    cols 5-11; PDB layout: 6-char field) with a new local serial."""
    if line.startswith("HETATM"):
        pre, rest = line[:6], line[13:]
    else:
        pre, rest = line[:4], line[11:]
    return pre + f"{serial:7d}" + rest


def _flex_residue_block(resname: str, chain: str, resseq: int,
                        root: dict, side: list[dict]):
    """Emit one AD4 flex-residue block (BEGIN_RES ... END_RES) for a side
    chain rooted at `root` (normally CA). The AD4 torsion graph is derived
    from 3D connectivity (bonded = < 1.9 A within the residue): every bond
    from the root becomes a rotatable branch, so all side-chain torsions
    are sampled. Returns (lines, n_local_atoms) or (None, 0)."""
    atoms = [root] + side
    n = len(atoms)
    adj: list[list[int]] = [[] for _ in range(n)]
    for i in range(n):
        xi, yi, zi = atoms[i]["xyz"]
        for j in range(i + 1, n):
            xj, yj, zj = atoms[j]["xyz"]
            d2 = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            if d2 < 1.9 ** 2:
                adj[i].append(j)
                adj[j].append(i)
    # BFS tree from the root (index 0)
    parent = {0: -1}
    queue = [0]
    while queue:
        i = queue.pop(0)
        for j in adj[i]:
            if j not in parent:
                parent[j] = i
                queue.append(j)
    for i in range(n):  # strays (defensive): attach to nearest tree atom
        if i not in parent:
            xi, yi, zi = atoms[i]["xyz"]
            best = min((t for t in parent),
                       key=lambda t: (xi - atoms[t]["xyz"][0]) ** 2
                       + (yi - atoms[t]["xyz"][1]) ** 2
                       + (zi - atoms[t]["xyz"][2]) ** 2)
            parent[i] = best
    children: dict[int, list[int]] = {}
    for j, p in parent.items():
        children.setdefault(p, []).append(j)
    lines: list[str] = []
    counter = [0]
    index_map: list[tuple[int, int]] = []

    def newnum(i: int) -> int:
        counter[0] += 1
        index_map.append((atoms[i]["serial"], counter[0]))
        return counter[0]

    def emit_branch(parent_n: int, i: int) -> None:
        cn = newnum(i)
        lines.append(f"BRANCH {parent_n} {cn}")
        lines.append(_renumber_pdbqt_line(atoms[i]["line"], cn))
        for g in children.get(i, []):
            if children.get(g):
                emit_branch(cn, g)
            else:  # leaf: listed directly in this branch (AD4 convention)
                gn = newnum(g)
                lines.append(_renumber_pdbqt_line(atoms[g]["line"], gn))
        lines.append(f"ENDBRANCH {parent_n} {cn}")

    root_n = newnum(0)
    lines.insert(0, f"BEGIN_RES {resname} {chain} {resseq}")
    lines.append("REMARK INDEX MAP "
                 + " ".join(f"{o} {l}" for o, l in index_map) + " 0")
    lines.append("ROOT")
    lines.append(_renumber_pdbqt_line(atoms[0]["line"], root_n))
    lines.append("ENDROOT")
    for c in children.get(0, []):
        emit_branch(root_n, c)  # every root child gets its own branch
    lines.append(f"END_RES {resname} {chain} {resseq}")
    return lines, counter[0]


def kabsch(moving: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Least-squares superposition of `moving` onto `target` (same atom order,
    n x 3). Returns (R, t) with aligned = moving @ R.T + t."""
    m = np.asarray(moving, dtype=float)
    t = np.asarray(target, dtype=float)
    cm, ct = m.mean(axis=0), t.mean(axis=0)
    h = (m - cm).T @ (t - ct)
    u, _, vt = np.linalg.svd(h)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:  # reflect -> rotate
        vt[-1, :] *= -1
        r = vt.T @ u.T
    return r, ct - r @ cm


def _collect_ca_pairs(entries, ref_by) -> list[tuple[int, tuple, tuple]]:
    """C-alpha candidate pairs via per-(model chain, ref chain) offset
    consensus with a bijective chain assignment.

    For every (model chain, ref chain) pair, each numbering offset is
    scored by how many resname identities it explains; the best offset
    (>= 3 matches) is kept, absorbing renumbered/merged chains (e.g. a
    binder appended to a target chain). Because homodimers let one model
    chain match several ref chains, the model chains are assigned to
    distinct ref chains (permutation search, few chains expected) so that
    symmetric copies cannot flood the global fit.
    """
    from collections import defaultdict
    import itertools
    counts: "defaultdict[tuple, dict]" = defaultdict(lambda: defaultdict(int))
    for idx, (l, a_) in enumerate(entries):
        if a_ is None or a_["name"] != "CA":
            continue
        for (rc, rn), seqs in ref_by.items():
            if rn != a_["resname"]:
                continue
            for rseq in seqs:
                counts[(a_["chain"], rc)][rseq - a_["resseq"]] += 1
    cand: dict[tuple, int] = {}
    for (mc, rc), offs in counts.items():
        off, n = max(offs.items(), key=lambda kv: kv[1])
        if n >= 3:
            cand[(mc, rc)] = off

    def pairs_for(mc: str, rc: str, off: int):
        ps = []
        for idx, (l, a_) in enumerate(entries):
            if a_ is None or a_["name"] != "CA" or a_["chain"] != mc:
                continue
            rxyz = ref_by.get((rc, a_["resname"]), {}).get(a_["resseq"] + off)
            if rxyz is not None:
                ps.append((idx, a_["xyz"], rxyz))
        return ps

    mc_list = sorted({c[0] for c in cand})
    rc_list = sorted({c[1] for c in cand})
    if not mc_list or not rc_list or len(rc_list) < len(mc_list):
        return []
    best: list[tuple[int, tuple, tuple]] = []
    best_rmsd = float("inf")
    for perm in itertools.permutations(rc_list, len(mc_list)):
        ps: list[tuple[int, tuple, tuple]] = []
        for mc, rc in zip(mc_list, perm):
            if (mc, rc) not in cand:
                ps = []
                break
            ps.extend(pairs_for(mc, rc, cand[(mc, rc)]))
        if len(ps) < 3:
            continue
        m = np.array([q[1] for q in ps])
        tg = np.array([q[2] for q in ps])
        R, tr = kabsch(m, tg)
        rmsd = float(np.sqrt(
            np.mean(((m @ R.T + tr - tg) ** 2).sum(axis=1))))
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best = ps
    return best


def _read_entries(pdb: Path) -> list[tuple[str, dict | None]]:
    entries = []
    for l in Path(pdb).read_text().splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            entries.append((l, None))
            continue
        a_ = {
            "name": l[12:16].strip(),
            "resname": l[17:20].strip(),
            "chain": l[21],
            "resseq": int(l[22:26]),
            "xyz": (float(l[30:38]), float(l[38:46]), float(l[46:54])),
        }
        entries.append((l, a_))
    return entries


def _ref_ca_index(ref_pdb: Path) -> "defaultdict[tuple, dict]":
    from collections import defaultdict
    ref_by: "defaultdict[tuple, dict]" = defaultdict(dict)
    for a_ in _atoms(ref_pdb):
        if a_["name"] == "CA":
            ref_by[(a_["chain"], a_["resname"])][a_["resseq"]] = a_["xyz"]
    return ref_by


def _cryst1_box(pdb: Path) -> np.ndarray | None:
    """Box (a, b, c) from a PDB CRYST1 line, or None."""
    for l in Path(pdb).read_text().splitlines():
        if l.startswith("CRYST1"):
            try:
                return np.array([float(l[6:15]), float(l[15:24]),
                                 float(l[24:33])])
            except ValueError:
                return None
    return None


def _mi_kabsch(moving: np.ndarray, target: np.ndarray,
               box: np.ndarray | None,
               iters: int = 10) -> tuple[np.ndarray, np.ndarray]:
    """R8: Kabsch with per-atom minimum-image refinement (PBC-aware).

    MD-rep PDBs (simulation-box frame) and crystal PDBs
    (crystallographic frame) differ by per-atom INTEGER box vectors;
    a plain least-squares Kabsch averages those vectors into a muddled
    transform (190 A residual on the 1HVI dimer). The image vectors
    live in the MODEL (simulation) frame, where the box lattice is
    axis-aligned, so the internal transform maps target -> model and
    alternates (1) per-atom image selection with the current transform
    and (2) a rigid refit against the selected images, accepting a step
    only while the model-frame score improves (atoms near a box
    boundary can flip images and make a naive loop oscillate). A
    compact core (consecutive residues share one box vector) provides
    an additional clean start. Returns the MODEL -> TARGET transform
    (the inverse), i.e. (R, t) with moving @ R.T + t ~= target."""
    moving = np.asarray(moving, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    if box is None:
        return kabsch(moving, target)
    box = np.asarray(box, dtype=np.float64)

    def score_model(r: np.ndarray, t: np.ndarray) -> float:
        # per-atom model-frame residual of the target under (r, t)
        gm = target @ r.T + t
        diff = gm - moving
        diff -= box * np.round(diff / box)
        return float(np.sqrt((diff * diff).sum(axis=1)).mean())

    def refine(r: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
        cur = score_model(r, t)
        for _ in range(iters):
            gm = target @ r.T + t
            diff = gm - moving
            k = np.round(diff / box)
            s = gm - k * box  # image of each target atom near its partner
            r2, t2 = kabsch(target, s)
            s2 = score_model(r2, t2)
            if s2 >= cur - 1e-6:
                break
            r, t, cur = r2, t2, s2
        return r, t, cur

    # starts: the raw full fit (target -> model) + compact cores
    starts = [kabsch(target, moving)]
    n = len(moving)
    if n >= 12:
        starts.append(kabsch(target[:10], moving[:10]))
        mid = n // 2
        starts.append(kabsch(target[mid:mid + 10], moving[mid:mid + 10]))
    best = None
    for st in starts:
        out = refine(*st)
        if best is None or out[2] < best[2]:
            best = out
    r, t = best[0], best[1]
    # invert: model -> target
    return r.T, -t @ r


def _ca_pair_list(mdl_pdb: Path, ref_pdb: Path):
    """C-alpha partner pairs (model index, model xyz, ref xyz) across all
    chain numbering offsets (see _collect_ca_pairs)."""
    ref_by = _ref_ca_index(ref_pdb)
    entries = _read_entries(mdl_pdb)
    return _collect_ca_pairs(entries, ref_by)


def kabsch_transform(mdl_pdb: Path, ref_pdb: Path):
    """Rigid transform (R, t) taking `mdl_pdb` into the frame of `ref_pdb`,
    from C-alpha pairs matched across all numbering offsets per (chain,
    resname), refined by an iterative 3 A filter. Returns (n_pairs, R, t);
    (0, None, None) if < 3 pairs survive."""
    from collections import defaultdict
    pairs = _ca_pair_list(mdl_pdb, ref_pdb)
    if len(pairs) < 3:
        return 0, None, None
    active = list(range(len(pairs)))
    box = _cryst1_box(mdl_pdb)
    m = np.array([pairs[i][1] for i in active])
    tg = np.array([pairs[i][2] for i in active])
    r, t = _mi_kabsch(m, tg, box)
    d = m @ r.T + t - tg
    if box is not None:
        br = np.diag(box) @ r
        d = d - np.round(d @ np.linalg.inv(br)) @ br
    d = np.sqrt((d * d).sum(axis=1))
    med = float(np.median(d))
    if med < 10.0:
        # consistent fit: trim outliers and refit once
        keep = [i for i, di in zip(active, d) if di < max(3.0, 2.0 * med)]
        if len(keep) >= 3:
            m = np.array([pairs[i][1] for i in keep])
            tg = np.array([pairs[i][2] for i in keep])
            r, t = _mi_kabsch(m, tg, box)
            active = keep
    return len(active), r, t


def _ca_by_residue(entries: list[tuple[str, dict | None]]) -> dict[tuple, tuple]:
    """(chain, resseq) -> C-alpha xyz of that residue."""
    out = {}
    for l, a_ in entries:
        if a_ is not None and a_["name"] == "CA":
            out[(a_["chain"], a_["resseq"])] = a_["xyz"]
    return out


def align_pdb_to_reference(mdl_pdb: Path, ref_pdb: Path, out: Path) -> int:
    """Superpose `mdl_pdb` onto `ref_pdb` (see kabsch_transform) and write
    the transformed model to `out`. Returns the matched C-alpha pair count.

    MD cluster representatives live in the simulation-box frame and often
    carry renumbered/merged chains; the docking grid is in the crystal
    frame, so MD conformers must be aligned (or the pocket transformed)
    before redocking.
    """
    n, r, t = kabsch_transform(mdl_pdb, ref_pdb)
    entries = _read_entries(mdl_pdb)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if n < 3 or r is None:
        out.write_text("\n".join(l for l, _ in entries) + "\n")
        return n
    # R8: PBC-frame compaction. The MI-Kabsch fit is exact only up to
    # PER-ATOM integer box vectors (the model lives in the
    # simulation-box frame, the reference in the crystallographic
    # frame), so a single rigid transform leaves the written atoms in
    # several box copies. Anchor instead on the C-alpha PARTNERS: each
    # model atom is written at its reference partner position plus the
    # local (intra-residue) offset measured in the model frame — no
    # box ambiguity (local offsets are << box/2), and the per-atom Cα
    # deviation carries the conformational change. Falls back to the
    # rigid transform for atoms without a usable partner.
    box = _cryst1_box(mdl_pdb)
    ca_res = _ca_by_residue(entries)
    anchor: dict[tuple, tuple] = {}
    if box is not None:
        for idx, mxyz_, rxyz_ in _ca_pair_list(mdl_pdb, ref_pdb):
            dev = np.array(mxyz_) - np.array(rxyz_)
            dev = dev - box * np.round(dev / box)  # model-frame MI
            # MI is valid for any offset < box/2; 10 A is a sanity
            # cap for genuinely displaced (unwrapped) atoms
            if np.linalg.norm(dev) < min(10.0, box.min() / 4):
                a_ = entries[idx][1]
                anchor[(a_["chain"], a_["resseq"])] = (np.array(rxyz_), dev)
    new_lines = []
    for l, a_ in entries:
        if a_ is None:
            new_lines.append(l)
            continue
        p_ = np.array(a_["xyz"])
        placed = False
        if (a_["chain"], a_["resseq"]) in anchor:
            rxyz_, dev = anchor[(a_["chain"], a_["resseq"])]
            # local offset of this atom from its Cα (intra-residue,
            # no box vector) + the anchored Cα position
            placed_xyz = rxyz_ + dev + (p_ - np.array(
                ca_res[(a_["chain"], a_["resseq"])]))
            new_lines.append(
                l[:30] + f"{placed_xyz[0]:8.3f}" + f"{placed_xyz[1]:8.3f}"
                + f"{placed_xyz[2]:8.3f}" + l[54:])
            placed = True
        if not placed:
            p_ = p_ @ r.T + t
            new_lines.append(
                l[:30] + f"{p_[0]:8.3f}" + f"{p_[1]:8.3f}" + f"{p_[2]:8.3f}" + l[54:])
    out.write_text("\n".join(new_lines) + "\n")
    return n


def flex_sidechain_pdbqt(rec_pdbqt: Path, lig_pdb: Path, out: Path,
                         cutoff: float = 5.0) -> dict:
    """Build the Vina --flex file for this build (AD4 lineage): one
    BEGIN_RES/END_RES block per selected residue with an AD4 torsion graph.

    This vina's flex parser (parse_pdbqt_flex) accepts ONLY BEGIN_RES blocks
    at top level — plain ATOM lists or MODEL tags throw 'Unknown or
    inappropriate tag'. Atom lines keep their PDBQT coordinates/charges and
    are renumbered locally within each residue. Side chains only: the
    backbone stays rigid; CA (or the closest backbone atom) is the immobile
    graph root. `rec_pdbqt` must be the obabel-converted PDBQT of the SAME
    receptor that will be docked (flex atoms anchor by coordinates)."""
    lig = [a for a in _atoms(lig_pdb) if a["name"] not in ("", "H")]
    if not lig:
        raise ValueError(f"no atoms found in ligand PDB {lig_pdb}")
    rec = [a for a in _atoms(rec_pdbqt)
           if a["resname"] not in SOLVENT_RES]
    sel: set[tuple[str, int, str]] = set()
    for a in rec:
        key = (a["chain"], a["resseq"], a["resname"])
        for b in lig:
            dx = a["xyz"][0] - b["xyz"][0]
            dy = a["xyz"][1] - b["xyz"][1]
            dz = a["xyz"][2] - b["xyz"][2]
            if dx * dx + dy * dy + dz * dz < cutoff * cutoff:
                sel.add(key)
                break
    out_lines: list[str] = []
    n_atoms = 0
    n_residues = 0
    for chain, resseq, resname in sorted(sel,
                                         key=lambda k: (k[0], k[1], k[2])):
        res_atoms = [a for a in rec
                     if (a["chain"], a["resseq"], a["resname"]) == (chain,
                                                                    resseq,
                                                                    resname)]
        side = [a for a in res_atoms
                if a["name"] not in BACKBONE_ATOMS]
        if not side:
            continue
        root = next((a for a in res_atoms if a["name"] == "CA"), side[0])
        block, n = _flex_residue_block(resname, chain, resseq, root, side)
        if not block:
            continue
        out_lines += block
        n_atoms += n
        n_residues += 1
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(out_lines) + "\n")
    return {"pdbqt": str(out), "n_residues": n_residues,
            "n_atoms": n_atoms, "cutoff": cutoff,
            "residues": sorted(f"{r} {c} {s}" for s, c, r in sel)}


def _rep_files(repdir: Path) -> tuple[Path, Path]:
    """(tpr, xtc) for representative extraction: prefer the R6 merged
    trajectory md_all.xtc — after auto-extension the clustering was run on
    the merged time axis, so cluster first-frame times must be dumped from
    it (they may lie beyond md.xtc's end)."""
    xtc = repdir / "md_all.xtc"
    if not xtc.is_file():
        xtc = repdir / "md.xtc"
    return repdir / "md.tpr", xtc


def _dump_group(tpr: Path, gmx: str, workdir: Path) -> str:
    """trjconv -dump group for the core biomolecule of this system
    (Backbone for proteins, DNA/RNA for nucleic systems, Protein/System as
    fallbacks) — a DNA-only tpr has no 'Protein' group, so the hardcoded
    prompt would fail every dump."""
    from .mdsim import _analysis_group, _index_groups
    try:
        return _analysis_group(
            _index_groups(tpr, gmx, Path(workdir) / "rep_groups.ndx"))
    except Exception:  # noqa: BLE001
        return "Protein"


def _ca_rmsd(pdb_a: Path | str, pdb_b: Path | str) -> float | None:
    """Cα Kabsch RMSD (nm) between two backbone PDBs; None when either has
    < 3 CA atoms.

    R9: PBC-aware — when `pdb_a` carries a CRYST1 box, the fit and the
    residual are done with per-atom minimum image (two MD cluster reps
    from the same system can sit in different box wraps; a plain Kabsch
    averages the integer box vectors into a muddled transform, and the
    dedup threshold in pool_representatives is then blind)."""
    pa = Path(pdb_a)
    a = [t["xyz"] for t in _atoms(pa) if t["name"] == "CA"]
    b = [t["xyz"] for t in _atoms(Path(pdb_b)) if t["name"] == "CA"]
    if len(a) < 3 or len(b) < 3 or len(a) != len(b):
        return None
    A = np.asarray(a)
    B = np.asarray(b)
    box = _cryst1_box(pa)
    R, t = _mi_kabsch(A, B, box)
    d = A @ R.T + t - B
    if box is not None:
        br = np.diag(box) @ R
        d = d - np.round(d @ np.linalg.inv(br)) @ br
    return float(np.sqrt(np.mean(np.sum(d ** 2, axis=1))))


def cluster_representatives(workdir: Path, *, rep: int = 1,
                            max_n: int = 3, min_pop: float = 0.05,
                            gmx: str | None = None) -> list[dict]:
    """Representative structures of the top GROMACS clusters.

    Preferred: full-atom (Protein group) extraction from the MD trajectory
    at each cluster's first frame via ``gmx trjconv -dump`` — required for
    side-chain flexible docking, since the ``gmx cluster`` representatives
    are backbone-only (N/CA/C). Fallback when no tpr/xtc (or gmx) is
    available: split the backbone-only ``clusters_r{rep}.pdb`` models.

    Keeps clusters with population >= ``min_pop``, top ``max_n`` by
    population. Returns [{cluster, population, pdb, full_atoms}].
    """
    from .mdsim import _parse_xvg
    anadir = Path(workdir) / "analysis"
    cl_pdb = anadir / f"clusters_r{rep}.pdb"
    idx_xvg = anadir / f"cluster_idx_r{rep}.xvg"
    if not cl_pdb.is_file() or not idx_xvg.is_file():
        return []
    t, y = _parse_xvg(idx_xvg)
    if not y:
        return []
    yr = np.round(np.asarray(y, dtype=float)).astype(int)
    pop: dict[int, float] = {}
    first_t: dict[int, float] = {}
    total = float(len(yr)) or 1.0
    for c in sorted(set(yr.tolist())):
        i = np.where(yr == c)[0]
        pop[int(c)] = float(len(i)) / total
        first_t[int(c)] = float(t[i[0]])
    selected = sorted(((c, p) for c, p in pop.items() if p >= min_pop),
                      key=lambda x: (-x[1], x[0]))[:max_n]
    outdir = anadir / "conformers"
    outdir.mkdir(parents=True, exist_ok=True)
    # full-atom path needs the replica's tpr/xtc and a gmx binary
    repdir = Path(workdir) / f"md_rep{rep}"
    tpr, xtc = _rep_files(repdir)
    full_ok = gmx is not None and tpr.is_file() and xtc.is_file()
    dump_group = _dump_group(tpr, gmx, repdir) if full_ok else "Protein"
    out: list[dict] = []
    for c, population in selected:
        pdb = outdir / f"reps_r{rep}_c{c}.pdb"
        full = False
        if full_ok and not pdb.is_file():
            # R8: mean over ALL cluster frames (MDAnalysis); fall back
            # to the legacy single-frame dump if MDAnalysis is missing
            i_frames = np.where(yr == c)[0]
            ok = _cluster_mean_pdb(tpr, xtc, i_frames, dump_group, pdb,
                                   gmx=gmx, workdir=repdir)
            if not ok:
                from ..utils import run_cmd
                try:
                    run_cmd([gmx, "trjconv", "-s", str(tpr),
                             "-f", str(xtc),
                             "-dump", f"{first_t[c]:.3f}", "-o", str(pdb)],
                            log_file=repdir / f"repdump_c{c}.log",
                            stdin=f"{dump_group}\n")
                except Exception:  # noqa: BLE001
                    pdb.unlink(missing_ok=True)
        if pdb.is_file():
            # backbone-only gmx cluster reps have exactly {N, CA, C};
            # full-atom extractions carry side chains (CB etc.)
            names = {l[12:16].strip() for l in pdb.read_text().splitlines()
                     if l.startswith("ATOM")}
            full = not names <= {"N", "CA", "C", "O"}
            out.append({"cluster": int(c), "population": round(population, 3),
                        "pdb": pdb, "full_atoms": full})
            continue
        # fallback: backbone-only model split
        _write_backbone_rep(cl_pdb, c, pdb)
        out.append({"cluster": int(c), "population": round(population, 3),
                    "pdb": pdb, "full_atoms": False})
    return out


def _cluster_mean_pdb(tpr: Path, xtc: Path, frame_idx: np.ndarray,
                      group: str, out: Path,
                      gmx: str | None = None, workdir: Path | None = None) -> bool:
    """R8: write the MEAN structure of `frame_idx` frames (all frames of
    one cluster) as a PDB — a cluster's first frame is an arbitrary
    boundary snapshot, the mean is its true centroid. MDAnalysis reads
    the GROMACS TPR as topology (residue identities survive) and the xtc
    is streamed frame-by-frame (no full-trajectory memory blowup).

    R8 pitfall: xtc coordinates are UNWRAPPED — the COM drifts across
    the run (nstcomm only removes COM velocity), so a multi-chain
    system's chains separate by tens of nm in unwrapped space (the
    1HVI dimer's monomers drift ~30 nm apart over 2 ns) and a raw
    frame-mean is a ghost structure that no Kabsch to the crystal can
    fix. Fit first: gmx trjconv -ur compact -center -fit rot+trans on
    the group (wrap into the box, COM-center, remove rigid drift),
    then mean the fitted trajectory. Without gmx, fall back to the raw
    mean (documented approximation). Returns True on success."""
    try:
        import MDAnalysis as mda
    except ImportError:
        return False
    sel = {"Protein": "protein", "DNA": "nucleic", "RNA": "nucleic",
           "System": "all"}.get(group, "protein")
    try:
        u = mda.Universe(str(tpr), str(xtc))
    except Exception:  # noqa: BLE001
        return False
    atoms = u.atoms.select_atoms(sel)
    if len(atoms) < 3:
        return False
    box = np.array(u.dimensions[:3], dtype=np.float64)
    if box.min() <= 0:
        box = np.array([100.0, 100.0, 100.0])  # no meaningful box
    first = None
    for k in frame_idx:
        try:
            u.trajectory[int(k)]
            first = int(k)
            break
        except (IndexError, ValueError):
            continue
    if first is None:
        return False
    # R8: per-atom PBC correction against the first cluster frame.
    # xtc coordinates are unwrapped: the COM (and individual chains of
    # a multi-chain system) drift tens of nm across a run, so a raw
    # frame average is a ghost. gmx rms/cluster handle this by
    # minimum-imaging every atom against the reference; we do the same
    # — valid while no atom diffuses more than half the box net
    # (true for bound biomolecules; waters would not qualify, but they
    # are not averaged here). No rigid transform needed: the per-atom
    # correction already places every frame in the reference frame.
    u.trajectory[first]
    ref = atoms.positions.copy()
    acc = np.zeros((len(atoms), 3), dtype=np.float64)
    n = 0
    for k in frame_idx:
        try:
            u.trajectory[int(k)]
        except (IndexError, ValueError):
            continue
        p = atoms.positions.copy()
        d = p - ref
        d -= box * np.round(d / box)
        acc += ref + d
        n += 1
    if n == 0:
        return False
    u.trajectory[first]
    # MA 2.x: the positions property returns a COPY (even for the full
    # universe on xtc universes) — item assignment is silently lost;
    # the setter persists
    new_pos = u.atoms.positions.copy()
    new_pos[atoms.indices] = acc / n
    u.atoms.positions = new_pos
    try:
        import warnings as _w
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            atoms.write(str(out))
    except Exception:  # noqa: BLE001
        out.unlink(missing_ok=True)
        return False
    return True


def pool_representatives(workdir: Path, *, max_n: int = 3,
                         min_pop: float = 0.05, gmx: str | None = None,
                         dedup_nm: float = 1.0) -> list[dict]:
    """R7: pool cluster representatives across ALL replicas.

    Candidates (per-replica top clusters, population >= min_pop) are sorted
    by population and selected greedily; a candidate whose Cα Kabsch RMSD
    to any already-selected representative is < `dedup_nm` is a redundant
    conformer (the same state seen in another replica) and is skipped.
    Returns [{rep, cluster, population, pdb, full_atoms}] — the ensemble
    for conformation-selection docking."""
    anadir = Path(workdir) / "analysis"
    avail = [i for i in range(1, 51)
             if (anadir / f"cluster_idx_r{i}.xvg").is_file()]
    cands: list[dict] = []
    for i in avail:
        for r in cluster_representatives(workdir, rep=i, max_n=max_n,
                                         min_pop=min_pop, gmx=gmx):
            r = dict(r)
            r["rep"] = i
            cands.append(r)
    cands.sort(key=lambda r: (-r["population"], r["rep"], r["cluster"]))
    chosen: list[dict] = []
    for c in cands:
        if len(chosen) >= max_n:
            break
        redundant = False
        for k in chosen:
            d = _ca_rmsd(c["pdb"], k["pdb"])
            if d is not None and d < dedup_nm:
                redundant = True
                break
        if not redundant:
            chosen.append(c)
    return chosen


def _write_backbone_rep(cl_pdb: Path, cluster: int, out: Path) -> None:
    """Write the `cluster`-th MODEL of a gmx clusters PDB (protein only)."""
    models: list[list[str]] = []
    cur: list[str] | None = None
    with open(cl_pdb) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                cur = []
                models.append(cur)
            elif line.startswith("ENDMDL"):
                cur = None
            elif cur is not None:
                cur.append(line)
    if not 1 <= cluster <= len(models):
        out.write_text("")
        return
    keep = [l for l in models[cluster - 1] if l.startswith("ATOM")
            and l[17:20].strip() not in SOLVENT_RES]
    out.write_text("".join(keep) + "END\n")


def decide_hits(scores: list[float], ref_score: float | None, n_wanted: int,
                brain: AgentBrain | None, context: str = "") -> dict:
    """Choose threshold + hit count.

    Rules: threshold = ref_score + 2.0 (if reference ligand available),
    else 5th percentile of the score distribution. Take at most n_wanted,
    at least 5 (if available). LLM can veto/adjust within sane bounds.
    """
    arr = np.array([s for s in scores if s is not None and not np.isnan(s)])
    if len(arr) == 0:
        return {"threshold": None, "n_hits": 0, "rationale": "no valid scores"}
    if ref_score is not None:
        thr = ref_score + 2.0
        basis = f"known ligand reference score {ref_score:.2f} + 2.0 kcal/mol"
    else:
        thr = float(np.percentile(arr, 5))
        basis = f"5th percentile of score distribution ({thr:.2f})"
    below = int((arr <= thr).sum())
    n_hits = max(min(n_wanted, below), min(5, len(arr)))
    decision = {"threshold": round(float(thr), 3), "n_hits": n_hits,
                "rationale": basis, "n_below_threshold": below}
    if brain is not None:
        dec = brain.decide(
            "screening",
            f"根据打分分布确定命中阈值与命中数量。当前建议 threshold={thr:.2f}, "
            f"n_hits={n_hits} (低于阈值 {below} 个)。请确认或给出新的 threshold/n_hits "
            f"(JSON: {{\"threshold\": num, \"n_hits\": int, \"rationale\": str}})",
            context=(f"scores: mean={arr.mean():.2f}, p5={np.percentile(arr,5):.2f}, "
                     f"p50={np.percentile(arr,50):.2f}, min={arr.min():.2f}, "
                     f"ref={ref_score}") + ("\n" + context if context else ""),
            expect="json",
        )
        ans = dec.answer if isinstance(dec.answer, dict) else {}
        try:
            new_thr = float(ans.get("threshold", thr))
            new_n = int(ans.get("n_hits", n_hits))
            # sanity bounds
            new_thr = float(np.clip(new_thr, arr.min(), np.percentile(arr, 50)))
            new_n = int(np.clip(new_n, 1, max(5, below)))
            decision.update({"threshold": round(new_thr, 3), "n_hits": new_n,
                             "rationale": dec.rationale or decision["rationale"],
                             "llm_adjusted": True})
        except (TypeError, ValueError):
            pass
    return decision


# --------------------------------------------------------------------------- #
# graph node
# --------------------------------------------------------------------------- #
def screen(state: dict) -> dict:
    workdir = Path(state["project_dir"]) / "02_screening"
    workdir.mkdir(parents=True, exist_ok=True)
    opts = state.get("options", {})
    prep = state["target_prep"]
    d = resolve_defaults(opts)
    brain = AgentBrain(project_dir=Path(state["project_dir"])) if not opts.get("no_llm") else None
    n_jobs = int(opts.get("n_jobs", 32))
    pocket = prep["pocket"]
    rec_pdbqt = prep["receptor_pdbqt"]

    # R10/G6: fallback-aware library resolution (flaky DTP mirror etc.)
    lib, lib_used = resolve_library(opts, d)

    # 1) standardize
    std = workdir / "std"
    std.mkdir(exist_ok=True)
    t0 = time.time()
    df = standardize_sdf(lib, std / "ok.sdf", n_jobs=min(n_jobs, 16))
    ok_df = df[df.status == "ok"].reset_index(drop=True)
    logger.info(f"standardized: {len(ok_df)}/{len(df)} ok in {time.time()-t0:.0f}s")

    # 2) features + filter
    fdf = compute_features(ok_df, n_jobs=n_jobs)
    fdf = physchem_filter(fdf)
    logger.info(f"physchem filter: {len(fdf)} remain")

    # 3) ML prefilter + adaptive N
    model_path = None
    from ..config import MODELS
    model_candidates = [Path(state["project_dir"]) / "data" / "prefilter_rf.pkl",
                        MODELS / "prefilter_rf.pkl"]
    for mc in model_candidates:
        if mc.exists():
            model_path = mc
            break
    fdf = ml_prefilter(fdf, model_path, n_keep=d.screen_max_ligands)

    # bench vina to estimate throughput, then cap N by time budget
    bench_prefix = std / "bench"
    bench_prefix.mkdir(exist_ok=True)
    sample = fdf.head(20)
    t0 = time.time()
    bench_args = []
    for k, (_, row) in enumerate(sample.iterrows()):
        lig = bench_prefix / f"b{k}.pdbqt"
        if write_ligand_pdbqt(row["smiles"], lig, seed=int(row["idx"])):
            bench_args.append((rec_pdbqt, str(lig), str(bench_prefix / f"b{k}"),
                               pocket, 4, 2))
    if bench_args:
        bench_res = pmap(dock_one, bench_args, n_jobs=min(8, len(bench_args)))
        bench_t = (time.time() - t0) / max(1, len(bench_args))
    else:
        bench_t = 2.0
    budget_s = d.screen_time_budget_h * 3600
    n_by_time = int(budget_s / max(bench_t, 0.05) * 0.9)
    n_dock = int(np.clip(n_by_time, d.screen_min_ligands, d.screen_max_ligands))
    fdf = fdf.head(n_dock)
    logger.info(f"docking {len(fdf)} ligands (bench {bench_t:.2f}s/ligand, "
                f"budget {d.screen_time_budget_h}h)")

    # 4) dock
    dockdir = workdir / "docks"
    dockdir.mkdir(exist_ok=True)
    cpu_per = max(1, os_cpu() // n_jobs)
    args = []
    for _, row in fdf.iterrows():
        lig = dockdir / f"{int(row['idx'])}.pdbqt"
        prefix = dockdir / f"{int(row['idx'])}"
        if write_ligand_pdbqt(row["smiles"], lig, seed=int(row["idx"])):
            args.append((rec_pdbqt, str(lig), str(prefix), pocket,
                         d.dock_exhaustiveness_fast, cpu_per))
    results = pmap(dock_one, args, n_jobs=n_jobs)
    rdf = pd.DataFrame(results)
    rdf = rdf[rdf.ok]
    rdf = rdf.merge(fdf[["idx", "smiles", "QED", "SA"]],
                    left_on="smiles", right_on="smiles", how="left") \
             if False else rdf
    # restore idx + real smiles (rdf["smiles"] carries the ligand idx string)
    idx_map = {str(int(row["idx"])): row for _, row in fdf.iterrows()}
    rdf["idx"] = rdf["smiles"].map(lambda s: idx_map.get(str(s), {}).get("idx", np.nan))
    rdf = rdf.dropna(subset=["idx"]).astype({"idx": int})
    rdf["smiles"] = rdf["idx"].map(lambda i: idx_map[str(i)]["smiles"])
    rdf = rdf.drop_duplicates(subset="smiles").sort_values("score")
    rdf.to_csv(workdir / "dock_results.csv", index=False)

    # 5) rescore top with GNINA (on the docked poses)
    top = rdf.head(d.dock_rescore_topn)
    gnina_rows = []
    for _, row in top.iterrows():
        pose = Path(row["top_pose_pdbqt"]) if row.get("top_pose_pdbqt") else None
        if pose is None or not pose.exists():
            continue
        g = gnina_rescore(rec_pdbqt, str(pose), pocket)
        gnina_rows.append({"idx": int(row["idx"]), **g})
    if gnina_rows:
        gdf = pd.DataFrame(gnina_rows)
        rdf = rdf.merge(gdf, on="idx", how="left")
        ok_g = rdf.gnina_score.notna()
        if ok_g.sum() > 0:
            rdf["final_score"] = np.where(
                ok_g,
                0.5 * rdf.score + 0.5 * rdf.gnina_score,
                rdf.score,
            )
            rdf = rdf.sort_values("final_score")
        else:
            rdf["final_score"] = rdf.score

    # 6) positive control + hit decision
    ref_score = None
    lig_pdbqt = prep.get("ligand_pdbqt")
    if lig_pdbqt and Path(lig_pdbqt).exists():
        prefix = std / "ref"
        r = dock_one((rec_pdbqt, lig_pdbqt, str(prefix), pocket,
                      d.dock_exhaustiveness_final, os_cpu() // 2))
        if r["ok"]:
            ref_score = r["score"]
            logger.info(f"reference ligand redocked: {ref_score}")

    hit = decide_hits(rdf["final_score"].tolist() if "final_score" in rdf else rdf["score"].tolist(),
                      ref_score, d.n_hits, brain)

    # 7) top hits with poses
    col = "final_score" if "final_score" in rdf.columns else "score"
    hits_df = rdf.head(hit["n_hits"]).copy()
    hits = []
    hits_dir = workdir / "hits"
    hits_dir.mkdir(exist_ok=True)
    for _, row in hits_df.iterrows():
        pose_src = dockdir / f"{int(row['idx'])}.pdbqt"
        pose_docked = row.get("top_pose_pdbqt")
        pose_path = (Path(pose_docked) if pose_docked and Path(pose_docked).is_file()
                     else pose_src)
        pose_out = hits_dir / f"hit_{len(hits)+1}_{int(row['idx'])}.pdb"
        if pose_path.is_file():
            with open(pose_path) as fi, open(pose_out, "w") as fo:
                for line in fi:
                    if line.startswith("ATOM"):
                        fo.write(line)
            pose_out.write_text(pose_out.read_text() + "END\n")
        hits.append({
            "rank": len(hits) + 1,
            "idx": int(row["idx"]),
            "smiles": row["smiles"],
            "vina_score": float(row["score"]),
            "gnina_score": float(row["gnina_score"]) if "gnina_score" in row and not np.isnan(row.get("gnina_score", np.nan)) else None,
            "final_score": float(row[col]),
            "pose_pdbqt": str(pose_path) if pose_path.is_file() else None,
            "pose_pdb": str(pose_out) if pose_out.is_file() else None,
            "qed": float(row.get("QED", np.nan)) if "QED" in row else None,
        })
    out = {
        "library": lib_used,
        "library_path_resolved": str(lib),
        "n_standardized": len(ok_df),
        "n_after_filter": int(len(fdf)),
        "n_docked": int(len(rdf)),
        "bench_seconds_per_ligand": round(float(bench_t), 3),
        "reference_ligand_score": ref_score,
        "hit_decision": hit,
        "hits": hits,
        "results_csv": str(workdir / "dock_results.csv"),
    }
    jsave(workdir / "screening.json", out)
    state_out = dict(state)
    state_out["screening"] = out
    return state_out


def os_cpu() -> int:
    import os
    return os.cpu_count() or 8
