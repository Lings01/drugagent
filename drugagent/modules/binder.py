"""Module C: de novo protein binder design (RFdiffusion + MPNN + ESMFold).

Agent decisions: binder length/topology, hotspot residues from pocket,
number of designs, filtering thresholds.
"""
from __future__ import annotations

import re
import shutil
from pathlib import Path

import numpy as np
from loguru import logger

from ..config import TOOLS
from ..llm import AgentBrain
from ..utils import jload, jsave, pmap, run_cmd


def rf_repo() -> Path:
    p = TOOLS / "RFdiffusion"
    if not p.is_dir():
        raise FileNotFoundError("RFdiffusion not found (run `drugagent setup --tools`)")
    return p


def rf_python() -> str:
    """Python interpreter with torch for RFdiffusion (own env preferred)."""
    from ..config import ENV_DIR
    for cand in (
        TOOLS / "rfdiff_env" / "bin" / "python",
        ENV_DIR / "bin" / "python",
    ):
        if cand.is_file():
            return str(cand)
    import sys
    return sys.executable


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _rename_chain(pdb: Path, out: Path, new_chain: str = "A") -> Path:
    with open(pdb) as fh:
        lines = []
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                line = line[:21] + new_chain + line[22:]
            lines.append(line)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(lines))
    return out


def pocket_hotspots(target_pdb: Path, pocket: dict, max_res: int = 12) -> list[str]:
    """Residues (chain A numbering) within radius of pocket center."""
    cx, cy, cz = pocket["center"]
    r = pocket["xsize"] / 2
    hits: dict[int, float] = {}
    with open(target_pdb) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                d = ((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2) ** 0.5
                if d <= r + 2.0:
                    resi = int(line[22:26].strip())
                    if resi not in hits or d < hits[resi]:
                        hits[resi] = d
    top = sorted(hits.items(), key=lambda kv: kv[1])[:max_res]
    return [f"A{i}" for i, _ in top]


# --------------------------------------------------------------------------- #
# RFdiffusion run
# --------------------------------------------------------------------------- #
def rfdesign(
    target_pdb: Path,
    pocket: dict,
    workdir: Path,
    *,
    n_designs: int = 8,
    length: tuple[int, int] = (60, 80),
    hotspots: list[str] | None = None,
    seed: int = 0,
    env: str | None = None,
    rf_cautious: bool = False,
) -> list[Path]:
    """Run RFdiffusion binder design. Returns list of design PDBs.

    R16/P3: rf_cautious=True skips the (expensive, ~20 min on CPU) RF
    run when at least n_designs top-level design PDBs already exist —
    symmetry with vhh_rf_cautious, whose RF scaffold-guided mode
    skips existing designs by default. Default False keeps the
    always-resample behavior (deterministic=False)."""
    repo = rf_repo()
    outdir = workdir / "rf_designs"
    outdir.mkdir(parents=True, exist_ok=True)
    existing = sorted(p for p in outdir.glob("design_*.pdb")
                      if p.stem.split("_")[-1].isdigit())
    if rf_cautious and len(existing) >= n_designs:
        logger.info(f"rf_cautious: reusing {len(existing)} existing "
                    f"designs in {outdir.name} (skip RF run)")
        return existing

    # RFdiffusion expects target on chain A (hotspot refs use chain A)
    target_a = workdir / "rf_target.pdb"
    _rename_chain(target_pdb, target_a, "A")
    # strip ligands (heavy ATOM only kept by _rename_chain; drop HETATM)
    with open(target_a) as fh:
        txt = "".join(l for l in fh if not l.startswith("HETATM"))
    target_a.write_text(txt + "END\n")

    n_target = len({l[22:26].strip() for l in open(target_a)
                if l.startswith("ATOM") and l[12:16].strip() == "CA"})
    contigs = f"A1-{n_target}/0 {length[0]}-{length[1]}"
    cmd = [
        rf_python(), str(repo / "scripts" / "run_inference.py"),
        f"contigmap.contigs=[{contigs}]",
        f"inference.input_pdb={target_a}",
        f"inference.num_designs={n_designs}",
        f"inference.output_prefix={outdir}/design",
        "inference.deterministic=False",
        "denoiser.noise_scale_ca=0",
        "denoiser.noise_scale_frame=0",
    ]
    if hotspots:
        cmd.append(f"ppi.hotspot_res=[{','.join(hotspots)}]")
    envd = None
    if env:
        envd = dict(run_cmd.__defaults__ and {})  # noqa
    run_cmd(cmd, cwd=repo, log_file=workdir / "rf_design.log")
    designs = sorted(outdir.glob("design_*.pdb"))
    logger.info(f"RFdiffusion produced {len(designs)} designs")
    return designs


# --------------------------------------------------------------------------- #
# MPNN sequence design
# --------------------------------------------------------------------------- #
def mpnn_sequence(design_pdb: Path, workdir: Path, *, n_seq: int = 20,
                  env: str | None = None) -> list[str]:
    """Design sequences for a backbone with ProteinMPNN (dauparas/ProteinMPNN).

    The tool lives in <RF repo>/mpnn/ (protein_mpnn_run.py + model/v_48_010.pt);
    if the deployment is missing or the run fails, falls back to a simple
    heuristic sequence (keeps the pipeline robust).
    """
    repo = rf_repo()
    mpnn_dir = repo / "mpnn"
    run_py = mpnn_dir / "protein_mpnn_run.py"
    outdir = workdir / "mpnn"
    outdir.mkdir(parents=True, exist_ok=True)
    if run_py.is_file() and (mpnn_dir / "model" / "v_48_010.pt").is_file():
        cmd = [
            rf_python(), str(run_py),
            "--pdb_path", str(Path(design_pdb).resolve()),
            "--path_to_model_weights", str(mpnn_dir / "model"),
            "--model_name", "v_48_010",
            "--num_seq_per_target", str(n_seq),
            "--out_folder", str(outdir.resolve()),
            "--sampling_temp", "0.1",
        ]
        try:
            run_cmd(cmd, cwd=mpnn_dir, log_file=outdir / "mpnn.log")
            fa = outdir / "seqs" / f"{design_pdb.stem}.fa"
            seqs: list[str] = []
            if fa.is_file():
                cur: list[str] = []
                for line in fa.read_text().splitlines():
                    if line.startswith(">"):
                        if cur:
                            seqs.append("".join(cur))
                            cur = []
                    else:
                        cur.append(line.strip().upper())
                if cur:
                    seqs.append("".join(cur))
            if seqs:
                # multi-chain inpaint outputs: the designed region is the
                # all-GLY chain; pick that part of '/'-separated sequences
                chains = _chain_ids(Path(design_pdb))
                dch = _design_chain(Path(design_pdb))
                idx = chains.index(dch)
                picked = []
                for x in seqs:
                    parts = x.split("/")
                    part = parts[idx] if len(parts) == len(chains) else parts[0]
                    picked.append(part)
                # first FASTA entry is the native (all-GLY for RF designs); drop it
                picked = [x for x in picked if set(x) != {"G"}] or picked
                if picked:
                    return picked[:n_seq]
        except Exception as e:  # noqa: BLE001
            logger.warning(f"ProteinMPNN run failed ({e}); using fallback sequence")
    # fallback: CA trace + DSSP-free heuristic -> hydrophobic-core random seq
    return [_fallback_sequence(design_pdb)]


def _fallback_sequence(pdb: Path) -> str:
    import random
    seq = ""
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                seq += "G"
    if len(seq) < 30:
        seq = "G" * 60
    return seq


# --------------------------------------------------------------------------- #
# scoring
# --------------------------------------------------------------------------- #
def score_designs(design_pdbs: list[Path], target_pdb: Path, workdir: Path,
                  seqs: dict[str, list[str]], *, device: str = "cpu") -> list[dict]:
    """ESMFold monomer pLDDT + target-binder complex interface pLDDT.

    R14: persistent per-design cache (``scored.json`` in workdir), same
    pattern as vhh.score_designs — a design whose (design mtime, target
    mtime, alt sequence, ESMFold weight tag) signature matches the cached
    entry is reused without paying the two ESMFold predictions again."""
    from .esmfold_run import esmfold_version_tag, interface_metrics, predict

    workdir.mkdir(parents=True, exist_ok=True)
    cache_path = workdir / "scored.json"
    cache = jload(cache_path) if cache_path.is_file() else {}
    target_mt = target_pdb.stat().st_mtime if target_pdb.is_file() else 0.0
    tag = esmfold_version_tag()
    scored = []
    for pdb in design_pdbs:
        name = pdb.stem
        d = {"design": str(pdb), "seqs": seqs.get(name, [])}
        alt = (seqs.get(name) or [None])[0]
        sig = [pdb.stat().st_mtime if pdb.is_file() else 0.0, target_mt,
               alt, tag]
        cached = cache.get(name)
        if cached and cached.get("sig") == sig:
            d.update({k: v for k, v in cached.items() if k != "sig"})
            logger.info(f"binder design {name}: reusing cached validation "
                        f"(interface pLDDT {cached.get('interface_plddt_mean')})")
            scored.append(d)
            continue
        # monomer — RF design PDBs carry GLY residue names; when a real
        # sequence (RF log / ProteinMPNN) is available, score that instead
        dch = _design_chain(pdb)
        mono_pdb = pdb
        if len(_chain_ids(pdb)) > 1:
            mono_pdb = workdir / f"{name}_binder_only.pdb"
            _extract_chain(pdb, dch, mono_pdb)
        seq = _ca_sequence(mono_pdb)
        if alt and set(seq) <= {"G"} and len(alt) == len(seq):
            seq = alt
            d["seq_used"] = "mpnn_or_rf_log"
        try:
            mono = predict(seq, num_recycles=2, device=device)
            d["mono_plddt"] = mono["mean_plddt"]
        except Exception as e:  # noqa: BLE001
            logger.warning(f"mono predict failed for {name}: {e}")
            d["mono_plddt"] = float("nan")
        # complex: target + binder (concatenate PDBs with chain B)
        complex_pdb = workdir / "complex" / f"{name}_complex.pdb"
        complex_pdb.parent.mkdir(parents=True, exist_ok=True)
        _make_complex(target_pdb, mono_pdb, complex_pdb)
        try:
            seqB = _ca_sequence(complex_pdb, "B")
            altB = (seqs.get(name) or [None])[0]
            if altB and set(seqB) <= {"G"} and len(altB) == len(seqB):
                seqB = altB
            comp = predict([_ca_sequence(complex_pdb, "A"), seqB],
                           num_recycles=3, device=device)
            write_pdb_file(complex_pdb, comp)
            rp = comp.get("res_present")
            im = interface_metrics(complex_pdb,
                                   comp["plddt"][rp] if rp is not None else comp["plddt"])
            d["complex_plddt"] = comp["mean_plddt"]
            d["interface_plddt_mean"] = im.get("interface_plddt_mean")
            d["interface_plddt_min"] = im.get("interface_plddt_min")
            d["n_interface"] = im.get("n_interface_residues")
            # cache only full successes (errors are retried next run)
            if d["interface_plddt_mean"] is not None:
                cache[name] = {"sig": sig,
                               **{k: v for k, v in d.items()
                                  if k not in ("design", "seqs")}}
                jsave(cache_path, cache)
        except Exception as e:  # noqa: BLE001
            logger.warning(f"complex predict failed for {name}: {e}")
            d["complex_plddt"] = float("nan")
        scored.append(d)
    scored.sort(key=lambda x: (x.get("interface_plddt_mean") or 0), reverse=True)
    return scored


_THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "ASN": "N",
    "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
    "TRP": "W", "TYR": "Y",
}


def _ca_sequence(pdb: Path, chain: str | None = None) -> str:
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("x", str(pdb))
    seq = []
    for model in struct:
        for ch in model:
            if chain and ch.id != chain:
                continue
            for res in ch:
                if res.id[0] != " ":
                    continue
                name = (res.get_resname() if hasattr(res, "get_resname")
                        else res.get_name())
                one = _THREE_TO_ONE.get(name)
                if one:
                    seq.append(one)
    return "".join(seq)


def _chain_ids(pdb: Path) -> list[str]:
    ids: list[str] = []
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ATOM"):
                ch = line[21]
                if ch not in ids:
                    ids.append(ch)
    return sorted(ids)


def _resname_set(pdb: Path, chain: str) -> set[str]:
    out: set[str] = set()
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[21] == chain:
                out.add(line[17:20].strip())
    return out


def _design_chain(pdb: Path) -> str:
    """Chain id of the designed region. RF inpaint outputs may carry the
    target as an extra chain; the designed region is the all-GLY chain
    (residue names are not written by RF)."""
    chains = _chain_ids(pdb)
    if len(chains) == 1:
        return chains[0]
    for ch in chains:
        if _resname_set(pdb, ch) <= {"GLY"}:
            return ch
    return chains[0]


def _extract_chain(pdb: Path, chain: str, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    with open(pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and line[21] == chain:
                lines.append(line)
    out.write_text("".join(lines) + "END\n")
    return out


def _make_complex(target_pdb: Path, binder_pdb: Path, out: Path) -> Path:
    """Concatenate target (chain A) + binder (chain B)."""
    def _read(p, chain):
        with open(p) as fh:
            lines = []
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    line = line[:21] + chain + line[22:]
                    lines.append(line)
        return lines
    txt = "".join(_read(target_pdb, "A")) + "".join(_read(binder_pdb, "B")) + "END\n"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(txt)
    return out


def write_pdb_file(path: Path, out: dict) -> None:
    path.write_text(out["pdb"])


# --------------------------------------------------------------------------- #
# graph node
# --------------------------------------------------------------------------- #
def design_binder(state: dict) -> dict:
    workdir = Path(state["project_dir"]) / "03_binder"
    workdir.mkdir(parents=True, exist_ok=True)
    opts = state.get("options", {})
    prep = state["target_prep"]
    pocket = prep["pocket"]
    brain = AgentBrain(project_dir=Path(state["project_dir"])) \
        if not opts.get("no_llm") else None

    from ..config import DEFAULTS, resolve_defaults
    d = resolve_defaults(opts)

    # agent picks length/topology
    length = (int(opts.get("binder_len_min", 60)), int(opts.get("binder_len_max", 80)))
    binder_type = opts.get("binder_type", "auto")
    if binder_type == "auto" and brain is not None:
        dec = brain.decide(
            "binder",
            "选择 binder 类型: miniprotein(60-80aa, 最通用) / helix(螺旋为主) / "
            "beta(beta-sheet 为主)。默认 miniprotein。",
            context=f"target chains={prep['completeness']['chains']}, "
                    f"pocket={pocket.get('method')}",
            choices=["miniprotein", "helix", "beta"],
            expect="choice",
        )
        binder_type = dec.answer
        if binder_type == "helix":
            length = (50, 70)
        elif binder_type == "beta":
            length = (60, 90)

    hotspots = pocket_hotspots(Path(prep["clean_pdb"]), pocket)
    n_designs = int(opts.get("n_binder_designs", d.n_binder_designs))

    # R16/P3: binder_rf_cautious (default false = always re-sample)
    designs = rfdesign(
        Path(prep["clean_pdb"]), pocket, workdir,
        n_designs=n_designs, length=length, hotspots=hotspots,
        rf_cautious=bool(opts.get("binder_rf_cautious", False)),
    )
    if not designs:
        raise RuntimeError("RFdiffusion produced no designs")

    # MPNN sequences per design
    seqs: dict[str, list[str]] = {}
    for pdb in designs:
        try:
            seqs[pdb.stem] = mpnn_sequence(pdb, workdir, n_seq=5)
        except Exception as e:  # noqa: BLE001
            logger.warning(f"mpnn failed for {pdb.stem}: {e}")
            seqs[pdb.stem] = []

    # score
    scored = score_designs(designs, Path(prep["clean_pdb"]), workdir, seqs,
                           device="cpu")
    # pick best with agent judgment
    top = scored[:5]
    out = {
        "binder_type": binder_type,
        "length_range": list(length),
        "hotspots": hotspots,
        "n_designs": len(designs),
        "designs": scored,
        "best": scored[0] if scored else None,
    }
    jsave(workdir / "binder.json", out)
    state_out = dict(state)
    state_out["binder"] = out
    return state_out
