"""Module D: nanobody (VHH) design & virtual screening (parallel tracks).

Track A: sequence library -> ESMFold modeling -> pLDDT filter -> Vina docking
Track B: RFdiffusion scaffold-guided de novo design on a VHH framework
Merge: composite affinity score (interface pLDDT + docking).
"""
from __future__ import annotations

import random
import re
from pathlib import Path

import numpy as np
from loguru import logger

from ..config import DEFAULTS, LIBRARIES, TOOLS, resolve_defaults
from ..llm import AgentBrain
from ..utils import jsave, pmap, run_cmd
from .screening import dock_one, os_cpu

# --------------------------------------------------------------------------- #
# synthetic VHH library
# --------------------------------------------------------------------------- #
# consensus human-like VHH framework (IMGT numbering), CDR gaps marked
FRAMEWORK = "EIYVASQGSSLVAPGQRFWMFWVRQAPGQNEKFRLYMYGPGQAFKYYGQWYWIGDTYNPSLRFSGSKSYNT"
# CDR1/CDR2 sampled from small germline-like repertoires
CDR1_POOL = [
    "GISGF", "GISGR", "GFTYD", "GITYD", "GFYTD", "GISGY", "GISTY",
]
CDR2_POOL = [
    "ISGNGGSY", "ISGSGGSL", "ISGNGGSN", "ISGKGGY", "ISGNGGSY",
]
CDR3_AA = "FYWGSTHILNDPCKAVMEQ"
# insertion points in the framework string (after these residues)
CDR1_POS = 12   # after "EIYVASQGSSLV"
CDR2_POS = 32   # after "APGQRFWMFWVRQAPGQNEK"


def generate_vhh_library(n: int, seed: int = 42) -> list[str]:
    """Generate a diverse synthetic VHH sequence library."""
    rng = random.Random(seed)
    seen: set[str] = set()
    out: list[str] = []
    fw = FRAMEWORK
    while len(out) < n:
        c1 = rng.choice(CDR1_POOL) + "".join(rng.choices("FYSDGN", k=rng.randint(2, 5)))
        c2 = rng.choice(CDR2_POOL)
        c3 = "".join(rng.choices(CDR3_AA, k=rng.randint(7, 17)))
        # avoid pathological CDR3 repeats
        if re.search(r"(.)\1{4,}", c3):
            continue
        seq = fw[:CDR1_POS] + c1 + fw[CDR1_POS:CDR2_POS] + c2 + fw[CDR2_POS:] + c3
        if seq in seen:
            continue
        seen.add(seq)
        out.append(seq)
    return out


def save_library(seqs: list[str], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for i, s in enumerate(seqs):
            fh.write(f">vhh_{i}\n{s}\n")
    return path


def load_fasta(path: Path) -> list[str]:
    seqs, cur = [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line.strip().upper())
    if cur:
        seqs.append("".join(cur))
    return seqs


# --------------------------------------------------------------------------- #
# ESMFold batch modeling
# --------------------------------------------------------------------------- #
def model_vhh_one(args: tuple) -> dict:
    """Model one VHH with ESMFold (separate process per item for RAM)."""
    idx, seq, outdir = args
    import os
    os.environ.setdefault("OMP_NUM_THREADS", "16")
    from .esmfold_run import predict, write_pdb
    res = {"idx": idx, "ok": False}
    pdbp = Path(outdir) / f"vhh_{idx}.pdb"
    # cache: reuse a previously modeled structure (pLDDT from B-factors)
    if pdbp.is_file():
        try:
            bf = []
            with open(pdbp) as fh:
                for line in fh:
                    if line.startswith("ATOM"):
                        bf.append(float(line[60:66]))
            if len(bf) > 10:
                import numpy as _np
                res.update(ok=True, plddt=float(_np.mean(bf)),
                           min_plddt=float(_np.min(bf)))
                return res
        except Exception:  # noqa: BLE001
            pass
    try:
        out = predict(seq, num_recycles=2, device="cpu")
        write_pdb(out, pdbp)
        res.update(ok=True, plddt=out["mean_plddt"], min_plddt=out["min_plddt"])
    except Exception as e:  # noqa: BLE001
        res["error"] = str(e)[:200]
    return res


# --------------------------------------------------------------------------- #
# track A: screening
# --------------------------------------------------------------------------- #
def dock_vhh_candidates(ok: list[dict], model_dir: Path, rec_pdbqt,
                        pocket, *, n_jobs: int = 8) -> list[dict]:
    """R10/G7: parallel VHH -> PDBQT -> Vina docking.

    Was a serial for-loop: with full-length VHHs (hundreds of rotatable
    bonds) a single vina run takes minutes, so 100 candidates meant hours.
    Now joblib-parallel; each worker gets os_cpu/n_jobs threads so total
    core usage stays ~constant. Idempotent: an existing <idx>.pdbqt is
    reused. Returns the candidate dicts updated with the docking result."""
    from ..modules.target_prep import to_pdbqt
    if not ok:
        return []
    n_jobs = max(1, min(n_jobs, len(ok), 16))
    cpu = max(1, os_cpu() // n_jobs)

    def _dock(r: dict) -> dict:
        pdb = model_dir / f"vhh_{r['idx']}.pdb"
        lig_pdbqt = model_dir / f"vhh_{r['idx']}.pdbqt"
        if not lig_pdbqt.is_file():
            to_pdbqt(pdb, lig_pdbqt)
        # full VHH is a huge "ligand" (100s of rotatable bonds); scale
        # exhaustiveness down or docking takes hours per candidate
        n_atoms = sum(1 for l in open(lig_pdbqt) if l.startswith("ATOM"))
        exh = 8 if n_atoms < 100 else 1
        prefix = model_dir / f"vhh_{r['idx']}_dock"
        d = dock_one((rec_pdbqt, str(lig_pdbqt), str(prefix), pocket, exh, cpu))
        return dict(r, **d)

    return pmap(_dock, ok, n_jobs=n_jobs)


def screen_vhh(state: dict, workdir: Path, *, n: int, n_jobs: int) -> dict:
    from ..modules.screening import dock_one
    prep = state["target_prep"]
    pocket = prep["pocket"]
    rec_pdbqt = prep["receptor_pdbqt"]

    lib_path = LIBRARIES / "vhh_library.fasta"
    if lib_path.exists():
        seqs = load_fasta(lib_path)
    else:
        seqs = generate_vhh_library(int(resolve_defaults(state.get("options") or {}).vhh_lib_size))
        save_library(seqs, lib_path)
    seqs = seqs[:n]
    logger.info(f"VHH track A: screening {len(seqs)} sequences")

    model_dir = workdir / "vhh_models"
    model_dir.mkdir(exist_ok=True)
    results = pmap(model_vhh_one,
                   [(i, s, str(model_dir)) for i, s in enumerate(seqs)],
                   n_jobs=n_jobs)
    ok = [r for r in results if r["ok"]]
    # pLDDT filter (relaxed in fast mode for synthetic libraries)
    plddt_min = float(state.get("options", {}).get(
        "vhh_plddt_min", 45.0 if state.get("options", {}).get("fast") else 70.0))
    ok = [r for r in ok if r["plddt"] > plddt_min]
    ok.sort(key=lambda r: r["plddt"], reverse=True)
    ok = ok[: int(state.get("options", {}).get("vhh_screen_n", 100))]
    logger.info(f"VHH modeling: {len(results)} total, {len(ok)} pass pLDDT>70 and top-n")

    # dock (R10/G7: parallel — the serial loop made Module D a 3-hour
    # bottleneck: a full VHH is a ~700-atom "ligand" and each vina call
    # takes minutes even at exhaustiveness 1)
    docked_all = dock_vhh_candidates(ok, model_dir, rec_pdbqt, pocket,
                                     n_jobs=n_jobs)
    docked = [r for r in docked_all if r.get("ok")]
    docked.sort(key=lambda r: r.get("score", 0))
    return {
        "track": "A_screening",
        "n_library": len(seqs),
        "n_modeled": len(results),
        "n_docked": len(docked),
        "results": docked[:50],
    }


# --------------------------------------------------------------------------- #
# track B: de novo design on VHH scaffold
# --------------------------------------------------------------------------- #
def design_vhh(state: dict, workdir: Path, *, n_designs: int) -> dict:
    from ..modules.binder import rf_repo, rf_python
    from ..modules.binder import pocket_hotspots

    prep = state["target_prep"]
    pocket = prep["pocket"]
    repo = rf_repo()
    outdir = workdir / "vhh_designs"
    outdir.mkdir(parents=True, exist_ok=True)

    scaffold_dir = TOOLS / "vhh_scaffolds"
    if not scaffold_dir.is_dir():
        logger.warning(f"vhh scaffold dir missing: {scaffold_dir}; track B skipped")
        return {"track": "B_de_novo", "n_designs": 0, "designs": []}

    # target pdb for scaffold-guided mode (chain A, no hetatm)
    target = workdir / "vhh_target.pdb"
    with open(prep["clean_pdb"]) as fh:
        txt = "".join(l for l in fh if not l.startswith("HETATM"))
    target.write_text(txt + "END\n")

    hotspots = pocket_hotspots(Path(prep["clean_pdb"]), pocket)
    cmd = [
        rf_python(), str(repo / "scripts" / "run_inference.py"),
        "scaffoldguided.scaffoldguided=True",
        f"scaffoldguided.scaffold_dir={scaffold_dir}",
        # use a .txt file (omegaconf ListConfig breaks negative slicing in RF)
        f"scaffoldguided.scaffold_list={scaffold_dir / 'scaffold_ids.txt'}",
        "scaffoldguided.target_pdb=True",
        f"scaffoldguided.target_path={target}",
        f"ppi.hotspot_res=[{','.join(hotspots)}]",
        f"inference.num_designs={n_designs}",
        f"inference.output_prefix={outdir}/vhh_design",
        "denoiser.noise_scale_ca=0",
        "denoiser.noise_scale_frame=0",
    ]
    # scaffoldguided wants Complex_Fold_base_ckpt.pt; fall back to the
    # plain complex model if the fold ckpt is not downloaded yet
    fold_ckpt = repo / "models" / "Complex_Fold_base_ckpt.pt"
    # treat a partial download (<400MB) as missing
    fold_ok = fold_ckpt.is_file() and fold_ckpt.stat().st_size > 400_000_000
    if not fold_ok:
        # need a ckpt trained with sec-struc/block-adjacency (d_t1d=28, d_t2d=47)
        fallback = repo / "models" / "InpaintSeq_Fold_ckpt.pt"
        if fallback.is_file():
            cmd.append(f"inference.ckpt_override_path={fallback}")
            logger.info("scaffoldguided: using InpaintSeq_Fold ckpt (fold ckpt missing)")
    run_cmd(cmd, cwd=repo, log_file=workdir / "vhh_design.log")
    designs = sorted(outdir.glob("vhh_design_*.pdb"))
    logger.info(f"VHH track B: {len(designs)} designs")

    # score each design with ESMFold complex
    from ..modules.binder import (_ca_sequence, _chain_ids, _design_chain,
                                  _extract_chain, _make_complex)
    from ..modules.esmfold_run import interface_metrics, predict

    scored = []
    for pdb in designs:
        name = pdb.stem
        entry = {"design": str(pdb)}
        try:
            # RF outputs may carry the target as an extra chain; score only
            # the designed (all-GLY) chain against the target
            binder_pdb = pdb
            if len(_chain_ids(pdb)) > 1:
                binder_pdb = outdir / f"{name}_binder.pdb"
                _extract_chain(pdb, _design_chain(pdb), binder_pdb)
            comp_pdb = outdir / f"{name}_complex.pdb"
            _make_complex(Path(prep["clean_pdb"]), binder_pdb, comp_pdb)
            out = predict([_ca_sequence(comp_pdb, "A"), _ca_sequence(comp_pdb, "B")],
                          num_recycles=3, device="cpu")
            comp_pdb.write_text(out["pdb"])
            rp = out.get("res_present")
            im = interface_metrics(comp_pdb,
                                   out["plddt"][rp] if rp is not None else out["plddt"])
            entry.update(
                complex_plddt=out["mean_plddt"],
                interface_plddt_mean=im.get("interface_plddt_mean"),
                interface_plddt_min=im.get("interface_plddt_min"),
                n_interface=im.get("n_interface_residues"),
            )
            # extract designed sequence (chain B)
            entry["sequence"] = _ca_sequence(comp_pdb, "B")
        except Exception as e:  # noqa: BLE001
            logger.warning(f"vhh design scoring failed {name}: {e}")
            entry["error"] = str(e)[:200]
        scored.append(entry)
    scored.sort(key=lambda x: (x.get("interface_plddt_mean") or 0), reverse=True)
    return {"track": "B_de_novo", "n_designs": len(designs), "designs": scored[:10]}


# --------------------------------------------------------------------------- #
# graph node
# --------------------------------------------------------------------------- #
def design_vhh_all(state: dict) -> dict:
    workdir = Path(state["project_dir"]) / "04_vhh"
    workdir.mkdir(parents=True, exist_ok=True)
    opts = state.get("options", {})
    d = resolve_defaults(opts)
    brain = AgentBrain(project_dir=Path(state["project_dir"])) if not opts.get("no_llm") else None
    n_jobs = int(opts.get("vhh_n_jobs", 4))

    track_a = screen_vhh(state, workdir, n=d.vhh_screen_n, n_jobs=n_jobs)
    track_b = design_vhh(state, workdir, n_designs=d.vhh_de_novo_designs)

    # merge with composite score
    candidates = []
    for r in track_a.get("results", []):
        candidates.append({
            "source": "screening",
            "plddt": r.get("plddt"),
            "docking_score": r.get("score"),
            "interface_plddt_mean": None,
            "idx": r.get("idx"),
        })
    for r in track_b.get("designs", []):
        candidates.append({
            "source": "de_novo",
            "plddt": r.get("complex_plddt"),
            "docking_score": None,
            "interface_plddt_mean": r.get("interface_plddt_mean"),
            "design": r.get("design"),
            "sequence": r.get("sequence"),
        })
    # normalize each metric to 0-1 and average
    def _norm(vals):
        arr = [v for v in vals if v is not None]
        if not arr:
            return {}
        lo, hi = min(arr), max(arr)
        span = (hi - lo) or 1.0
        return {i: (v - lo) / span for i, v in enumerate(vals) if v is not None}

    pl = _norm([c["plddt"] for c in candidates])
    dk = _norm([-(c["docking_score"] or 0) for c in candidates if c["docking_score"] is not None])
    im = _norm([c["interface_plddt_mean"] for c in candidates])
    for i, c in enumerate(candidates):
        parts = []
        if i in pl:
            parts.append(pl[i])
        if i in dk:
            parts.append(dk[i])
        if i in im:
            parts.append(im[i])
        c["composite_score"] = float(np.mean(parts)) if parts else 0.0
    candidates.sort(key=lambda c: c["composite_score"], reverse=True)

    out = {
        "track_a": track_a,
        "track_b": track_b,
        "ranked": candidates[:20],
        "best": candidates[0] if candidates else None,
    }
    if brain is not None and candidates:
        ctx = "\n".join(
            f"{i+1}. {c['source']} plddt={c['plddt']} dock={c['docking_score']} "
            f"ifpLDDT={c['interface_plddt_mean']} score={c['composite_score']}"
            for i, c in enumerate(candidates[:10]))
        dec = brain.decide("vhh", "从候选中选出最可能高亲和力的纳米抗体（编号1-10）",
                           context=ctx, choices=[str(i + 1) for i in range(min(10, len(candidates)))],
                           expect="choice")
        out["llm_pick"] = int(dec.answer) - 1
        out["llm_rationale"] = dec.rationale
    jsave(workdir / "vhh.json", out)
    state_out = dict(state)
    state_out["vhh"] = out
    return state_out
