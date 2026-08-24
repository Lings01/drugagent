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
from .target_prep import make_rigid_pdbqt, to_pdbqt

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
def vhh_cdr_fragments(pdb: Path, *, plddt_cutoff: float = 50.0, pad: int = 2,
                      min_res: int = 4, max_frags: int = 3) -> list[list[tuple[str, int]]]:
    """R11/G10-v2: CDR/loop-candidate fragments from a modeled VHH.

    Residues whose mean pLDDT (PDB B-factor) is below `plddt_cutoff` are
    loop/CDR-like; take maximal contiguous runs (resSeq consecutive within
    a chain), pad each side by `pad` residues, keep runs >= `min_res`.
    Returns up to `max_frags` fragments (largest first) as
    [(chain, resSeq), ...]. Returns [] when nothing qualifies or a single
    run covers most of the structure (caller falls back to full-VHH
    docking). Benchmark basis (scripts/bench_vhh_dock.py): docking cost
    scales ~O(n^1.9) with ligand atom count, so a 100-200-atom fragment
    docks in ~3-5 min vs ~80 min for the full 773-atom VHH."""
    per_res: dict[tuple[str, int], list[float]] = {}
    try:
        with open(pdb) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")) and line[17:20] != "WAT":
                    try:
                        res = (line[21:22].strip() or "A", int(line[22:26]))
                        per_res.setdefault(res, []).append(float(line[60:66]))
                    except (ValueError, IndexError):
                        continue
    except OSError:
        return []
    if not per_res:
        return []
    # group by chain, sorted by resSeq
    chains: dict[str, list[tuple[int, float]]] = {}
    for (ch, n), bs in per_res.items():
        chains.setdefault(ch, []).append((n, float(np.mean(bs))))
    frags: list[list[tuple[str, int]]] = []
    for ch, items in chains.items():
        items.sort(key=lambda x: x[0])
        nums = [n for n, _ in items]
        lo = {n: v for n, v in items}
        # maximal contiguous low-pLDDT runs
        run: list[int] = []
        def _flush():
            if len(run) >= min_res:
                a = max(min(nums), run[0] - pad)
                b = min(max(nums), run[-1] + pad)
                frags.append([(ch, k) for k in range(a, b + 1) if k in lo])
        for n in nums:
            if lo[n] < plddt_cutoff:
                run.append(n)
            else:
                _flush()
                run = []
        _flush()
    # drop runs that are basically the whole structure
    n_res = len(per_res)
    frags = [f for f in frags if len(f) < 0.6 * n_res]
    frags.sort(key=len, reverse=True)
    return frags[:max_frags]


def _write_fragment_pdb(pdb: Path, frag: list[tuple[str, int]], out: Path) -> Path:
    """Write a PDB containing only the fragment's residues (all atoms)."""
    keep = set(frag)
    with open(pdb) as fh:
        src = fh.read().splitlines()
    lines = [l for l in src
             if not l.startswith(("ATOM", "HETATM"))
             or (l[21:22].strip() or "A", int(l[22:26])) in keep]
    out.write_text("".join(l + "\n" for l in lines))
    return out


def _pdbqt_is_flex(pdbqt: Path) -> bool:
    """True if the ligand PDBQT has ACTIVE torsions (TORSDOF n>0).

    This vina build needs the ROOT/ENDROOT graph on every ligand, so a
    RIGID-body PDBQT is still graph-wrapped — with TORSDOF 0. Distinguish
    by the torsion count, not by the presence of the graph."""
    try:
        with open(pdbqt) as fh:
            for l in fh:
                if l.startswith("TORSDOF"):
                    parts = l.split()
                    return len(parts) > 1 and int(parts[1]) > 0
    except (OSError, ValueError):
        pass
    return False


def _dock_ligand(model_dir: Path, src_pdb: Path, idx: str, rec_pdbqt,
                 pocket, cpu: int, flex: bool) -> dict:
    """Convert src_pdb -> pdbqt (cached, flex/rigid consistent) and dock."""
    lig_pdbqt = model_dir / f"vhh_{idx}.pdbqt"
    if lig_pdbqt.is_file() and _pdbqt_is_flex(lig_pdbqt) != flex:
        lig_pdbqt.unlink()
    if not lig_pdbqt.is_file():
        to_pdbqt(src_pdb, lig_pdbqt)          # graph + element fixes
        if not flex:
            make_rigid_pdbqt(lig_pdbqt)       # TORSDOF 0 -> rigid body
    # full VHH is a huge "ligand"; scale exhaustiveness down or docking
    # takes hours per candidate (even rigid, >100 atoms is a big grid)
    n_atoms = sum(1 for l in open(lig_pdbqt) if l.startswith("ATOM"))
    exh = 8 if n_atoms < 100 else 1
    prefix = model_dir / f"vhh_{idx}_dock"
    return dock_one((rec_pdbqt, str(lig_pdbqt), str(prefix), pocket, exh, cpu))


def dock_vhh_candidates(ok: list[dict], model_dir: Path, rec_pdbqt,
                        pocket, *, n_jobs: int = 8,
                        flex: bool = False, cdr_only: bool = False) -> list[dict]:
    """R10/G7 (+R11/G10 +G10-v2): parallel VHH -> PDBQT -> Vina docking.

    Was a serial for-loop: with full-length VHHs (hundreds of rotatable
    bonds) a single vina run takes minutes, so 100 candidates meant hours.
    Now joblib-parallel; each worker gets os_cpu/n_jobs threads so total
    core usage stays ~constant. Idempotent: an existing <idx>.pdbqt is
    reused — but a stale file whose flex/rigid mode does not match the
    requested `flex` is reconverted.

    R11/G10: `flex=False` (default) docks each VHH as a RIGID body. A
    folded domain from a single ESMFold model carries one conformation, so
    the torsional search space buys little while dominating the runtime
    (this vina build also runs single-core on large ligands — `--cpu 64`
    measured 1 core, see ROUNDLOG R10).

    R11/G10-v2: `cdr_only=True` (fast default) docks the CDR/loop
    FRAGMENTS instead of the full VHH. Benchmark basis: cost ~O(n^1.9) in
    ligand atom count, and the full-VHH score is clash-dominated (773
    atoms in a 25.84 Å box) — fragments (100-200 atoms) dock in ~3-5 min
    and their scores are interpretable. The composite candidate score is
    the best (lowest) fragment score; per-fragment scores are kept in
    `fragment_scores`. Falls back to full-VHH docking when no fragment
    qualifies. Returns the candidate dicts updated with the docking
    result."""
    if not ok:
        return []
    n_jobs = max(1, min(n_jobs, len(ok), 16))
    cpu = max(1, os_cpu() // n_jobs)

    def _dock(r: dict) -> dict:
        pdb = model_dir / f"vhh_{r['idx']}.pdb"
        if cdr_only:
            frags = vhh_cdr_fragments(pdb)
            if frags:
                frag_scores = []
                for fi, frag in enumerate(frags):
                    tag = f"{r['idx']}_frag{fi}"
                    frag_pdb = model_dir / f"vhh_{tag}.pdb"
                    if not frag_pdb.is_file():
                        _write_fragment_pdb(pdb, frag, frag_pdb)
                    d = _dock_ligand(model_dir, frag_pdb, tag, rec_pdbqt,
                                     pocket, cpu, flex)
                    if d.get("ok") and d.get("score") == d.get("score"):
                        frag_scores.append(d["score"])
                best = min(frag_scores) if frag_scores else None
                return dict(r, ok=bool(frag_scores),
                            score=best if best is not None else np.nan,
                            rmsd_lb=np.nan, rmsd_ub=np.nan,
                            top_pose_pdbqt=None,
                            n_fragments=len(frags),
                            fragment_scores=frag_scores)
        # full-VHH path (also the cdr_only fallback)
        d = _dock_ligand(model_dir, pdb, str(r["idx"]), rec_pdbqt, pocket,
                         cpu, flex)
        return dict(r, **d)

    return pmap(_dock, ok, n_jobs=n_jobs)


def screen_vhh(state: dict, workdir: Path, *, n: int, n_jobs: int) -> dict:
    from ..modules.screening import dock_one
    prep = state["target_prep"]
    pocket = prep["pocket"]
    rec_pdbqt = prep["receptor_pdbqt"]
    d = resolve_defaults(state.get("options") or {})

    lib_path = LIBRARIES / "vhh_library.fasta"
    if lib_path.exists():
        seqs = load_fasta(lib_path)
    else:
        seqs = generate_vhh_library(int(d.vhh_lib_size))
        save_library(seqs, lib_path)
    seqs = seqs[:n]
    logger.info(f"VHH track A: screening {len(seqs)} sequences")

    model_dir = workdir / "vhh_models"
    model_dir.mkdir(exist_ok=True)
    results = pmap(model_vhh_one,
                   [(i, s, str(model_dir)) for i, s in enumerate(seqs)],
                   n_jobs=n_jobs)
    ok = [r for r in results if r["ok"]]
    # pLDDT filter — R11/G9: default comes from config (fast 35 / full 50),
    # overridable via options.vhh_plddt_min
    plddt_min = float(d.vhh_plddt_min)
    ok = [r for r in ok if r["plddt"] > plddt_min]
    ok.sort(key=lambda r: r["plddt"], reverse=True)
    ok = ok[: int(d.vhh_screen_n)]
    logger.info(f"VHH modeling: {len(results)} total, {len(ok)} pass "
                f"pLDDT>{plddt_min:g} and top-{int(d.vhh_screen_n)}")

    # dock (R10/G7: parallel; R11/G10: rigid bodies by default; R11/G10-v2:
    # CDR-fragment docking in fast mode — full-VHH dock is ~80-100 min and
    # clash-dominated; fragments ~3-5 min with interpretable scores)
    flex = bool(d.vhh_dock_flex)
    cdr_only = bool(d.vhh_dock_cdr_only)
    docked_all = dock_vhh_candidates(ok, model_dir, rec_pdbqt, pocket,
                                     n_jobs=n_jobs, flex=flex,
                                     cdr_only=cdr_only)
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
