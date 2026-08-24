"""GROMACS / MD tools: the agent owns force field + MDP + protocol."""
from __future__ import annotations

import os
import re
import time
from pathlib import Path

import numpy as np

from ..config import DEFAULTS, resolve_defaults, n_cores
from ..modules import mdsim as md
from ..utils import jsave, run_cmd
from .loop import Ctx, Tool, ToolError

SOLVENT_RES = {"HOH", "WAT", "SOL", "NA", "CL", "K", "MG", "ZN", "CA", "NA+", "CL-"}


def _env() -> dict:
    return md.gromacs()


def _is_ligand_pdb(pdb: Path) -> bool:
    het = set()
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("HETATM"):
                het.add(line[17:20].strip())
    return any(r.upper() not in SOLVENT_RES for r in het)


def _parse_mdrun_log(log: Path) -> dict:
    info = {"finished": False, "final_step": None, "final_pe": None,
            "performance_ns_per_day": None, "runtime": None}
    if not log.is_file():
        return info
    txt = log.read_text(errors="ignore")
    info["finished"] = "Finished mdrun" in txt
    m = re.search(r"Finished mdrun in (.+)", txt)
    if m:
        info["runtime"] = m.group(1).strip()
    m = re.search(r"Performance:\s+([\d.]+)\s+ns/day", txt)
    if m:
        info["performance_ns_per_day"] = float(m.group(1))
    steps = re.findall(r"^\s*(\d{4,})\s+(-?[\d.eE+]+)\s", txt, re.M)
    if steps:
        info["final_step"] = int(steps[-1][0])
        info["final_pe"] = float(steps[-1][1])
    return info


def _gro_atoms(gro: Path) -> int:
    if not gro.is_file():
        return 0
    # GRO dialects differ: GROMACS 2023 writes the title first, older builds
    # the atom count first. Accept either order.
    with open(gro) as fh:
        lines = fh.readlines(256)
    for l in lines[:2]:
        l = l.strip()
        if l.isdigit():
            return int(l)
    return 0


# --------------------------------------------------------------------------- #
def gmx_env(ctx: Ctx) -> dict:
    """GROMACS build + available force fields (the agent chooses the ff)."""
    env = _env()
    top = Path(env["gmxdata"]) / "top"
    ffs = sorted(p.name[:-3] for p in top.glob("*.ff"))
    return {"gmx": env["gmx"], "version": env["ver"],
            "detected_ff": env["ff"], "ff_ligand": env["ff_ligand"],
            "available_force_fields": ffs,
            "gmxdata": env["gmxdata"]}


def md_prepare(ctx: Ctx, complex_pdb: str, ff: str | None = None,
               salt: float | None = None, box_margin: float | None = None,
               water: str = "spce") -> dict:
    """Build the solvated ionized system: pdb2gmx (+ACPYPE ligand) -> box ->
    genion -> energy minimization.  Idempotent (reuses a finished EM)."""
    env = _env()
    ff = ff or env["ff"]
    salt = salt if salt is not None else \
        resolve_defaults(ctx.options).md_salt_m
    box_margin = box_margin if box_margin is not None else 1.0
    p = Path(complex_pdb)
    if not p.is_file():
        raise ToolError(f"not a file: {p}")
    workdir = ctx.stage_dir("05_md")
    is_lig = _is_ligand_pdb(p)
    sysinfo = md.build_system(p, workdir, env, is_ligand=is_lig,
                              salt=salt, box_margin=box_margin)
    build = Path(sysinfo["build_dir"])
    n_atoms = _gro_atoms(Path(sysinfo["gro"]))
    st = ctx.state()
    m = st.get("md") or {}
    m.update({
        "complex_pdb": str(p),
        "is_ligand": is_lig,
        "gromacs": {"binary": env["gmx"], "ff": ff, "version": env["ver"]},
        "build_dir": sysinfo["build_dir"],
    })
    st["md"] = m
    ctx.save_state()
    return {"ok": True, "ff": ff, "is_ligand": is_lig, "salt_m": salt,
            "box_margin_nm": box_margin, "n_atoms": n_atoms,
            "gro": sysinfo["gro"], "top": sysinfo["top"],
            "em_tpr": sysinfo["em_tpr"], "build_dir": sysinfo["build_dir"],
            "note": "EM 已完成; 生产 MD 的 mdp 用 mdp_template + write_file 自定, "
                    "grompp 时 -c 用 build/em.gro"}


def mdp_template(ctx: Ctx, name: str = "md", ns: float = 100.0) -> dict:
    """Default MDP template text (a starting point — edit as you see fit)."""
    if name == "em":
        text = md._mdp_em()
    elif name == "md":
        text = md._mdp_md(ns)
    else:
        raise ToolError("name must be 'em' or 'md'")
    return {"name": name, "md": name, "ns": ns if name == "md" else None,
            "text": text,
            "note": "默认模板: Verlet/PME/v-rescale/Berendsen; 你可用 write_file "
                    "写自己的 MDP (如 c-rescale、Lindman 盒子、不同盐浓度)"}


def grompp(ctx: Ctx, mdp: str, gro: str | None = None,
           workdir: str | None = None, maxwarn: int = 1) -> dict:
    """Validate a custom MDP and build a TPR.

    gro defaults: EM-style mdp -> ions.gro; production -> em.gro (start from
    the EM-relaxed coordinates; keep gen_vel=yes in the mdp).
    """
    env = _env()
    mdp_p = Path(mdp)
    if not mdp_p.is_absolute():
        mdp_p = ctx.project_dir / mdp_p
    if not mdp_p.is_file():
        raise ToolError(f"mdp not found: {mdp_p}")
    build = Path(workdir) if workdir else ctx.stage_dir("05_md") / "build"
    if not build.is_absolute():
        build = ctx.project_dir / build
    top = build / "solvated.top"
    if not top.is_file():
        raise ToolError(f"topology not found: {top} (先 md_prepare)")
    if gro:
        gro_p = Path(gro)
        if not gro_p.is_absolute():
            gro_p = ctx.project_dir / gro_p
    else:
        base = mdp_p.stem
        gro_p = build / ("em.gro" if "em" not in base else "ions.gro")
        if not gro_p.is_file():
            gro_p = build / ("em.gro" if "em" in base else "em.gro")
    if not gro_p.is_file():
        raise ToolError(f"gro not found: {gro_p} (传 gro 参数指定)")
    out_tpr = build / f"{mdp_p.stem}.tpr"
    rdir = build if mdp_p.parent == build else mdp_p.parent
    out_tpr = rdir / f"{mdp_p.stem}.tpr"
    run_cmd([env["gmx"], "grompp", "-f", str(mdp_p),
             "-c", str(gro_p), "-p", str(top),
             "-o", str(out_tpr), "-maxwarn", str(maxwarn)],
            log_file=rdir / f"grompp_{mdp_p.stem}.log")
    warnings = []
    glog = rdir / f"grompp_{mdp_p.stem}.log"
    if glog.is_file():
        for line in glog.read_text(errors="ignore").splitlines():
            if "WARNING" in line:
                warnings.append(line.strip())
    return {"ok": True, "tpr": str(out_tpr), "gro_used": str(gro_p),
            "n_warnings": len(warnings), "warnings": warnings[:20]}


def mdrun(ctx: Ctx, tpr: str, workdir: str,
          max_hours: float = 24.0, ntomp: int | None = None) -> dict:
    """Run one production/EM mdrun (idempotent; reuses a finished trajectory)."""
    env = _env()
    tpr_p = Path(tpr)
    if not tpr_p.is_absolute():
        tpr_p = ctx.project_dir / tpr_p
    if not tpr_p.is_file():
        raise ToolError(f"tpr not found: {tpr_p}")
    rdir = Path(workdir)
    if not rdir.is_absolute():
        rdir = ctx.project_dir / rdir
    rdir.mkdir(parents=True, exist_ok=True)
    xtc, mlog = rdir / "md.xtc", rdir / "md.log"
    if xtc.is_file() and mlog.is_file() and \
            "Finished mdrun" in mlog.read_text(errors="ignore"):
        info = _parse_mdrun_log(mlog)
        info.update({"ok": True, "reused": True, "dir": str(rdir)})
        return info
    ntomp = ntomp or max(4, n_cores() // 4)
    t0 = time.time()
    run_cmd([env["gmx"], "mdrun", "-deffnm", str(rdir / "md"),
             "-ntmpi", "1", "-ntomp", str(ntomp)],
            log_file=rdir / "mdrun_run.log",
            timeout=int(max_hours * 3600),
            env=dict(os.environ, OMP_NUM_THREADS=str(ntomp)))
    info = _parse_mdrun_log(mlog)
    if not info["finished"]:
        run_log = rdir / "mdrun_run.log"
        tail = run_log.read_text(errors="ignore")[-1500:] \
            if run_log.is_file() else ""
        raise ToolError(f"mdrun 未完成 (step={info['final_step']}); "
                        f"日志尾部:\n{tail}")
    info.update({"ok": True, "reused": False, "dir": str(rdir),
                 "wall_h": round((time.time() - t0) / 3600, 2)})
    return info


def gmx_analyze(ctx: Ctx, kind: str = "all", workdir: str = "05_md",
                replica: int = 1) -> dict:
    """Single-replica analysis: rmsd/rmsf/gyrate/ligand_rmsd/cluster."""
    env = _env()
    base = Path(workdir)
    if not base.is_absolute():
        base = ctx.project_dir / base
    rdir = base / f"md_rep{replica}"
    anadir = base / "analysis"
    anadir.mkdir(parents=True, exist_ok=True)
    gmx = env["gmx"]
    out: dict = {"rep": replica, "kind": kind}
    is_lig = bool((ctx.state().get("md") or {}).get("is_ligand"))

    def _xvg_sum(path: Path) -> dict:
        t, y = md._parse_xvg(path)
        if not y:
            return {"exists": False}
        return {"exists": True, "xvg": str(path),
                "final": round(float(y[-1]), 4),
                "mean": round(float(np.mean(y)), 4),
                "max": round(float(np.max(y)), 4),
                "n_frames": len(y)}

    if kind in ("rmsd", "all"):
        o = anadir / f"rmsd_r{replica}.xvg"
        if not o.is_file():
            run_cmd([gmx, "rms", "-s", str(rdir / "md.tpr"),
                     "-f", str(rdir / "md.xtc"), "-o", str(o),
                     "-fit", "rot+trans"],
                    log_file=rdir / "rms.log", stdin="Backbone\nBackbone\n")
        out["rmsd"] = _xvg_sum(o)
    if kind in ("gyrate", "all"):
        o = anadir / f"rg_r{replica}.xvg"
        if not o.is_file():
            run_cmd([gmx, "gyrate", "-s", str(rdir / "md.tpr"),
                     "-f", str(rdir / "md.xtc"), "-o", str(o)],
                    log_file=rdir / "gyrate.log", stdin="Protein\n")
        out["gyrate"] = _xvg_sum(o)
    if kind in ("rmsf", "all"):
        o = anadir / f"rmsf_r{replica}.xvg"
        if not o.is_file():
            run_cmd([gmx, "rmsf", "-s", str(rdir / "md.tpr"),
                     "-f", str(rdir / "md.xtc"), "-o", str(o), "-res"],
                    log_file=rdir / "rmsf.log", stdin="Backbone\n")
        t, y = md._parse_xvg(o)
        out["rmsf"] = {"exists": bool(y), "xvg": str(o),
                       "mean": round(float(np.mean(y)), 4) if y else None,
                       "max": round(float(np.max(y)), 4) if y else None,
                       "max_residue_index": int(np.argmax(y)) if y else None}
    if kind == "ligand_rmsd" or (kind == "all" and is_lig):
        o = anadir / f"lig_rmsd_r{replica}.xvg"
        if not o.is_file():
            run_cmd([gmx, "rms", "-s", str(rdir / "md.tpr"),
                     "-f", str(rdir / "md.xtc"), "-o", str(o),
                     "-fit", "rot+trans"],
                    log_file=rdir / "ligrms.log",
                    stdin="Protein\nLigand\n")
        out["ligand_rmsd"] = _xvg_sum(o)
    if kind in ("cluster", "all"):
        o = anadir / f"cluster_idx_r{replica}.xvg"
        if not o.is_file():
            run_cmd([gmx, "cluster", "-s", str(rdir / "md.tpr"),
                     "-f", str(rdir / "md.xtc"), "-method", "gromos",
                     "-cutoff", "1.5",
                     "-cl", str(anadir / f"clusters_r{replica}"),
                     "-clid", str(o)],
                    log_file=rdir / "cluster.log",
                    stdin="Backbone\nBackbone\n")
        t, y = md._parse_xvg(o)
        clusters = {}
        if y:
            uniq, counts = np.unique(np.round(y), return_counts=True)
            total = float(counts.sum()) or 1.0
            clusters = {int(k): round(float(c) / total, 3)
                        for k, c in zip(uniq, counts)}
        out["clusters"] = clusters
    if kind in ("chain", "all"):
        # per-chain RMSD (needs chain.ndx built by analyze_replicas; fit on
        # the largest chain, measure each other chain)
        chain_ndx = base / "chain.ndx"
        if chain_ndx.is_file():
            groups = {n: c for n, c in md._parse_ndx(chain_ndx)}
            # R10: DNA systems name their chains DNA_chainN (splitch),
            # not Protein_chainN — accept both so per-chain RMSD works
            # for nucleic-acid targets too
            chains = [{"name": n, "atoms": c} for n, c in groups.items()
                      if n.startswith(("Protein_chain", "DNA_chain"))
                      and c >= 50]
            if len(chains) >= 2:
                ref = max(chains, key=lambda c: c["atoms"])
                for c in chains:
                    if c["name"] == ref["name"]:
                        continue
                    if c["name"].startswith("DNA_chain"):
                        short = "dna_" + c["name"].replace(
                            "DNA_chain", "chain")
                    else:
                        short = c["name"].replace("Protein_chain", "chain")
                    o = anadir / f"rmsd_{short}_r{replica}.xvg"
                    if not o.is_file():
                        run_cmd([gmx, "rms", "-s", str(rdir / "md.tpr"),
                                 "-f", str(rdir / "md.xtc"),
                                 "-n", str(chain_ndx), "-o", str(o),
                                 "-fit", "rot+trans"],
                                log_file=rdir / f"rms_{short}.log",
                                stdin=f"{ref['name']}\n{c['name']}\n")
                    out[short] = _xvg_sum(o)
        else:
            out["chain"] = {"exists": False,
                            "note": "chain.ndx not built yet (run analyze_replicas / run_md first)"}
    if kind in ("ss", "all"):
        # DSSP-like secondary structure (MDAnalysis; optional dep)
        try:
            ss = md.analyze_ss(rdir / "md.tpr", rdir / "md.xtc")
            frac = ss["ss_frac"]
            out["ss"] = {"exists": True, "n_frames": ss["n_frames"],
                         "n_residues": ss["n_residues"],
                         "structured_initial": round(frac[0], 3),
                         "structured_final": round(frac[-1], 3),
                         "structured_min": round(min(frac), 3)}
        except Exception as e:  # noqa: BLE001
            out["ss"] = {"exists": False, "error": str(e)}
    out["units"] = ("gmx rms/gyrate/rmsf xvg are all in nm (report as Å = nm*10); "
                    "cluster cutoff 1.5 is in nm (GromOS standard); "
                    "ss structured fraction = H/G/I/E/B residues per frame")
    return out


def md_summary(ctx: Ctx, label: str = "", workdir: str = "05_md") -> dict:
    """Assemble state.md from the replica artifacts (report-compatible)."""
    env = _env()
    base = Path(workdir)
    if not base.is_absolute():
        base = ctx.project_dir / base
    replicas = []
    per_rep = []
    for i, rdir in enumerate(sorted(base.glob("md_rep*")), start=1):
        mlog = rdir / "md.log"
        info = _parse_mdrun_log(mlog)
        # R6: after auto-extension the last step lives in the extension
        # run logs (mdrun_ext{round}.log) — take the highest final_step
        ext_logs = sorted(rdir.glob("mdrun_ext*.log"),
                          key=lambda f: int("".join(c for c in f.stem
                                                    if c.isdigit()) or 0))
        for el in ext_logs:
            ei = _parse_mdrun_log(el)
            if ei["final_step"] and (not info["final_step"] or
                                     ei["final_step"] > info["final_step"]):
                info["final_step"] = ei["final_step"]
        rep = {"rep": i, "dir": str(rdir),
               "final_step": info.get("final_step")}
        anadir = base / "analysis"
        t, y = md._parse_xvg(anadir / f"rmsd_r{i}.xvg")
        if y:
            rep["rmsd_final"] = round(float(y[-1]), 4)
            per_rep.append({"rep": i, "rmsd": y,
                            "rg": md._parse_xvg(anadir / f"rg_r{i}.xvg")[1],
                            "rmsf_profile":
                                md._parse_xvg(anadir / f"rmsf_r{i}.xvg")[1]})
            # R1: per-chain RMSD artifacts (relative + self-fit)
            for f in sorted(anadir.glob(f"rmsd_chain*_r{i}.xvg")):
                key = f.stem[: -len(f"_r{i}")]
                y2 = md._parse_xvg(f)[1]
                if y2:
                    per_rep[-1][key] = y2
        t, y = md._parse_xvg(anadir / f"cluster_idx_r{i}.xvg")
        if y:
            uniq, counts = np.unique(np.round(y), return_counts=True)
            total = float(counts.sum()) or 1.0
            rep["clusters"] = {int(k): round(float(c) / total, 3)
                               for k, c in zip(uniq, counts)}
        replicas.append(rep)
    summary: dict = {}
    if per_rep:
        def _avg(key: str) -> dict:
            mats = [np.array(r[key]) for r in per_rep if r.get(key)]
            if not mats:
                return {}
            L = min(len(m) for m in mats)
            arr = np.stack([m[:L] for m in mats])
            return {"mean": arr.mean(axis=0).tolist(),
                    "std": arr.std(axis=0).tolist()}
        summary = {
            "rmsd": _avg("rmsd"), "rg": _avg("rg"),
            "rmsf_profile_mean": (np.mean([np.array(r["rmsf_profile"])
                                           for r in per_rep
                                           if r.get("rmsf_profile")], axis=0)
                                  .tolist()
                                  if any(r.get("rmsf_profile") for r in per_rep)
                                  else []),
            "final_rmsd_mean": float(np.mean([r["rmsd"][-1]
                                              for r in per_rep])),
            "final_rg_mean": float(np.mean([r["rg"][-1]
                                            for r in per_rep])),
            "replicas": per_rep,
        }
        # R1: aggregated cluster populations (int keys)
        cl_all = [{int(k): float(v) for k, v in
                   (r.get("clusters") or {}).items()}
                  for r in per_rep if r.get("clusters")]
        cl_keys = {k for c in cl_all for k in c}
        summary["clusters"] = {int(k): round(float(np.mean(
            [c.get(k, 0.0) for c in cl_all])), 3) for k in sorted(cl_keys)}
        # R1: per-chain RMSD (relative + self-fit), if xvg artifacts exist
        chain_keys = sorted({k for r in per_rep for k in r
                             if k.startswith("rmsd_chain")})
        for key in chain_keys:
            vals = [np.array(r[key]) for r in per_rep if r.get(key)]
            if vals:
                L = min(len(v) for v in vals)
                arr = np.stack([v[:L] for v in vals])
                summary[key] = {
                    "mean": arr.mean(axis=0).tolist(),
                    "final": round(float(np.mean([v[-1] for v in vals])), 4)}
        # R1: secondary-structure persistence (per replica, then average)
        ss_mats, ss_st_mats = [], []
        for r in per_rep:
            try:
                ss = md.analyze_ss(Path(r["dir"]) / "md.tpr",
                                   Path(r["dir"]) / "md.xtc")
                ss_mats.append(np.array(ss["ss_frac"]))
                ss_st_mats.append(np.array(ss["ss_stable"]))
            except Exception:  # noqa: BLE001
                pass
        if ss_mats:
            L = min(len(m) for m in ss_mats)
            arr = np.stack([m[:L] for m in ss_mats])
            summary["ss_frac_mean"] = arr.mean(axis=0).tolist()
            summary["initial_ss_mean"] = float(arr[:, 0].mean())
            summary["final_ss_mean"] = float(arr[:, -1].mean())
        if ss_st_mats:
            L = min(len(m) for m in ss_st_mats)
            summary["ss_stable_mean"] = (
                np.stack([m[:L] for m in ss_st_mats]).mean(axis=0).tolist())
        # R10/G3 (R1 收尾): region-level flexibility from the mean RMSF
        # profile, so the report/interpretation can cite flexible loops
        summary["flexible_regions"] = md.flexible_regions(
            summary["rmsf_profile_mean"])
        summary["interpretation"] = md.interpret_stability(summary)
    ns = None
    if replicas and replicas[0].get("final_step"):
        ns = round(replicas[0]["final_step"] * 0.002 / 1000.0, 3)
    st = ctx.state()
    m = st.get("md") or {}
    m.update({
        "system": {"label": label or m.get("system", {}).get("label", ""),
                   "type": "ligand" if m.get("is_ligand") else "protein"},
        "gromacs": m.get("gromacs") or {"binary": env["gmx"], "ff": env["ff"],
                                        "version": env["ver"]},
        "ns": ns, "reps": len(replicas),
        "replicas": replicas, "summary": summary,
    })
    st["md"] = m
    jsave(base / "md.json", m)
    ctx.save_state(md=m)
    return {"ns": ns, "reps": len(replicas),
            "final_rmsd_mean": summary.get("final_rmsd_mean"),
            "final_rg_mean": summary.get("final_rg_mean"),
            "clusters": {r["rep"]: r.get("clusters") for r in replicas},
            "md_json": str(base / "md.json")}


def run_md(ctx: Ctx) -> dict:
    """Full deterministic MD stage (1.0 pipeline as one tool)."""
    st = ctx.state()
    state_out = md.run_md(st)           # returns full state
    out = state_out["md"]
    ctx.save_state(md=out)
    return {"system": out.get("system"), "ns": out.get("ns"),
            "reps": out.get("reps"),
            "final_rmsd_mean": (out.get("summary") or {}).get("final_rmsd_mean"),
            "final_rg_mean": (out.get("summary") or {}).get("final_rg_mean")}


# --------------------------------------------------------------------------- #
def build() -> list[Tool]:
    return [
        Tool("gmx_env",
             "GROMACS 信息: 版本 + 可用力场列表 (你选力场)。",
             {"type": "object", "properties": {}, "required": []}, gmx_env),
        Tool("md_prepare",
             "建体系: pdb2gmx(+配体 ACPYPE)→盒子→加离子→能量最小化。"
             "ff/salt/box_margin 由你指定 (缺省用构建自带力场/0.15M/1.0nm)。幂等。",
             {"type": "object",
              "properties": {"complex_pdb": {"type": "string"},
                             "ff": {"type": "string"},
                             "salt": {"type": "number"},
                             "box_margin": {"type": "number"},
                             "water": {"type": "string"}},
              "required": ["complex_pdb"]}, md_prepare, long_running=True),
        Tool("mdp_template",
             "默认 MDP 模板文本 (em/md; 作为起点, 建议按需修改后 write_file)。",
             {"type": "object",
              "properties": {"name": {"type": "string",
                                      "enum": ["em", "md"]},
                             "ns": {"type": "number"}},
              "required": ["name"]}, mdp_template),
        Tool("grompp",
             "用你的 MDP 校验 + 生成 TPR (返回 WARNING 列表)。生产 mdp 默认 "
             "-c build/em.gro (需保留 gen_vel=yes)。",
             {"type": "object",
              "properties": {"mdp": {"type": "string"},
                             "gro": {"type": "string"},
                             "workdir": {"type": "string"},
                             "maxwarn": {"type": "integer"}},
              "required": ["mdp"]}, grompp),
        Tool("mdrun",
             "跑一次 mdrun (幂等; 完成后返回 final_step/PE/性能)。长任务, "
             "max_hours 设大些。",
             {"type": "object",
              "properties": {"tpr": {"type": "string"},
                             "workdir": {"type": "string",
                                         "description": "如 05_md/md_rep1"},
                             "max_hours": {"type": "number"},
                             "ntomp": {"type": "integer"}},
              "required": ["tpr", "workdir"]}, mdrun, long_running=True),
        Tool("gmx_analyze",
             "单副本分析: rmsd/rmsf/gyrate/ligand_rmsd/cluster/chain(分链RMSD)/ss(二级结构) (幂等, 复用已有 xvg)。",
             {"type": "object",
              "properties": {"kind": {"type": "string",
                                      "enum": ["rmsd", "rmsf", "gyrate",
                                               "ligand_rmsd", "cluster",
                                               "chain", "ss", "all"]},
                             "workdir": {"type": "string"},
                             "replica": {"type": "integer"}},
              "required": []}, gmx_analyze),
        Tool("md_summary",
             "汇总所有副本 → state.md (报告格式; label 描述模拟的体系)。",
             {"type": "object",
              "properties": {"label": {"type": "string"},
                             "workdir": {"type": "string"}},
              "required": []}, md_summary),
        Tool("run_md",
             "整段 MD (确定性标准流程: 选体系→建体系→固定模板→3 副本→分析)。",
             {"type": "object", "properties": {}, "required": []},
         run_md, long_running=True),
    ]
