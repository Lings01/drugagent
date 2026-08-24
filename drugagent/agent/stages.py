"""R11/G8: stage-level idempotency.

The ``run_*`` tools each execute a whole deterministic stage.  Without any
stage-level reuse, restarting a crashed e2e re-runs already finished stages
from scratch (round 10 lost ~35 min re-docking the screening stage).  These
helpers let a stage tool skip its work when the stage's ``state.json``
section already carries its completion marker, and expose ``force=True`` to
bypass.

Marker semantics (intentionally strict — a missing marker means "re-run",
which is always safe because the stage functions themselves reuse their
low-level artifacts: PDBQT caches, mdrun trajectories, xvg files, ...):

================  ==========================================================
stage             marker in state.json
================  ==========================================================
target_prep       receptor_pdbqt + clean_pdb + pocket all set
screening         hit_decision (written by decide_hits / sc.screen)
binder            designs non-empty and best is not None
vhh               track_a and track_b both present
md                summary and replicas present (run_md or an equivalent
                  md_prepare + mdrun + md_summary agent path finished)
report            report.html exists on disk
================  ==========================================================
"""
from __future__ import annotations

from pathlib import Path

from loguru import logger

from .loop import Ctx


# --------------------------------------------------------------------------- #
# completion checks
# --------------------------------------------------------------------------- #
def _check_target_prep(state: dict, project_dir: Path | None) -> tuple[bool, str]:
    tp = state.get("target_prep") or {}
    missing = [k for k in ("receptor_pdbqt", "clean_pdb", "pocket") if not tp.get(k)]
    if missing:
        return False, f"target_prep 缺 {missing}"
    return True, "ok"


def _check_screening(state: dict, project_dir: Path | None) -> tuple[bool, str]:
    scr = state.get("screening") or {}
    if not scr.get("hit_decision"):
        return False, "screening 无 hit_decision (未到命中判定)"
    return True, "ok"


def _check_binder(state: dict, project_dir: Path | None) -> tuple[bool, str]:
    b = state.get("binder") or {}
    if not b.get("designs"):
        return False, "binder.designs 为空"
    if b.get("best") is None:
        return False, "binder.best 为 None (打分未完成)"
    return True, "ok"


def _check_vhh(state: dict, project_dir: Path | None) -> tuple[bool, str]:
    v = state.get("vhh") or {}
    if "track_a" not in v or "track_b" not in v:
        return False, "vhh 缺 track_a/track_b (整段未跑完)"
    return True, "ok"


def _check_md(state: dict, project_dir: Path | None) -> tuple[bool, str]:
    m = state.get("md") or {}
    if not m.get("summary") or not m.get("replicas"):
        return False, "md 无 summary/replicas (run_md 未完成)"
    return True, "ok"


def _check_report(state: dict, project_dir: Path | None) -> tuple[bool, str]:
    r = state.get("report") or {}
    html = r.get("html")
    if not html or not Path(html).is_file():
        return False, "report.html 不存在"
    return True, "ok"


_CHECKS = {
    "target_prep": _check_target_prep,
    "screening": _check_screening,
    "binder": _check_binder,
    "vhh": _check_vhh,
    "md": _check_md,
    "report": _check_report,
}


def stage_complete(state: dict, stage: str,
                   project_dir: Path | None = None) -> tuple[bool, str]:
    """(complete, reason) for a stage from the state.json section."""
    check = _CHECKS.get(stage)
    if check is None:
        return False, f"未知阶段 {stage!r}"
    return check(state, project_dir)


# --------------------------------------------------------------------------- #
# cached summaries (same shapes the run_* tools return after a real run)
# --------------------------------------------------------------------------- #
def summarize_target_prep(state: dict) -> dict:
    out = state.get("target_prep") or {}
    return {k: out.get(k) for k in
            ("raw_pdb", "clean_pdb", "pocket", "receptor_pdb",
             "receptor_pdbqt", "ligand_pdbqt", "ligand_resnames",
             "judgment", "completeness")}


def summarize_screening(state: dict) -> dict:
    out = state.get("screening") or {}
    return {"n_docked": out.get("n_docked"),
            "reference_ligand_score": out.get("reference_ligand_score"),
            "hit_decision": out.get("hit_decision"),
            "n_hits": len(out.get("hits", [])),
            "top3": [{"rank": h["rank"], "smiles": h["smiles"],
                      "final_score": h.get("final_score")}
                     for h in out.get("hits", [])[:3]]}


def summarize_binder(state: dict) -> dict:
    out = state.get("binder") or {}
    return {"n_designs": out.get("n_designs"),
            "binder_type": out.get("binder_type"),
            "best": out.get("best")}


def summarize_vhh(state: dict) -> dict:
    out = state.get("vhh") or {}
    return {"pick": out.get("pick") or out.get("best"),
            "n_track_a": len(out.get("track_a_results",
                                    (out.get("track_a") or {}).get("results", []))
                              or []),
            "n_track_b": len((out.get("track_b") or {}).get("designs", []) or [])}


def summarize_md(state: dict) -> dict:
    out = state.get("md") or {}
    return {"system": out.get("system"), "ns": out.get("ns"),
            "reps": out.get("reps"),
            "final_rmsd_mean": (out.get("summary") or {}).get("final_rmsd_mean"),
            "final_rg_mean": (out.get("summary") or {}).get("final_rg_mean")}


def summarize_report(state: dict) -> dict:
    out = state.get("report") or {}
    return {"html": out.get("html"), "pdf": out.get("pdf")}


_SUMMARIES = {
    "target_prep": summarize_target_prep,
    "screening": summarize_screening,
    "binder": summarize_binder,
    "vhh": summarize_vhh,
    "md": summarize_md,
    "report": summarize_report,
}


def maybe_reuse(ctx: Ctx, stage: str, force: bool = False) -> dict | None:
    """Stage-level idempotency gate.

    Returns the cached summary (tagged ``reused=True``) when the stage's
    state section is complete and ``force`` is False; otherwise None, and
    the caller runs the stage for real."""
    if force:
        return None
    ok, why = stage_complete(ctx.state(), stage, ctx.project_dir)
    if not ok:
        logger.info(f"[stage:{stage}] 未完成 ({why}) → 重跑")
        return None
    summary = _SUMMARIES[stage](ctx.state())
    summary["reused"] = True
    logger.info(f"[stage:{stage}] 已完成 → 复用 state.{stage} "
                f"(force=true 强制重跑)")
    return summary


# --------------------------------------------------------------------------- #
# compact per-stage status lines (for `drugagent status`)
# --------------------------------------------------------------------------- #
def status_lines(state: dict, project_dir: Path | None = None) -> list[dict]:
    """One compact line per stage: done + key numbers or what's missing."""
    out = []
    tp = state.get("target_prep") or {}
    if tp:
        pocket = tp.get("pocket") or {}
        detail = f"pocket={pocket.get('method', '?')}"
        if tp.get("ligand_resnames"):
            detail += f", 配体={','.join(tp['ligand_resnames'][:3])}"
    else:
        detail = "未开始"
    out.append({"stage": "target_prep",
                "done": stage_complete(state, "target_prep", project_dir)[0],
                "detail": detail})

    scr = state.get("screening") or {}
    if scr:
        detail = (f"{scr.get('n_docked', '?')} 对接, "
                  f"{len(scr.get('hits', []))} 命中")
        if scr.get("library"):
            # keep the full label: "chembl35_small (fallback for dtp)" is
            # the substitution the user needs to see
            detail += f" [{scr['library']}]"
    else:
        detail = "未开始"
    out.append({"stage": "screening",
                "done": stage_complete(state, "screening", project_dir)[0],
                "detail": detail})

    b = state.get("binder") or {}
    if b:
        best = b.get("best") or {}
        detail = (f"{b.get('n_designs', '?')} 设计, "
                  f"best interface pLDDT="
                  f"{best.get('interface_plddt_mean', '?')}")
    else:
        detail = "未开始"
    out.append({"stage": "binder",
                "done": stage_complete(state, "binder", project_dir)[0],
                "detail": detail})

    v = state.get("vhh") or {}
    if v:
        ta = v.get("track_a") or {}
        tb = v.get("track_b") or {}
        detail = (f"trackA {ta.get('n_docked', '?')} 对接 / "
                  f"trackB {len(tb.get('designs', []) or [])} 设计")
        best = v.get("best") or {}
        if best:
            detail += f", 最优={best.get('source', '?')}"
    else:
        detail = "未开始"
    out.append({"stage": "vhh",
                "done": stage_complete(state, "vhh", project_dir)[0],
                "detail": detail})

    m = state.get("md") or {}
    if m:
        s = m.get("summary") or {}
        detail = (f"{m.get('ns', '?')} ns x {m.get('reps', '?')} 副本, "
                  f"final RMSD={s.get('final_rmsd_mean', '?')}")
    else:
        detail = "未开始"
    out.append({"stage": "md",
                "done": stage_complete(state, "md", project_dir)[0],
                "detail": detail})

    r = state.get("report") or {}
    done_r, why_r = stage_complete(state, "report", project_dir)
    if done_r:
        detail = str(r.get("html", ""))
    else:
        detail = why_r if not r else "未生成"
    out.append({"stage": "report", "done": done_r, "detail": detail})
    return out
