"""LangGraph state machine for the DrugAgent pipeline.

init -> target_prep -> [cp1] -> screening? -> [cp2] -> binder? -> vhh?
      -> [cp3] -> md? -> [cp4] -> report -> END
Human checkpoints use langgraph interrupt(); --auto skips them.
"""
from __future__ import annotations

import json
from pathlib import Path

from langgraph.checkpoint.sqlite import SqliteSaver
from langgraph.graph import END, START, StateGraph
from langgraph.types import Command, interrupt
from loguru import logger

from .config import PROJECTS
from .modules import binder, mdsim, screening, target_prep, vhh
from .report import report as report_mod
from .state import AgentState

ALL_MODULES = ["screen", "binder", "vhh", "md"]


# --------------------------------------------------------------------------- #
# checkpoint helper
# --------------------------------------------------------------------------- #
def _make_checkpoint(name: str, summarize):
    def node(state: dict) -> dict:
        opts = state.get("options", {})
        summary = summarize(state)
        if opts.get("auto"):
            logger.info(f"[checkpoint:{name}] auto-approved")
            return {}
        payload = interrupt({
            "checkpoint": name,
            "summary": summary,
            "choices": ["approve", "modify", "abort"],
        })
        payload = payload or {"action": "approve"}
        action = payload.get("action", "approve")
        if action == "abort":
            raise RunAborted(name)
        if action == "modify" and payload.get("params"):
            merged = dict(state.get("options", {}))
            merged.update(payload["params"])
            return {"options": merged}
        return {}
    node.__name__ = f"cp_{name}"
    return node


class RunAborted(Exception):
    pass


# --------------------------------------------------------------------------- #
# summarizers for humans
# --------------------------------------------------------------------------- #
def _sum_target(state: dict) -> str:
    t = state.get("target_prep", {})
    c = t.get("completeness", {})
    p = t.get("pocket", {})
    return (
        f"靶点: chains={c.get('chains')} 配体={c.get('ligands')}\n"
        f"判断: {t.get('judgment', {}).get('action')} - {t.get('judgment', {}).get('rationale', '')}\n"
        f"口袋: {p.get('method')} center={p.get('center')} size={p.get('xsize')}x{p.get('ysize')}x{p.get('zsize')}"
    )


def _sum_screen(state: dict) -> str:
    s = state.get("screening", {})
    h = s.get("hit_decision", {})
    hits = s.get("hits", [])[:5]
    lines = [
        f"筛选: 标准化 {s.get('n_standardized')} -> 过滤后 {s.get('n_after_filter')} -> 对接 {s.get('n_docked')}",
        f"参考配体打分: {s.get('reference_ligand_score')}",
        f"命中判定: threshold={h.get('threshold')} n_hits={h.get('n_hits')} ({h.get('rationale')})",
    ]
    for hit in hits:
        lines.append(f"  top{hit['rank']}: {hit['smiles'][:40]} score={hit.get('final_score')}")
    return "\n".join(lines)


def _sum_design(state: dict) -> str:
    b = state.get("binder", {})
    v = state.get("vhh", {})
    lines = []
    if b:
        best = b.get("best") or {}
        lines.append(f"Binder: {b.get('binder_type')} x{b.get('n_designs')}, "
                      f"best if-pLDDT={best.get('interface_plddt_mean')}")
    if v:
        best = v.get("best") or {}
        lines.append(f"VHH: 筛选{v.get('track_a', {}).get('n_docked', 0)}个 + "
                      f"de novo {v.get('track_b', {}).get('n_designs', 0)}个, "
                      f"best composite={best.get('composite_score')}")
    return "\n".join(lines)


def _sum_md(state: dict) -> str:
    m = state.get("md", {})
    s = m.get("summary", {})
    return (
        f"MD: {m.get('system', {}).get('label')} ns={m.get('ns')} reps={m.get('reps')}\n"
        f"RMSD final mean = {s.get('final_rmsd_mean')} A, Rg final = {s.get('final_rg_mean')} A\n"
        f"cluster populations: { {r['rep']: r.get('clusters') for r in m.get('replicas', [])} }"
    )


# --------------------------------------------------------------------------- #
# nodes
# --------------------------------------------------------------------------- #
def init_node(state: dict) -> dict:
    opts = state.setdefault("options", {})
    modules = set(opts.get("modules", ALL_MODULES))
    unknown = modules - set(ALL_MODULES)
    if unknown:
        raise ValueError(f"unknown modules: {unknown}")
    pdir = Path(state["project_dir"])
    pdir.mkdir(parents=True, exist_ok=True)
    (pdir / "run_config.json").write_text(json.dumps(
        {"target": state["target"], "options": opts}, indent=2, default=str))
    return {"status": "running", "modules": sorted(modules)}


def route_after_target(state: dict) -> str:
    mods = set(state.get("options", {}).get("modules", ALL_MODULES))
    if "screen" in mods:
        return "screening"
    if "binder" in mods:
        return "binder"
    if "vhh" in mods:
        return "vhh"
    if "md" in mods:
        return "md"
    return "report"


def route_after_screen(state: dict) -> str:
    mods = set(state.get("options", {}).get("modules", ALL_MODULES))
    if "binder" in mods:
        return "binder"
    if "vhh" in mods:
        return "vhh"
    if "md" in mods:
        return "md"
    return "report"


def route_after_binder(state: dict) -> str:
    mods = set(state.get("options", {}).get("modules", ALL_MODULES))
    if "vhh" in mods:
        return "vhh"
    return "cp_design"


def route_after_vhh(state: dict) -> str:
    mods = set(state.get("options", {}).get("modules", ALL_MODULES))
    if "md" in mods:
        return "md"
    return "report"


def report_node(state: dict) -> dict:
    try:
        out = report_mod.build_report(state)
        state_out = dict(state)
        state_out["report"] = out
        return state_out
    except Exception as e:  # noqa: BLE001
        logger.exception(f"report generation failed: {e}")
        state_out = dict(state)
        state_out["report"] = {"error": str(e)}
        state_out["errors"] = state.get("errors", []) + [{"node": "report", "error": str(e)}]
        return state_out


# --------------------------------------------------------------------------- #
# graph construction
# --------------------------------------------------------------------------- #
def build_graph():
    g = StateGraph(AgentState)
    g.add_node("init", init_node)
    g.add_node("target_prep", target_prep.prepare_target)
    g.add_node("cp_target", _make_checkpoint("target", _sum_target))
    g.add_node("screening", screening.screen)
    g.add_node("cp_screen", _make_checkpoint("screening", _sum_screen))
    g.add_node("binder", binder.design_binder)
    g.add_node("vhh", vhh.design_vhh_all)
    g.add_node("cp_design", _make_checkpoint("design", _sum_design))
    g.add_node("md", mdsim.run_md)
    g.add_node("cp_md", _make_checkpoint("md", _sum_md))
    g.add_node("report", report_node)

    g.add_edge(START, "init")
    g.add_edge("init", "target_prep")
    g.add_edge("target_prep", "cp_target")
    g.add_conditional_edges("cp_target", route_after_target,
                            {"screening": "screening", "binder": "binder",
                             "vhh": "vhh", "md": "md", "report": "report"})
    g.add_edge("screening", "cp_screen")
    g.add_conditional_edges("cp_screen", route_after_screen,
                            {"binder": "binder", "vhh": "vhh", "md": "md",
                             "report": "report"})
    g.add_conditional_edges("binder", route_after_binder,
                            {"vhh": "vhh", "cp_design": "cp_design"})
    g.add_edge("vhh", "cp_design")
    g.add_conditional_edges("cp_design", lambda s: "md"
                            if "md" in set(s.get("options", {}).get("modules", []))
                            else "report",
                            {"md": "md", "report": "report"})
    g.add_edge("md", "cp_md")
    g.add_edge("cp_md", "report")
    g.add_edge("report", END)
    return g


def compile_graph(project_dir: str | Path):
    import sqlite3
    pdir = Path(project_dir)
    pdir.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(pdir / "graph.sqlite"),
                           check_same_thread=False)
    checkpointer = SqliteSaver(conn)
    if hasattr(checkpointer, "__enter__"):
        checkpointer.__enter__()
    if hasattr(checkpointer, "setup"):
        checkpointer.setup()
    return build_graph().compile(checkpointer=checkpointer), conn
