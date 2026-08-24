"""VHH (nanobody) tools: library screening + de novo design."""
from __future__ import annotations

from pathlib import Path

from ..config import DEFAULTS
from ..modules import vhh as vh
from .loop import Ctx, Tool, ToolError
from .stages import maybe_reuse, summarize_vhh


def _prep(ctx: Ctx) -> dict:
    prep = ctx.state().get("target_prep") or {}
    if not prep.get("receptor_pdbqt"):
        raise ToolError("state.target_prep 未就绪 (先做靶点准备)")
    return prep


def vhh_screen(ctx: Ctx, n: int | None = None) -> dict:
    """Track A: screen a generated VHH library against the target."""
    _prep(ctx)
    from ..config import resolve_defaults
    d = resolve_defaults(ctx.options)
    n = n or d.vhh_screen_n
    n_jobs = int(ctx.options.get("vhh_n_jobs", 4))
    workdir = ctx.stage_dir("04_vhh")
    st = ctx.state()
    out = vh.screen_vhh(st, workdir, n=n, n_jobs=n_jobs)
    v = st.get("vhh") or {}
    v["track_a"] = out
    st["vhh"] = v
    ctx.save_state()
    res = out.get("results", [])
    return {"n_screened": len(res),
            "top5": [{"idx": r.get("idx"), "score": r.get("score"),
                      "plddt": r.get("plddt")} for r in res[:5]],
            "n_docked": sum(1 for r in res if r.get("score") is not None)}


def vhh_design(ctx: Ctx, n_designs: int = 2) -> dict:
    """Track B: de novo nanobody design on a VHH scaffold (RF scaffold-guided)."""
    _prep(ctx)
    workdir = ctx.stage_dir("04_vhh")
    st = ctx.state()
    out = vh.design_vhh(st, workdir, n_designs=n_designs)
    v = st.get("vhh") or {}
    v["track_b"] = out
    st["vhh"] = v
    ctx.save_state()
    designs = out.get("designs", [])
    return {"n_designs": len(designs),
            "top": [{"complex_plddt": d.get("complex_plddt"),
                     "interface_plddt_mean": d.get("interface_plddt_mean"),
                     "design": d.get("design")} for d in designs[:3]]}


def run_vhh(ctx: Ctx, force: bool = False) -> dict:
    """Full deterministic VHH stage (track A + B + merged pick).

    R11/G8: reuses state.vhh when both tracks are present; force=true to
    re-run."""
    cached = maybe_reuse(ctx, "vhh", force)
    if cached is not None:
        return cached
    st = ctx.state()
    state_out = vh.design_vhh_all(st)   # returns full state
    ctx.save_state(vhh=state_out["vhh"])
    return summarize_vhh(ctx.state())


# --------------------------------------------------------------------------- #
def build() -> list[Tool]:
    return [
        Tool("vhh_screen",
             "VHH 轨道A: 生成库筛选 (建模+打分+对接)。",
             {"type": "object",
              "properties": {"n": {"type": "integer"}}, "required": []},
         vhh_screen),
        Tool("vhh_design",
             "VHH 轨道B: 基于 VHH scaffold 的 RF 从头设计 (scaffold-guided)。",
             {"type": "object",
              "properties": {"n_designs": {"type": "integer"}}, "required": []},
         vhh_design),
        Tool("run_vhh",
             "整段 VHH (轨道A+轨道B+合并择优, 确定性标准流程)。"
             "已完成阶段自动复用, force=true 强制重跑。",
             {"type": "object",
              "properties": {"force": {"type": "boolean",
                                       "description": "强制重跑 (默认复用已完成阶段)"}},
              "required": []}, run_vhh),
    ]
