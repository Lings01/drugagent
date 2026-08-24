"""Report tools."""
from __future__ import annotations

from .loop import Ctx, Tool
from .stages import maybe_reuse, summarize_report


def build_report(ctx: Ctx, force: bool = False) -> dict:
    """R11/G8: reuses an existing report (state.report + report.html on
    disk); force=true to regenerate."""
    cached = maybe_reuse(ctx, "report", force)
    if cached is not None:
        return cached
    from ..report import report as rep
    st = ctx.state()
    out = rep.build_report(st)
    st["report"] = out
    ctx.save_state(report=out)
    return summarize_report(ctx.state())


def build() -> list[Tool]:
    return [
        Tool("build_report",
             "生成最终 HTML + PDF 报告 (汇总所有阶段 + 决策日志)。"
             "任务收尾前调用。已存在时自动复用, force=true 重新生成。",
             {"type": "object",
              "properties": {"force": {"type": "boolean",
                                       "description": "强制重新生成"}},
              "required": []},
         build_report),
    ]
