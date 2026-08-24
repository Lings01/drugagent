"""Report tools."""
from __future__ import annotations

from .loop import Ctx, Tool


def build_report(ctx: Ctx) -> dict:
    from ..report import report as rep
    st = ctx.state()
    out = rep.build_report(st)
    st["report"] = out
    ctx.save_state(report=out)
    return {"html": out.get("html"), "pdf": out.get("pdf")}


def build() -> list[Tool]:
    return [
        Tool("build_report",
             "生成最终 HTML + PDF 报告 (汇总所有阶段 + 决策日志)。"
             "任务收尾前调用。",
             {"type": "object", "properties": {}, "required": []},
         build_report),
    ]
