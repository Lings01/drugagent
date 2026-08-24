"""Meta tools: file access, shell, decisions, human interaction, finish."""
from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path

from loguru import logger

from .loop import Ctx, Tool, ToolError

FIXED_CHECKPOINTS = ("target", "screening", "design", "md")


def _resolve(ctx: Ctx, path: str) -> Path:
    p = Path(path)
    if p.is_absolute():
        return p
    return ctx.project_dir / p


# --------------------------------------------------------------------------- #
# filesystem
# --------------------------------------------------------------------------- #
def fs_list(ctx: Ctx, path: str = ".") -> dict:
    p = _resolve(ctx, path)
    if not p.exists():
        raise ToolError(f"path not found: {p}")
    if p.is_file():
        return {"path": str(p), "is_file": True,
                "size": p.stat().st_size}
    entries = []
    for e in sorted(p.iterdir()):
        entries.append({"name": e.name, "dir": e.is_dir(),
                        "size": 0 if e.is_dir() else e.stat().st_size})
    return {"path": str(p), "entries": entries[:300]}


def fs_read(ctx: Ctx, path: str, offset: int = 1, limit: int = 400) -> dict:
    p = _resolve(ctx, path)
    if not p.is_file():
        raise ToolError(f"not a file: {p}")
    lines = p.read_text(errors="ignore").splitlines()
    chunk = lines[max(0, offset - 1): max(0, offset - 1) + limit]
    text = "\n".join(f"{i + offset}: {l}" for i, l in enumerate(chunk))
    return {"path": str(p), "n_lines_total": len(lines),
            "showing": f"{max(1, offset)}-{offset + len(chunk) - 1}",
            "content": text[: 8000]}


def fs_write(ctx: Ctx, path: str, content: str) -> dict:
    p = _resolve(ctx, path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(content)
    return {"path": str(p), "bytes": len(content.encode())}


def fs_edit(ctx: Ctx, path: str, old: str, new: str,
            replace_all: bool = False) -> dict:
    p = _resolve(ctx, path)
    if not p.is_file():
        raise ToolError(f"not a file: {p}")
    text = p.read_text()
    n = text.count(old)
    if n == 0:
        raise ToolError("old_string not found in file")
    if n > 1 and not replace_all:
        raise ToolError(f"old_string appears {n} times; give a more specific "
                        "string or set replace_all=true")
    p.write_text(text.replace(old, new) if replace_all else text.replace(old, new, 1))
    return {"path": str(p), "replaced": n if replace_all else 1}


# --------------------------------------------------------------------------- #
# shell
# --------------------------------------------------------------------------- #
def shell(ctx: Ctx, cmd: str, timeout: int = 600, cwd: str | None = None) -> dict:
    workdir = _resolve(ctx, cwd) if cwd else ctx.project_dir
    workdir.mkdir(parents=True, exist_ok=True)
    log_path = ctx.out_dir / f"step{ctx._step:03d}_shell_{int(time.time())}.log"
    t0 = time.time()
    try:
        r = subprocess.run(cmd, shell=True, cwd=workdir, capture_output=True,
                           text=True, timeout=timeout)
        out, err, code = r.stdout, r.stderr, r.returncode
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or b"").decode(errors="ignore") if isinstance(e.stdout, bytes) else (e.stdout or "")
        err, code = "timeout", 124
    with open(log_path, "w") as fh:
        fh.write(f"$ {cmd}\n--- stdout ---\n{out}\n--- stderr ---\n{err}\n")
    return {"ok": code == 0, "exit": code,
            "seconds": round(time.time() - t0, 1),
            "stdout_tail": out[-4000:],
            "stderr_tail": err[-2000:],
            "log": str(log_path)}


# --------------------------------------------------------------------------- #
# state / decisions
# --------------------------------------------------------------------------- #
def state_get(ctx: Ctx, key: str | None = None) -> dict:
    st = ctx.state()
    if key is None:
        # compact view: stage keys + status, not full payloads
        view = {k: st.get(k) for k in ("status", "options")}
        for k in Ctx.STATE_KEYS:
            v = st.get(k)
            view[k] = v if isinstance(v, dict) else v
            if isinstance(v, dict):
                view[k] = {kk: vv for kk, vv in v.items()
                           if not isinstance(vv, (list, dict)) or kk in (
                               "pocket", "hit_decision", "system", "gromacs")}
        return view
    return {"key": key, "value": st.get(key)}


def record_decision(ctx: Ctx, stage: str, choice, rationale: str = "") -> dict:
    return {"ok": True,
            "decision": ctx.record_decision(stage, choice, rationale)}


# --------------------------------------------------------------------------- #
# human interaction
# --------------------------------------------------------------------------- #
def _prompt(ctx: Ctx, title: str, question: str, options: list[str] | None,
            context: str = "") -> str:
    if ctx.auto:
        return ""
    if ctx.prompt_fn is not None:
        return ctx.prompt_fn(title, question, options, context)
    # last resort: plain stdin
    print("\n" + "=" * 72)
    print(f"[{title}] {question}")
    if context:
        print(context[:1500])
    if options:
        for i, o in enumerate(options):
            print(f"  {i + 1}. {o}")
    return input("回答/选择 > ").strip()


def ask_human(ctx: Ctx, question: str, options: list[str] | None = None,
              context: str = "") -> dict:
    """Dynamic human checkpoint, callable at any moment."""
    if ctx.auto:
        return {"ok": True, "auto": True, "answer": None,
                "hint": "auto 模式：请自行做出你认为最优的选择，"
                        "并用 record_decision 记录"}
    answer = _prompt(ctx, "动态确认", question, options, context)
    return {"ok": True, "auto": False, "answer": answer}


def checkpoint(ctx: Ctx, stage: str, summary: str,
               options: list[str] | None = None) -> dict:
    """Fixed milestone checkpoint (target/screening/design/md)."""
    if stage not in FIXED_CHECKPOINTS:
        raise ToolError(f"stage must be one of {FIXED_CHECKPOINTS}")
    options = options or ["批准并继续"]
    if ctx.auto:
        choice = options[0]
        ctx.record_decision(f"checkpoint_{stage}", choice,
                            f"auto 通过: {summary[:200]}")
        return {"ok": True, "auto": True, "approved": True, "choice": choice}
    answer = _prompt(ctx, f"里程碑检查点: {stage}",
                     "请审阅后选择", options, summary)
    approved = answer.lower() not in ("x", "abort", "中止", "0")
    ctx.record_decision(f"checkpoint_{stage}",
                        f"{'批准' if approved else '中止'}: {answer}",
                        summary[:300])
    return {"ok": True, "auto": False, "approved": approved, "choice": answer}


# --------------------------------------------------------------------------- #
# finish
# --------------------------------------------------------------------------- #
def finish(ctx: Ctx, status: str, summary: str) -> dict:
    if status not in ("success", "needs_human", "failed"):
        raise ToolError("status must be success|needs_human|failed")
    return {"ok": True, "status": status, "summary": summary}


# --------------------------------------------------------------------------- #
def build() -> list[Tool]:
    return [
        Tool("list_dir",
             "列出目录内容（路径相对项目目录，也接受绝对路径）。",
             {"type": "object", "properties": {"path": {"type": "string"}},
              "required": []}, fs_list),
        Tool("read_file",
             "读取文本文件（PDB/MDP/JSON/日志），返回带行号的内容（截断到 8000 字符）。",
             {"type": "object",
              "properties": {"path": {"type": "string"},
                             "offset": {"type": "integer"},
                             "limit": {"type": "integer"}},
              "required": ["path"]}, fs_read),
        Tool("write_file",
             "写文件（创建父目录）。用于写 MDP、JSON 等；覆盖已有文件。",
             {"type": "object",
              "properties": {"path": {"type": "string"},
                             "content": {"type": "string"}},
              "required": ["path", "content"]}, fs_write),
        Tool("edit_file",
             "对文件做精确文本替换（old_string 必须唯一，除非 replace_all）。",
             {"type": "object",
              "properties": {"path": {"type": "string"}, "old": {"type": "string"},
                             "new": {"type": "string"},
                             "replace_all": {"type": "boolean"}},
              "required": ["path", "old", "new"]}, fs_edit),
        Tool("shell",
             "运行 shell 命令（逃生舱；优先用专用工具）。返回退出码 + 输出尾部 + 完整日志路径。",
             {"type": "object",
              "properties": {"cmd": {"type": "string"},
                             "timeout": {"type": "integer",
                                         "description": "秒, 默认 600"},
                             "cwd": {"type": "string"}},
              "required": ["cmd"]}, shell),
        Tool("state_get",
             "查看项目状态（state.json 摘要；key 可取 target_prep/screening/"
             "binder/vhh/md/report 取完整阶段结果）。",
             {"type": "object",
              "properties": {"key": {"type": "string"}}, "required": []},
         state_get),
        Tool("record_decision",
             "记录一个关键判断（进入 decisions.json，报告会展示）。stage 用 "
             "target/screening/design/vhh/md 等。",
             {"type": "object",
              "properties": {"stage": {"type": "string"},
                             "choice": {},
                             "rationale": {"type": "string"}},
              "required": ["stage", "choice"]}, record_decision),
        Tool("ask_human",
             "动态人工确认（任何时刻可用；auto 模式返回提示让你自行决策）。",
             {"type": "object",
              "properties": {"question": {"type": "string"},
                             "options": {"type": "array", "items": {"type": "string"}},
                             "context": {"type": "string"}},
              "required": ["question"]}, ask_human),
        Tool("checkpoint",
             "里程碑检查点（target/screening/design/md）。auto 模式自动通过并记录。",
             {"type": "object",
              "properties": {"stage": {"type": "string",
                                       "enum": list(FIXED_CHECKPOINTS)},
                             "summary": {"type": "string",
                                         "description": "给人看的阶段总结"},
                             "options": {"type": "array", "items": {"type": "string"}}},
              "required": ["stage", "summary"]}, checkpoint),
        Tool("finish",
             "结束任务。status: success(全部完成)/needs_human(需人介入)/failed。",
             {"type": "object",
              "properties": {"status": {"type": "string",
                                        "enum": ["success", "needs_human", "failed"]},
                             "summary": {"type": "string"}},
              "required": ["status", "summary"]}, finish),
    ]
