"""ReAct main loop: the LLM drives the drug-discovery toolkit.

- ``Ctx``        shared project context passed to every tool.
- ``Tool``       one callable tool (name + JSON schema + handler).
- ``AgentLoop``  the ReAct loop with transcript, budget and checkpoints.
"""
from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from loguru import logger

from ..utils import jload, jsave


class ToolError(Exception):
    """Expected tool failure (the agent can react to the message)."""


# --------------------------------------------------------------------------- #
# project context
# --------------------------------------------------------------------------- #
class Ctx:
    """Everything a tool may need. The project dir is the single source of
    truth; ``state.json`` keeps the 1.0 state shape so reports/status work
    unchanged."""

    STATE_KEYS = ("target_prep", "screening", "binder", "vhh", "md", "report")

    def __init__(self, project_dir: str | Path, brain, options: dict, *,
                 auto: bool = False, prompt_fn: Callable | None = None,
                 max_result_chars: int = 6000, target: dict | None = None):
        # absolute: tools resolve relative paths against this, and resume may
        # be invoked with a relative --project path
        self.project_dir = Path(project_dir).expanduser().resolve()
        self.project_dir.mkdir(parents=True, exist_ok=True)
        self.brain = brain
        self.options = options or {}
        self.auto = auto
        self.prompt_fn = prompt_fn
        self.max_result_chars = max_result_chars
        self.agent_dir = self.project_dir / "agent"
        self.out_dir = self.agent_dir / "out"
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self._step = 0
        self._state: dict | None = None
        state_path = self.project_dir / "state.json"
        if state_path.exists():
            self._state = jload(state_path)
        if self._state is None:
            self._state = {
                "project_dir": str(self.project_dir),
                "target": target or {},
                "options": self.options,
                "status": "running",
                "decisions": [],
                "errors": [],
            }
            jsave(state_path, self._state)
        else:
            self._state["options"] = self.options

    # -- state --------------------------------------------------------------
    def state(self) -> dict:
        if self._state is None:
            self._state = jload(self.project_dir / "state.json")
        return self._state

    def save_state(self, **updates) -> None:
        st = self.state()
        st.update(updates)
        jsave(self.project_dir / "state.json", st)
        self._state = st

    def stage_dir(self, name: str) -> Path:
        p = self.project_dir / name
        p.mkdir(parents=True, exist_ok=True)
        return p

    # -- decisions ----------------------------------------------------------
    def record_decision(self, stage: str, choice: Any, rationale: str = "") -> dict:
        dec = {"node": stage, "question": "agent decision",
               "answer": choice, "rationale": rationale,
               "ts": time.strftime("%Y-%m-%dT%H:%M:%S")}
        st = self.state()
        st.setdefault("decisions", []).append(dec)
        self.save_state()
        # decisions.json (same format as 1.0; report reads it)
        log_path = self.project_dir / "decisions.json"
        log = jload(log_path) if log_path.exists() else []
        log.append(dec)
        jsave(log_path, log)
        logger.info(f"[decision@{stage}] {choice}")
        return dec

    # -- output helpers -------------------------------------------------------
    def dump_output(self, name: str, obj: Any) -> Path:
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name)[:80]
        p = self.out_dir / f"step{self._step:03d}_{safe}.json"
        jsave(p, obj)
        return p

    def serialize(self, obj: Any) -> str:
        """JSON-serialize a tool result; truncate large payloads and point at
        the full file on disk."""
        try:
            s = json.dumps(obj, ensure_ascii=False, default=str)
        except Exception:  # noqa: BLE001
            s = str(obj)
        if len(s) <= self.max_result_chars:
            return s
        p = self.dump_output("result", obj)
        return s[: self.max_result_chars] + \
            f' …[截断, 完整 {len(s)} 字符见 {p}]'


# --------------------------------------------------------------------------- #
# tool registry
# --------------------------------------------------------------------------- #
@dataclass
class Tool:
    name: str
    description: str
    parameters: dict                      # JSON Schema (object root)
    fn: Callable[[Ctx], Any] = field(repr=False)
    long_running: bool = False            # hint for the system prompt

    def schema(self) -> dict:
        return {"type": "function",
                "function": {"name": self.name,
                             "description": self.description,
                             "parameters": self.parameters}}

    def call(self, ctx: Ctx, **args: Any) -> dict:
        try:
            res = self.fn(ctx, **_resolve_ctx_paths(ctx, args))
        except ToolError as e:
            return {"ok": False, "error": str(e)}
        except Exception as e:  # noqa: BLE001
            logger.exception(f"tool {self.name} failed")
            return {"ok": False, "error": f"{type(e).__name__}: {e}"}
        if not isinstance(res, dict):
            res = {"result": res}
        res.setdefault("ok", True)
        return res


def _resolve_ctx_paths(ctx: "Ctx", args: dict) -> dict:
    """R18: LLM/scripted steps pass stage-relative paths
    ('01_target/receptor.pdb'). Tools already re-root OUTPUT paths via
    _p(), but INPUT paths used to stay CWD-relative — fatal for
    `run --root .` from any folder. Central fix: a relative arg that
    EXISTS under the project dir and NOT relative to CWD is re-rooted.
    Existence-based, so SMILES/resname-like strings are never touched."""
    for k, v in list(args.items()):
        if not (isinstance(v, str) and "/" in v and len(v) < 400):
            continue
        p = Path(v)
        if p.is_absolute():
            continue
        proj = ctx.project_dir / p
        try:
            if proj.exists() and not p.exists():
                args[k] = str(proj)
        except OSError:
            pass
    return args


def build_tools() -> list[Tool]:
    from .tools_core import build as core
    from .tools_target import build as target
    from .tools_screen import build as screen
    from .tools_design import build as design
    from .tools_vhh import build as vhh
    from .tools_md import build as md
    from .tools_report import build as report
    return core() + target() + screen() + design() + vhh() + md() + report()


# --------------------------------------------------------------------------- #
# ReAct loop
# --------------------------------------------------------------------------- #
class AgentLoop:
    """LLM main loop. One iteration = one chat completion plus the tool calls
    it requested. Everything is appended to ``agent/transcript.jsonl`` so a
    crashed/paused run can be resumed verbatim."""

    def __init__(self, ctx: Ctx, tools: list[Tool], system: str, goal: str,
                 *, max_steps: int = 300, max_tokens: int = 8192):
        self.ctx = ctx
        self.tools = {t.name: t for t in tools}
        self.system = system
        self.goal = goal
        self.max_steps = max_steps
        # thinking models put long reasoning into content; 3000 tokens often
        # truncates the turn before the tool call appears
        self.max_tokens = max_tokens
        self.transcript_path = ctx.agent_dir / "transcript.jsonl"
        self.finished: dict | None = None
        self.step = 0
        self.messages: list[dict] = []
        self._rebuild_messages()

    # -- transcript ---------------------------------------------------------
    def _rebuild_messages(self) -> None:
        msgs: list[dict] = [
            {"role": "system", "content": self.system},
            {"role": "user", "content": self.goal},
        ]
        if self.transcript_path.exists():
            for line in self.transcript_path.read_text().splitlines():
                if not line.strip():
                    continue
                try:
                    obj = json.loads(line)
                except json.JSONDecodeError:
                    continue
                self.step = max(self.step, int(obj.get("step", 0)))
                role = obj.get("role")
                if role == "assistant":
                    m: dict = {"role": "assistant",
                               "content": obj.get("content") or ""}
                    tcs = obj.get("tool_calls") or []
                    if tcs:
                        m["tool_calls"] = [
                            {"id": tc["id"], "type": "function",
                             "function": {"name": tc["name"],
                                          "arguments": json.dumps(tc.get("args", {}),
                                                                 ensure_ascii=False)}}
                            for tc in tcs]
                    msgs.append(m)
                elif role == "tool":
                    msgs.append({"role": "tool",
                                 "tool_call_id": obj.get("tool_call_id", ""),
                                 "content": obj.get("content", "")})
                elif role in ("user", "system"):
                    msgs.append({"role": role, "content": obj.get("content", "")})
        self.messages = msgs

    def _append(self, obj: dict) -> None:
        with self.transcript_path.open("a") as fh:
            fh.write(json.dumps(obj, ensure_ascii=False) + "\n")

    # -- main ----------------------------------------------------------------
    def run(self) -> dict:
        nudges = 0
        while self.finished is None and self.step < self.max_steps:
            self.step += 1
            self.ctx._step = self.step
            logger.info(f"--- agent step {self.step}/{self.max_steps} ---")
            try:
                msg = self.ctx.brain.chat_tools(
                    self.messages, [t.schema() for t in self.tools.values()],
                    max_tokens=self.max_tokens)
            except Exception as e:  # noqa: BLE001
                logger.exception("LLM endpoint failed")
                self._finish({"status": "needs_human",
                              "summary": f"LLM 端点调用失败: {e} "
                                         f"(step {self.step})"})
                break

            tcs = msg.get("tool_calls") or []
            if not tcs:
                content = (msg.get("content") or "").strip()
                self._append({"step": self.step, "role": "assistant",
                              "content": content})
                nudges += 1
                if nudges >= 3 or re.search(r"\bfinish\b", content, re.I) \
                        and nudges >= 1:
                    self._finish({"status": "needs_human",
                                  "summary": content[:400] or "agent 未调用工具"})
                    break
                nudge = ("请调用一个工具继续工作；若任务已全部完成，请调用 "
                         "finish 工具并给出 summary。")
                self.messages.append({"role": "user", "content": nudge})
                self._append({"step": self.step, "role": "user", "content": nudge})
                continue
            nudges = 0

            self._append({"step": self.step, "role": "assistant",
                          "content": msg.get("content") or "",
                          "tool_calls": [
                              {"id": tc["id"], "name": tc["name"],
                               "args": tc.get("args", {})} for tc in tcs]})
            for tc in tcs:
                name = tc["name"]
                args = tc.get("args") or {}
                tool = self.tools.get(name)
                if tool is None:
                    result = {"ok": False,
                              "error": f"未知工具 {name!r}; 可用: {sorted(self.tools)}"}
                else:
                    result = tool.call(self.ctx, **args)
                res_text = self.ctx.serialize(result)
                self.messages.append({"role": "tool",
                                      "tool_call_id": tc.get("id", ""),
                                      "content": res_text})
                self._append({"step": self.step, "role": "tool",
                              "tool_call_id": tc.get("id", ""), "name": name,
                              "args": args, "content": res_text})
                logger.info(f"  tool {name} -> ok={result.get('ok')}")
                if name == "finish" and result.get("ok"):
                    self._finish(result)
                    break
            if self.finished is not None:
                break

        if self.finished is None:
            self._finish({"status": "needs_human",
                          "summary": f"预算耗尽 (max_steps={self.max_steps})，"
                                     f"停在 step {self.step}"})
        return self.finished

    def _finish(self, result: dict) -> None:
        status = result.get("status", "needs_human")
        self.ctx.save_state(status=status)
        self.finished = {"status": status,
                         "summary": result.get("summary", ""),
                         "step": self.step}
        logger.info(f"agent finished: {status} @ step {self.step}")


# --------------------------------------------------------------------------- #
# no-LLM mode: deterministic fallback through the stage tools
# --------------------------------------------------------------------------- #
def scripted_run(ctx: Ctx, tools: list[Tool]) -> dict:
    """--no-llm: run the standard stage tools in a fixed order.

    The agent concept degrades to a script; same state/report contract.
    """
    name2 = {t.name: t for t in tools}
    opts = ctx.options
    mods = opts.get("modules", ["screen", "binder", "vhh", "md"])
    seq: list[tuple[str, dict]] = [("run_target_prep", {}),
                                   ("checkpoint", {"stage": "target",
                                                   "summary": "靶点准备完成"})]
    if "screen" in mods:
        seq += [("run_screening", {}),
                ("checkpoint", {"stage": "screening",
                                "summary": "小分子筛选完成"})]
    if "binder" in mods:
        seq += [("run_design", {}),
                ("checkpoint", {"stage": "design",
                                "summary": "binder 设计完成"})]
    if "vhh" in mods:
        seq.append(("run_vhh", {}))
    if "md" in mods:
        seq += [("checkpoint", {"stage": "md", "summary": "MD 参数确认"}),
                ("run_md", {})]
    seq.append(("build_report", {}))

    transcript = ctx.agent_dir / "transcript.jsonl"
    step = 0
    for name, args in seq:
        tool = name2.get(name)
        if tool is None:
            continue
        step += 1
        ctx._step = step
        result = tool.call(ctx, **args)
        with transcript.open("a") as fh:
            fh.write(json.dumps({"step": step, "role": "tool",
                                 "name": name, "args": args,
                                 "content": ctx.serialize(result)},
                                ensure_ascii=False) + "\n")
        logger.info(f"scripted {name} -> ok={result.get('ok')}")
        if name == "build_report" and not result.get("ok"):
            return {"status": "needs_human",
                    "summary": f"build_report 失败: {result.get('error')}",
                    "step": step}
    ctx.save_state(status="success")
    return {"status": "success",
            "summary": "脚本模式 (no-llm) 完成全部阶段", "step": step}
