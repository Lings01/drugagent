"""Tests for the 2.0 agent core (ReAct loop + tools) with a fake LLM."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from drugagent.agent import build_tools, Ctx
from drugagent.agent.loop import AgentLoop, ToolError


class FakeBrain:
    """Scripted chat_tools: pops canned responses, one per call."""

    def __init__(self, responses: list[dict]):
        self.responses = list(responses)
        self.calls = 0

    def chat_tools(self, messages, tools=None, **kw):
        self.calls += 1
        if not self.responses:
            return {"content": "done?", "tool_calls": None}
        r = self.responses.pop(0)
        # count tool messages in the history for state checks
        return r


def _mk_ctx(tmp_path: Path, auto: bool = True) -> Ctx:
    return Ctx(tmp_path / "proj", None,
               {"modules": ["md"], "fast": True, "auto": auto,
                "no_llm": False},
               auto=auto, target={"kind": "pdb_id", "value": "1HVI"})


def _finish_call(status="success", summary="ok"):
    return {"content": "", "tool_calls": [
        {"id": "c1", "name": "finish",
         "args": {"status": status, "summary": summary}}]}


# --------------------------------------------------------------------------- #
def test_tool_schemas():
    tools = build_tools()
    assert len(tools) >= 35
    names = {t.name for t in tools}
    for need in ("read_file", "write_file", "find_pockets", "dock_screen",
                 "rf_design", "mdp_template", "grompp", "mdrun",
                 "gmx_analyze", "checkpoint", "ask_human", "finish",
                 "build_report", "run_md", "gmx_env"):
        assert need in names
    for t in tools:
        s = t.schema()
        assert s["function"]["name"] == t.name
        assert s["function"]["parameters"]["type"] == "object"
        json.dumps(s, ensure_ascii=False)


def test_ctx_state_and_decisions(tmp_path):
    ctx = _mk_ctx(tmp_path)
    assert (ctx.project_dir / "state.json").exists()
    ctx.save_state(target_prep={"pocket": {"center": [1, 2, 3]}})
    ctx2 = Ctx(ctx.project_dir, None, ctx.options, auto=True)
    assert ctx2.state()["target_prep"]["pocket"]["center"] == [1, 2, 3]
    ctx.record_decision("md", "amber19sb", "test")
    log = json.loads((ctx.project_dir / "decisions.json").read_text())
    assert log[-1]["node"] == "md"
    assert log[-1]["answer"] == "amber19sb"


def test_fs_tools(tmp_path):
    ctx = _mk_ctx(tmp_path)
    tools = {t.name: t for t in build_tools()}
    r = tools["write_file"].call(ctx, path="a/b.txt", content="line1\nline2")
    assert r["ok"]
    r = tools["read_file"].call(ctx, path="a/b.txt")
    assert "line2" in r["content"]
    r = tools["edit_file"].call(ctx, path="a/b.txt", old="line1", new="L1")
    assert r["ok"]
    r = tools["read_file"].call(ctx, path="a/b.txt")
    assert "L1" in r["content"]
    r = tools["list_dir"].call(ctx, path="a")
    assert any(e["name"] == "b.txt" for e in r["entries"])
    r = tools["edit_file"].call(ctx, path="a/b.txt", old="nope", new="x")
    assert not r["ok"]


def test_checkpoint_auto_and_finish(tmp_path):
    ctx = _mk_ctx(tmp_path, auto=True)
    tools = {t.name: t for t in build_tools()}
    r = tools["checkpoint"].call(ctx, stage="target", summary="summary")
    assert r["ok"] and r["auto"] and r["approved"]
    log = json.loads((ctx.project_dir / "decisions.json").read_text())
    assert any(d["node"] == "checkpoint_target" for d in log)
    r = tools["checkpoint"].call(ctx, stage="bogus", summary="x")
    assert not r["ok"]
    r = tools["finish"].call(ctx, status="success", summary="done")
    assert r["ok"] and r["status"] == "success"


def test_loop_runs_and_finishes(tmp_path):
    ctx = _mk_ctx(tmp_path)
    brain = FakeBrain([
        # step 1: read a file we write first
        {"content": "先看状态", "tool_calls": [
            {"id": "c1", "name": "state_get", "args": {}}]},
        # step 2: finish
        _finish_call(summary="test done"),
    ])
    ctx.brain = brain
    loop = AgentLoop(ctx, build_tools(), system="S", goal="G", max_steps=10)
    result = loop.run()
    assert result["status"] == "success"
    assert brain.calls == 2
    # transcript has assistant + tool + assistant
    lines = (ctx.project_dir / "agent" / "transcript.jsonl").read_text().splitlines()
    roles = [json.loads(l)["role"] for l in lines]
    assert roles == ["assistant", "tool", "assistant", "tool"]
    assert ctx.state()["status"] == "success"


def test_loop_nudges_then_needs_human(tmp_path):
    ctx = _mk_ctx(tmp_path)
    brain = FakeBrain([
        {"content": "让我想想...", "tool_calls": None},
        {"content": "还在想", "tool_calls": None},
        {"content": "想不出来", "tool_calls": None},
    ])
    ctx.brain = brain
    loop = AgentLoop(ctx, build_tools(), system="S", goal="G", max_steps=10)
    result = loop.run()
    assert result["status"] == "needs_human"
    assert ctx.state()["status"] == "needs_human"


def test_loop_budget(tmp_path):
    ctx = _mk_ctx(tmp_path)
    brain = FakeBrain([
        {"content": "", "tool_calls": [
            {"id": f"c{i}", "name": "state_get", "args": {}}]}
        for i in range(10)
    ])
    ctx.brain = brain
    loop = AgentLoop(ctx, build_tools(), system="S", goal="G", max_steps=3)
    result = loop.run()
    assert result["status"] == "needs_human"
    assert "预算" in result["summary"]


def test_loop_unknown_tool(tmp_path):
    ctx = _mk_ctx(tmp_path)
    brain = FakeBrain([
        {"content": "", "tool_calls": [
            {"id": "c1", "name": "no_such_tool", "args": {}}]},
        _finish_call(),
    ])
    ctx.brain = brain
    loop = AgentLoop(ctx, build_tools(), system="S", goal="G", max_steps=10)
    result = loop.run()
    assert result["status"] == "success"
    # the tool error was surfaced
    lines = (ctx.project_dir / "agent" / "transcript.jsonl").read_text().splitlines()
    tool_msg = [l for l in lines if json.loads(l).get("role") == "tool"]
    assert "未知工具" in tool_msg[0]


def test_loop_resume_replays_transcript(tmp_path):
    ctx = _mk_ctx(tmp_path)
    brain = FakeBrain([
        {"content": "step1", "tool_calls": [
            {"id": "c1", "name": "record_decision",
             "args": {"stage": "md", "choice": "ff=x", "rationale": "t"}}]},
        _finish_call(),
    ])
    ctx.brain = brain
    loop = AgentLoop(ctx, build_tools(), system="S", goal="G", max_steps=10)
    loop.run()
    # now resume with a fresh loop object; it must rebuild messages
    brain2 = FakeBrain([_finish_call()])
    ctx2 = Ctx(ctx.project_dir, brain2, ctx.options, auto=True)
    loop2 = AgentLoop(ctx2, build_tools(), system="S", goal="G", max_steps=10)
    assert loop2.step >= 2                      # replayed steps counted
    n_assistant = sum(1 for m in loop2.messages
                      if m.get("role") == "assistant" and m.get("tool_calls"))
    assert n_assistant == 2                     # both replayed tool turns
    assert len(loop2.messages) == 6             # sys+goal+2*(assistant+tool)
    result = loop2.run()
    assert result["status"] == "success"
    assert brain2.calls == 1


def test_scripted_run(tmp_path, monkeypatch):
    """no-LLM mode: stage tools run in order without an LLM."""
    ctx = _mk_ctx(tmp_path)
    tools = {t.name: t for t in build_tools()}
    calls = []

    def fake(tool_name):
        def fn(ctx_, **kw):
            calls.append(tool_name)
            return {"ok": True}
        return fn

    for name in ("run_target_prep", "run_screening", "run_design",
                 "run_vhh", "run_md", "build_report",
                 "checkpoint", "record_decision"):
        tools[name].fn = fake(name)
    from drugagent.agent.loop import scripted_run
    ctx.options["modules"] = ["screen", "binder", "vhh", "md"]
    result = scripted_run(ctx, list(tools.values()))
    assert result["status"] == "success"
    assert "run_target_prep" in calls
    assert "build_report" in calls
    assert ctx.state()["status"] == "success"
    t = ctx.project_dir / "agent" / "transcript.jsonl"
    assert t.exists()


def test_gmx_env_and_mdp_template(tmp_path):
    ctx = _mk_ctx(tmp_path)
    tools = {t.name: t for t in build_tools()}
    r = tools["gmx_env"].call(ctx)
    assert r["ok"], r
    assert r["version"] in ("2023.1", "4.6.5")
    assert r["detected_ff"] in r["available_force_fields"]
    r = tools["mdp_template"].call(ctx, name="md", ns=5.0)
    assert r["ok"] and "nsteps" in r["text"]
    # R8: production barostat is C-rescale (Berendsen segfaulted the
    # 1HVI smoke system ~20 ps into production)
    assert "C-rescale" in r["text"]
    assert "2500000" in r["text"]            # 5 ns / 2 fs


def test_prompt_flexible_target_workflow():
    """R10/G4: the system prompt must steer the agent to the R5 flexible
    target workflow (conformational selection + flexible docking) and name
    the tools involved."""
    from drugagent.agent import system_prompt
    p = system_prompt({"fast": False, "modules": ["screen", "md"]})
    assert "柔性靶点" in p
    assert "dock_conformer_set" in p
    assert "make_flex_receptor" in p
    assert "consensus" in p
