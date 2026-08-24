import pytest

from drugagent.graph import build_graph, compile_graph


def test_graph_builds():
    g = build_graph()
    compiled = g.compile()
    nodes = set(compiled.get_graph().nodes)
    for expected in ("init", "target_prep", "cp_target", "screening",
                     "cp_screen", "binder", "vhh", "cp_design", "md",
                     "report"):
        assert expected in nodes


def test_graph_mini_run_auto(tmp_path, monkeypatch, hivp_pdb):
    """Full graph with stubbed heavy nodes; target_prep is real but light."""
    from drugagent.modules import screening, binder, vhh, mdsim

    calls = []

    def stub(name, out_key):
        def _node(state):
            calls.append(name)
            state_out = dict(state)
            state_out[out_key] = {"stub": name}
            return state_out
        return _node

    monkeypatch.setattr(screening, "screen", stub("screen", "screening"))
    monkeypatch.setattr(binder, "design_binder", stub("binder", "binder"))
    monkeypatch.setattr(vhh, "design_vhh_all", stub("vhh", "vhh"))
    monkeypatch.setattr(mdsim, "run_md", stub("md", "md"))

    app, conn = compile_graph(tmp_path / "proj")
    try:
        state = {
            "project_dir": str(tmp_path / "proj"),
            "target": {"kind": "pdb_file", "value": str(hivp_pdb)},
            "options": {"modules": ["screen", "binder", "vhh", "md"],
                        "auto": True, "no_llm": True, "fast": True,
                        "n_jobs": 2},
        }
        final = app.invoke(state, {"configurable": {"thread_id": "t1"}})
        assert "target_prep" in final
        assert final["target_prep"]["pocket"]["ligand"] == "A77"
        for c in ("screen", "binder", "vhh", "md"):
            assert c in calls
        assert "report" in final
    finally:
        conn.close()


def test_graph_screen_only(tmp_path, monkeypatch, hivp_pdb):
    from drugagent.modules import screening, binder, vhh, mdsim

    calls = []
    monkeypatch.setattr(screening, "screen",
                        lambda s: calls.append("screen") or {**s, "screening": {}})
    monkeypatch.setattr(binder, "design_binder",
                        lambda s: calls.append("binder") or {**s, "binder": {}})
    monkeypatch.setattr(vhh, "design_vhh_all",
                        lambda s: calls.append("vhh") or {**s, "vhh": {}})
    monkeypatch.setattr(mdsim, "run_md",
                        lambda s: calls.append("md") or {**s, "md": {}})
    app, conn = compile_graph(tmp_path / "proj2")
    try:
        state = {
            "project_dir": str(tmp_path / "proj2"),
            "target": {"kind": "pdb_file", "value": str(hivp_pdb)},
            "options": {"modules": ["screen"], "auto": True,
                        "no_llm": True, "n_jobs": 2},
        }
        app.invoke(state, {"configurable": {"thread_id": "t1"}})
    finally:
        conn.close()
    assert calls == ["screen"]


def test_checkpoint_interrupts(tmp_path, monkeypatch, hivp_pdb):
    """Without --auto, the graph must pause at the first checkpoint."""
    from drugagent.modules import screening
    monkeypatch.setattr(screening, "screen",
                        lambda s: {**s, "screening": {}})
    app, conn = compile_graph(tmp_path / "proj3")
    try:
        state = {
            "project_dir": str(tmp_path / "proj3"),
            "target": {"kind": "pdb_file", "value": str(hivp_pdb)},
            "options": {"modules": ["screen"], "auto": False,
                        "no_llm": True, "n_jobs": 2},
        }
        app.invoke(state, {"configurable": {"thread_id": "t1"}})
        snap = app.get_state({"configurable": {"thread_id": "t1"}})
        assert snap.next  # paused
        assert snap.tasks[0].interrupts[0].value["checkpoint"] == "target"
        # resume
        from langgraph.types import Command
        app.invoke(Command(resume={"action": "approve"}),
                   {"configurable": {"thread_id": "t1"}})
        snap = app.get_state({"configurable": {"thread_id": "t1"}})
        assert "screening" in snap.values
    finally:
        conn.close()
