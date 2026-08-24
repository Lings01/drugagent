"""R11/G8: stage-level idempotency (run_* reuse complete state sections)."""
import json
from pathlib import Path

import pytest

from drugagent.agent.loop import Ctx
from drugagent.agent import stages as stg
from drugagent.agent.tools_screen import run_screening
from drugagent.agent.tools_vhh import run_vhh
from drugagent.agent.tools_design import run_design
from drugagent.agent.tools_report import build_report


def _ctx(tmp_path: Path, state: dict) -> Ctx:
    """Ctx whose state.json is exactly `state` (project files included)."""
    pdir = tmp_path / "proj"
    pdir.mkdir(parents=True, exist_ok=True)
    (pdir / "state.json").write_text(json.dumps(state, ensure_ascii=False))
    ctx = Ctx(pdir, brain=None, options=state.get("options", {}))
    return ctx


def _complete_state(tmp_path: Path) -> dict:
    """All stage sections carry their completion markers."""
    report_html = tmp_path / "proj" / "reports" / "report.html"
    report_html.parent.mkdir(parents=True, exist_ok=True)
    report_html.write_text("<html></html>")
    return {
        "status": "running",
        "project_dir": str(tmp_path / "proj"),
        "options": {"no_llm": True},
        "target_prep": {
            "raw_pdb": "x.pdb", "clean_pdb": "x_clean.pdb",
            "receptor_pdbqt": "rec.pdbqt",
            "pocket": {"center": [0, 0, 0], "xsize": 12,
                       "ysize": 12, "zsize": 12},
        },
        "screening": {
            "n_docked": 10, "hit_decision": {"n_hits": 2, "threshold": -8.0},
            "hits": [{"rank": 1, "smiles": "c1ccccc1", "final_score": -9.0},
                     {"rank": 2, "smiles": "CCO", "final_score": -8.5}],
            "library": "chembl35_small (fallback for dtp)",
        },
        "binder": {
            "n_designs": 2, "binder_type": "miniprotein",
            "designs": [{"design": "d0.pdb", "interface_plddt_mean": 80.0}],
            "best": {"design": "d0.pdb", "interface_plddt_mean": 80.0},
        },
        "vhh": {
            "track_a": {"results": [{"idx": 0, "score": -5.0}]},
            "track_b": {"designs": [{"design": "v0.pdb"}]},
            "best": {"source": "screening", "idx": 0},
        },
        "md": {
            "ns": 5.0, "reps": 3,
            "summary": {"final_rmsd_mean": 0.2, "final_rg_mean": 1.2},
            "replicas": [{"rep": 1, "final_step": 2500}],
        },
        "report": {"html": str(report_html), "pdf": str(report_html) + "f"},
    }


# --------------------------------------------------------------------------- #
# completion markers
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("stage", ["target_prep", "screening", "binder",
                                   "vhh", "md", "report"])
def test_empty_state_incomplete(tmp_path, stage):
    ok, why = stg.stage_complete({}, stage, tmp_path)
    assert not ok
    assert why


def test_complete_state_all_stages(tmp_path):
    state = _complete_state(tmp_path)
    for stage in ("target_prep", "screening", "binder", "vhh", "md", "report"):
        ok, why = stg.stage_complete(state, stage, tmp_path)
        assert ok, (stage, why)


@pytest.mark.parametrize("stage,mutate", [
    ("target_prep", lambda s: s["target_prep"].pop("pocket")),
    ("screening", lambda s: s["screening"].pop("hit_decision")),
    ("binder", lambda s: s["binder"].update(best=None)),
    ("vhh", lambda s: s["vhh"].pop("track_b")),
    ("md", lambda s: s["md"].pop("summary")),
    ("report", lambda s: None),  # file existence checked separately
])
def test_partial_state_incomplete(tmp_path, stage, mutate):
    state = _complete_state(tmp_path)
    if stage == "report":
        Path(state["report"]["html"]).unlink()
    else:
        mutate(state)
    ok, _ = stg.stage_complete(state, stage, tmp_path)
    assert not ok


# --------------------------------------------------------------------------- #
# maybe_reuse gate
# --------------------------------------------------------------------------- #
def test_maybe_reuse_returns_cached_summary(tmp_path, monkeypatch):
    state = _complete_state(tmp_path)
    ctx = _ctx(tmp_path, state)
    called = {"n": 0}

    def fake_screen(st):
        called["n"] += 1
        return st

    import drugagent.modules.screening as sc
    monkeypatch.setattr(sc, "screen", fake_screen)
    res = run_screening(ctx)
    assert res["reused"] is True
    assert called["n"] == 0
    assert res["n_hits"] == 2
    assert res["top3"][0]["smiles"] == "c1ccccc1"


def test_maybe_reuse_force_reruns(tmp_path, monkeypatch):
    state = _complete_state(tmp_path)
    ctx = _ctx(tmp_path, state)
    called = {"n": 0}

    def fake_screen(st):
        called["n"] += 1
        return dict(st, screening=dict(
            st["screening"], hit_decision={"n_hits": 1, "threshold": -7.0},
            hits=[{"rank": 1, "smiles": "CCO", "final_score": -7.5}]))

    import drugagent.modules.screening as sc
    monkeypatch.setattr(sc, "screen", fake_screen)
    res = run_screening(ctx, force=True)
    assert called["n"] == 1
    assert "reused" not in res
    assert res["n_hits"] == 1


def test_maybe_reuse_incomplete_runs(tmp_path, monkeypatch):
    state = _complete_state(tmp_path)
    state["screening"].pop("hit_decision")  # incomplete
    ctx = _ctx(tmp_path, state)
    called = {"n": 0}

    def fake_screen(st):
        called["n"] += 1
        return dict(st, screening={"n_docked": 1,
                                   "hit_decision": {"n_hits": 1},
                                   "hits": [{"rank": 1, "smiles": "c",
                                             "final_score": -1.0}]})

    import drugagent.modules.screening as sc
    monkeypatch.setattr(sc, "screen", fake_screen)
    res = run_screening(ctx)
    assert called["n"] == 1
    assert res["n_hits"] == 1


def test_run_vhh_reuse(tmp_path, monkeypatch):
    state = _complete_state(tmp_path)
    ctx = _ctx(tmp_path, state)
    import drugagent.modules.vhh as vh
    called = {"n": 0}

    def fake_all(st):
        called["n"] += 1
        return st

    monkeypatch.setattr(vh, "design_vhh_all", fake_all)
    res = run_vhh(ctx)
    assert res["reused"] is True and called["n"] == 0
    assert res["n_track_a"] == 1 and res["n_track_b"] == 1


def test_run_design_reuse(tmp_path, monkeypatch):
    state = _complete_state(tmp_path)
    ctx = _ctx(tmp_path, state)
    import drugagent.modules.binder as bd
    called = {"n": 0}

    def fake_all(st):
        called["n"] += 1
        return st

    monkeypatch.setattr(bd, "design_binder", fake_all)
    res = run_design(ctx)
    assert res["reused"] is True and called["n"] == 0
    assert res["n_designs"] == 2


def test_build_report_reuse(tmp_path):
    state = _complete_state(tmp_path)
    ctx = _ctx(tmp_path, state)
    res = build_report(ctx)
    assert res["reused"] is True
    assert Path(res["html"]).is_file()
    res2 = build_report(ctx, force=True)
    assert "reused" not in res2  # actually regenerated


# --------------------------------------------------------------------------- #
# scripted run skips finished stages end-to-end (no-llm)
# --------------------------------------------------------------------------- #
def test_scripted_run_skips_completed_stages(tmp_path, monkeypatch):
    from drugagent.agent.loop import scripted_run
    from drugagent.agent import build_tools

    state = _complete_state(tmp_path)
    # drop report so build_report still runs once (it's cheap)
    state.pop("report")
    (tmp_path / "proj" / "state.json").write_text(
        json.dumps(state, ensure_ascii=False))
    ctx = Ctx(tmp_path / "proj", brain=None,
              options={"modules": ["screen", "binder", "vhh", "md"],
                       "fast": True, "no_llm": True},
              auto=True)
    tools = build_tools()

    calls = {"target": 0, "screen": 0, "binder": 0, "vhh": 0, "md": 0}

    def _mk(key):
        def f(st):
            calls[key] += 1
            return st
        return f

    import drugagent.modules.target_prep as tp
    import drugagent.modules.screening as sc
    import drugagent.modules.binder as bd
    import drugagent.modules.vhh as vh
    import drugagent.modules.mdsim as mds
    monkeypatch.setattr(tp, "prepare_target", _mk("target"))
    monkeypatch.setattr(sc, "screen", _mk("screen"))
    monkeypatch.setattr(bd, "design_binder", _mk("binder"))
    monkeypatch.setattr(vh, "design_vhh_all", _mk("vhh"))
    monkeypatch.setattr(mds, "run_md", _mk("md"))

    result = scripted_run(ctx, tools)
    assert result["status"] == "success"
    # every finished stage was reused, none re-executed
    assert calls == {"target": 0, "screen": 0, "binder": 0,
                     "vhh": 0, "md": 0}, calls
    # report was (re)generated from the cached stages
    assert (tmp_path / "proj" / "reports" / "report.html").is_file()
