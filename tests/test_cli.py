import pytest
from typer.testing import CliRunner

from drugagent.cli import _parse_target, app

runner = CliRunner()


def test_parse_target_pdb_file(hivp_pdb):
    t = _parse_target(str(hivp_pdb))
    assert t["kind"] == "pdb_file"


def test_parse_target_id():
    t = _parse_target("1hvi")
    assert t["kind"] == "pdb_id"
    assert t["value"] == "1HVI"


def test_parse_target_fasta(tmp_path):
    p = tmp_path / "x.fasta"
    p.write_text(">seq\nACDEFGHIKLMNPQRSTVWYACDEFGHIK\n")
    t = _parse_target(str(p))
    assert t["kind"] == "fasta"


def test_parse_target_sequence():
    seq = "ACDEFGHIKLMNPQRSTVWY" * 5
    t = _parse_target(seq)
    assert t["kind"] == "sequence"


def test_parse_target_bad():
    with pytest.raises(Exception):
        _parse_target("hello world")


def test_status_no_project():
    res = runner.invoke(app, ["status"])
    assert res.exit_code == 0


def test_status_stage_detail_and_errors(tmp_path):
    """R11: status shows per-stage completion + key numbers, artifacts and
    recent tool failures (not just ✓/✗)."""
    import json
    pdir = tmp_path / "proj"
    (pdir / "agent").mkdir(parents=True)
    (pdir / "02_screening").mkdir(parents=True)
    state = {
        "project_dir": str(pdir),
        "status": "running",
        "options": {},
        "target_prep": {
            "raw_pdb": "r.pdb", "clean_pdb": "c.pdb",
            "receptor_pdbqt": "rec.pdbqt",
            "pocket": {"method": "ligand_centroid(A77)",
                       "center": [0, 0, 0], "xsize": 12,
                       "ysize": 12, "zsize": 12},
            "ligand_resnames": ["A77"],
        },
        "screening": {
            "n_docked": 42, "hit_decision": {"n_hits": 2, "threshold": -8.0},
            "hits": [{"rank": 1, "smiles": "c1ccccc1", "final_score": -9.0},
                     {"rank": 2, "smiles": "CCO", "final_score": -8.5}],
            "library": "chembl35_small (fallback for dtp)",
        },
    }
    (pdir / "state.json").write_text(json.dumps(state, ensure_ascii=False))
    (pdir / "02_screening" / "screening.json").write_text("{}")
    # transcript with one failed tool call
    (pdir / "agent" / "transcript.jsonl").write_text(
        json.dumps({"step": 2, "role": "tool", "name": "dock_screen",
                    "content": json.dumps({"ok": False,
                                           "error": "vina exploded"})}) + "\n"
        + json.dumps({"step": 3, "role": "tool", "name": "checkpoint",
                      "content": json.dumps({"ok": True})}) + "\n")
    res = runner.invoke(app, ["status", "--project", str(pdir)])
    assert res.exit_code == 0
    out = res.output
    assert "[x] target_prep" in out and "[x] screening" in out
    assert "[ ] binder" in out and "[ ] vhh" in out and "[ ] md" in out
    assert "[ ] report" in out
    assert "42 对接, 2 命中" in out
    assert "fallback for dtp" in out
    assert "02_screening/screening.json" in out  # artifacts listed
    assert "最近错误" in out and "vina exploded" in out
    assert "dock_screen" in out


# ---------------------------------------------------------------- R9


def test_resolve_defaults_options_override():
    """R9: options layer on top of resolved Defaults — known fields win,
    unknown keys are ignored, fast still applies first."""
    from drugagent.config import resolve_defaults
    d = resolve_defaults({"fast": True,
                          "dock_exhaustiveness_final": 32,
                          "n_hits": 3,
                          "name": "x", "max_steps": 50, "auto": True})
    assert d.dock_exhaustiveness_final == 32
    assert d.n_hits == 3
    assert d.md_ns == 5.0  # fast value survived
    d2 = resolve_defaults({"dock_md_rep_exhaustiveness": 16})
    assert d2.dock_md_rep_exhaustiveness == 16
    assert d2.dock_exhaustiveness_final == 16  # untouched
    d3 = resolve_defaults(None)
    assert d3.flex_cutoff_ang == 5.0


def test_resolve_defaults_none_values_ignored():
    from drugagent.config import resolve_defaults
    d = resolve_defaults({"n_hits": None, "fast": False})
    assert d.n_hits == 20


def test_run_exposes_md_knobs():
    """R10/G5: the key MD knobs are user-facing CLI options."""
    import os
    os.environ["COLUMNS"] = "200"  # keep rich from eliding long option names
    res = runner.invoke(app, ["run", "--help"])
    assert res.exit_code == 0
    for opt in ("--md-salt", "--md-divalent", "--md-divalent-m",
                "--md-extend-ns", "--md-max-extensions", "--md-burn-in-ps"):
        assert opt in res.output, opt


def test_run_exposes_vhh_knobs():
    """R11: the VHH knobs are user-facing CLI options."""
    import os
    os.environ["COLUMNS"] = "200"
    res = runner.invoke(app, ["run", "--help"])
    assert res.exit_code == 0
    assert "--vhh-plddt-min" in res.output
    assert "--vhh-dock-flex" in res.output
    assert "--vhh-dock-cdr-only" in res.output
