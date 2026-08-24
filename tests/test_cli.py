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
    res = runner.invoke(app, ["run", "--help"])
    assert res.exit_code == 0
    for opt in ("--md-salt", "--md-divalent", "--md-divalent-m",
                "--md-extend-ns", "--md-max-extensions", "--md-burn-in-ps"):
        assert opt in res.output, opt
