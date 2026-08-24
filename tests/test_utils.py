import pytest
from drugagent.utils import (centroid_from_pdb, is_pdb_id, pdb_chains,
                             pdb_extent)


def test_is_pdb_id():
    assert is_pdb_id("1HVI")
    assert is_pdb_id("1hvi")
    assert not is_pdb_id("1HV")
    assert not is_pdb_id("/home/a.pdb")
    assert not is_pdb_id("1HVIA")


def test_pdb_chains(hivp_pdb):
    ch = pdb_chains(hivp_pdb)
    assert set(ch) == {"A", "B"}
    assert ch["A"] > 800


def test_centroid_and_extent(hivp_pdb):
    c = centroid_from_pdb(hivp_pdb, resnames=["A77"])
    assert all(abs(v) < 50 for v in c)
    ext = pdb_extent(hivp_pdb, c, resnames=["A77"])
    assert 3 < ext < 25


# --------------------------------------------------------------------------- #
# R10: error surfacing — a failed *logged* command must carry its log tail in
# the exception, otherwise the agent/user cannot see WHY (e.g. gmx solvate
# "File does not exist") without a second read_file round-trip
# --------------------------------------------------------------------------- #
def test_run_cmd_failure_includes_log_tail(tmp_path):
    from drugagent.utils import run_cmd
    log = tmp_path / "cmd.log"
    with pytest.raises(RuntimeError) as ei:
        run_cmd(["/bin/sh", "-c",
                 "echo line1; echo line2; echo FATAL: boom >&2; exit 3"],
                log_file=log, check=True)
    msg = str(ei.value)
    assert "FATAL: boom" in msg          # stderr is in the log -> in message
    assert "line2" in msg                # stdout lines too
    assert "exit" not in msg[:0]         # (no-op, keep simple)
    # and the log path is mentioned so the full context is findable
    assert "cmd.log" in msg or str(log) in msg


def test_run_cmd_success_not_affected(tmp_path):
    from drugagent.utils import run_cmd
    log = tmp_path / "ok.log"
    p = run_cmd(["/bin/sh", "-c", "echo hello"], log_file=log, check=True)
    assert p.returncode == 0


def test_run_cmd_failure_no_log_still_raises(tmp_path):
    from drugagent.utils import run_cmd
    with pytest.raises(RuntimeError):
        run_cmd(["/bin/false"], check=True)
