"""LangGraph state definition."""
from __future__ import annotations

from typing import Any, TypedDict


class AgentState(TypedDict, total=False):
    # inputs
    project_dir: str            # output directory for this run
    target: dict                # {"kind": pdb_file|pdb_id|fasta|sequence, "value": str}
    options: dict               # modules, fast, auto, n_jobs, md overrides, llm overrides
    # module outputs (dicts, see each module for schema)
    target_prep: dict
    screening: dict
    binder: dict
    vhh: dict
    md: dict
    report: dict
    # bookkeeping
    decisions: list[dict]
    errors: list[dict]
    status: str                 # running|paused|done|failed
