#!/usr/bin/env python
"""R8 e2e: pool-dock the smoke project (crystal + MD cluster reps,
flexible side chains) after the v5 smoke MD finishes.

Usage: PYTHONNOUSERSITE=1 env/bin/python scripts/r8_pool_dock.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from drugagent.agent.loop import Ctx  # noqa: E402
from drugagent.agent import tools_screen  # noqa: E402

PROJ = ROOT / "projects" / "agent_smoke_0821_0404"


def main() -> int:
    ctx = Ctx(PROJ, brain=None, options={
        "dock_exhaustiveness_final": 32,
        "n_hits": 3,
    })
    out = tools_screen.dock_conformer_set(ctx, max_conformers=3,
                                          min_pop=0.05)
    print(json.dumps(out, indent=2, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
