"""R11/G10 benchmark: rigid vs flexible full-VHH docking (vhh_30, 1HVI).

Usage: env/bin/python scripts/bench_vhh_dock.py
Writes timing to logs/bench_vhh_dock.json (and stdout).
"""
import json
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
import sys
sys.path.insert(0, str(ROOT))

from drugagent.config import TOOLS  # noqa: E402
from drugagent.modules.target_prep import make_rigid_pdbqt, to_pdbqt  # noqa: E402
from drugagent.utils import run_cmd  # noqa: E402

MODEL = ROOT / "projects/r10_e2e/04_vhh/vhh_models/vhh_30.pdb"
REC = ROOT / "projects/r10_e2e/01_target/receptor.pdbqt"
OUT = ROOT / "data/fixtures/bench_vhh_dock"
OUT.mkdir(parents=True, exist_ok=True)
POCKET = json.loads((ROOT / "projects/r10_e2e/state.json").read_text()) \
    ["target_prep"]["pocket"]
VINA = TOOLS / "vina" / "bin" / "vina"

n_atoms = sum(1 for l in open(MODEL) if l.startswith("ATOM"))
print(f"model: {MODEL.name} atoms={n_atoms}")

results = {}
for mode in ("rigid", "flex"):
    pdbqt = OUT / f"vhh_30_{mode}.pdbqt"
    if not pdbqt.is_file():
        t0 = time.time()
        to_pdbqt(MODEL, pdbqt)  # graph + element fixes (required by this build)
        if mode == "rigid":
            make_rigid_pdbqt(pdbqt)  # TORSDOF 0 -> rigid body
        print(f"{mode}: pdbqt in {time.time() - t0:.1f}s")
    prefix = OUT / f"vhh_30_{mode}_dock"
    cmd = [str(VINA), "--receptor", str(REC), "--ligand", str(pdbqt),
           "--center_x", f"{POCKET['center'][0]:.4f}",
           "--center_y", f"{POCKET['center'][1]:.4f}",
           "--center_z", f"{POCKET['center'][2]:.4f}",
           "--size_x", f"{POCKET['xsize']:.2f}",
           "--size_y", f"{POCKET['ysize']:.2f}",
           "--size_z", f"{POCKET['zsize']:.2f}",
           "--exhaustiveness", "1", "--cpu", "8", "--seed", "42",
           "--out", str(prefix)]
    t0 = time.time()
    run_cmd(cmd, log_file=OUT / f"vina_{mode}.log", check=True)
    dt = time.time() - t0
    score = None
    with open(prefix) as fh:
        for line in fh:
            if line.startswith("REMARK") and "VINA RESULT" in line:
                score = float(line.split()[-1])
                break
    results[mode] = {"seconds": round(dt, 1), "min": round(dt / 60, 2),
                     "best_score": score}
    print(f"{mode}: {dt / 60:.1f} min, best score={score}")

ratio = (results["flex"]["seconds"] / results["rigid"]["seconds"]
         if results["rigid"]["seconds"] else None)
print(f"flex/rigid ratio: {ratio:.1f}x" if ratio else "n/a")
(OUT / "bench.json").write_text(json.dumps(
    {"n_atoms": n_atoms, "results": results,
     "flex_over_rigid": round(ratio, 1) if ratio else None}, indent=2))
print("wrote", OUT / "bench.json")
