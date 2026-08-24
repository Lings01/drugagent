"""De novo binder design tools."""
from __future__ import annotations

import re
from pathlib import Path

from ..config import DEFAULTS
from ..llm import AgentBrain
from ..modules import binder as bd
from ..modules.esmfold_run import interface_metrics, predict
from ..utils import jsave
from .loop import Ctx, Tool, ToolError
from .stages import maybe_reuse, summarize_binder

AA = "ACDEFGHIKLMNPQRSTVWY"


def _rf_final_seq(log_path: Path) -> str | None:
    """The designed sequence only appears in the RF denoising log
    ('input to next step' lines); the last such line is the final design."""
    if not log_path.is_file():
        return None
    lines = log_path.read_text(errors="ignore").splitlines()
    for line in reversed(lines):
        if "input to next step" in line:
            m = re.search(r"[" + AA + r"]{15,}", line)
            if m:
                return m.group(0)
    return None


def _prep(ctx: Ctx) -> dict:
    prep = ctx.state().get("target_prep") or {}
    if not prep.get("clean_pdb"):
        raise ToolError("state.target_prep 未就绪 (先做靶点准备)")
    return prep


# --------------------------------------------------------------------------- #
def rf_design(ctx: Ctx, target_pdb: str | None = None,
              n_designs: int = 2, length_min: int = 60, length_max: int = 80,
              hotspots: list[str] | None = None) -> dict:
    """RFdiffusion inpaint binder design (backbone only).

    NOTE: design PDBs carry GLY residue names; the real sequence is extracted
    from the RF denoising log (returned as `seqs`).
    """
    from ..config import resolve_defaults
    d = resolve_defaults(ctx.options)
    prep = _prep(ctx)
    workdir = ctx.stage_dir("03_binder")
    target = Path(target_pdb) if target_pdb else Path(prep["clean_pdb"])
    if not target.exists():
        raise ToolError(f"target pdb not found: {target}")
    pocket = prep.get("pocket") or {}
    hotspots = hotspots or bd.pocket_hotspots(target, pocket)
    designs = bd.rfdesign(target, pocket, workdir,
                          n_designs=n_designs,
                          length=(length_min, length_max),
                          hotspots=hotspots)
    if not designs:
        raise ToolError("RFdiffusion produced no designs (看 03_binder/rf_design.log)")
    info = []
    for pdb in designs:
        seq = _rf_final_seq(workdir / "rf_design.log")
        info.append({"design": str(pdb), "seq": seq, "seqs": [seq] if seq else []})
    st = ctx.state()
    binder = st.get("binder") or {}
    binder.update({
        "n_designs": len(designs),
        "length_range": [length_min, length_max],
        "hotspots": hotspots,
        "designs": info,
        "best": None,
    })
    st["binder"] = binder
    ctx.save_state()
    return {"n_designs": len(designs),
            "designs": [{"pdb": i["design"], "seq": i["seq"]} for i in info],
            "log": str(workdir / "rf_design.log"),
            "note": "PDB 残基名为 GLY (RF 输出特性); 真实序列在 seqs/log 里"}


def mpnn_seq(ctx: Ctx, design_pdb: str, n_seq: int = 5) -> dict:
    """MPNN sequence design for a designed backbone (fallback if env missing)."""
    p = Path(design_pdb)
    if not p.is_file():
        raise ToolError(f"not a file: {p}")
    workdir = ctx.stage_dir("03_binder")
    seqs = bd.mpnn_sequence(p, workdir, n_seq=n_seq)
    return {"design": str(p), "seqs": seqs}


def make_complex(ctx: Ctx, target_pdb: str, design_pdb: str,
                 out: str | None = None) -> dict:
    a = Path(target_pdb)
    b = Path(design_pdb)
    if not a.is_file() or not b.is_file():
        raise ToolError(f"missing file: {a if not a.is_file() else b}")
    prep = _prep(ctx)
    out_p = (Path(out) if out else
             ctx.stage_dir("03_binder") / "complex" / f"{b.stem}_complex.pdb")
    bd._make_complex(a, b, out_p)
    return {"complex_pdb": str(out_p)}


def esmfold_score(ctx: Ctx, complex_pdb: str, num_recycles: int = 3) -> dict:
    """ESMFold-refold a complex PDB; returns pLDDT + interface metrics.

    NOTE: this scores ESMFold's re-folding confidence, not how well RF
    reproduced the designed conformation.
    """
    p = Path(complex_pdb)
    if not p.is_file():
        raise ToolError(f"not a file: {p}")
    seqs = {}
    # per-chain CA sequences
    from ..modules.binder import _ca_sequence
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("x", str(p))
    chains = [c.id for c in struct[0]]
    for ch in chains:
        s = _ca_sequence(p, ch)
        if s:
            seqs[ch] = s
    if len(seqs) < 2:
        raise ToolError(f"complex needs >=2 protein chains, got {list(seqs)}")
    out = predict([seqs[c] for c in sorted(seqs)], num_recycles=num_recycles,
                  device="cpu")
    rp = out.get("res_present")
    im = interface_metrics(p, out["plddt"][rp] if rp is not None else out["plddt"])
    res = {
        "complex_pdb": str(p),
        "plddt_mean": round(float(out["mean_plddt"]), 2),
        "plddt_min": round(float(out["min_plddt"]), 2),
        "interface_plddt_mean": im.get("interface_plddt_mean"),
        "interface_plddt_min": im.get("interface_plddt_min"),
        "n_interface": im.get("n_interface_residues"),
    }
    # write the ESMFold-refined complex next to the input
    refined = p.with_name(p.stem + "_esmfold.pdb")
    from ..modules.binder import write_pdb_file
    write_pdb_file(refined, out)
    res["refined_pdb"] = str(refined)
    return res


def score_designs(ctx: Ctx) -> dict:
    """ESMFold-score every design in state.binder (mono + complex interface)."""
    st = ctx.state()
    binder = st.get("binder") or {}
    prep = _prep(ctx)
    designs = [Path(d["design"]) for d in binder.get("designs", [])
               if Path(d["design"]).is_file()]
    if not designs:
        raise ToolError("state.binder.designs 为空 (先 rf_design)")
    seqs = {Path(d["design"]).stem: d.get("seqs", [])
            for d in binder.get("designs", [])}
    workdir = ctx.stage_dir("03_binder")
    scored = bd.score_designs(designs, Path(prep["clean_pdb"]), workdir, seqs,
                              device="cpu")
    binder["designs"] = scored
    binder["best"] = scored[0] if scored else None
    st["binder"] = binder
    jsave(workdir / "binder.json", binder)
    return {"n_scored": len(scored),
            "top3": [{"design": s["design"],
                      "mono_plddt": s.get("mono_plddt"),
                      "interface_plddt_mean": s.get("interface_plddt_mean")}
                     for s in scored[:3]]}


def run_design(ctx: Ctx, force: bool = False) -> dict:
    """Full deterministic binder stage (1.0 pipeline as one tool).

    R11/G8: reuses state.binder when complete (designs + best); force=true
    to re-run."""
    cached = maybe_reuse(ctx, "binder", force)
    if cached is not None:
        return cached
    st = ctx.state()
    state_out = bd.design_binder(st)    # returns full state
    ctx.save_state(binder=state_out["binder"])
    return summarize_binder(ctx.state())


# --------------------------------------------------------------------------- #
def build() -> list[Tool]:
    return [
        Tool("rf_design",
             "RFdiffusion 从头设计 binder 骨架 (inpaint; 注意 PDB 残基名为 GLY, "
             "真实序列从去噪日志提取并返回)。",
             {"type": "object",
              "properties": {"target_pdb": {"type": "string"},
                             "n_designs": {"type": "integer"},
                             "length_min": {"type": "integer"},
                             "length_max": {"type": "integer"},
                             "hotspots": {"type": "array",
                                           "items": {"type": "string"}}},
              "required": []}, rf_design),
        Tool("mpnn_seq",
             "MPNN 序列设计 (给定骨架; 环境缺失时退化为回退序列)。",
             {"type": "object",
              "properties": {"design_pdb": {"type": "string"},
                             "n_seq": {"type": "integer"}},
              "required": ["design_pdb"]}, mpnn_seq),
        Tool("make_complex",
             "拼接 靶点(链A) + 设计物(链B) 复合物 PDB。",
             {"type": "object",
              "properties": {"target_pdb": {"type": "string"},
                             "design_pdb": {"type": "string"},
                             "out": {"type": "string"}},
              "required": ["target_pdb", "design_pdb"]}, make_complex),
        Tool("esmfold_score",
             "ESMFold 复折叠打分: pLDDT + 界面指标 (注意: 评的是 ESMFold "
             "重折叠置信度, 不是 RF 构象复现度)。",
             {"type": "object",
              "properties": {"complex_pdb": {"type": "string"},
                             "num_recycles": {"type": "integer"}},
              "required": ["complex_pdb"]}, esmfold_score),
        Tool("score_designs",
             "对 state.binder 里所有设计做 ESMFold 单体+复合物打分并排序。",
             {"type": "object", "properties": {}, "required": []},
         score_designs),
        Tool("run_design",
             "整段 binder 设计 (确定性标准流程: RF设计→MPNN→ESMFold打分)。"
             "已完成阶段自动复用, force=true 强制重跑。",
             {"type": "object",
              "properties": {"force": {"type": "boolean",
                                       "description": "强制重跑 (默认复用已完成阶段)"}},
              "required": []}, run_design),
    ]
