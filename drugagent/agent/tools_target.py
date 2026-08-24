"""Target preparation tools."""
from __future__ import annotations

from pathlib import Path

from ..modules import target_prep as tp
from ..modules.screening import library_path
from ..utils import pdb_ligands
from .loop import Ctx, Tool, ToolError


# --------------------------------------------------------------------------- #
def resolve_target(ctx: Ctx, kind: str | None = None,
                   value: str | None = None) -> dict:
    target = {"kind": kind, "value": value} if (kind and value) else \
        (ctx.state().get("target") or {})
    if not target.get("kind"):
        raise ToolError("需要 kind + value (或 state.target 已设置)")
    out = tp.resolve_target(target, ctx.stage_dir("01_target"),
                            modeler=ctx.options.get("modeler", "esmfold"))
    st = ctx.state()
    st.setdefault("target_prep", {})
    st["target_prep"]["raw_pdb"] = str(out["pdb_path"])
    st["target_prep"]["resolved"] = {k: (str(v) if isinstance(v, Path) else v)
                                     for k, v in out.items()}
    ctx.save_state()
    return out


def analyze_pdb(ctx: Ctx, path: str | None = None) -> dict:
    p = Path(path) if path else (ctx.state().get("target_prep", {}) or {}).get("raw_pdb")
    if not p or not Path(p).exists():
        raise ToolError(f"PDB not found: {p!r} (先 resolve_target 或给 path)")
    return tp.analyze_completeness(Path(p))


def repair_structure(ctx: Ctx, path: str, actions: list[str],
                     out: str = "01_target/target_repaired.pdb") -> dict:
    src = Path(path)
    dst = _p(ctx, out)
    st = tp.repair_structure(src, dst, actions=actions)
    tp_st = ctx.state().setdefault("target_prep", {})
    tp_st["repaired_pdb"] = str(dst)
    ctx.save_state()
    return st


def clean_pdb(ctx: Ctx, path: str, out: str = "01_target/target_clean.pdb",
              keep_resnames: list[str] | None = None,
              keep_waters: bool = True,
              keep_chain: str | None = None) -> dict:
    return tp.clean_pdb(Path(path), _p(ctx, out),
                        keep_resnames=keep_resnames,
                        keep_waters=keep_waters, keep_chain=keep_chain)


def _p(ctx: Ctx, s: str) -> Path:
    p = Path(s)
    return p if p.is_absolute() else ctx.project_dir / p


def find_pockets(ctx: Ctx, pdb: str | None = None,
                 lig_resname: str | None = None,
                 grid: float = 0.5, pad: float = 5.0,
                 top_n: int = 5) -> dict:
    """Rank docking cavities. Known ligand wins. The result (and candidates)
    are stored into state.target_prep.pocket for downstream tools."""
    st = ctx.state()
    prep = st.get("target_prep") or {}
    p = Path(pdb) if pdb else Path(prep.get("clean_pdb") or prep.get("raw_pdb") or "")
    if not p.exists():
        raise ToolError(f"PDB not found: {p!r}")
    if lig_resname is None:
        lig_resname = (prep.get("ligand_resnames") or [None])[0]

    if lig_resname:
        pocket = tp.pocket_from_ligand(p, lig_resname)
        pocket["rationale"] = f"已知配体 {lig_resname} 所在位点，最可靠的对接位点"
        cands = [pocket]
    else:
        cands = tp.grid_pockets(p, grid=grid, pad=pad)
        if not cands:
            from ..utils import centroid_from_pdb
            c = centroid_from_pdb(p)
            pocket = {"center": list(c), "xsize": 24, "ysize": 24,
                      "zsize": 24, "method": "protein_centroid_fallback",
                      "rationale": "未检测到空腔，用蛋白中心兜底"}
            cands = [pocket]
        else:
            # the agent (main LLM) makes the final pick; here we only record
            # the candidates and preselect the top one as default
            pocket = cands[0]
            pocket["rationale"] = "默认取体积最大空腔 (agent 可改选其它候选)"
    cands_summary = [{
        "rank": i + 1,
        "center": c["center"],
        "volume_A3": c.get("volume_A3"),
        "xsize": c.get("xsize"),
    } for i, c in enumerate(cands[:top_n])]
    prep["pocket"] = pocket
    prep["pocket_candidates"] = cands_summary
    st["target_prep"] = prep
    ctx.save_state()
    return {"pocket": pocket, "candidates": cands_summary,
            "note": "candidates 中可改选: 用 record_decision 记录后, "
                    "通过 state 或直接以 pocket 参数传给对接工具"}


def pdb_to_pdbqt(ctx: Ctx, pdb: str, out: str | None = None,
                 keep_resnames: list[str] | None = None, flex: bool = True) -> dict:
    src = Path(pdb)
    dst = _p(ctx, out) if out else src.with_suffix(".pdbqt")
    tp.to_pdbqt(src, dst, keep_resnames=keep_resnames, flex=flex)
    return {"pdbqt": str(dst)}


def run_target_prep(ctx: Ctx) -> dict:
    """Full deterministic target stage (1.0 pipeline as one tool)."""
    st = ctx.state()
    state_out = tp.prepare_target(st)   # returns full state
    out = state_out["target_prep"]
    ctx.save_state(target_prep=out)
    return {k: out.get(k) for k in
            ("raw_pdb", "clean_pdb", "pocket", "receptor_pdb",
             "receptor_pdbqt", "ligand_pdbqt", "ligand_resnames",
             "judgment", "completeness")}


# --------------------------------------------------------------------------- #
def build() -> list[Tool]:
    return [
        Tool("resolve_target",
             "把靶点输入 (PDB文件/PDB ID/FASTA/裸序列) 解析成 PDB 文件。"
             "裸序列会用 ESMFold 建模 (较慢)。",
             {"type": "object",
              "properties": {"kind": {"type": "string",
                                      "enum": ["pdb_file", "pdb_id", "fasta",
                                               "sequence"]},
                             "value": {"type": "string"}},
              "required": []}, resolve_target),
        Tool("analyze_pdb",
             "结构完整性+坑检测: 链/配体/残基/多MODEL/altloc/缺失残基/金属/核酸/无序末端, 返回 issues[] 及修复建议。",
             {"type": "object", "properties": {"path": {"type": "string"}},
              "required": []}, analyze_pdb),
        Tool("repair_structure",
             "按 actions 修复 PDB 并写出新文件: dedupe_models(只保留第一个 MODEL, "
             "修 NMR 系综/多构象叠加), keep_altloc_a(去备选位置), "
             "trim_disordered(裁末端无序 backbone-only 残基), "
             "drop_metals(去金属离子), drop_hetatm(去全部 HETATM)。"
             "analyze_pdb 返回的 issues[].fix 给出建议 action。",
             {"type": "object",
              "properties": {"path": {"type": "string"},
                             "actions": {"type": "array",
                                         "items": {"type": "string",
                                                   "enum": ["dedupe_models",
                                                            "keep_altloc_a",
                                                            "trim_disordered",
                                                            "drop_metals",
                                                            "drop_hetatm"]}},
                             "out": {"type": "string"}},
              "required": ["path", "actions"]}, repair_structure),
        Tool("clean_pdb",
             "清洗 PDB (去水/去无关 HETATM, 可选保留配体与指定链)。",
             {"type": "object",
              "properties": {"path": {"type": "string"},
                             "out": {"type": "string"},
                             "keep_resnames": {"type": "array",
                                                "items": {"type": "string"}},
                             "keep_waters": {"type": "boolean"},
                             "keep_chain": {"type": "string"}},
              "required": ["path"]}, clean_pdb),
        Tool("find_pockets",
             "对接口袋检测: 有已知配体时用配体位点; 否则 0.5A 网格空腔按体积排序"
             "(返回 top 候选, 结果存入 state.target_prep.pocket)。",
             {"type": "object",
              "properties": {"pdb": {"type": "string"},
                             "lig_resname": {"type": "string"},
                             "grid": {"type": "number"},
                             "pad": {"type": "number"},
                             "top_n": {"type": "integer"}},
              "required": []}, find_pockets),
        Tool("pdb_to_pdbqt",
             "PDB -> PDBQT (obabel)。flex=True (默认, 配体/VHH 等柔性配体, "
             "自动补 ROOT/TORSDOF); flex=False (刚性受体)。",
             {"type": "object",
              "properties": {"pdb": {"type": "string"},
                             "out": {"type": "string"},
                             "keep_resnames": {"type": "array",
                                                "items": {"type": "string"}},
                             "flex": {"type": "boolean"}},
              "required": ["pdb"]}, pdb_to_pdbqt),
        Tool("run_target_prep",
             "整段靶点准备 (确定性标准流程: 解析→分析→LLM判断→清洗→口袋→pdbqt)。"
             "需要定制时改用上面的细粒度工具。",
             {"type": "object", "properties": {}, "required": []},
         run_target_prep),
    ]
