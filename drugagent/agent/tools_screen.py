"""Small-molecule screening tools (fine-grained + stage fallback)."""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import DEFAULTS, MODELS
from ..llm import AgentBrain
from ..modules import screening as sc
from ..utils import jload, jsave, pmap
from .loop import Ctx, Tool, ToolError
from .stages import maybe_reuse, summarize_screening


def _brain(ctx: Ctx) -> AgentBrain | None:
    if ctx.options.get("no_llm") or ctx.brain is None:
        return None
    return ctx.brain


def _d(ctx: Ctx):
    from ..config import resolve_defaults
    return resolve_defaults(ctx.options)


def _write_sdf(blocks: list[tuple[str, dict, str]], out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        for title, props, molblock in blocks:
            fh.write(title + "\n")
            for k, v in props.items():
                fh.write(f">  <{k}>\n{v}\n")
            fh.write(molblock)
            if not molblock.endswith("\n"):
                fh.write("\n")
            fh.write("$$$$\n")
    return out


def _blocks(sdf: Path) -> list[tuple[str, dict, str, str]]:
    """blocks + smiles (canonical_smiles prop preferred)."""
    out = []
    for i, (title, props, molblock) in enumerate(sc._parse_sdf_blocks(sdf)):
        smi = props.get("canonical_smiles") or props.get("smiles") or ""
        out.append((i, title, props, molblock, smi))
    return [(i, t, p, m, s) for i, t, p, m, s in out]


# --------------------------------------------------------------------------- #
def standardize_library(ctx: Ctx, sdf: str | None = None,
                        out: str = "02_screening/std/ok.sdf") -> dict:
    """Sanitize + 3D-embed a vendor SDF; write the ok subset."""
    d = _d(ctx)
    workdir = ctx.stage_dir("02_screening") / "std"
    workdir.mkdir(parents=True, exist_ok=True)
    if sdf:
        lib = Path(sdf)
        if not lib.exists():
            raise ToolError(f"library not found: {lib}")
        lib_used = "custom"
    else:
        # R10/G6: fallback-aware (flaky DTP mirror -> ChEMBL subsample)
        lib, lib_used = sc.resolve_library(ctx.options, d)
    n_jobs = min(int(ctx.options.get("n_jobs", 32)), 16)
    t0 = time.time()
    df = sc.standardize_sdf(lib, workdir / "ok.sdf", n_jobs=n_jobs)
    ok = df[df.status == "ok"].reset_index(drop=True)
    # write ok sdf
    blocks = _blocks(lib)
    ok_idx = set(int(i) for i in ok["idx"])
    _write_sdf([(t, p, m) for i, t, p, m, s in blocks if i in ok_idx],
               workdir / Path(out).name)
    stats = jload(workdir / "ok.sdf.stats.json") \
        if (workdir / "ok.sdf.stats.json").exists() else {}
    return {"library": lib_used, "library_path": str(lib),
            "n_total": int(stats.get("n_total", len(df))),
            "n_ok": int(len(ok)),
            "ok_sdf": str(workdir / Path(out).name),
            "seconds": round(time.time() - t0, 1)}


def prefilter_ligands(ctx: Ctx, sdf: str | None = None,
                      n_keep: int | None = None) -> dict:
    """Features + physchem filter + ML ranking; write kept.sdf + kept.csv."""
    d = _d(ctx)
    n_jobs = int(ctx.options.get("n_jobs", 32))
    workdir = ctx.stage_dir("02_screening")
    if sdf:
        p = Path(sdf)
    elif (workdir / "std" / "ok.sdf").exists():
        p = workdir / "std" / "ok.sdf"
    else:
        raise ToolError("sdf not given and no 02_screening/std/ok.sdf (先 standardize_library)")
    ok_df = pd.DataFrame(
        [{"idx": i, "smiles": s} for i, t, pr, m, s in _blocks(p) if s])
    if not len(ok_df):
        raise ToolError("sdf 里没有带 SMILES 的分子")
    t0 = time.time()
    fdf = sc.compute_features(ok_df, n_jobs=min(n_jobs, 16))
    n_after_feat = len(fdf)
    fdf = sc.physchem_filter(fdf)
    n_after_phys = len(fdf)
    model_path = None
    for mc in (ctx.project_dir / "data" / "prefilter_rf.pkl", MODELS / "prefilter_rf.pkl"):
        if mc.exists():
            model_path = str(mc)
            break
    n_keep = n_keep or d.screen_max_ligands
    fdf = sc.ml_prefilter(fdf, model_path, n_keep=n_keep)
    kept_csv = workdir / "kept.csv"
    fdf.to_csv(kept_csv, index=False)
    keep_idx = set(int(i) for i in fdf["idx"])
    blocks = [(t, pr, m) for i, t, pr, m, s in _blocks(p) if i in keep_idx]
    kept_sdf = _write_sdf(blocks, workdir / "kept.sdf")
    return {"n_input": int(len(ok_df)), "n_after_features": n_after_feat,
            "n_after_physchem": n_after_phys, "n_kept": int(len(fdf)),
            "model": model_path, "kept_sdf": str(kept_sdf),
            "kept_csv": str(kept_csv),
            "seconds": round(time.time() - t0, 1),
            "top5": fdf.head(5)[["idx", "smiles", "ml_score", "QED", "SA"]]
            .to_dict("records")}


def dock_screen(ctx: Ctx, sdf: str | None = None,
                pocket: dict | None = None,
                n_dock: int | None = None,
                n_jobs: int | None = None,
                rescore_topn: int | None = None) -> dict:
    """Dock the (prefiltered) ligands into the pocket.

    Writes 02_screening/dock_results.csv (same schema as the 1.0 stage) and
    a partial screening.json; call redock_ligand + decide_hits to finish.
    """
    d = _d(ctx)
    prep = ctx.state().get("target_prep") or {}
    rec_pdbqt = prep.get("receptor_pdbqt")
    if not rec_pdbqt or not Path(rec_pdbqt).exists():
        raise ToolError("receptor_pdbqt 不存在 (先做靶点准备 / run_target_prep)")
    pocket = pocket or prep.get("pocket")
    if not pocket:
        raise ToolError("pocket 未设置 (先 find_pockets)")
    if sdf:
        p = Path(sdf)
    elif (ctx.project_dir / "02_screening" / "kept.sdf").exists():
        p = ctx.project_dir / "02_screening" / "kept.sdf"
    elif (ctx.project_dir / "02_screening" / "std" / "ok.sdf").exists():
        p = ctx.project_dir / "02_screening" / "std" / "ok.sdf"
    else:
        raise ToolError("sdf 未指定且没有 kept.sdf/ok.sdf")
    rows = [(i, s) for i, t, pr, m, s in _blocks(p) if s]
    n_dock = n_dock or len(rows)
    rows = rows[:n_dock]

    dockdir = ctx.stage_dir("02_screening") / "docks" / "agent"
    dockdir.mkdir(parents=True, exist_ok=True)
    n_jobs = n_jobs or int(ctx.options.get("n_jobs", 32))
    cpu_per = max(1, sc.os_cpu() // n_jobs)
    exh = d.dock_exhaustiveness_fast
    t0 = time.time()

    def _one(item):
        i, smi = item
        lig = dockdir / f"{i}.pdbqt"
        prefix = dockdir / f"{i}"
        if not sc.write_ligand_pdbqt(smi, lig, seed=42 + i % 1000):
            return {"idx": i, "smiles": smi, "ok": False, "score": np.nan,
                    "top_pose_pdbqt": None}
        r = sc.dock_one((rec_pdbqt, str(lig), str(prefix), pocket, exh, cpu_per))
        return {"idx": i, "smiles": smi, "ok": r["ok"],
                "score": r.get("score", np.nan),
                "rmsd_lb": r.get("rmsd_lb"), "rmsd_ub": r.get("rmsd_ub"),
                "top_pose_pdbqt": r.get("top_pose_pdbqt"),
                "error": r.get("error")}

    results = pmap(_one, rows, n_jobs=n_jobs)
    rdf = pd.DataFrame(results)
    rdf_ok = rdf[rdf.ok].sort_values("score")
    # GNINA rescore of the top poses
    rescore_topn = rescore_topn or d.dock_rescore_topn
    top = rdf_ok.head(rescore_topn)
    gmap = {}
    for _, row in top.iterrows():
        pose = row.get("top_pose_pdbqt")
        if pose and Path(pose).exists():
            g = sc.gnina_rescore(rec_pdbqt, str(pose), pocket)
            if g.get("ok"):
                gmap[int(row["idx"])] = g
    if gmap:
        rdf["gnina_score"] = rdf["idx"].map(
            lambda i: gmap.get(int(i), {}).get("gnina_score"))
        ok_g = rdf.gnina_score.notna()
        rdf["final_score"] = np.where(ok_g,
                                      0.5 * rdf.score + 0.5 * rdf.gnina_score,
                                      rdf.score)
        rdf = rdf.sort_values("final_score")
        score_col = "final_score"
    else:
        score_col = "score"
    rdf.to_csv(ctx.stage_dir("02_screening") / "dock_results.csv", index=False)
    scores = [float(s) for s in rdf[score_col].dropna().tolist()]
    st = ctx.state()
    scr = st.get("screening") or {}
    scr.update({
        "n_docked": int(len(rdf_ok)),
        "n_failed": int((~rdf.ok).sum()),
        "results_csv": str(ctx.project_dir / "02_screening" / "dock_results.csv"),
        "score_stats": {"mean": round(float(np.mean(scores)), 3) if scores else None,
                        "p5": round(float(np.percentile(scores, 5)), 3) if scores else None,
                        "p50": round(float(np.percentile(scores, 50)), 3) if scores else None,
                        "min": round(float(min(scores)), 3) if scores else None,
                        "score_column": score_col},
    })
    st["screening"] = scr
    ctx.save_state()
    return {"n_docked": int(len(rdf_ok)), "n_failed": int((~rdf.ok).sum()),
            "score_stats": scr["score_stats"],
            "top10": rdf.head(10)[["idx", "smiles", score_col]].to_dict("records"),
            "results_csv": scr["results_csv"],
            "seconds": round(time.time() - t0, 1)}


def redock_ligand(ctx: Ctx) -> dict:
    """Positive control: redock the crystallographic ligand at the pocket."""
    d = _d(ctx)
    prep = ctx.state().get("target_prep") or {}
    lig = prep.get("ligand_pdbqt")
    rec = prep.get("receptor_pdbqt")
    pocket = prep.get("pocket")
    if not (lig and rec and pocket and Path(lig).exists()):
        raise ToolError("需要 state.target_prep.ligand_pdbqt + receptor_pdbqt + pocket")
    prefix = ctx.stage_dir("02_screening") / "ref_redock"
    r = sc.dock_one((rec, lig, str(prefix), pocket,
                     d.dock_exhaustiveness_final, sc.os_cpu() // 2))
    if not r["ok"]:
        return {"ok": False, "error": r.get("error", "redock failed")}
    st = ctx.state()
    scr = st.get("screening") or {}
    scr["reference_ligand_score"] = float(r["score"])
    st["screening"] = scr
    ctx.save_state()
    return {"ref_score": float(r["score"]), "rmsd_lb": r.get("rmsd_lb"),
            "rmsd_ub": r.get("rmsd_ub")}


def make_flex_receptor(ctx: Ctx, rec_pdb: str | None = None,
                       lig_pdb: str | None = None,
                       cutoff: float = 5.0) -> dict:
    """Build the Vina --flex file: side chains of receptor residues near the
    ligand (default: crystal ligand / target ligand at the pocket)."""
    prep = ctx.state().get("target_prep") or {}
    rec_pdb = rec_pdb or prep.get("raw_pdb")
    lig_pdb = lig_pdb or prep.get("ligand_pdb")
    if not lig_pdb and prep.get("ligand_pdbqt"):
        # state 里没存 ligand_pdb 时用 ligand.pdbqt 的同名 .pdb
        cand = Path(prep["ligand_pdbqt"]).with_suffix(".pdb")
        if cand.exists():
            lig_pdb = str(cand)
    if not (rec_pdb and lig_pdb and Path(rec_pdb).exists()
            and Path(lig_pdb).exists()):
        raise ToolError(f"需要受体+配体 PDB (got {rec_pdb} / {lig_pdb})")
    # this vina build (AD4 lineage) parses the AD4 PDBQT layout, not plain
    # PDB columns — filter the obabel-converted PDBQT, not the raw PDB
    rec_pdbqt = ctx.stage_dir("02_screening") / "receptor_for_flex.pdbqt"
    if not rec_pdbqt.exists():
        from ..modules.target_prep import to_pdbqt
        to_pdbqt(Path(rec_pdb), rec_pdbqt, flex=False)
    out = ctx.stage_dir("02_screening") / "flex_receptor.pdbqt"
    r = sc.flex_sidechain_pdbqt(rec_pdbqt, Path(lig_pdb), out,
                                cutoff=cutoff)
    if r["n_residues"] == 0:
        raise ToolError(f"cutoff {cutoff} A 内没有受体残基 — 增大 cutoff 或检查坐标")
    st = ctx.state()
    scr = st.get("screening") or {}
    scr["flex_receptor_pdbqt"] = str(out)
    scr["flex_residues"] = r["residues"]
    st["screening"] = scr
    ctx.save_state()
    return r


def dock_conformer_set(ctx: Ctx, lig_pdbqt: str | None = None,
                       rec_pdbs: list[str] | None = None,
                       flex_pdbqt: str | None = None,
                       pocket: dict | None = None,
                       exhaustiveness: int | None = None,
                       max_conformers: int = 3, min_pop: float = 0.05) -> dict:
    """Dock one ligand against a SET of receptor conformations (crystal
    + MD cluster representatives) with optional flexible side chains.
    Returns per-conformer scores + consensus (mean = consensus value).

    This is the 柔性靶点工作流: 构象选择 (cluster reps) + 柔性对接 (--flex).
    """
    d = _d(ctx)
    prep = ctx.state().get("target_prep") or {}
    scr = ctx.state().get("screening") or {}
    lig_pdbqt = lig_pdbqt or prep.get("ligand_pdbqt")
    if not lig_pdbqt or not Path(lig_pdbqt).exists():
        raise ToolError("需要 ligand PDBQT (state.target_prep.ligand_pdbqt)")
    pocket = pocket or prep.get("pocket")
    if not pocket:
        raise ToolError("需要 pocket (state.target_prep.pocket)")
    flex_pdbqt = flex_pdbqt or scr.get("flex_receptor_pdbqt")
    exh = exhaustiveness or d.dock_exhaustiveness_final
    # R9: MD reps are consensus conformers — lower exhaustiveness unless
    # the caller pinned one explicitly
    md_exh = exhaustiveness or d.dock_md_rep_exhaustiveness

    # receptor conformations: explicit list, else crystal + cluster reps
    recs: list[tuple[str, Path]] = []
    if rec_pdbs:
        recs = [("conformer%d" % (i + 1), Path(p))
                for i, p in enumerate(rec_pdbs[:max_conformers])]
    else:
        # ligand-free receptor (receptor_pdb) if available, else raw
        raw = prep.get("receptor_pdb") or prep.get("raw_pdb")
        if raw and Path(raw).exists():
            recs.append(("crystal", Path(raw)))
        md_work = ctx.project_dir / "05_md"
        if (md_work / "analysis" / "clusters_r1.pdb").is_file():
            gmx = ((ctx.state().get("md") or {}).get("gromacs") or {}).get(
                "binary")
            # R7: pool cluster representatives across ALL replicas
            # (population-ordered, Cα-RMSD deduplicated) — the ensemble
            # for conformation-selection docking
            reps = sc.pool_representatives(md_work,
                                           max_n=max_conformers,
                                           min_pop=min_pop, gmx=gmx)
            for r in reps:
                recs.append((f"md_r{r['rep']}_c{r['cluster']}", r["pdb"]))
    if not recs:
        raise ToolError("没有可用的受体构象 (raw_pdb 与 MD 聚类代表结构都缺失)")

    base = ctx.stage_dir("02_screening") / "conformer_dock"
    base.mkdir(parents=True, exist_ok=True)
    # flex atoms anchor by coordinates -> each conformer needs its OWN flex
    # file built from its own PDBQT (crystal ligand position selects the
    # residues; the pocket region is what matters)
    prep = ctx.state().get("target_prep") or {}
    lig_pdb = prep.get("ligand_pdb")
    if not lig_pdb and prep.get("ligand_pdbqt"):
        cand = Path(prep["ligand_pdbqt"]).with_suffix(".pdb")
        if cand.exists():
            lig_pdb = str(cand)
    results = []
    for name, pdb in recs:
        from ..modules.target_prep import to_pdbqt
        work_pdb = pdb
        # MD cluster representatives live in the simulation-box frame; the
        # docking grid is in the crystal frame -> align onto the receptor
        if name.startswith("md_"):
            ref = prep.get("receptor_pdb") or prep.get("raw_pdb")
            if ref and Path(ref).exists():
                aligned = base / f"{name}_aligned.pdb"
                n_ca = sc.align_pdb_to_reference(pdb, Path(ref), aligned)
                if n_ca >= 3:
                    work_pdb = aligned
        rec_pdbqt = base / f"rec_{name}.pdbqt"
        to_pdbqt(work_pdb, rec_pdbqt, flex=False)
        flex_file = flex_pdbqt
        if flex_file and flex_file == (ctx.stage_dir("02_screening")
                                       / "flex_receptor.pdbqt"):
            flex_file = None  # crystal default: rebuild per conformer below
        if flex_file is not None and not str(flex_file).startswith(str(base)):
            # external flex file: only valid for the crystal conformer
            flex_file = flex_file if name == "crystal" else None
        if flex_file is None and lig_pdb:
            per = base / f"flex_{name}.pdbqt"
            sc.flex_sidechain_pdbqt(rec_pdbqt, Path(lig_pdb), per,
                                    cutoff=float(d.flex_cutoff_ang))
            flex_file = str(per)
        prefix = base / f"pose_{name}"
        exh_i = int(md_exh) if name.startswith("md_") else int(exh)
        args = (str(rec_pdbqt), lig_pdbqt, str(prefix), pocket,
                exh_i, sc.os_cpu() // max(1, len(recs)))
        if flex_file:
            args = args + (flex_file,)
        r = sc.dock_one(args)
        results.append({"conformer": name, "receptor_pdb": str(pdb),
                        "ok": r["ok"], "score": r.get("score"),
                        "rmsd_lb": r.get("rmsd_lb"),
                        "rmsd_ub": r.get("rmsd_ub"),
                        "error": r.get("error")})
    ok = [r for r in results if r["ok"]]
    out = {"n_conformers": len(results), "results": results}
    if ok:
        stats = sc.consensus_stats([r["score"] for r in ok])
        out["consensus"] = stats
        out["consensus_score"] = stats["mean"]
    st = ctx.state()
    s2 = st.get("screening") or {}
    s2["conformer_dock"] = out
    st["screening"] = s2
    ctx.save_state()
    return out


def decide_hits(ctx: Ctx, n_wanted: int | None = None) -> dict:
    """Choose threshold + hits from the docked scores; write final screening.json."""
    d = _d(ctx)
    csv = ctx.project_dir / "02_screening" / "dock_results.csv"
    if not csv.exists():
        raise ToolError("dock_results.csv 不存在 (先 dock_screen)")
    rdf = pd.read_csv(csv)
    col = "final_score" if "final_score" in rdf.columns else "score"
    scores = [float(s) for s in rdf[col].dropna().tolist()]
    ref = (ctx.state().get("screening") or {}).get("reference_ligand_score")
    n_wanted = n_wanted or d.n_hits
    hit = sc.decide_hits(scores, ref, n_wanted, _brain(ctx))
    # poses for the top hits
    hits_dir = ctx.stage_dir("02_screening") / "hits"
    hits_dir.mkdir(exist_ok=True)
    hits = []
    for _, row in rdf.head(hit["n_hits"]).iterrows():
        pose_src = row.get("top_pose_pdbqt")
        pose_path = Path(pose_src) if pose_src and Path(pose_src).is_file() else None
        pose_out = hits_dir / f"hit_{len(hits)+1}_{int(row['idx'])}.pdb"
        if pose_path is not None:
            with open(pose_path) as fi, open(pose_out, "w") as fo:
                for line in fi:
                    if line.startswith("ATOM"):
                        fo.write(line)
            pose_out.write_text(pose_out.read_text() + "END\n")
        hits.append({
            "rank": len(hits) + 1,
            "idx": int(row["idx"]),
            "smiles": str(row["smiles"]),
            "vina_score": float(row["score"]) if pd.notna(row.get("score")) else None,
            "gnina_score": float(row["gnina_score"]) if "gnina_score" in row and pd.notna(row.get("gnina_score")) else None,
            "final_score": float(row[col]),
            "pose_pdbqt": str(pose_path) if pose_path is not None else None,
            "pose_pdb": str(pose_out) if pose_out.is_file() else None,
        })
    st = ctx.state()
    scr = st.get("screening") or {}
    scr.update({"hit_decision": hit, "hits": hits,
                "library": scr.get("library", ""),
                "n_docked": scr.get("n_docked", int(len(rdf)))})
    st["screening"] = scr
    jsave(ctx.stage_dir("02_screening") / "screening.json", scr)
    return {"decision": hit, "n_hits": len(hits),
            "top3": [{"rank": h["rank"], "smiles": h["smiles"],
                      "final_score": h["final_score"]} for h in hits[:3]]}


def run_screening(ctx: Ctx, force: bool = False) -> dict:
    """Full deterministic screening stage (1.0 pipeline as one tool).

    R11/G8: reuses state.screening when it is already complete
    (hit_decision present); force=true to re-run."""
    cached = maybe_reuse(ctx, "screening", force)
    if cached is not None:
        return cached
    st = ctx.state()
    if "target_prep" not in st:
        raise ToolError("state.target_prep 未就绪 (先 run_target_prep)")
    state_out = sc.screen(st)           # returns full state
    ctx.save_state(screening=state_out["screening"])
    return summarize_screening(ctx.state())


# --------------------------------------------------------------------------- #
def build() -> list[Tool]:
    return [
        Tool("standardize_library",
             "小分子库标准化 + 3D 构象 (RDKit); 写 ok.sdf。sdf 缺省用 options.library。",
             {"type": "object",
              "properties": {"sdf": {"type": "string"},
                             "out": {"type": "string"}},
              "required": []}, standardize_library),
        Tool("prefilter_ligands",
             "特征计算 + 物理化学过滤 + ML 排序, 保留 top N (默认 fast:50/full:200), "
             "写 kept.sdf + kept.csv。",
             {"type": "object",
              "properties": {"sdf": {"type": "string"},
                             "n_keep": {"type": "integer"}},
              "required": []}, prefilter_ligands),
        Tool("dock_screen",
             "批量对接 kept 分子 (vina 并行 + 顶部 GNINA 重打分)。写 "
             "dock_results.csv。pocket 可显式传入覆盖 state。",
             {"type": "object",
              "properties": {"sdf": {"type": "string"},
                             "pocket": {"type": "object"},
                             "n_dock": {"type": "integer"},
                             "n_jobs": {"type": "integer"},
                             "rescore_topn": {"type": "integer"}},
              "required": []}, dock_screen),
        Tool("make_flex_receptor",
             "生成 Vina --flex 文件: 口袋附近受体残基的侧链 (柔性对接用)。",
             {"type": "object",
              "properties": {"rec_pdb": {"type": "string"},
                             "lig_pdb": {"type": "string"},
                             "cutoff": {"type": "number"}},
              "required": []}, make_flex_receptor),
        Tool("dock_conformer_set",
             "多构象+柔性对接: 配体对接到 crystal + MD 聚类代表构象 "
             "(--flex 侧链), 返回每个构象的分数 + consensus (均值)。"
             "柔性靶点的构象选择工作流。",
             {"type": "object",
              "properties": {"lig_pdbqt": {"type": "string"},
                             "rec_pdbs": {"type": "array",
                                          "items": {"type": "string"}},
                             "flex_pdbqt": {"type": "string"},
                             "pocket": {"type": "object"},
                             "exhaustiveness": {"type": "integer"},
                             "max_conformers": {"type": "integer"},
                             "min_pop": {"type": "number"}},
              "required": []}, dock_conformer_set),
        Tool("redock_ligand",
             "阳性对照: 已知配体重对接, 得到参考分。",
             {"type": "object", "properties": {}, "required": []},
         redock_ligand),
        Tool("decide_hits",
             "根据打分分布 + 参考分确定阈值与命中数 (LLM 参与), 写最终 hits 与 "
             "screening.json。",
             {"type": "object",
              "properties": {"n_wanted": {"type": "integer"}},
              "required": []}, decide_hits),
        Tool("run_screening",
             "整段筛选 (确定性标准流程: 标准化→过滤→对接→GNINA→命中判定)。" "已完成阶段自动复用, force=true 强制重跑。",
             {"type": "object",
               "properties": {"force": {"type": "boolean",
                                        "description": "强制重跑 (默认复用已完成阶段)"}},
               "required": []},
         run_screening),
    ]
