"""DrugAgent CLI (typer) — 2.0 agent semantics.

Commands:
  setup   - one-click environment configuration
  run     - run the agent (the LLM drives the pipeline toolkit)
  resume  - resume a crashed/paused agent run (replays the transcript)
  status  - show run status
  report  - regenerate the report for a project
"""
from __future__ import annotations

import json
import os
import sys
from datetime import datetime
from pathlib import Path

# hide a read-only user site-packages that shadows our env
_us = os.path.join(os.path.expanduser("~"), ".local", "lib")
if os.path.isdir(_us):
    os.environ.setdefault("PYTHONNOUSERSITE", "1")
    sys.path = [p for p in sys.path
                if not (p.startswith(_us) or ".local/lib" in p)]

import typer
from loguru import logger

from .config import DEFAULTS, LLMConfig, PROJECTS
from .utils import is_pdb_id

app = typer.Typer(help="DrugAgent - LLM 驱动的药物筛选 agent (工具: 靶点/筛选/设计/MD/报告)",
                  no_args_is_help=True)


DEFAULT_MODULES = ["screen", "binder", "vhh", "md"]


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _parse_target(value: str) -> dict:
    v = value.strip()
    if v.endswith((".pdb", ".PDB")) or (Path(v).exists() and Path(v).suffix.lower() == ".pdb"):
        return {"kind": "pdb_file", "value": v}
    if v.endswith((".fasta", ".fa", ".fasta.gz")):
        return {"kind": "fasta", "value": v}
    if is_pdb_id(v):
        return {"kind": "pdb_id", "value": v.upper()}
    aa = set("ACDEFGHIKLMNPQRSTVWY")
    if len(v) >= 30 and set(v.upper()) <= aa:
        return {"kind": "sequence", "value": v.upper()}
    raise typer.BadParameter(
        f"无法识别的靶点输入: {v!r} (支持: PDB文件路径 / PDB ID / FASTA / 裸序列)")


def _project_dir(name: str | None) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    slug = name or f"run_{stamp}"
    p = PROJECTS / slug
    p.mkdir(parents=True, exist_ok=True)
    return p


def _load_state(project: Path) -> dict:
    """State from state.json (2.0), falling back to 1.0 layout."""
    st_path = project / "state.json"
    if st_path.exists():
        import json as _json
        state = _json.loads(st_path.read_text())
        state["project_dir"] = str(project)
        return state
    state: dict = {"project_dir": str(project)}
    cfg_path = project / "run_config.json"
    if cfg_path.exists():
        cfg = json.loads(cfg_path.read_text())
        state["target"] = cfg["target"]
        state["options"] = cfg["options"]
    for stage, key in (("01_target", "target_prep"), ("02_screening", "screening"),
                       ("03_binder", "binder"), ("04_vhh", "vhh"),
                       ("05_md", "md")):
        jp = project / stage / f"{key}.json"
        if jp.exists():
            state[key] = json.loads(jp.read_text())
    return state


def _make_prompt_fn() -> callable:
    def prompt_fn(title: str, question: str,
                  options: list[str] | None, context: str = "") -> str:
        print("\n" + "=" * 72)
        print(f"[{title}]")
        print(question)
        if context:
            print(context[:2000])
        if options:
            for i, o in enumerate(options):
                print(f"  {i + 1}. {o}")
        while True:
            try:
                ans = input("回答/选择 > ").strip()
            except EOFError:
                print("(EOF -> 默认批准)")
                return options[0] if options else "approve"
            if options and ans.isdigit() and 1 <= int(ans) <= len(options):
                return options[int(ans) - 1]
            if ans or not options:
                return ans
            print("  (输入编号或文字)")
    return prompt_fn


# --------------------------------------------------------------------------- #
# setup
# --------------------------------------------------------------------------- #
@app.command()
def setup(
    libraries: str = typer.Option("dtp", help="下载哪些库: dtp,chembl,pdbbind,vhh,all"),
    gromacs: bool = typer.Option(True, help="构建/复用 GROMACS 2023"),
    tools: bool = typer.Option(True, help="安装 vina/gnina/fpocket/3Dmol"),
    rfdiffusion: bool = typer.Option(True, help="部署 RFdiffusion + 权重"),
):
    """一键配置环境 (幂等, 可重复执行)."""
    from .setup import run_setup
    run_setup(libraries=libraries, gromacs=gromacs, tools=tools,
              rfdiffusion=rfdiffusion)


# --------------------------------------------------------------------------- #
# run
# --------------------------------------------------------------------------- #
@app.command()
def run(
    target: str = typer.Option(..., help="PDB文件 / PDB ID / FASTA / 序列"),
    modules: str = typer.Option("all", help="逗号分隔子集: screen,binder,vhh,md"),
    fast: bool = typer.Option(False, help="快速模式 (小规模验证)"),
    auto: bool = typer.Option(False, help="检查点自动通过 (不等人)"),
    name: str = typer.Option(None, help="项目名 (默认时间戳)"),
    library: str = typer.Option("dtp", help="小分子库: dtp/chembl/pdbbind 或 SDF 路径"),
    n_jobs: int = typer.Option(32, help="并行度"),
    md_ns: float = typer.Option(None, help="MD 时长 (ns), 覆盖默认"),
    md_reps: int = typer.Option(None, help="MD 重复次数, 覆盖默认"),
    md_salt: float = typer.Option(None, help="MD 离子浓度 (M, 默认 0.15)"),
    md_divalent: str = typer.Option(None, help="MD 二价抗衡离子 (如 MG/CA, 默认无)"),
    md_divalent_m: float = typer.Option(None, help="二价离子浓度 (M, 默认 0.01)"),
    md_extend_ns: float = typer.Option(None, help="MD 自动延长步长 (ns, 0=用 md_ns)"),
    md_max_extensions: int = typer.Option(None, help="MD 收敛前自动延长最大轮数"),
    md_burn_in_ps: float = typer.Option(None, help="生产 MD 烧入段剔除 (ps)"),
    max_steps: int = typer.Option(300, help="agent 步数预算 (LLM 轮数)"),
    llm_base: str = typer.Option(None, help="LLM base url 覆盖"),
    llm_model: str = typer.Option(None, help="LLM model 覆盖"),
    no_llm: bool = typer.Option(False, help="不用 LLM (确定性脚本模式)"),
    vhh_plddt_min: float = typer.Option(None, help="VHH pLDDT 门槛 (默认 fast 35 / full 50)"),
    vhh_dock_flex: bool = typer.Option(None, help="VHH 柔性对接 (默认刚性, R11/G10)"),
    vhh_dock_cdr_only: bool = typer.Option(None, help="VHH CDR 片段对接 (fast 默认开, R11/G10-v2)"),
):
    """运行 agent (LLM 主循环驱动工具; --no-llm 退化为脚本模式)."""
    from .agent import (AgentLoop, Ctx, build_tools, goal_text,
                        system_prompt)
    from .agent.loop import scripted_run
    from .llm import AgentBrain

    tinfo = _parse_target(target)
    # absolute target paths keep the agent from confusing project dir with repo root
    if tinfo["kind"] in ("pdb_file", "fasta") and not Path(tinfo["value"]).is_absolute():
        tinfo["value"] = str(Path(tinfo["value"]).resolve())
    mods = (sorted(set(m.strip() for m in modules.split(",")))
            if modules != "all" else DEFAULT_MODULES)
    pdir = _project_dir(name)
    logger.info(f"project: {pdir}")

    d = DEFAULTS.resolved(fast)
    options = {
        "modules": mods,
        "fast": fast,
        "auto": auto,
        "library": library if not Path(library).exists() else Path(library).name,
        "library_path": str(Path(library)) if Path(library).exists() else None,
        "n_jobs": n_jobs,
        "md_ns": md_ns or d.md_ns,
        "md_reps": md_reps or d.md_reps,
        "no_llm": no_llm,
        "llm_model": llm_model or LLMConfig.from_env().model,
        "max_steps": max_steps,
    }
    # MD 旋钮透传 (resolve_defaults 会把已知字段覆盖到 Defaults 上)
    for key, val in (("md_salt_m", md_salt), ("md_divalent", md_divalent),
                     ("md_divalent_m", md_divalent_m),
                     ("md_extend_ns", md_extend_ns),
                     ("md_max_extensions", md_max_extensions),
                     ("md_burn_in_ps", md_burn_in_ps)):
        if val is not None:
            options[key] = val
    if llm_base:
        os.environ["DRUGAGENT_LLM_BASE_URL"] = llm_base
    if llm_model:
        os.environ["DRUGAGENT_LLM_MODEL"] = llm_model
    # R11: VHH 旋钮透传 (resolve_defaults 会把已知 Defaults 字段叠上去)
    if vhh_plddt_min is not None:
        options["vhh_plddt_min"] = vhh_plddt_min
    if vhh_dock_flex is not None:
        options["vhh_dock_flex"] = vhh_dock_flex
    if vhh_dock_cdr_only is not None:
        options["vhh_dock_cdr_only"] = vhh_dock_cdr_only

    brain = None if no_llm else AgentBrain(project_dir=pdir)
    ctx = Ctx(pdir, brain, options, auto=auto, target=tinfo,
              prompt_fn=None if auto else _make_prompt_fn())
    tools = build_tools()

    if no_llm:
        result = scripted_run(ctx, tools)
    else:
        loop = AgentLoop(ctx, tools,
                         system=system_prompt(options),
                         goal=goal_text(tinfo, options),
                         max_steps=max_steps)
        result = loop.run()

    status = result.get("status")
    icon = {"success": "✔", "needs_human": "?", "failed": "✘"}.get(status, "?")
    typer.secho(f"\n=== agent 结束: {status} @ step {result.get('step')} ===",
                fg=typer.colors.GREEN if status == "success" else
                typer.colors.YELLOW)
    typer.echo(result.get("summary", ""))
    r = pdir / "reports" / "report.html"
    if r.exists():
        typer.secho(f"报告: {r}", fg=typer.colors.CYAN)
    typer.echo(f"transcript: {pdir / 'agent' / 'transcript.jsonl'}")
    if status != "success":
        typer.echo(f"恢复: drugagent resume --project {pdir}")


# --------------------------------------------------------------------------- #
# resume
# --------------------------------------------------------------------------- #
@app.command()
def resume(
    project: str = typer.Option(..., help="项目目录"),
    auto: bool = typer.Option(False, help="检查点自动通过"),
):
    """恢复 agent 运行 (重放 transcript, 从断点继续)."""
    from .agent import AgentLoop, Ctx, build_tools, goal_text, system_prompt
    from .agent.loop import scripted_run
    from .llm import AgentBrain

    pdir = Path(project)
    state = _load_state(pdir)
    options = state.get("options", {})
    tinfo = state.get("target", {})
    brain = None if options.get("no_llm") else AgentBrain(project_dir=pdir)
    ctx = Ctx(pdir, brain, options, auto=auto, target=tinfo,
              prompt_fn=None if auto else _make_prompt_fn())
    tools = build_tools()
    if options.get("no_llm"):
        result = scripted_run(ctx, tools)
    else:
        loop = AgentLoop(ctx, tools,
                         system=system_prompt(options),
                         goal=goal_text(tinfo, options),
                         max_steps=int(options.get("max_steps", 300)))
        loop.messages.append({"role": "user",
                              "content": "继续 (运行已恢复; 先 state_get 查看进度, "
                                         "已完成的产物会被幂等复用)。"})
        with (ctx.agent_dir / "transcript.jsonl").open("a") as fh:
            fh.write(json.dumps({"step": loop.step, "role": "user",
                                 "content": "继续 (运行已恢复)"},
                                ensure_ascii=False) + "\n")
        result = loop.run()
    typer.echo(f"status={result.get('status')} step={result.get('step')} "
               f"{result.get('summary', '')}")


# --------------------------------------------------------------------------- #
# rerun (R12/G8: force the CLI-ification)
# --------------------------------------------------------------------------- #
_STAGE_TOOL = {"target_prep": "run_target_prep", "screening": "run_screening",
               "binder": "run_design", "vhh": "run_vhh", "md": "run_md",
               "report": "build_report"}


@app.command()
def rerun(project: str = typer.Option(None, help="项目目录 (默认最新)"),
          stage: str = typer.Option(..., help="要重跑的阶段: "
          "target_prep|screening|binder|vhh|md|report"),
          with_report: bool = typer.Option(True, help="阶段成功后重建报告")):
    """强制重跑单个阶段 (G8: force=true 绕过 stage 复用; 其余阶段产物不动)."""
    from .agent import Ctx, build_tools
    pdir = Path(project) if project else _latest_project()
    if pdir is None:
        typer.echo("没有项目")
        raise typer.Exit(1)
    if not (pdir / "state.json").is_file():
        typer.echo("state.json 不存在 — 先跑 `run`")
        raise typer.Exit(1)
    tool_name = _STAGE_TOOL.get(stage)
    if tool_name is None:
        typer.echo(f"未知阶段: {stage} (可选 {', '.join(_STAGE_TOOL)})")
        raise typer.Exit(1)
    options = dict(_load_state(pdir).get("options") or {})
    brain = None  # stage tools are deterministic; no LLM needed
    ctx = Ctx(pdir, brain, options, auto=True)
    tools = {t.name: t for t in build_tools()}
    ctx._step = 1
    result = tools[tool_name].call(ctx, force=True)
    ok = bool(result.get("ok"))
    ctx.save_state(status="success" if ok else "failed")
    icon = "✔" if ok else "✘"
    typer.secho(f"=== rerun {stage}: {'成功' if ok else '失败'} ===",
                fg=typer.colors.GREEN if ok else typer.colors.RED)
    summary = str(result.get("summary", ""))
    if summary:
        typer.echo(summary[:400])
    if not ok:
        err = result.get("error")
        if err:
            typer.echo(f"error: {str(err)[:300]}")
        raise typer.Exit(1)
    if with_report and stage != "report":
        typer.echo("重建报告...")
        r = tools["build_report"].call(ctx, force=True)
        if r.get("ok"):
            typer.echo(f"报告: {r.get('report_html', pdir / 'reports' / 'report.html')}")
        else:
            typer.echo(f"报告重建失败: {str(r.get('error'))[:200]}")


# --------------------------------------------------------------------------- #
# status
# --------------------------------------------------------------------------- #
_STAGE_JSON = {"target_prep": "01_target", "screening": "02_screening",
               "binder": "03_binder", "vhh": "04_vhh", "md": "05_md"}


def _transcript_failures(tpath: Path, n: int = 3) -> list[str]:
    """Last n failed tool calls from the transcript (content carries the
    serialized tool result, which has ok=false on failure)."""
    if not tpath.is_file():
        return []
    fails = []
    for line in tpath.read_text(errors="ignore").splitlines():
        if not line.strip():
            continue
        try:
            obj = json.loads(line)
        except json.JSONDecodeError:
            continue
        if obj.get("role") != "tool":
            continue
        try:
            res = json.loads(obj.get("content") or "{}")
        except (json.JSONDecodeError, TypeError):
            continue
        if isinstance(res, dict) and res.get("ok") is False:
            fails.append(f"step {obj.get('step')} {obj.get('name')}: "
                         f"{str(res.get('error'))[:160]}")
    return fails[-n:]


def _fmt_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit != "B" else f"{n} B"
        n /= 1024
    return f"{n:.0f} GB"


@app.command()
def status(project: str = typer.Option(None, help="项目目录 (默认最新)")):
    """查看运行状态 (阶段完成度 + 关键数字 + 产物 + 最近错误)."""
    from .agent.stages import status_lines

    pdir = Path(project) if project else _latest_project()
    if pdir is None:
        typer.echo("没有项目")
        return
    state = _load_state(pdir)
    typer.echo(f"项目: {pdir}")
    st = pdir / "state.json"
    if st.exists():
        d = json.loads(st.read_text())
        typer.echo(f"状态: {d.get('status')}")
    # per-stage completion + key numbers (R11: stage-level detail)
    for line in status_lines(state, pdir):
        mark = "x" if line["done"] else " "
        typer.echo(f"  [{mark}] {line['stage']:<12} {line['detail']}")
    # stage json artifacts on disk
    artifacts = []
    for key, subdir in _STAGE_JSON.items():
        jp = pdir / subdir / f"{key}.json"
        if jp.is_file():
            artifacts.append((str(jp.relative_to(pdir)), jp.stat().st_size))
    report_html = pdir / "reports" / "report.html"
    if report_html.is_file():
        artifacts.append(("reports/report.html", report_html.stat().st_size))
    if artifacts:
        typer.echo("产物:")
        for name, size in artifacts:
            typer.echo(f"  {name} ({_fmt_size(size)})")
    # recent tool failures (state.errors + transcript ok=false)
    fails = [str(e) for e in state.get("errors", [])][-3:]
    fails += _transcript_failures(pdir / "agent" / "transcript.jsonl")
    if fails:
        typer.echo("最近错误:")
        for f in fails[-3:]:
            typer.echo(f"  {f}")
    t = pdir / "agent" / "transcript.jsonl"
    if t.exists():
        lines = t.read_text(errors="ignore").splitlines()
        typer.echo(f"transcript: {len(lines)} 条消息")
        if lines:
            try:
                last = json.loads(lines[-1])
                typer.echo(f"最后: step {last.get('step')} {last.get('role')} "
                           f"{str(last.get('name', last.get('content', '')))[:80]}")
            except json.JSONDecodeError:
                pass


def _latest_project() -> Path | None:
    if not PROJECTS.exists():
        return None
    runs = sorted([p for p in PROJECTS.iterdir() if p.is_dir()],
                  key=lambda p: p.stat().st_mtime)
    return runs[-1] if runs else None


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #
@app.command()
def report(project: str = typer.Option(..., help="项目目录")):
    """(重新)生成 HTML + PDF 报告."""
    from .report import report as report_mod
    pdir = Path(project)
    state = _load_state(pdir)
    out = report_mod.build_report(state)
    typer.echo(f"HTML: {out['html']}")
    typer.echo(f"PDF:  {out['pdf']}")


if __name__ == "__main__":
    app()
