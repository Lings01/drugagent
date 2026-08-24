"""Build the project report: interactive HTML + PDF.

HTML: Plotly interactive charts + 3Dmol.js structure viewers + decision log.
PDF:  same content with static PNG charts (WeasyPrint).
"""
from __future__ import annotations

import base64
import html
import json
from datetime import datetime
from pathlib import Path

from loguru import logger

from ..config import TOOLS


CSS = """
body { font-family: 'Helvetica Neue', Arial, 'Noto Sans CJK SC', sans-serif;
       margin: 2em auto; max-width: 1100px; color: #222; line-height: 1.5; }
h1 { border-bottom: 3px solid #2c7fb8; padding-bottom: .3em; }
h2 { color: #2c7fb8; margin-top: 1.8em; }
table { border-collapse: collapse; margin: 1em 0; font-size: .92em; }
th, td { border: 1px solid #ccc; padding: 5px 9px; text-align: left; }
th { background: #eef4fa; }
tr:nth-child(even) { background: #f8fafc; }
.card { background: #f6f9fc; border-left: 4px solid #2c7fb8;
        padding: .8em 1.2em; margin: 1em 0; border-radius: 4px; }
.small { color: #666; font-size: .85em; }
code { background: #eef; padding: 1px 4px; border-radius: 3px; }
.viewer { width: 100%; height: 420px; border: 1px solid #ddd; border-radius: 6px;
          margin: .6em 0; background: #fff; }
.badge { display: inline-block; background: #2c7fb8; color: #fff;
         border-radius: 10px; padding: 1px 10px; font-size: .8em; }
"""


def _esc(x) -> str:
    return html.escape(str(x))


def _fmt(v, nd: int = 3) -> str:
    try:
        return f"{float(v):.{nd}f}"
    except (TypeError, ValueError):
        return "?"


def _table(headers: list[str], rows: list[list]) -> str:
    if not rows:
        return '<p class="small">(空)</p>'
    th = "".join(f"<th>{_esc(h)}</th>" for h in headers)
    trs = []
    for r in rows:
        tds = "".join(f"<td>{_esc(v)}</td>" for v in r)
        trs.append(f"<tr>{tds}</tr>")
    return f"<table><tr>{th}</tr>{''.join(trs)}</table>"


def _plot_div(fig, plotlyjs: str = "cdn") -> str:
    try:
        return fig.to_html(full_html=False, include_plotlyjs=plotlyjs,
                           config={"displayModeBar": False})
    except TypeError:
        # plotly >= 6 removed include_plotlyjs kw
        return fig.to_html(full_html=False,
                           config={"displayModeBar": False})


def _plot_png_b64(fig) -> str:
    png = fig.to_image(format="png", width=900, height=320, scale=1)
    return base64.b64encode(png).decode()


def _viewer_pdb(pdb_path: str | Path, name: str, ligand_resnames: list[str] | None = None) -> str:
    """3Dmol.js viewer with embedded PDB."""
    if not pdb_path:
        return '<p class="small">(无结构)</p>'
    p = Path(pdb_path)
    if not p.is_file():
        return f'<p class="small">结构文件缺失: {_esc(p)}</p>'
    pdb_text = p.read_text()
    lig_sel = ",".join(f'elem == "{c[0]}" and resname == "{r}"'
                       for r in (ligand_resnames or []) for c in [r[:1]])
    script = f"""
(function(){{
  if (typeof $3Dmol === 'undefined') {{
    document.getElementById('viewer_{name}').innerHTML = '<p class="small">3Dmol.js 未加载</p>';
    return;
  }}
  var v = $3Dmol.createViewer('viewer_{name}');
  var addw = v.addModel({json.dumps(pdb_text)}, 'pdb');
  v.setStyle({{}});
  if ({json.dumps(bool(ligand_resnames))}) {{
    var lig = {json.dumps(ligand_resnames)};
    v.setStyle({{ resn: lig }}, {{stick: {{radius: 0.15, color: 'red'}}}});
  }}
  v.zoomTo(); v.render();
}})();
"""
    return (f'<div id="viewer_{name}" class="viewer"></div>'
            f'<script type="text/plain" id="pdb_{name}">{html.escape(pdb_text)}</script>'
            f'<script>{script}</script>')


def _3dmol_tag() -> str:
    js = TOOLS / "3Dmol" / "3Dmol-min.js"
    if js.is_file():
        return f'<script src="3Dmol/3Dmol-min.js"></script>'
    return '<script src="https://3dmol.org/build/3Dmol-min.js"></script>'


# --------------------------------------------------------------------------- #
def _sec_target(state: dict) -> str:
    t = state.get("target_prep")
    if not t:
        return ""
    c = t.get("completeness", {})
    p = t.get("pocket", {})
    j = t.get("judgment", {})
    viewer = _viewer_pdb(t.get("clean_pdb", ""), "target",
                         t.get("ligand_resnames"))
    sev_icon = {"error": "🔴", "warn": "🟠", "info": "⚪"}
    rows = []
    for i in c.get("issues", []):
        flag = ""
        if i.get("auto_fixed"):
            flag = " <span class='small'>(已自动修复)</span>"
        elif i.get("fix"):
            flag = f" <span class='small'>(可用 repair_structure: {_esc(i.get('fix'))})</span>"
        rows.append([f"{sev_icon.get(i.get('severity'), '•')} {_esc(i.get('type'))}{flag}",
                     f"{_esc(i.get('detail'))}<br>"
                     f"<span class='small'>建议: {_esc(i.get('suggestion', ''))}</span>"])
    repaired = c.get("repaired")
    rep_note = ""
    if repaired:
        rep_note = (f"<div class='card'><b>已自动修复:</b> "
                    f"{_esc(', '.join(repaired.get('applied', [])))} "
                    f"(移除 {sum(repaired.get('removed_atoms', {}).values())} 个原子)</div>")
    precheck = (f"<h3>结构预检</h3>"
                + (_table(['问题', '说明'], rows) if rows
                   else "<p class='small'>未检测到结构坑。</p>")
                + rep_note) if c.get("issues") is not None else ""
    return f"""
<h2>1. 靶点准备 <span class="badge">Module A</span></h2>
<div class="card">
<b>判断:</b> {_esc(j.get('action'))}<br>
<span class="small">{_esc(j.get('rationale', ''))}</span>
</div>
{_table(['指标', '值'], [
    ['链', c.get('chains')],
    ['配体', c.get('ligands')],
    ['原子数', c.get('n_atoms')],
    ['分辨率(来自PDB头, 若缺失)', '见原始文件'],
    ['多聚体', c.get('multimer')],
])}
{precheck}
<h3>对接位点</h3>
{_table(['项', '值'], [
    ['方法', p.get('method')],
    ['中心 (A)', p.get('center')],
    ['盒子 (A)', f"{p.get('xsize')} x {p.get('ysize')} x {p.get('zsize')}"],
    ['依据', p.get('rationale', '')],
])}
{viewer}
"""


def _sec_screen(state: dict, static: bool) -> str:
    s = state.get("screening")
    if not s:
        return ""
    import plotly.graph_objects as go
    from ..utils import jload
    csv_path = s.get("results_csv")
    scores = []
    if csv_path and Path(csv_path).exists():
        import pandas as pd
        df = pd.read_csv(csv_path)
        col = "final_score" if "final_score" in df.columns else "score"
        scores = df[col].dropna().tolist()
    fig = go.Figure()
    if scores:
        fig.add_trace(go.Histogram(x=scores, nbinsx=80,
                                   name="Vina 打分分布"))
    fig.update_layout(title="对接打分分布", height=320, template="plotly_white")
    chart = _plot_png_b64(fig) if static else _plot_div(fig)

    h = s.get("hit_decision", {})
    rows = []
    for hit in s.get("hits", [])[:20]:
        rows.append([hit["rank"], hit["smiles"][:50],
                     f"{hit.get('vina_score', float('nan')):.2f}",
                     f"{hit.get('gnina_score', float('nan')):.2f}"
                     if hit.get("gnina_score") is not None else "-",
                     f"{hit.get('final_score', float('nan')):.2f}"])
    return f"""
<h2>2. 小分子虚拟筛选 <span class="badge">Module B</span></h2>
<div class="card">
库: <code>{_esc(s.get('library'))}</code><br>
标准化 {s.get('n_standardized')} → 过滤后 {s.get('n_after_filter')} → 对接 {s.get('n_docked')}<br>
参考配体重对接打分: <b>{_esc(s.get('reference_ligand_score'))}</b><br>
命中标准: <b>threshold = {_esc(h.get('threshold'))}</b>, n_hits = {_esc(h.get('n_hits'))}
<span class="small">({_esc(h.get('rationale', ''))})</span>
</div>
{f'<img src="data:image/png;base64,{chart}" style="width:100%">' if static else chart}
<h3>Top 命中</h3>
{_table(['排名', 'SMILES', 'Vina', 'GNINA', 'final'], rows)}
"""


def _sec_binder(state: dict, static: bool) -> str:
    b = state.get("binder")
    if not b:
        return ""
    import plotly.graph_objects as go
    names = [d["design"].split("/")[-1][:14] for d in b.get("designs", [])]
    ifv = [d.get("interface_plddt_mean") or 0 for d in b.get("designs", [])]
    mfv = [d.get("mono_plddt") or 0 for d in b.get("designs", [])]
    fig = go.Figure()
    fig.add_trace(go.Bar(x=names, y=mfv, name="单体 pLDDT"))
    fig.add_trace(go.Bar(x=names, y=ifv, name="界面 pLDDT(均值)"))
    fig.update_layout(height=320, template="plotly_white",
                      title="binder 设计 pLDDT")
    chart = _plot_png_b64(fig) if static else _plot_div(fig)
    rows = []
    for i, d in enumerate(b.get("designs", [])[:10]):
        rows.append([
            i + 1, d["design"].split("/")[-1],
            f"{d.get('mono_plddt', float('nan')):.1f}",
            f"{d.get('complex_plddt', float('nan')):.1f}",
            f"{(d.get('interface_plddt_mean') or float('nan')):.1f}",
            f"{(d.get('interface_plddt_min') or float('nan')):.1f}",
            (d.get("seqs") or ["-"])[0][:40] if d.get("seqs") else "-",
        ])
    viewer = ""
    if b.get("designs"):
        viewer = _viewer_pdb(b["designs"][0]["design"], "binder_best")
    return f"""
<h2>3. 蛋白 Binder 设计 <span class="badge">Module C</span></h2>
<div class="card">
类型: <b>{_esc(b.get('binder_type'))}</b> 长度范围 {b.get('length_range')} 共 {b.get('n_designs')} 个设计<br>
热点残基: <code>{_esc(b.get('hotspots'))}</code>
</div>
{f'<img src="data:image/png;base64,{chart}" style="width:100%">' if static else chart}
{_table(['#', '设计', '单体pLDDT', '复合物pLDDT', '界面pLDDT均值', '界面pLDDT最小', 'MPNN序列(top1)'], rows)}
{viewer}
"""


def _sec_vhh(state: dict, static: bool) -> str:
    v = state.get("vhh")
    if not v:
        return ""
    a = v.get("track_a", {})
    b = v.get("track_b", {})
    rows = []
    for i, c in enumerate(v.get("ranked", [])[:15]):
        # R12/G10-v3: CDR-fragment docking details (fast mode)
        fr = ""
        if c.get("n_fragments"):
            fs = c.get("fragment_scores") or []
            best = min(fs) if fs else None
            fr = f"{c['n_fragments']} 片段" + (
                f", 最佳 {best:.2g}" if best is not None else "")
        rows.append([
            i + 1, c["source"],
            f"{c.get('plddt', float('nan')):.1f}" if c.get("plddt") else "-",
            f"{c.get('docking_score', float('nan')):.2f}"
            if c.get("docking_score") is not None else "-",
            fr,
            f"{c.get('interface_plddt_mean', float('nan')):.1f}"
            if c.get("interface_plddt_mean") else "-",
            f"{c.get('composite_score', float('nan')):.3f}",
        ])
    llm = v.get("llm_pick")
    llm_note = f"<div class='card'>LLM 推荐: 候选 <b>{llm+1}</b> — {_esc(v.get('llm_rationale', ''))}</div>" if llm is not None else ""
    viewer = ""
    de_novo = [c for c in v.get("ranked", []) if c["source"] == "de_novo"]
    if de_novo and de_novo[0].get("design"):
        viewer = _viewer_pdb(de_novo[0]["design"], "vhh_best")
    # R13: pLDDT distribution of the modeled library (G9 threshold
    # visualization) — the pass rate at the gate is visible at a glance
    hist = ""
    vals = a.get("plddt_all") or []
    if vals:
        from ..config import resolve_defaults
        import plotly.graph_objects as go
        d = resolve_defaults(state.get("options") or {})
        thr = float(d.vhh_plddt_min)
        fig = go.Figure()
        fig.add_trace(go.Histogram(x=vals, nbinsx=20,
                                   marker_color="#4C78A8", name="pLDDT"))
        fig.add_vline(x=thr, line_dash="dash", line_color="#E45756",
                      annotation_text=f"门槛 {thr:g}")
        fig.update_layout(title=(f"VHH 文库 pLDDT 分布 (n={len(vals)}, "
                                 f"fast/full 门槛 {thr:g})"),
                          height=300, template="plotly_white",
                          xaxis_title="pLDDT", yaxis_title="条数")
        hist = _plot_png_b64(fig) if static else _plot_div(fig)
    return f"""
<h2>4. 纳米抗体 (VHH) 筛选 + 设计 <span class="badge">Module D</span></h2>
<div class="card">
Track A (文库筛选): 文库 {a.get('n_library')} → 建模 {a.get('n_modeled')} → 对接 {a.get('n_docked')}<br>
Track B (de novo 设计): {b.get('n_designs')} 个 (VHH scaffold 约束)
</div>
{hist}
{llm_note}
{_table(['#', '来源', 'pLDDT', '对接打分', 'CDR片段', '界面pLDDT', '综合分'], rows)}
{viewer}
"""


def _sec_md(state: dict, static: bool) -> str:
    m = state.get("md")
    if not m:
        return ""
    import plotly.graph_objects as go
    s = m.get("summary", {})
    fig = go.Figure()
    # gmx rms / gyrate write nm; display RMSD in Å, Rg in nm
    if s.get("rmsd", {}).get("mean"):
        x = list(range(len(s["rmsd"]["mean"])))
        fig.add_trace(go.Scatter(x=x, y=[v * 10.0 for v in s["rmsd"]["mean"]],
                                 name="RMSD 均值 (Å)", mode="lines"))
    if s.get("rg", {}).get("mean"):
        x = list(range(len(s["rg"]["mean"])))
        fig.add_trace(go.Scatter(x=x, y=s["rg"]["mean"], name="Rg 均值 (nm)", mode="lines"))
    # R1: per-chain RMSD (nm -> A)
    for k in sorted(s):
        if k.startswith("rmsd_chain") and isinstance(s.get(k), dict) \
                and s[k].get("mean"):
            x = list(range(len(s[k]["mean"])))
            fig.add_trace(go.Scatter(x=x,
                                     y=[v * 10.0 for v in s[k]["mean"]],
                                     name=f"{k.replace('rmsd_', '')} RMSD (Å, 相对最大链)",
                                     mode="lines"))
    # R12/R1: per-domain RMSD (each domain self-fit; nm -> A) +
    # R13/R1-v2: domain vs rest of protein (hinge/allosteric signal)
    dom_rmsd = s.get("domain_rmsd") or {}
    if dom_rmsd:
        for name in sorted(dom_rmsd):
            ser = dom_rmsd[name].get("series") or []
            if ser:
                fig.add_trace(go.Scatter(
                    x=list(range(len(ser))), y=[v * 10.0 for v in ser],
                    name=f"结构域 {name} 自拟合 RMSD (Å)", mode="lines"))
    dom_vs = s.get("domain_rmsd_vs_rest") or {}
    if dom_vs:
        for name in sorted(dom_vs):
            ser = dom_vs[name].get("series") or []
            if ser:
                fig.add_trace(go.Scatter(
                    x=list(range(len(ser))), y=[v * 10.0 for v in ser],
                    name=f"结构域 {name} 相对其余 (Å)", mode="lines",
                    line={"dash": "dot"}))
    fig.update_layout(title=f"MD 稳定性 ({m.get('reps', 1)} 次重复平均)",
                      height=320, template="plotly_white")
    chart = _plot_png_b64(fig) if static else _plot_div(fig)

    fig2 = go.Figure()
    if s.get("rmsf_profile_mean"):
        # gmx rmsf -res writes nm; display in Å (R10: was mislabeled)
        rmsf = [v * 10.0 for v in s["rmsf_profile_mean"]]
        fig2.add_trace(go.Scatter(x=list(range(len(rmsf))), y=rmsf,
                                  mode="lines", name="RMSF (Å)"))
        fig2.update_layout(title="每残基 RMSF (均值)", height=320,
                           template="plotly_white")
        chart2 = f'<img src="data:image/png;base64,{_plot_png_b64(fig2)}" style="width:100%">' if static else _plot_div(fig2)
    else:
        chart2 = ""

    # R1: secondary-structure persistence figure
    chart3 = ""
    ss_text = ""
    if s.get("ss_frac_mean"):
        fig3 = go.Figure()
        x3 = list(range(len(s["ss_frac_mean"])))
        fig3.add_trace(go.Scatter(x=x3,
                                  y=[v * 100.0 for v in s["ss_frac_mean"]],
                                  mode="lines", name="结构化残基占比 (%)"))
        fig3.update_layout(title="二级结构稳定性 (DSSP 式氢键分类, H/G/I/E/B)",
                           height=280, template="plotly_white")
        chart3 = (f'<img src="data:image/png;base64,{_plot_png_b64(fig3)}" '
                  f'style="width:100%">' if static else _plot_div(fig3))
        if s.get("initial_ss_mean") is not None:
            ss_text = (f"二级结构占比: 起始 <b>{_fmt(s['initial_ss_mean'] * 100, 1)}%</b> → "
                       f"最终 <b>{_fmt(s.get('final_ss_mean', 0) * 100, 1)}%</b> "
                       f"(最低 {_fmt(min(s['ss_frac_mean']) * 100, 1)}%)<br>")
    fr = s.get("flexible_regions") or []
    fr_text = ""
    if fr:
        top = sorted(fr, key=lambda r: -r["mean_rmsf_nm"])[:3]
        fr_text = ("柔性区: " + "; ".join(
            f"残基 {r['res_start']}-{r['res_end']} "
            f"(平均 {_fmt(r['mean_rmsf_nm'] * 10, 1)} Å)" for r in top)
            + "<br>")
    interp = s.get("interpretation") or []
    interp_html = ""
    if interp:
        lis = "".join(f"<li>{_esc(t)}</li>" for t in interp)
        interp_html = (f"<div class='card'><b>柔性解读 (规则引擎自动判断)</b>"
                       f"<ul>{lis}</ul></div>")
    sysinfo = m.get("system", {})
    # R11: per-replica cluster populations live in summary.replicas (the
    # run records in m.replicas carry only rep/dir/wall_h)
    clus = []
    per_rep = s.get("replicas") or []
    for r in per_rep:
        c = r.get("clusters") or {}
        parts = ", ".join(f"cluster {k}: {v * 100:.1f}%"
                          for k, v in sorted(c.items(), key=lambda kv: -kv[1]))
        clus.append(f"rep{r.get('rep', '?')}: " + (parts or "无聚类数据"))
    if not per_rep and s.get("clusters"):
        c = s["clusters"]
        parts = ", ".join(f"cluster {k}: {v * 100:.1f}%"
                          for k, v in sorted(c.items(), key=lambda kv: -kv[1]))
        clus.append("均值: " + parts)
    clus_text = chr(10).join(clus) if clus else "无聚类数据"
    # R12/R1: structural-domain table (residue ranges + RMSD in A);
    # R13: + vs-rest (hinge) column
    dom_rows = []
    dom_vs = s.get("domain_rmsd_vs_rest") or {}
    for d in s.get("domains") or []:
        st = (s.get("domain_rmsd") or {}).get(d["name"]) or {}
        vs = dom_vs.get(d["name"]) or {}
        dom_rows.append([
            d.get("name", "?"),
            f"{d.get('res_start', '?')}-{d.get('res_end', '?')}",
            d.get("n_res", "?"),
            f"{st['final'] * 10:.1f}" if st.get("final") is not None else "-",
            f"{st['mean'] * 10:.1f}" if st.get("mean") is not None else "-",
            f"{vs['final'] * 10:.1f}" if vs.get("final") is not None else "-",
        ])
    dom_html = (_table(["结构域", "残基", "残基数", "末端自拟合 (Å)",
                        "均值自拟合 (Å)", "末端相对其余 (Å)"], dom_rows)
                if dom_rows else "")
    return f"""
<h2>5. MD 模拟 <span class="badge">Module E</span></h2>
<div class="card">
体系: <b>{_esc(sysinfo.get('label'))}</b><br>
引擎: GROMACS {_esc(m.get('gromacs', {}).get('version'))} 力场 {_esc(m.get('gromacs', {}).get('ff'))}<br>
时长: {m.get('ns')} ns x {m.get('reps')} 次重复<br>
最终 RMSD (均值) = <b>{_fmt(s.get('final_rmsd_mean') * 10, 2) if s.get('final_rmsd_mean') else '?'} Å</b>, 最终 Rg = <b>{_fmt(s.get('final_rg_mean'), 3) if s.get('final_rg_mean') else '?'} nm</b><br>
{ss_text}{fr_text}</div>
{f'<img src="data:image/png;base64,{chart}" style="width:100%">' if static else chart}
{chart2}
{chart3}
{interp_html}
{dom_html}
<h3>聚类 (gmx cluster / GromOS, 截断 1.5 nm)</h3>
<pre>{_esc(clus_text)}</pre>
"""


def _sec_decisions(state: dict) -> str:
    pdir = Path(state.get("project_dir", "."))
    log = pdir / "decisions.json"
    if not log.exists():
        return ""
    rows = []
    for d in json.loads(log.read_text()):
        rows.append([d["node"], d["question"][:90],
                     str(d["answer"])[:120], d.get("rationale", "")[:160]])
    if not rows:
        return ""
    return f"""
<h2>6. Agent 判断日志</h2>
<p class="small">LLM: {_esc(state.get('options', {}).get('llm_model', 'qwen'))} — 每次关键判断的依据都记录如下</p>
{_table(['节点', '问题', '决定', '依据'], rows)}
"""


def build_report(state: dict) -> dict:
    pdir = Path(state["project_dir"])
    rdir = pdir / "reports"
    rdir.mkdir(parents=True, exist_ok=True)
    (rdir / "3Dmol").mkdir(exist_ok=True)
    # bundle 3dmol if present
    src3d = TOOLS / "3Dmol" / "3Dmol-min.js"
    if src3d.is_file():
        (rdir / "3Dmol" / "3Dmol-min.js").write_bytes(src3d.read_bytes())

    opts = state.get("options", {})
    parts = [
        "<!DOCTYPE html><html><head><meta charset='utf-8'>",
        f"<title>DrugAgent 报告 - {pdir.name}</title>",
        f"<style>{CSS}</style>",
        _3dmol_tag(),
        "</head><body>",
        "<h1>DrugAgent 药物筛选报告</h1>",
        f"<p class='small'>生成时间: {datetime.now().isoformat(timespec='seconds')} "
        f"| 模块: {','.join(opts.get('modules', []))} "
        f"| LLM: {opts.get('llm_model', 'qwen3.8-27b')}"
        f" | 快速模式: {opts.get('fast', False)}</p>",
    ]
    parts.append(_sec_target(state))
    parts.append(_sec_screen(state, static=False))
    parts.append(_sec_binder(state, static=False))
    parts.append(_sec_vhh(state, static=False))
    parts.append(_sec_md(state, static=False))
    parts.append(_sec_decisions(state))
    parts.append("</body></html>")
    html_path = rdir / "report.html"
    html_path.write_text("".join(parts))

    # PDF via WeasyPrint (static charts)
    pdf_path = rdir / "report.pdf"
    try:
        parts2 = [
            "<!DOCTYPE html><html><head><meta charset='utf-8'>",
            f"<style>{CSS}</style></head><body>",
            "<h1>DrugAgent 药物筛选报告 (PDF)</h1>",
        ]
        parts2.append(_sec_target(state))
        parts2.append(_sec_screen(state, static=True))
        parts2.append(_sec_binder(state, static=True))
        parts2.append(_sec_vhh(state, static=True))
        parts2.append(_sec_md(state, static=True))
        parts2.append(_sec_decisions(state))
        parts2.append("</body></html>")
        from weasyprint import HTML
        HTML(string="".join(parts2), base_url=str(rdir)).write_pdf(str(pdf_path))
    except Exception as e:  # noqa: BLE001
        logger.warning(f"PDF generation failed: {e}")
        pdf_path = None

    out = {"html": str(html_path), "pdf": str(pdf_path) if pdf_path else None}
    return out
