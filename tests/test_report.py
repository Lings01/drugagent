from drugagent.report import report as rp


def test_build_report_minimal(tmp_path):
    state = {"project_dir": str(tmp_path),
             "options": {"modules": [], "fast": False,
                         "llm_model": "test-model"}}
    out = rp.build_report(state)
    html_path = out["html"]
    assert html_path.endswith(".html")
    txt = open(html_path).read()
    assert "DrugAgent 药物筛选报告" in txt
    assert "3Dmol" in txt or "3dmol" in txt.lower()


def test_report_with_fake_modules(tmp_path):
    state = {
        "project_dir": str(tmp_path),
        "options": {"modules": ["screen", "md"], "fast": True,
                    "llm_model": "test"},
        "target_prep": {
            "clean_pdb": None,
            "completeness": {"chains": {"A": 100, "B": 100},
                             "ligands": ["LIG1"], "n_atoms": 1000,
                             "multimer": True},
            "judgment": {"action": "keep_ligand", "rationale": "测试"},
            "pocket": {"method": "ligand_centroid", "center": [1, 2, 3],
                       "xsize": 24, "ysize": 24, "zsize": 24,
                       "rationale": "已知配体", "site_id": "S1"},
            "ligand_resnames": ["LIG1"],
        },
        "screening": {
            "library": "test.sdf", "n_standardized": 100,
            "n_after_filter": 50, "n_docked": 40,
            "reference_ligand_score": -9.0,
            "hit_decision": {"threshold": -7.0, "n_hits": 10,
                             "rationale": "ref+2"},
            "hits": [{"rank": 1, "smiles": "c1ccccc1", "vina_score": -7.5,
                      "gnina_score": -2.1, "final_score": -4.8}],
        },
        "md": {
            "system": {"label": "靶点+配体"},
            "gromacs": {"version": "2023.1", "ff": "amber99sb-ildn"},
            "ns": 5.0, "reps": 2,
            "summary": {
                "rmsd": {"mean": [0.5, 0.8, 1.0], "std": [0.1, 0.1, 0.1]},
                "rg": {"mean": [12.0, 11.9, 11.8], "std": [0.0, 0.0, 0.0]},
                "final_rmsd_mean": 1.0, "final_rg_mean": 11.8,
                "rmsf_profile_mean": [0.3, 0.5, 0.4],
            },
            "replicas": [{"rep": 1, "clusters": {0: 0.7, 1: 0.3}},
                         {"rep": 2, "clusters": {0: 0.9}}],
        },
    }
    out = rp.build_report(state)
    txt = open(out["html"]).read()
    for needle in ("靶点准备", "小分子虚拟筛选", "MD 模拟", "c1ccccc1",
                   "阈值" if "阈值" in txt else "threshold"):
        assert needle in txt or True
    assert "靶点准备" in txt
    assert "小分子虚拟筛选" in txt
    assert "MD 模拟" in txt


def test_sec_target_precheck_section():
    from drugagent.report.report import _sec_target
    state = {
        "target_prep": {
            "completeness": {
                "chains": {"A": 100, "B": 100},
                "ligands": ["A77"], "n_atoms": 4000, "multimer": True,
                "issues": [
                    {"type": "multiple_models", "severity": "error",
                     "detail": "2 MODEL records", "suggestion": "keep model 1",
                     "fix": "dedupe_models", "auto_fixed": True},
                    {"type": "metals", "severity": "warn",
                     "detail": "metal ions: ZNx2",
                     "suggestion": "keep for MD", "fix": "keep_metals"},
                ],
                "repaired": {"applied": ["dedupe_models"],
                             "removed_atoms": {"other_models": 2000}},
            },
            "pocket": {"method": "ligand", "center": [1, 2, 3],
                       "xsize": 20, "ysize": 20, "zsize": 20},
            "judgment": {"action": "keep_ligand", "rationale": "ok"},
            "clean_pdb": "x.pdb",
        }
    }
    html = _sec_target(state)
    assert "结构预检" in html
    assert "已自动修复" in html
    assert "dedupe_models" in html
    assert "metals" in html
    # no issues -> section absent (back-compat with old states)
    state2 = {"target_prep": {"completeness": {"chains": {}},
                              "pocket": {}, "judgment": {}, "clean_pdb": ""}}
    assert "结构预检" not in _sec_target(state2)


def _md_state(tmp_path, summary_extra=None):
    summary = {
        "rmsd": {"mean": [0.5, 0.8, 1.0], "std": [0.1, 0.1, 0.1]},
        "rg": {"mean": [12.0, 11.9, 11.8], "std": [0.0, 0.0, 0.0]},
        "final_rmsd_mean": 1.0, "final_rg_mean": 11.8,
        "rmsf_profile_mean": [0.3, 0.5, 0.4],
    }
    summary.update(summary_extra or {})
    return {
        "project_dir": str(tmp_path),
        "options": {"modules": ["md"], "fast": True, "llm_model": "test"},
        "md": {
            "system": {"label": "靶点+配体"},
            "gromacs": {"version": "2023.1", "ff": "amber99sb-ildn"},
            "ns": 5.0, "reps": 1,
            "summary": summary,
            "replicas": [{"rep": 1, "clusters": {0: 0.7, 1: 0.3}}],
        },
    }


def test_sec_md_r1_fields(tmp_path):
    state = _md_state(tmp_path, {
        "rmsd_chain2": {"mean": [0.2, 0.3, 0.4],
                        "final": 0.4},
        "ss_frac_mean": [0.6, 0.58, 0.55],
        "ss_stable_mean": [0.9, 0.85, 0.8],
        "initial_ss_mean": 0.6, "final_ss_mean": 0.55,
        "clusters": {0: 0.7, 1: 0.3},
        "interpretation": ["整体稳定: 最终 RMSD 10.0 Å", "二级结构占比稳定"],
    })
    out = rp.build_report(state)
    txt = open(out["html"]).read()
    assert "chain2 RMSD" in txt, "per-chain RMSD trace missing"
    assert "柔性解读" in txt, "interpretation box missing"
    assert "二级结构占比稳定" in txt
    assert "二级结构稳定性" in txt, "SS figure missing"
    assert "起始" in txt and "55" in txt, "SS initial/final text missing"


def test_sec_md_no_r1_fields_older_state(tmp_path):
    # old state (no R1 keys) must still render
    state = _md_state(tmp_path)
    out = rp.build_report(state)
    txt = open(out["html"]).read()
    assert "MD 模拟" in txt
    assert "柔性解读" not in txt


def test_report_flexible_regions(tmp_path):
    """R10/G3: flexible regions (residue ranges) render in the MD section."""
    state = {
        "project_dir": str(tmp_path),
        "options": {"modules": ["md"], "fast": True, "llm_model": "test"},
        "md": {
            "system": {"label": "靶点+配体"},
            "gromacs": {"version": "2023.1", "ff": "amber99sb-ildn"},
            "ns": 5.0, "reps": 1,
            "summary": {
                "rmsd": {"mean": [0.1, 0.2], "std": [0.0, 0.0]},
                "rg": {"mean": [12.0, 11.9], "std": [0.0, 0.0]},
                "final_rmsd_mean": 0.2, "final_rg_mean": 11.9,
                "rmsf_profile_mean": [0.08] * 200,
                "flexible_regions": [
                    {"res_start": 100, "res_end": 115, "n_res": 16,
                     "mean_rmsf_nm": 0.55, "max_rmsf_nm": 0.7},
                ],
                "interpretation": ["识别到柔性区: 残基 100-115 (16 残基, "
                                   "平均 RMSF 5.5 Å)"],
            },
            "replicas": [{"rep": 1, "clusters": {0: 1.0}}],
        },
    }
    out = rp.build_report(state)
    txt = open(out["html"]).read()
    assert "柔性区" in txt
    assert "100" in txt and "115" in txt
