# DrugAgent 迭代轮次记录（循环：计划 → 检索 → TDD → e2e → 文档 → 反思）

> 每轮一节。缺口编号沿用用户优先级：R1 MD 柔性诊断 / R2 受体构象选择+柔性对接 /
> R3 MD 金属离子 / R4 核酸力场 / R5 柔性靶点系综。代码内另有 R6-R10 为历史轮次
> 引入的编号（收敛自动延长 / 构象池 / 辅因子 / MDP 指纹 / 离子），勿混淆。

## 轮次状态总览

| 缺口 | 状态 | 位置 |
|---|---|---|
| R1 分链 RMSD（相对+自拟合） | 已实现 | mdsim.analyze_replicas / report |
| R1 二级结构稳定性 + 自动解读 | 已实现 | mdsim.analyze_ss / interpret_stability |
| R1 区域/域级柔性定位 | **本轮 G3** | mdsim.flexible_regions |
| R2 多构象选择 + 侧链柔性 + consensus | 已实现 | screening.pool_representatives / flex_sidechain_pdbqt / consensus_stats |
| R2 策略引导（何时走柔性靶点工作流） | **本轮 G4** | agent/prompts.py |
| R3 金属离子配位（pair restraints） | 已实现 | mdsim.find_metals / _integrate_metals |
| R4 核酸原生参数化（pdb2gmx） | 已实现 | mdsim（NA 链分流） |
| R5 平衡段（NVT→NPT）+ 收敛自动延长 | 已实现 | mdsim.run_equilibration / _md_converged / extend_replicas |
| R5 多构象系综 → 对接 | 已实现（工具层） | tools_screen.dock_conformer_set |

## 第 10 轮（本轮）

### 计划

基线：快速套件仅 1 个失败 = `test_md_auto_extend`（docstring 标 @slow 但缺
`pytest.mark.slow`，进了快速套件；失败本身是 /tmp tmpfs 后台清理导致的偶发，
本地盘重跑通过）。

| # | 目标 | 做法 |
|---|---|---|
| G1 | run_cmd 错误可见性（易用性+agent 自排错） | 带 log_file 的命令失败时，RuntimeError 附日志尾部（~15 行），agent 不用二次 read_file 就能看到 gmx 报错 |
| G2 | 测试卫生 | test_md_auto_extend 补 @pytest.mark.slow |
| G3 | R1 收尾：区域级柔性定位 | `flexible_regions(rmsf_profile)`：滑窗找连续高 RMSF 区段（≥5 残基、> 全局均值 2×）→ 残基区间+均值；进 summary/interpretation/md_summary/report；顺修 gmx_analyze 只认 Protein_chain 的 DNA 盲区 |
| G4 | R5 策略进系统提示词 | "柔性靶点工作流"段：短 MD 系综 → make_flex_receptor → dock_conformer_set（consensus 打分）→ 设计/全长 MD |
| G5 | CLI 易用性 | `run` 暴露关键 MD 旋钮：--md-salt / --md-divalent / --md-divalent-m / --md-extend-ns / --md-max-extensions / --md-burn-in-ps |
| G6 | 小分子库回退（e2e 时发现） | dtp.sdf.gz 是 0 字节（DTP 镜像不稳，历史问题）。`resolve_library()`：缺失/<1MB → 自动回退 chembl35_small (50k) 并在 state/报告里标明 "fallback for dtp"；自定义 SDF 路径不回退。顺修 RMSF 单位（gmx rmsf 是 nm，报告原标 Å）与 gmx_analyze 只认 Protein_chain 的 DNA 盲区 |

### 验证
- [ ] 快速套件全绿（--basetemp 本地盘）
- [ ] e2e：`run --target 1HVI --modules all --fast --auto --no-llm`（脚本模式确定性）
- [ ] 文档：README/DESIGN 更新 + 本节结果

### 结果
（待填）

### 反思 / 下轮缺口
（待填）
