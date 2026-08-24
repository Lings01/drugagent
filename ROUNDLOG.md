# DrugAgent 迭代轮次记录（循环：计划 → 检索 → TDD → e2e → 文档 → 反思）

> 每轮一节。缺口编号沿用用户优先级：R1 MD 柔性诊断 / R2 受体构象选择+柔性对接 /
> R3 MD 金属离子 / R4 核酸力场 / R5 柔性靶点系综。代码内另有 R6-R10 为历史轮次
> 引入的编号（收敛自动延长 / 构象池 / 辅因子 / MDP 指纹 / 离子），勿混淆。

## 轮次状态总览

| 缺口 | 状态 | 位置 |
|---|---|---|
| R1 分链 RMSD（相对+自拟合） | 已实现 | mdsim.analyze_replicas / report |
| R1 二级结构稳定性 + 自动解读 | 已实现 | mdsim.analyze_ss / interpret_stability |
| R1 区域级柔性定位（残基区间） | 已实现（第 10 轮 G3） | mdsim.flexible_regions / interpret / report |
| R2 多构象选择 + 侧链柔性 + consensus | 已实现 | screening.pool_representatives / flex_sidechain_pdbqt / consensus_stats |
| R2 策略引导（何时走柔性靶点工作流） | 已实现（第 10 轮 G4） | agent/prompts.py |
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
| G7 | VHH 对接并行化（e2e 时发现） | Module D 的对接是串行 for 循环：全长 VHH 当 ~700 原子"配体"，单个 vina 即使 exh=1 也要分钟级 → 100 个候选 ≈ 数小时，历史上 VHH 阶段从未完整跑完。`dock_vhh_candidates()` 用 pmap 并行（16 路 × os_cpu/16 核），fast 模式 vhh_screen_n 100→40 |

### 验证
- [x] 快速套件全绿（--basetemp 本地盘）：147 通过（162 含 15 slow）
- [x] e2e：`run --target 1HVI --modules all --fast --auto --no-llm`（projects/r10_e2e，
  收尾时进行到 VHH 对接阶段；target/screening/binder 均 ok，G6 回退与 G7 并行化
  已在真实流程中触发验证）
- [x] 文档：README（模块表/参数/排错提示/测试数）+ 本文件 + git 仓库初始化

### 结果
1. **G1 错误可见性**：带 log_file 的命令失败时 RuntimeError 附日志尾部（40 行/2KB 上限）。
   基线失败（gmx solvate "File does not exist"）正是这种"看不见原因"的典型。
2. **G2**：test_md_auto_extend 补 @pytest.mark.slow（之前混进快速套件且偶发受
   /tmp tmpfs 清理影响）。
3. **G3 R1 收尾**：`flexible_regions()` 从均值 RMSF 剖面找连续柔性段
   （阈值 max(2×mean, mean+2σ, 0.1nm)，≥5 残基），进 summary/interpretation/
   md_summary/报告（"柔性区: 残基 X-Y (平均 Z Å)"）；interpretation 提示
   "若在口袋内建议柔性靶点工作流"。顺修：报告 RMSF 图单位（nm→Å 显示）、
   gmx_analyze 的 DNA_chain 盲区（tools_md）。
4. **G4 R5 策略**：系统提示词新增"柔性靶点策略"段（短 MD 系综 →
   make_flex_receptor → dock_conformer_set → consensus 判命中），并修正
   rmsf 单位表述（xvg 都是 nm）。
5. **G5 CLI**：--md-salt / --md-divalent / --md-divalent-m / --md-extend-ns /
   --md-max-extensions / --md-burn-in-ps 透传到 resolve_defaults。
6. **G6 库回退**：dtp.sdf.gz 是 0 字节（dtpbase.org 镜像超时，实测 curl 000）→
   `resolve_library()` 回退 chembl35_small.sdf (50k, 132MB) 并在 state/报告标明
   "fallback for dtp"；e2e 中真实触发并跑通 50k 标准化 + ~3600 对接（13 min）。
7. **G7 VHH 并行**：`dock_vhh_candidates()` pmap 并行（≤16 路，os_cpu/路），
   fast vhh_screen_n 100→40；单测验证 4 worker 真并行 + pdbqt 幂等复用。
   注意：loky 子进程不继承父进程 monkeypatch（单测用 n_jobs=1 走进程内路径）。
8. **仓库**：git 初始化 + .gitignore（env/data/projects/日志/gmx 备份文件）。

### 反思 / 下轮缺口（按优先级）
1. **G8 阶段级幂等**：`run_screening`/`run_vhh`/`run_design` 等 run_* 工具
   总是整段重跑（本轮重启 e2e 重跑了 35 min 筛选）。应在 stage json 完整时
   直接复用（加 force 参数），mdrun 级幂等已有。
2. **G9 VHH 命中面太窄**：40 条合成序列仅 1 条过 pLDDT>45（fast 阈值），
   对接样本=1，Module D 统计意义弱。考虑：fast 阈值再放宽 + 提高 n，或对
   合成库先按序列多样性分层再建模 top-N。
3. **G10 单 VHH 对接仍慢**：全长 VHH（~700 原子、数百扭转）exh=1 也要
   ~15-30 min；实测 `vina --cpu 64` 只跑 1 核（该 AD4 版 vina 对大配体
   并行失效）→ 并行化靠"多配体"而非"单配体多线程"。可试：exh 降到 0.5/
   截断框架只柔 CDR、或 GNINA 替代复打分。
4. **R1 域级 RMSD**：当前"分链 + 柔性区"已覆盖大部分需求；若需真正的
   结构域划分（DSSP 域边界/动态域检测 NMDYN 式），下轮可加
   gmx_analyze(kind="domain")。
5. **R5 多构象系综的自动触发**：目前靠提示词引导 agent 选择柔性工作流；
   可加自动判据（靶点无配体 + RMSF 高 → 建议短 MD 系综）进
   analyze_pdb/repair 的 issues。
6. **status 命令增强**：显示各阶段产物与错误摘要（目前只有 ✓/✗ + 最后一条
   transcript）。
7. DTP 下载重试（dtpbase.org 偶活）；pdbbind.tar.gz 138 字节也是坏的。
