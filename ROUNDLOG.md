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
| R5 柔性工作流自动判据（apo+高RMSF→建议系综） | 已实现（第 11 轮） | target_prep.analyze_completeness (apo_target issue) / mdsim.interpret_stability |
| G8 阶段级幂等（run_* 完成即跳过 + force） | 已实现（第 11 轮） | agent/stages.py |
| G9 VHH 命中面（pLDDT 门槛 35/50 + fast n=80） | 已实现（第 11 轮） | config.vhh_plddt_min / modules.vhh.screen_vhh |
| G10 VHH 刚性对接（TORSDOF 0 + 缓存一致性） | 已实现（第 11 轮） | target_prep.make_rigid_pdbqt / modules.vhh.dock_vhh_candidates |
| G10-v2 VHH CDR 片段对接（fast 默认，综合分=最佳片段） | 已实现（第 11 轮） | modules.vhh.vhh_cdr_fragments / dock_vhh_candidates(cdr_only) |
| G10-v3 片段拆分 + 自适应盒子（撞墙罚分消失） | 已实现（第 12 轮） | modules.vhh._frag_box / vhh_cdr_fragments(max_res) |
| R1 真·结构域 RMSD（SS 域提取 + 域自拟合 Kabsch） | 已实现（第 12 轮） | mdsim.find_ss_domains / domain_rmsd_series / analyze_ss / interpret 注 8 |
| G8 收尾：rerun --stage X（force CLI 化） | 已实现（第 12 轮） | cli.rerun |
| R1 纯蛋白靶点 apo MD 通路（select_system/build apo 分支） | 已实现（第 16 轮） | mdsim.select_system / build_complex_pdb |
| R1 scaffold 保真度字段 scaffold_rmsd_a（Kabsch Å） | 已实现（第 16 轮） | vhh.scaffold_fidelity / design_vhh |
| G8 对称：binder_rf_cautious 选项 | 已实现（第 16 轮） | binder.rfdesign(rf_cautious) |
| R1 单位哨兵（MDA vs gmx self-fit 比值带） | 已实现（第 16 轮） | tests/test_mdsim.py::test_md_unit_sentinel_mda_vs_gmx |
| R1 跨靶点刚性基线（ubq apo 2ns×3） | 已实现（第 16 轮，数据） | projects/rigid_ubq + ROUNDLOG 标定表 |
| R1 gmx 分析 compact-unwrap（紧凑盒 flapping 修复） | 已实现（第 17 轮） | mdsim._ensure_unwrapped / analyze_replicas |
| R1 刚性基线判据（√t 折算 + 域级 norm 0.08） | 已实现（第 17 轮） | mdsim.interpret_stability 注 1b |
| R4 修饰残基过滤（apo 构建 + 系统选择偏好） | 已实现（第 17 轮） | mdsim._keep_standard_res / _MODIFIED_RES |

## 第 12 轮（本轮）

### 计划

起点：第 11 轮反思优先级：G10 v3 片段打磨 → R1 真结构域 RMSD → rerun CLI →
pLDDT 直方图/片段分进报告 → DTP 重试。bench.json flex 收尾数字（第 11 轮
遗留）本轮先收掉。

| # | 目标 | 做法 |
|---|---|---|
| P5 | bench 收尾 | 修分数解析 bug（误取 VINA RESULT 行尾 rmsd_ub=0.000 → 应取 score 字段）；同条件数字入档 |
| P1 | G10 v3 片段打磨 | (a) 长片段拆分：>20 残基的 run 切成 ≤20 残基块（frag 成本超线性）；(b) **片段自适应盒子**：中心=pocket 中心，边长=片段直径+8 Å（clip 12-30 Å）——片段装进自己盒子，撞墙罚分消失；(c) `n_fragments`/`fragment_scores` 进 merge → 报告新列 |
| P3 | `rerun --stage X` | G8 force 的 CLI 化：`drugagent rerun --project X --stage {target_prep,screening,binder,vhh,md,report} [--no-with-report]`，成功后自动重建报告 |
| P2 | R1 真结构域 RMSD | 此 GROMACS build 无 `doomain` 工具 → 用已有 DSSP-like SS 分类（frame 0）提取**结构域**（≥8 连续结构化残基），每域 CA 原子 Kabsch 自拟合（对 frame 0 自身）→ 每帧域 RMSD 序列。与全局 `rms`（整体一次拟合）不同，域自拟合隔离**域刚体运动**。进 summary（domain_rmsd final/mean/series）+ 解读注 8（域末端 RMSD>4 Å → 铰链/变构域；整体 RMSD 高但各域稳定 → 域间相对运动）+ 报告表格+曲线 |
| P6 | DTP 重试 | dtpbase.org 再试（第 11 轮起偶活）|

P4 pLDDT 直方图：若时间不够留第 13 轮（报告侧纯展示，风险低）。

### 验证
- [x] 快速套件全绿：191（+7：CDR chunk/box 2、no_frag skip 1、domain 3、kabsch 约定 1）
- [x] R1：r10_e2e md_rep1 真实轨迹 domain 分析——5 结构域（8-17 残基），
      末端自拟合 RMSD 3.6-6.7 Å（apo 5ns 合理量级）；排查出 **3 个 PBC 坑**：
      (a) 残基跨盒边界（res24 N 在 x≈0.7、CA/C/O 在 x≈58.6）毒化中心法
      → 两级 make-whole（残基内一致 + 链走）；(b) 残基沿边界 flapping
      （res95 x 在 58.6↔0.1 间往返）MDAnalysis 2.x `bb.unwrap()`
      （compound/COM 法）抓不到 → 手做逐原子 min-image 累积；(c) frame 0
      本身 95/96 跨盒 → 参考帧必须 make-whole。另修 **Kabsch 约定 bug**：
      H=AᵀB 时 R=U D Vᵀ（原 (U D Vᵀ)ᵀ 在准线性点云上碰巧对、良态 3D 云
      全错）——brute-force SO(3) 验证 + 回归测试
- [x] e2e：`rerun --stage vhh`（rerun CLI + G10 v3 CDR 全链路）——9 片段
      全 dock 成功，vhh_30 frag1 = -8.66 kcal/mol，vhh_30 综合分 1.0 居首，
      报告 CDR 片段列渲染 ✓（详见结果 6）
- [x] `rerun --stage md` 重分析（3 副本 ~3 min，gmx 产物全缓存）→ state
      5 结构域 + domain_rmsd（末端 4.5-5.9 Å）+ 解读注 8 触发（dom3 构象
      漂移 5.9 Å）+ 报告域表格渲染 ✓（详见结果 4）
- [x] 文档（README/HANDOFF 已改）

### 结果
1. **P5 bench 收尾**：修分数解析 bug（误取 VINA RESULT 行尾 rmsd_ub=0.000；
   此构建 pose 写到 --out 原路径无 .pdbqt 后缀）→ bench.json 最终：
   rigid 100.96 min / flex 109.73 min（1.10×），**同条件最优分完全相同**
   （214373513.9，差 <0.001 kcal/mol）——全长 VHH 刚柔之分无意义，搜索
   收敛到同一个撞墙 pose。
2. **P1 G10 v3 片段打磨**（三项全上）：
   - **长片段切块**：>20 残基的 low-pLDDT run 自动切 ≤20 残基块
     （vhh_30 的 30 残基 run → 20+10 两块；dock 成本超线性）。
   - **片段自适应盒子**：中心=pocket 中心，边长=片段直径+8 Å（clip
     12-30 Å）。实测 vhh_30 frag0（257 原子、旧 25.84 Å 盒）：**6.25 min
     → 1.3 min，分数 2.3e7 → 145 kcal/mol（~5 个数量级）**——片段装进
     自己盒子后撞墙罚分消失，分数进入可解读量级。
   - **无片段跳过**：`vhh_dock_full_fallback`（fast 默认 False）——CDR 全
     无低 pLDDT 区时不再付 ~100 min 全长 dock（撞墙分无信息量），候选
     score=NaN、pLDDT 仍参与综合排名。full 模式默认 True（保守）。
   - `n_fragments`/`fragment_scores` 进 merge → 报告 VHH 表新列
     "CDR片段"（n 片段 + 最佳分）。
3. **P3 `rerun --stage X`**：`drugagent rerun --project X --stage
   {target_prep,screening,binder,vhh,md,report} [--no-with-report]`——
   G8 force 的 CLI 化，成功后默认重建报告。
4. **P2 R1 真结构域 RMSD**：frame 0 DSSP-like SS 分类提取 ≥8 连续结构化
   残基的结构域（1HVI 二聚体 → 5 域 8-17 残基）；每域 CA Kabsch 自拟合到
   自身 frame 0 → 逐帧域 RMSD（隔离域内部形变）。**r10_e2e md_rep1 实测
   末端自拟合 RMSD 3.6-6.7 Å**（apo 5 ns 合理量级）。排查出 **4 个坑**
   （全部修复+回归测试）：
   (a) **残基跨盒**（res24：N 在 x≈0.7、CA/C/O 在 x≈58.6）毒化中心法
   → frame 0 两级 make-whole（残基内 min-image 一致 + 链走）；
   (b) **边界 flapping**（res95 x 在 58.6↔0.1 往返）——MDAnalysis 2.x
   `bb.unwrap()` 是 compound/COM 法，大蛋白边缘原子抓不到 → 手做逐原子
   min-image 累积；
   (c) **frame 0 本身跨盒**（95/96 键跨 x 边界）参考帧必须 make-whole；
   (d) **Kabsch 约定 bug**：H=AᵀB 时 R=U D Vᵀ（原 (U D Vᵀ)ᵀ 在准线性
   点云上碰巧对、良态 3D 云全错，brute-force SO(3) 验证）。
   解读注 8：域末端自拟合 RMSD>4 Å → "结构域构象漂移"注；整体 RMSD 高
   但各域稳定 → "域间相对运动而非域内展开"。报告：域表格（残基/末端/
   均值 Å）+ 每域 RMSD 曲线。
5. **P6 DTP**：再试仍 SSL EOF（dtpbase.org 持续不稳）——fallback
   chembl35_small 兜底，`setup --libraries dtp` 软失败不阻断。
6. **e2e（`rerun --stage vhh`，一次验证 rerun CLI + G10 v3 全链路）**：
   9 个 CDR 片段全部成功对接（30 Å 自适应盒），**vhh_30 frag1 = -8.66
   kcal/mol**（真正有利的片段结合分！），vhh_30 综合分 1.0 居首（pLDDT
   46.8 + 最佳片段 -8.66）；vhh_60/22/78 最佳片段 126-187；vhh_79 无
   片段 → 旧代码路径走了全长 fallback（4.48e8 撞墙分，新代码会跳过）。
   报告重建：CDR 片段列 ✓。总时长 ~3h（大头 = vhh_79 全长 dock ~2.5h +
   track B 6 设计 × 6 min 无缓存重验证——见反思 2）。
7. 快速套件 191/191 绿（+7：CDR 切块/盒子 2、no_frag 1、domain 3、
   kabsch 约定 1）。

### 反思 / 下轮缺口
1. ✅ G10 v3 验证完成：自适应盒子是本轮最大单项收益（5 个数量级分数
   修正 + 1.3 min/片段）；fragment_scores 进报告后可直接横向比较候选。
2. **track B 设计验证无缓存**：`design_vhh` 对 outdir 里**全部**
   `vhh_design_*.pdb`（含历史轮的）逐个 ESMFold 复杂验证（~6 min/个、
   无缓存判据）——rerun 时 6 个设计 = 36 min 纯重复。下轮：设计验证
   结果（complex_plddt/interface）落 `vhh_designs/scored.json`，存在即
   跳过（G8 同款思路下沉到设计粒度）。
3. **R1-v2 铰链检测**：当前域自拟合 RMSD 度量的是域**内部形变**；
   "域相对其余蛋白的刚体运动"（真铰链/变构）需要"fit 其余蛋白、测该
   域"的相对 RMSD——下轮补 `domain_vs_rest` 指标 + 解读注。
4. **pLDDT 直方图（P4 顺延）**：报告 VHH 节加 pLDDT 分布直方图（G9
   门槛 35/50 可视化）——纯展示，风险低。
5. **MD 分析 rerun 成本**：`rerun --stage md` 会重算全部副本分析
   （gmx 产物有缓存但 SS/域分析 + 聚类每副本 ~1-2 min）——可接受，但
   若下轮 R1-v2 加域-其余相对 RMSD，注意别在 analyze 里再加载一遍
   轨迹（与 SS 共用一次 MDAnalysis 加载）。
6. bench.json 的 flex 数字（109.73 min）是在 8 核 + 其他 7 个并行
   vina 争抢下测的；同条件分数一致已足够支撑结论，计时不必复测。

## 第 11 轮（本轮）

### 计划

起点：第 10 轮反思优先级 G8→G9→G10→R1 域级→R5 自动判据→status 增强→
pdbbind 修复。本轮先确认上一轮 e2e（projects/r10_e2e）：仍在运行
（target/screening/binder 完成，VHH 对接中），收尾检查留到本轮 e2e 验证时做。

| # | 目标 | 做法 |
|---|---|---|
| G8 | 阶段级幂等 | 新 `drugagent/agent/stages.py`：每阶段的"完成判据"（state.json 标记位）+ 缓存摘要；6 个 `run_*`/`build_report` 工具加 `force` 参数，stage json 完整 → 整体跳过。低层产物幂等（pdbqt/mdrun/xvg 缓存）不动 |
| G9 | VHH 命中面 | 实测 100 条已建模 VHH 的 pLDDT 分布：p50=31.5 / p90=34.5 / >45 仅 1 条 → 旧 45 门槛 fast 模式只剩 1 个对接样本。fast 门槛 45→35、fast 屏 40→80、full 70→50；`vhh_plddt_min` 进 Defaults（options/CLI 可覆盖）；修 screen_vhh 写死 "pLDDT>70" 的误导日志 |
| G10 | VHH 对接提速 | 默认**刚性**对接：全长 VHH 是 ESMFold 单模型单构象，200+ 扭转只拖慢不增益；且该 vina 构建大配体单核（R10 证据）。发现坑：此构建**所有**配体（含小分子）都要 ROOT/ENDROOT 分子图，无图 ligand 报 "Unknown or inappropriate tag" → 刚性 = `TORSDOF 0`（图保留、扭转清零），`make_rigid_pdbqt()` + PDBQT 缓存带 flex/rigid 一致性检查。`vhh_dock_flex` 可开柔性。基准：`scripts/bench_vhh_dock.py`（rigid vs flex 同条件）。**G10-v2 CDR 片段对接**（基准的直接推论，见结果 3）：fast 默认 `vhh_dock_cdr_only`，按 pLDDT 低值连续区切 1-3 个片段分别对接，取最佳分 |
| R5 | 柔性工作流自动判据 | 前半：`analyze_completeness` 对无配体/无金属结构加 info 级 `apo_target` issue（提示后续 MD 判据）；后半：`interpret_stability` 加判据——apo + 平均 RMSF > 2.5 Å → 自动建议"短 MD 系综 + 柔性靶点工作流"（md_summary 把 is_ligand 注入 summary） |
| status | 命令增强 | 每阶段完成度+关键数字（对接数/命中/库名含回退标注/设计数/MD ns+final RMSD）、磁盘阶段 JSON 产物清单、最近 3 条工具失败（state.errors + transcript ok=false）、transcript 尾行 |
| pdbbind | 库修复 | 根因查实：pdbbind.org.cn 两个 tar.gz URL 已 404（下载改为注册制 download.php），138 字节文件就是 nginx 404 页。`_setup_pdbbind` 校验大小+tar 可读、丢坏件、软失败（resolve_library 回退 chembl35_small 兜底），不再打断整体 setup |

R1 真·结构域 RMSD（DSSP 域边界/NMDYN）留第 12 轮：需要域分割算法设计 +
独立 e2e 验证，单列一轮更稳。

### 验证
- [x] 快速套件全绿：182 通过（新增 35：stages 20 / vhh 10 / r5 1+1 / cli 3 / report 1）
- [x] 文档：README（模块 D 行/两层幂等/status/VHH 对接与 pLDDT 门槛/库回退）
- [x] G10 基准（scripts/bench_vhh_dock.py，vhh_30 773 原子 @1HVI，exh=1）
- [x] e2e：projects/r10_e2e 收尾验证（上轮遗留）——13:26 success，6 阶段
  全绿（track A 1 柔性对接 76 min → track B 2 设计 → MD 5ns×3 → 报告）。
  报告核对：柔性区 ✓ / 解读 ✓ / 回退库标注 "chembl35_small (fallback for
  dtp)" ✓ / 无 apo 误报 ✓
- [x] G8 幂等 e2e：同命令重跑 → 6 阶段全部 maybe_reuse 复用，~2 s 完成
  （原跑 5.3 h）

### 结果
1. **G8 阶段级幂等**：`agent/stages.py` 定义 6 阶段完成判据 + 缓存摘要；
   6 个整段工具加 `force`。单测覆盖：完成→复用（模块函数零调用）、
   force→真跑、未完成→真跑、scripted no-llm 全阶段跳过端到端。
2. **G9 VHH 命中面**：100 条已建模 VHH（真实库前 100）pLDDT 实测
   p50=31.5 / p75=33.2 / p90=34.5 / >45 仅 1 条 → fast 门槛 45→35
   （~6-10 条过阈 vs 旧 1 条）、fast n 40→80、full 70→50；
   `vhh_plddt_min`/`vhh_dock_flex` 进 Defaults + CLI。
3. **G10 刚性对接 + 基准（本轮最重要的测量）**：
   - 发现此 vina（AD4 血缘 f458505-mod）**所有** ligand 必须带 ROOT/ENDROOT
     分子图（小分子筛选一直是这么跑的），无图 ligand 报
     "Unknown or inappropriate tag" → 刚性 = 图保留 + `TORSDOF 0`。
   - 计时（exh=1，同网格；bench 同条件复跑）：200 原子片段 3.4 min；
     **200 原子 TORSDOF 0 vs TORSDOF 20 = 3.4 vs 3.6 min（无差异 → 该构建
     忽略/自推扭转数，成本由原子数主导，~O(n^1.9-2.2)）**；全长 773 原子
     **rigid=100.96 min / flex=109.73 min（1.10×）**（e2e 早期 flex 实测
     76 min 是轻负载下数据，量级一致）。
   - **同条件（同 seed/盒/核数）rigid vs flex 最优 pose 分数完全一致**：
     均 214373513.9（~2.1e8 kcal/mol，差 <0.001）→ 搜索收敛到同一个
     "撞墙" pose（773 原子塞进 25.84 Å 盒，绝大部分原子在盒外吃罚分）。
     修了 bench 脚本分数解析 bug（误取 VINA RESULT 行尾 rmsd_ub）。
   - 结论：刚性默认 = 语义正确 + pose 质量不劣 + 不更慢；**但绝对分
     目前由撞墙主导、意义有限** → 真正杠杆是缩小配体/盒子（CDR 片段
     对接：~100-200 原子 ≈ 3-5 min/候选，~15-20×）。
3b. **G10 v2 CDR 片段对接（本轮补做，基准的直接推论）**：
   - `vhh_cdr_fragments(pdb)`：按残基平均 pLDDT < 50 的连续区（=CDR/
     loop 候选）切片段，两侧各 pad 2 残基，≥4 残基，取最大 3 个；整结构
     都低（单 run 覆盖 >60%）则回退全长对接。
   - `dock_vhh_candidates(cdr_only=True)`：每个候选对接其片段，综合分 =
     最佳（最低）片段分，`fragment_scores` 保留各片段分。
   - fast 默认开（`vhh_dock_cdr_only`，Defaults+CLI 可覆盖），full 默认关。
   - **实测**（vhh_30 真实模型）：frag0=30 残基/257 原子对接 **6.25 min**
     （vs 全长 ~80-100 min，~13-16×）；frag1/2 ≈100 原子更快。
   - 修了 screen_vhh 一个潜伏 NameError（日志行 `plddtt_min`，无测试
     覆盖路径 → 新增 screen_vhh 假路径冒烟测试捕获）。
4. **R5 自动判据**：apo 结构在 analyze_completeness 出 info 级
   `apo_target` issue；interpret_stability 判据（apo + 平均 RMSF>2.5 Å
   → 短 MD 系综 + 柔性靶点工作流建议），md_summary 注入 is_ligand。
5. **status 增强**：阶段 ✓/✗+关键数字、产物清单、最近 3 条工具失败。
   实测 r10_e2e 输出正确（含 dtp 回退标注 + 早期 dtp 失败记录）。
6. **pdbbind**：138 字节文件 = nginx 404 页（两 URL 均 404，下载改注册制
   download.php）。setup 校验 + 软失败；坏件已删。
7. **r10_e2e 收尾验证**：13:26 success。MD 5ns×3（final RMSD 0.9 nm、
   自拟合 RMSD 各链 < 1.8 Å、二级结构 66%→69%、cluster 1 占 97.8-100%）。
   报告核对通过（柔性区/解读/回退库标注/无 apo 误报）。**G8 幂等 e2e**：
   同命令重跑 ~2 s 完成，6 阶段全复用。
8. **报告两处修复（收尾验证时发现）**：
   - 聚类段读错数据源（`m.replicas` 运行记录无 clusters 键 → 渲染
     `rep1: {}`）→ 改读 `summary.replicas`，百分比化展示（rep1: cluster 1:
     97.8%…）。回归测试 test_md_cluster_text_uses_summary_replicas。
   - 解读措辞矛盾（"平均 RMSF 偏高" 与 "各残基 RMSF 均低于阈值" 同屏）
     → 柔性区判据是"≥5 连续残基 ≥2×均值"，个别高值残基不构成"区"；
     文案改为 "无明显连续柔性区…个别高 RMSF 残基见上条"。

### 反思 / 下轮缺口（按优先级）
1. **G10 v3 片段对接打磨**：frag0（257 原子）仍偏大（~6 min），可加
   片段上限（如 ≤200 原子时截断到 CDR 核心区）或把盒子随片段自适应
   （当前用 pocket 盒，片段在盒外仍吃罚分，2.3e7 量级仍偏撞墙）。
   e2e 验证：r10_e2e 重跑时 vhh 阶段已完成（G8 会跳过）→ 需
   `run_vhh(force=true)` 或新项目验证 CDR 对接全链路。
2. **R1 真·结构域 RMSD**：DSSP 域边界/NMDYN 式动态域，`gmx_analyze
   kind="domain"`（分链+柔性区已覆盖大部分需求，收尾项）。
3. **G8 收尾验证 ✅（本轮完成）**：同命令重跑 ~2 s、6 阶段全复用。
   剩余：`drugagent rerun --stage X`（force 的 CLI 化）——易用性收尾。
4. DTP 下载重试（dtpbase.org 偶活）；pdbbind 注册后手动补库。
5. 易用性：vhh_screen 的 pLDDT 分布直方图进报告（G9 门槛选择可视化）；
   CDR 片段对接在报告里展示各片段分（fragment_scores 目前只在 state）。
6. **报告检查 ✅（本轮完成）**：柔性区/interpretation/回退库标注全部
   核对通过（并修了聚类数据源 + 措辞矛盾两处）。

## 第 10 轮

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

---

# 第 13 轮

## 计划
| # | 优先级 | 任务 | 依据（第 12 轮反思） |
|---|--------|------|----------------------|
| P1 | 高 | track B 设计验证缓存（scored.json，G8 下沉到设计粒度） | 反思 2：rerun 时 6 设计 × 6 min 无缓存重验证 |
| P2 | 高 | R1-v2 域-其余蛋白相对 RMSD（铰链/变构检测）+ 解读注 + 报告列 | 反思 3：自拟合 RMSD 度量域内部形变，真铰链需相对 RMSD |
| P3 | 中 | pLDDT 直方图进报告（G9 门槛 35/50 可视化） | 反思 4（原 P4 顺延） |

### 验证
- [x] 快速套件全绿：197（+5：scored.json 缓存、glob 副产物、hinge 签名、
      直方图 ×2）
- [x] r10_e2e md_rep1 真实轨迹 vs-rest（单副本）：5 域 vs-rest 末端
      7.5-14.7 Å，全部 > 自拟合（3.6-6.7 Å）→ 域间相对运动签名
- [x] `rerun --stage md` → state domain_rmsd_vs_rest + 解读注 9 触发
      （dom3：相对其余末端 14.1 Å、自身 5.9 Å → "构象重排为主"分支）
      + 报告域表 "末端相对其余 (Å)" 列 + 点线曲线 ✓
- [x] `rerun --stage vhh` → scored.json 缓存（2 条目）+ 缓存复用日志
      （"reusing cached validation"）+ plddt_all (n=80) + 报告 pLDDT
      直方图 + 门槛线 ✓；总时长 ~3 h → **7.5 min**
- [x] 文档 + git 提交

### 结果
1. **P1 track B 设计验证缓存**：`score_designs()` 抽出 + `scored.json`
   持久缓存（按设计 PDB mtime 匹配；只缓存成功、错误下轮重试；每个
   成功后增量写盘抗崩溃）。附带修掉一个 R12 遗留 bug：`vhh_design_*.pdb`
   glob 把打分副产物（`_binder.pdb`/`_complex.pdb` 等）也当成设计验证
   ——R12 rerun 的 "6 designs" 实为 2 设计 × 3 文件，且验证副产物又
   生成新副产物（12 个）。现只验证 `vhh_design_N.pdb` 顶层设计。
   **实测**：rerun vhh 从 ~3 h（R12）降到 **7.5 min**（track A 全缓存
   复用 + 1 设计缓存复用 + 1 个被 RF 覆盖的设计重验证 6 min）。
2. **P2 R1-v2 域-其余蛋白相对 RMSD**：`_kabsch_transform()` 抽出
   （行向量约定 mobile @ M + t，M = U D Vt；顺带修掉上一轮重构引入的
   转置 bug——准线性点云上 SVD 解的转置变体不是最优解，brute-force
   验证 tr(R*H)=34.7 < λ=42）；`domain_vs_rest_rmsd_series()`：每帧
   把其余蛋白（所有其他 CA）fit 到 frame 0，用同一变换测该域 →
   域刚体相对运动（自拟合序列按构造消掉的铰链信号）。summary 新增
   `domain_rmsd_vs_rest`（final/mean/series，多副本均值）；解读注 9
   两分支（内部稳定 <3 Å → "铰链/变构结构域特征"；否则 → "域运动以
   构象重排为主"）；报告域表新列 + 点线曲线。
   **实测（r10_e2e 3 副本）**：vs-rest 末端 8.4-14.1 Å（自拟合
   4.5-5.9 Å）——HIV-PR β-三明治 apo 5 ns 呈半刚性域重排，dom3
   （82-98）最大：vs-rest 14.1 Å + 自身 5.9 Å → 注 9 柔性分支触发。
3. **P3 pLDDT 直方图**：screen_vhh 返回 `plddt_all`（全建模库分布，
   n=80 浮点）；报告 VHH 节渲染直方图 + `vhh_plddt_min` 门槛虚线
   （fast 35 / full 50 随 options 自适应）。
4. **e2e 汇总（r10_e2e 最新状态）**：vhh_30 综合分 1.0（pLDDT 46.8 +
   frag1 -8.66 kcal/mol）居首；track B 2 设计（interface pLDDT 46.4，
   综合 0.31）；MD 域表 6 列 + 注 9；报告 pLDDT 直方图 n=80。

### 反思 / 下轮缺口
1. **binder 模块同款缓存缺失**：`binder.score_designs`（小分子/肽
   设计验证，tools_design.py）与 VHH track B 同构——无 per-design
   缓存、glob 可能同样吃到副产物。下轮把 `score_designs` 的缓存
   模式复用过去（或抽公共 helper）。
2. **vs-rest 阈值 4 Å 偏严/偏松未标定**：HIV-PR 5 ns 全部域超阈值
   （8-17 残基短域在 5 ns 尺度本就有大位移）。需要"域尺寸归一"
   （如 vs-rest RMSD / 域直径，或按域残基数分档）才能跨靶点可比；
   下轮在 2-3 个不同靶点（大/小、刚/柔）上标定阈值。
3. **scored.json 无失效判据（除 mtime）**：ESMFold 权重/参数变化后
   缓存应失效（加版本号或 ckpt hash 进缓存键）。
4. **直方图只有门槛线**：可加 p50/p90 标注 + 过门槛计数文字
   （"80 条中 6 条过 35"）——纯展示增强。
5. **track B 2 设计 interface pLDDT 相同（46.4）**：可能是 RF
   scaffold-guided 的噪声尺度问题或验证复用了同一复杂——下轮抽查
   两个设计的序列/构象是否实质不同（若 RF 每次生成近似设计，n_designs
   >1 的价值存疑）。

---

# 第 14 轮

## 计划
| # | 优先级 | 任务 | 依据（第 13 轮反思） |
|---|--------|------|----------------------|
| P1 | 高 | binder.score_designs 同款缓存（scored.json） | 反思 1：binder 设计验证同构无缓存 |
| P2 | 高 | track B 两设计 interface pLDDT 相同根因核查 + 几何接口度量 | 反思 5：46.4 完全相同，RF 设计是否实质不同 |
| P3 | 中 | scored.json 缓存键加 ESMFold 权重版本；直方图 p50/p90 + 过门槛计数 | 反思 3/4 |
| P4 | 中 | vs-rest 域尺寸归一（diameter，报告列） | 反思 2：4 Å 阈值跨域尺寸不可比 |

### 验证
- [x] 快速套件全绿：202（+5：binder 缓存、几何接口、版本失效、
      直径归一、报告统计）
- [x] P2 根因：两设计 PDB 逐原子对比（几何不同：目标相对帧位姿差
      2.2 nm；**均在接触**：min CA 距离 3.4/3.44 Å、<6 Å CA 对 30/31）
- [x] `rerun --stage md` → 域直径 11.4-24.4 nm + final_norm
      （0.041-0.109，dom1 相对摆动最大）进 state/report ✓
- [x] `rerun --stage vhh` → min_dist_a/contact/n_contacts 进
      track_b + ranked，报告 "目标距离 (Å)" 列（✓/⚠未接触 标记）✓；
      版本标签使旧缓存整体失效重验证（2×6 min，符合设计）
- [x] 文档 + git 提交

### 结果
1. **P2 根因（track B 两设计 interface pLDDT 完全相同 = 46.428…）**：
   逐原子对比确认**两个设计几何实质不同**（目标相对帧 VHH 位姿差
   2.2 nm，min CA 距离 3.4 Å、接触 CA 对 30/31——设计是好的、在
   结合面上）。根因在验证：RF 不写残基名（设计链全 GLY），
   `score_designs` 用 `_ca_sequence` 取序列后 ESMFold **只收序列**
   → 两设计输入序列完全相同（target + 200×GLY）→ pLDDT 位点盲。
   附带发现：`_complex.pdb` 落盘的是 ESMFold 输出（两文件 md5 相同
   因此），输入差异被覆盖。修复：新增 `design_interface_geom()`
   （设计链 vs 目标链 min CA 距离 / <6 Å 接触对 / contact<8 Å 标志）
   进每个设计条目 + ranked 候选 + 报告 "目标距离 (Å)" 列
   （⚠未接触 标记）。**插曲**：中途误判"VHH 漂离靶点 34 Å"——实为
   手测脚本把 Å 当 nm ×10（与 `design_interface_geom` 首版同型
   bug，已修+测试锁定）。
2. **P1 binder 缓存**：`binder.score_designs` 同款 `scored.json`
   缓存，签名 = [设计 mtime, 靶点 mtime, alt 序列, ESMFold 权重
   标签]；只缓存完整成功。与 vhh 版共享 `esmfold_version_tag()`。
3. **P3 缓存版本 + 直方图**：`esmfold_run.esmfold_version_tag()`
   （ckpt 名+大小）进 vhh/binder 缓存签名——权重更新自动整体失效
   （实测 R13 缓存无 `__version__` → R14 重验证 2 设计）。直方图
   标题加 p50/p90 + 过门槛条数（注意 plotly JSON 把 `/` 转义成
   \u002f，标题里避开斜杠）。
4. **P4 vs-rest 归一**：`domain_diameters()`（frame 0 CA 最大
   成对距离）进域条目 `diameter_nm`；summary vs-rest 加
   `final_norm`/`mean_norm`（÷直径）；报告域表 "归一化 (÷直径)"
   列。**实测 r10_e2e**：直径 11.4-24.4 nm，final_norm 0.041-0.109
   ——dom1（11.4 nm 小域）12.4 Å 位移 = 相对摆动最大（0.109），
   dom2（20.3 nm）8.4 Å 只有 0.041。绝对 4 Å 阈值与归一值并存，
   解读注暂用绝对值（归一阈值待标定，见反思 2）。
5. **e2e 汇总**：202/202 快速测试绿；r10_e2e 报告含目标距离列
   （两设计 3.4/3.4 Å ✓ 接触）、归一化列、直方图 p50/p90。

### 反思 / 下轮缺口
1. **R14-v2 设计验证升级**：当前 interface pLDDT 是序列级（全 GLY
   位点盲），几何接触是位姿级但只测"接没接触"。下一步：RF 输出
   里恢复真实序列（scaffold-guided 的 1EWN 序列是已知的——scaffold
   库 `data/tools/vhh_scaffolds/vhh1ewn.pdb`；把 scaffold 序列作为
   seqs 传入 ESMFold，pLDDT 才是真序列打分），或 ESMFold 带结构
   输入（pseudo-contacts）验证位姿。
2. **归一化阈值标定**：final_norm 0.04-0.11 在 5 ns apo 单靶点上
   没有"正常/异常"参照。下轮跑 1-2 个已知刚性的复合物靶点
   （如带配体的 HIV-PR、或 GFP 单体）建立基线分布，再定铰链判据
   （如 norm>0.2 且绝对>4 Å → 强铰链信号）。
3. **binder 缓存 e2e 未走真实路径**：单元测试覆盖，但 03_binder
   的 scored.json 要等下次 `rerun --stage binder` 才生成（binder
   阶段有 G8 幂等，正常流程不重跑设计验证——缓存主要惠及
   force-rerun 场景，符合预期）。
4. **RF 谨慎模式（cautious）**：`vhh_design_N.pdb` 已存在则跳过
   设计——`rerun --stage vhh` 永远复用首跑设计（n_designs 的
   多样性只在首跑生效）。若要每次 rerun 重新采样，需
   `--force-designs`（先挪走旧 PDB）；下轮可加 options.vhh_rf_cautious。

---

# 第 15 轮

## 计划
| # | 优先级 | 任务 | 依据（第 14 轮反思） |
|---|--------|------|----------------------|
| P1 | 高 | R14-v2：scaffold 真实序列进 track B ESMFold 验证（全 GLY 位点盲 → 真序列打分） | 反思 1 |
| P2 | 高 | vhh_rf_cautious 选项（cautious 默认导致 rerun 永远复用首跑设计） | 反思 4 |
| P3 | 中 | 铰链归一阈值标定：用带配体 HIV-PR（agent_smoke，2 ns×2）做刚性参照 | 反思 2 |
| P4 | 中 | binder 缓存 e2e（rerun --stage binder 走真实路径生成 scored.json） | 反思 3 |

### 验证
- [x] 快速套件全绿：205（+3：scaffold 序列、rf_cautious、轨迹单位回归）
- [x] P4：`rerun --stage binder` → 03_binder/scored.json 生成（2 设计，
      签名含 MPNN 真实序列 + 权重标签），新设计 iface 52.5/39.6（不同
      序列 → 不同分数，与全 GLY 时代对比鲜明）
- [x] **P0（本轮新增，最高优先）**：轨迹单位 Å→nm 修复（见下）
- [x] P3：带配体参照重算（见标定数据）
- [x] `rerun --stage vhh` → scaffold 序列生效 + 报告（进行中→已完成）
- [x] 文档 + git 提交

### 结果
1. **P0 轨迹单位 bug（本轮最大发现）**：`_ss_backbone_trajectory` 返回
   **MDAnalysis 原始单位 Å**，但 R12-R14 所有下游代码按 nm 处理
   （阈值 nm、报告 ×10 显示 Å）→ 所有域 RMSD 放大 10 倍，解读注 8/9
   的"构象漂移 5.9 Å"等全是假阳性（真值 0.59 Å）。修复：函数出口
   `/10.0`（unwrapping/make-whole 全程 Å 自洽，只在出口换算）+ 回归
   测试（真实轨迹：相邻 CA 间距 0.25-0.50 nm、蛋白尺度 < 8 nm）。
   **连带发现**：`HBOND_DIST = 4.5`（注释 Å、按 nm 坐标比较）→
   4.5 nm 匹配所有 O-N 对，ss_frac 0.99、"域"=整条链（1-98/100-197）。
   修为 0.45 nm 后 ss_frac 回到 0.75（与 DSSP ~73% 校准一致）、5 域
   恢复。合成 SS 测试夹具同步 ×0.1。
2. **P1 scaffold 序列（R14-v2）**：`design_vhh` 读 scaffold PDB（1EWN，
   200 残基 HLTRLGLEFF…）序列传入 `score_designs(seqs=…)`；复杂打分
   时若 PDB 序列全 GLY 且长度匹配则替换为 scaffold 序列（`seq_used =
   "scaffold"`）。缓存签名加 alt 序列。**语义**：pLDDT 现在回答
   "1EWN 序列作为此靶点 binder 的可折叠性"（序列级），设计间的
   判别仍靠位姿级几何（min_dist_a/接触对）——两者正交、互补。
3. **P2 vhh_rf_cautious**：选项（默认 True=旧行为）；False 时把顶层
   vhh_design_*.pdb 归档到 `archive_<ts>/` 再跑 RF → rerun 真正重采样。
4. **P3 标定（修复后真值）**：
   | 体系 | 时长 | 域 self-fit 末端 (Å) | 域 vs-rest 末端 (Å) | final_norm |
   |------|------|---------------------|--------------------|-----------|
   | apo HIV-PR dimer（r10_e2e） | 5 ns×3 | 0.5-0.6 | 0.8-1.4 | 0.041-0.109 |
   | 带配体 HIV-PR（smoke，参照） | 2 ns×2 | 0.5-1.7 | 1.7-2.4 | 0.068-0.098 |
   结论：同靶点两状态在 2-5 ns 内 final_norm 区间**重叠**（0.04-0.11），
   norm 不足以区分"配体稳定 vs apo 柔性"；绝对值上带配体参照 vs-rest
   (1.7-2.4 Å) 反而更高（2 ns 尚处松弛段）。**判定**：当前 note 8/9
   的 4 Å 绝对阈值在真值下两体系均不触发 → 无误报；norm 列保留作
   跨尺寸参照，铰链判据维持绝对阈值（跨靶点标定留给 R16+，需 ≥50
   残基差异的刚性参照系综）。解读注从"dom3 构象漂移 5.9 Å"（假阳性）
   变为"各结构域内部均稳定（最大 0.6 Å < 3 Å）— 高整体 RMSD 主要来自
   域间相对运动"（正确且更准）。
5. **e2e 汇总**：205/205 快测绿；r10_e2e md 重算（5 域真值 + 归一列）、
   vhh 重验证（scaffold 序列）、报告含目标距离列/归一化列/直方图统计。

### 反思 / 下轮缺口
1. **P0 的长尾**：`_ss_backbone_trajectory` 只被 analyze_ss 用（域分析
   链路），但 gmx 系（chain RMSD/RMSF/聚类）一直用 xvg（nm 正确）——
   本轮修复只影响域链路。R16 加一个"单位哨兵"测试：对同一轨迹
   比较 MDA self-fit 与 gmx rms 同帧值（±20% 容差），防止单位回归。
2. **scaffold 序列假设待证**：设计链 200 残基 = scaffold 200 残基，
   长度匹配支持"序列保留"，但未读 RF 输出 log 确认（RF 的
   scaffold-guided 不保证全序列保留——若 RF 重设计了 scaffold 区，
   pLDDT 是"1EWN 序列"而非"实际设计序列"的打分）。R16：解析 RF log
   或跑一次 ProteinMPNN 验证。
3. **norm 标定样本不足**：目前 1 靶点 2 状态。跨靶点标定需要一个
   已知铰链蛋白（如带门控回路的激酶）+ 一个刚性参照（如 GFP 单体
   5 ns）——需要新的 e2e 项目（~2×2 h MD）。
4. **binder 设计多样性**：`deterministic=False` 但每次 rerun 都重跑 RF
   （03_binder 无设计缓存——与 vhh 的 cautious 模式不对称）。R16 可给
   binder 加同款 cautious/归档选项统一行为。

# 第 16 轮

## 计划
| # | 优先级 | 任务 | 依据（第 15 轮反思） |
|---|--------|------|----------------------|
| P1 | 高 | scaffold 序列保留验证 + 保真度守卫（fidelity 字段进 state/报告） | 反思 2 |
| P2 | 高 | 单位哨兵测试：同一轨迹 MDA self-fit vs gmx rms（防 10× 单位回归） | 反思 1 |
| P3 | 中 | binder-cautious 对称（rerun --stage binder 不再无脑重跑 RF） | 反思 4 |
| P4 | 中 | 跨靶点铰链标定：刚性参照蛋白新 e2e（apo MD） | 反思 3 |
| — | — | 每轮互联网检索 + 文档 + git 提交 | 循环要求 |

### 验证
- [x] 快速套件全绿：212（+7：fidelity、sentinel、binder cautious、
      纯蛋白 target_prep、apo select/build）
- [x] P0a/P0b（本轮新增）：纯蛋白靶点 e2e 两处崩溃修复（见下）
- [x] P1：scaffold_rmsd_a 字段 + 真实设计 15.08 Å 发现 + RF 官方文档佐证
- [x] P2：单位哨兵测试（比预期复杂：三处 PBC/参考伪影先于单位问题）
- [x] P3：binder_rf_cautious（默认 False 保行为）
- [x] P4：ubq 刚性参照 e2e（2 ns×3，apo 通路）+ 标定表
- [x] 文档 + git 提交

### 结果
1. **P0a 纯蛋白 target_prep 崩溃（校准的直接拦路虎）**：靶点无配体时
   `lig_pdb` 未初始化 → `UnboundLocalError`（rigid_ubq/rigid_gfl 两次
   e2e 均死在此）。修复：`lig_pdb = None` 初始化 + 测试（纯蛋白 PDB
   走完全程，ligand_* 为 None，pocket 走表面回退）。
2. **P0b MD 无系统可选**：`select_system` 只认 ligand/screening/binder/
   vhh 四类复合物，纯蛋白靶点 → "no MD system candidates"。修复：无
   配体时追加 `apo` 选项（"apo 靶点，蛋白单独 MD — 柔性/铰链基线"），
   `build_complex_pdb` 加 apo 分支（仅靶点去水、保留原生链号）。两个
   e2e 测试。apo 通路即本轮 P4 的运行通道。
3. **P1 scaffold 保真度守卫**：`scaffold_fidelity()`（设计链 vs scaffold
   CA 的 Kabsch RMSD，Å，scale-invariant 不乘 10）+ `scaffold_rmsd_a`
   字段挂进每个设计。**真实 r10_e2e 发现**：设计链对 1EWN scaffold
   骨架漂移 **15.08 Å**（noise_scale_ca=0 并未刚性钉住 scaffold），但
   局部结构保留（40-mer 滑窗最小 5.9 Å、长度 200=200 精确匹配、SS
   70% 链状一致）。**身份发现**：`vhh1ewn.pdb` 实为 1EWN = 人 AAG
   DNA 修复糖基化酶核心（非 VHH，目录名误导）。**检索佐证**：RF
   scaffold-guided 官方文档明确其用"SS 模式 + 块邻接"粗粒度条件、
   "允许序列与精细结构变化"（sequence and fine structural variation）
   → 15 Å 漂移是**预期行为而非 bug**；scaffold 序列作为打分先验仍合理
   （近似），fidelity 字段现在把这个近似程度显式化。
4. **P2 单位哨兵（比"单位"复杂）**：对 r10_e2e rep1 轨迹逐层剥伪影——
   (a) `rmsd_r1.xvg` 的 gmx 参考是**过期 md.tpr**（与 xtc frame 0 相差
   ~21-23 Å，全 8847 原子）→ 核心 RMSD 2.083 与 per-chain 漂移系列
   系统性偏高（**R17 独立修复项**）；(b) 紧凑盒子（蛋白几乎充满盒子）
   下 wrapped 坐标在边界 flapping，gmx wrapped 视图 0→1.36 nm 而 MDA
   逐原子 unwrap 仅 0.13 nm；(c) `trjconv -ur compact` 自身成像选择
   会把 frame 0 平移整盒子向量（flat 2.54 nm 伪影）。哨兵最终形态：
   trjconv compact 展开 + **其自身 frame 0 作参考** + gmx rms vs MDA
   CA self-fit，三副本基线比 0.13-0.40（CA vs 全 backbone + 展开策略
   之差，非单位），断言带 [0.05, 1.0]——10× 单位错误（双向）都出带。
5. **P3 binder_rf_cautious**：`rfdesign(rf_cautious=True)` 在顶层
   design_*.pdb ≥ n_designs 时跳过 RF 复用现有设计（默认 False = 永远
   重采样，行为不变）；`run_design` 从 options 读 `binder_rf_cautious`。
   测试：cautious 不触发 RF、非 cautious 必触发、设计不足时仍触发。
6. **P4 跨靶点标定（刚性参照 = 泛素 1UBQ，apo，2 ns×3，~18 min）**：
   | 指标 | 刚性 ubq（2 ns×3） | 柔性 HIV-PR（5 ns×3, r10_e2e） |
   |------|-------------------|-------------------------------|
   | 整体 self-fit 末端 (nm) | 0.110-0.208（中位 0.154） | 0.289-2.083（中位 0.317，2.083 为过期 TPR 副本） |
   | 域 vs-rest 末端 (nm) | 0.076-0.183（中位 0.107） | 0.075-0.160（中位 0.115） |
   | final_norm（每域聚合） | 0.040-0.080（中位 0.056） | 0.041-0.109（中位 0.056） |
   **结论**：① 域级 vs-rest/final_norm 在刚性与柔性靶点间**几乎完全
   重叠** → 域级 norm 不足以单独分类"铰链 vs 刚性"，刚性蛋白基线
   ~0.04-0.08 必须先扣；② **整体 self-fit 末端能干净分开**（ubq 上限
   0.21 nm < HIV-PR 下限 0.29 nm，无重叠）→ 建议把它纳入铰链/柔性
   判据（当前 interpret 仅用域级）；③ 绝对 4 Å 铰链阈值远高于刚性
   基线（vs-rest 上限 1.83 Å）→ 对刚性蛋白无误报，阈值保守合理。
   GFP（1GFL）尝试作第二刚性参照失败：色原体 CRO 残基（SER 65 缺 O）
   使 pdb2gmx 报错 → 修饰残基处理是 R17 缺口。
7. **e2e 汇总**：212/212 快测绿；新项目 rigid_ubq（apo 通路首跑成功）、
   rigid_gfl（ pdb2gmx 修饰残基失败，留档）。r10_e2e 状态未动（R15 值）。

### 反思 / 下轮缺口（R17）
1. **过期 TPR 参考（本轮最大未修 bug）**：gmx 链级 RMSD 系列
   （`rmsd_rN.xvg` + per-chain drift）以 md.tpr 为参考，而 tpr 是
   2 ns 段产物、xtc 是延长段轨迹（起点已漂移 ~21 Å）→ rep1 "final
   2.083 nm" 应读作 ~0.3。修复候选：分析前 `trjconv -b 0 -dt 大值`
   抽 xtc frame 0 作参考，或延长前重新 grompp。影响：所有项目的
   链级漂移系列与"最终 RMSD 0.896"均值。
2. **铰链判据升级**：把整体 self-fit 末端（本轮证明可分刚性/柔性）
   并入 interpret 注 8/9，形成"整体 + 域级"两级判据；刚性基线
   （ubq 0.04-0.08 norm / ≤0.21 nm self-fit）写进注释。
3. **第二刚性参照 + 修饰残基**：修 pdb2gmx 对 CRO 类修饰残基的
   处理（或选无修饰残基的刚性蛋白如 crambin 1CPS 去交叉残基），
   补第二刚性点 + 一个已知铰链体系（无辅因子的腺苷酸激酶），
   把标定从"1 刚性 + 1 柔性"扩成可发表级小矩阵。
4. **scaffold 内容**：1EWN 是 AAG 酶核心非 VHH——要么换真 VHH
   scaffold（如 1HYW 类），要么把目录/注释改名（vhh_scaffolds/
   aag_core_1ewn）并文档化"fold 条件而非序列条件"。
5. **展开策略一致性**：gmx rep2 展开视图 0.984 nm vs MDA 0.130 nm
   差异提示 compact-com 与逐原子 unwrap 在紧凑盒子下不等价；域分析
   链路（MDA 逐原子）正确，但 gmx 链级 xvg 仍用 wrapped——过期 TPR
   修复时顺带评估是否统一 unwrap。

# 第 17 轮

## 计划
| # | 优先级 | 任务 | 依据（第 16 轮反思） |
|---|--------|------|----------------------|
| P1 | 高 | 修 gmx 链级 RMSD "过期 TPR 参考"（rmsd_rN.xvg + per-chain drift 系统性偏高） | 反思 1 |
| P2 | 高 | 铰链判据升级：整体 self-fit 并入 interpret + 刚性基线注释 | 反思 2 |
| P3 | 中 | 修饰残基处理（第二刚性参照 crambin + 柔性参照 CDK2 启动） | 反思 3 |
| P4 | 低 | scaffold 内容文档化（1EWN=AAG 非 VHH） | 反思 4 |
| — | — | 每轮互联网检索 + 文档 + git 提交 | 循环要求 |

### 验证
- [x] 快速套件全绿：218（+6：unwrap 幂等/接线、修饰残基过滤×2、
      修饰残基偏好、刚性基线注释）
- [x] P1：r10_e2e md 重算，rep1 2.083 → 0.328（触发 R6 自动延长，见下）
- [x] P2：interpret 注 1b（刚性基线 √t 折算 + 域级 norm 基线 0.08）
- [x] P3：crambin（CPM 交叉链接）走 apo 通路启动；CDK2(1M17+AQ4)
      配体通路启动（两作业后台运行中）
- [x] P4：scaffold 目录 NOTE.md
- [x] 文档 + git 提交

### 结果
1. **P1 "过期 TPR"真相 + 修复（本轮核心）**：重新测量发现 R16 的
   "tpr 与 xtc frame 0 差 21-23 Å"是 **MDA 读 xtc 给 Å、读 tpr 给 nm 的
   单位假象**——单位一致后 tpr==xtc f0（Kabsch 0.05 Å，三副本全同）。
   真正病因：**紧凑盒子下 wrapped 坐标在盒边界 flapping**（蛋白几乎充满
   盒子；整体平动 > box 后单帧原子跨盒断裂，单帧 rot+trans fit +
   min-image 距离被破坏）。GROMACS 论坛同款案例（"reference broken
   across PBCs → 突刺"）。修复：`_ensure_unwrapped()`（trjconv
   `-ur compact -pbc mol`，G8 幂等：源 xtc 更新才重展开）→ analyze 全部
   per-frame-fit 调用（rms/rmsf/gyrate/cluster/链级 rms）改用展开轨迹；
   tpr 参考保留（每帧独立 fit 吸收恒定盒向量偏移，实测 tpr 参考与展开
   f0 参考结果恒等 0.328）。**重算结果**：rep1 2.083→**0.328**、rep2
   0.289→**0.984**、rep3 0.317→**0.317**——新画像：链自拟合全副本
   0.05-0.14 nm（内部刚性），链相对 1.15-3.88 nm（**二聚体相对运动**
   主导，rep2 最显著）。另触发 R6 自动延长：修正后序列尾-头漂移 0.7 Å
   > 0.5 Å 判"未收敛"→ 3 副本各延长 5 ns（后台运行中）。
2. **P2 刚性基线判据（interpret 注 1b）**：泛素（76 残基单体，2 ns×3）
   整体 RMSD 末端 ≤0.21 nm、域级 final_norm ≤0.08 作为刚性基线；
   RMSD √t 标度 → 基线上限 = 0.21×√(ns/2) nm。summary 新增
   `ns_run`（时间轴）。实测两分支：rigid-like 合成例 → "整体处于刚性
   蛋白水平"；HIV-PR 修正值（均值 0.543 > 0.332 上限）→ "超出刚性
   蛋白水平"+ 域级 norm 0.109 > 0.08 → "真实域运动/铰链信号"。
3. **P3 修饰残基**：(a) `_keep_standard_res()`（20 AA + 常见变体
   MSE/HID/ACE/NME + 核酸，HETATM→ATOM 规范化）进 apo 构建分支；
   (b) `select_system` 新增偏好：靶点"配体"全部属修饰残基集合
   （CPM/CRO/金属/辅因子，64 名）时 apo 优先（1CPS 的 CPM+ZN 不再走
   acpype 配体通路；A77 等真配体不受影响；`--md-system` 仍可强制）。
   crambin（48 残基）apo 2 ns×3 与 CDK2 1M17（290 残基+AQ4 配体通路）
   2 ns×2 已后台启动——第二刚性点 + 已知柔性参照（gate loop）。
   GFP 1GFL 的 CRO 色原体（SER 65 缺 O）仍需残基级修复（R18）。
4. **P4 scaffold NOTE**：`data/tools/vhh_scaffolds/NOTE.md` 文档化
   1EWN=人 AAG 糖基化酶核心（非 VHH）、scaffold-guided 语义（SS+邻接
   条件、序列/精细结构可变）、替换真 VHH scaffold 的方法。
5. **检索佐证**：GROMACS 论坛 "PBC effect on RMSD calculation" 与本轮
   诊断完全一致（wrapped 参考跨盒断裂 → RMSD 突刺；合并/延长轨迹用
   错 tpr 参考是经典坑）；推荐做法 -center -pbc mol -ur compact +
   正确 frame-0 参考，与本实现一致。
6. **e2e 汇总**：218/218 快测绿。r10_e2e md 重算中（3 副本延长 5 ns，
   完成后 state/报告自动更新）。

### 反思 / 下轮缺口（R18）
1. **收口后台作业**：crambin/CDK2 完成后提取 final_norm/self-fit 进
   标定表（目标：刚性 ≥2 点 + 柔性 ≥2 点）；r10_e2e 延长完成后确认
   收敛判据与报告。
2. **收敛判据与 P1 的耦合**：R6 自动延长现在基于修正后序列（更真实），
   二聚体相对运动 5 ns 未平台 → 可能需 10-15 ns 才"收敛"；考虑收敛
   判据区分"内部柔性漂移"与"寡聚体相对运动"（用链自拟合而非整体
   RMSD 判漂移）。
3. **GFP/CRO 类修饰残基**：残基级修复（补缺的骨架 O 或删受影响残基）
   才能作第二刚性参照候选；或改用无修饰残基的刚性蛋白。
4. **unwrap 策略一致性（R16 遗留）**：compact-com（gmx 链）与逐原子
   （MDA 域）在紧凑盒子下不等价（rep2 0.984 vs 0.130）；域分析链路
   正确，但两类指标并列展示时需在报告注明展开策略差异。
