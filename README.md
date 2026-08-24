# DrugAgent — 药物发现一体化 Agent

**流水线是工具箱，LLM 是主程序。** DrugAgent 2.0 用一个 ReAct 主循环
（LLM + 原生 function calling）驱动整个药物发现流程：LLM 自己规划、调用
~41 个细粒度工具、自己写 MDP/选力场、失败时自己读日志排错重试；
4 个固定里程碑检查点 + LLM 随时发起的动态确认，`--auto` 时自动通过。

架构细节见 [DESIGN.md](DESIGN.md)。

## 模块（工具箱覆盖的功能）

| 模块 | 内容 | 关键工具 |
|---|---|---|
| A 靶点准备 | PDB 文件 / PDB ID / FASTA / 裸序列输入；完整性分析 + 结构坑预检（多 MODEL/altloc/缺失残基/金属/核酸/无序末端）；agent 判定；自动+手动修复 PDB；清洗；口袋检测；PDBQT 转换 | RCSB, obabel, ESMFold（仅序列输入） |
| B 小分子筛选 | 大库（DTP/ChEMBL/PDBBind 或自定义 SDF）→ 标准化 → 理化/ML 预过滤 → Vina 并行对接 → GNINA 复打分 → agent 定命中标准（参考配体重对接做阳性对照）；**柔性靶点工作流（R2/R5：MD 构象系综 → 多构象选择 + 侧链 `--flex`，consensus 平均分）** | RDKit, Vina, GNINA |
| C binder 设计 | RFdiffusion 从头设计 + ProteinMPNN 序列 + ESMFold 单体/复合物打分（界面 pLDDT） | RFdiffusion, ProteinMPNN, ESMFold |
| D 纳米抗体 | 轨道A：VHH 库（data/libraries/vhh_library.fasta 或合成）→ ESMFold 建模 → pLDDT 过滤（R11/G9：fast 35 / full 50，可 `vhh_plddt_min` 覆盖）→ **刚性对接**（R11/G10；fast 默认 **CDR 片段对接**（R11/G10-v2：pLDDT 低值区切 1-3 片段，~15× 提速，`vhh_dock_cdr_only` 可关；`vhh_dock_flex` 可开柔性））→ 并行筛选；轨道B：RFdiffusion scaffold-guided（1EWN）从头设计；综合评分 + agent 选择 | ESMFold, RFdiffusion, Vina |
| E MD 模拟 | 体系选择 → pdb2gmx + ACPYPE (GAFF2) → 溶胀/加离子/EM → **NVT→NPT 平衡段（R5，带位置约束；R8 改 C-rescale）+ 生产 MD 烧入段剔除** → N ns × R 副本（自平衡态分叉）→ **收敛判定 + 自动延长（R6：RMSD 平台 + 主导簇，未收敛自动续跑合并）** → RMSD/RMSF/Rg/聚类 + **柔性诊断（分链/自拟合 RMSD + DSSP 式二级结构 + 柔性区定位（连续高 RMSF 段→残基区间）+ 规则化柔性解读）** + **金属离子协调（R3）** + **核酸链原生参数化（R4）** + **辅因子/血红素自动参数化（R8：ACPYPE + 嵌入金属独立离子）** | GROMACS 2023.1, ACPYPE, MDAnalysis |

## 快速开始

```bash
cd /home/data/lrs/drug/drugagent

# 1) 一键环境配置（幂等，可重复执行）
env/bin/python -m drugagent.setup

# 2) 快速端到端验证（1HVI 靶点，全模块，小规模，无人值守）
env/bin/python -m drugagent.cli run --target 1HVI --modules all --fast --auto

# 3) 生产运行（检查点交互：批准/修改/中止）
env/bin/python -m drugagent.cli run --target 1HVI --modules screen,binder,vhh,md
env/bin/python -m drugagent.cli resume --project projects/<项目>

# 4) 状态 / 报告
env/bin/python -m drugagent.cli status
env/bin/python -m drugagent.cli report --project projects/<项目>
```

常用 `run` 参数：`--library dtp|chembl35|pdbbind|<SDF路径>`，`--n-jobs 32`，
`--md-ns 100 --md-reps 3`，`--max-steps 300`（agent 步数预算），
`--no-llm`（确定性脚本模式），`--llm-base-url/--llm-model`（覆盖 LLM）。
MD 细调：`--md-salt`（离子浓度 M）、`--md-divalent MG --md-divalent-m 0.01`
（二价抗衡离子）、`--md-extend-ns`（自动延长步长 ns）、
`--md-max-extensions`（自动延长轮数）、`--md-burn-in-ps`（烧入段剔除）。

排错提示：外部命令（gmx/vina/...）失败时异常信息自带日志尾部（最近 40 行），
直接看报错即可，不用手动翻日志文件。

## Agent 架构（2.0）

- **主循环** `drugagent/agent/loop.py`：LLM（OpenAI 兼容端点）每轮返回
  tool_calls → 执行 → 结果回灌 → 循环；终止于 `finish` 工具 / 预算耗尽 / 人工介入。
- **工具** `drugagent/agent/tools_*.py`：41 个工具 = 元工具（文件/shell/决策/确认）
  + 五个阶段的细粒度工具 + 五个 `run_*` 整段确定性兜底工具（1.0 流水线整体保留）。
- **参数主权**：力场由 `gmx_env` 列出、agent 选择；MDP 由 agent 自己
  `write_file`（模板只是起点）；盒子/盐浓度/barostat/dt 均可改；
  每个关键决定经 `record_decision` 留痕（`decisions.json`，报告中展示）。
- **排错**：工具失败 → agent 读日志（`read_file`）→ 诊断 → `edit_file` 打补丁 →
  重试；同一问题修 3 次仍失败则 `ask_human`。
- **检查点**：固定 4 个（target / screening / design / md）+ 动态 `ask_human`；
  `--auto` 自动通过并记录。
- **状态**：项目目录是唯一事实源。`state.json` 保存阶段索引；
  `agent/transcript.jsonl` 记录完整对话+工具调用（resume = 重放 transcript 继续）；
  两层幂等：工具产物级（已完成的 mdrun/对接/PDBQT 自动复用）+
  **阶段级（R11/G8：`run_*` 工具在 state.json 阶段段完整时整体跳过，
  `force=true` 强制重跑；重启崩溃的 e2e 不再整段重跑已完成的 35 分钟筛选）**。

## LLM

默认 `http://127.0.0.1:18080/v1` / `qwen3.8-27b-uncensored`（本地 llama.cpp，
需支持 function calling）。可用环境变量覆盖：
`DRUGAGENT_LLM_BASE_URL` / `DRUGAGENT_LLM_MODEL` / `DRUGAGENT_LLM_API_KEY`。

## 目录结构

```
/home/data/lrs/drug/drugagent/        # 项目根（持久盘）
├── env/                              # python 环境 (py3.12, conda)
├── data/
│   ├── libraries/                    # 小分子库 SDF / VHH 库 fasta
│   ├── weights/                      # ESMFold/ESM2 权重、torch/hf 缓存
│   └── tools/
│       ├── gromacs/                  # GROMACS 2023.1 (自编译, amber99sb-ildn)
│       ├── vina/                     # autodock-vina
│       ├── gnina/                    # GNINA 1.3.1 ELF (CPU)
│       ├── RFdiffusion/              # RFdiffusion (hydra) + models/ 权重
│       ├── vhh_scaffolds/            # 1EWN VHH 支架 + secstruc 邻接
│       └── 3Dmol/                    # 3Dmol-min.js
├── drugagent/                        # 主包
│   ├── agent/                        # 2.0 核心: loop/提示词/40 个工具
│   ├── vendor/openfold/              # aqlaboratory/openfold + CPU 补丁
│   ├── modules/                      # A-E 五个模块（工具后端）
│   ├── graph.py                      # 1.0 LangGraph 状态机（兜底工具复用）
│   ├── cli.py                        # typer CLI (run/resume/status/report/setup)
│   └── report/                       # 交互式 HTML + PDF (WeasyPrint)
├── DESIGN.md                         # 2.0 架构设计
├── projects/                         # 每次运行一个目录 (01_target…05_md, reports/, agent/)
├── tests/                            # pytest 套件 (176 快测 + 15 slow e2e)
└── logs/                             # 构建/运行日志
```

## 运行状态与断点

- `status`（R11 增强）：每个阶段的完成状态 + 关键数字（对接数/命中数/
  设计数/MD ns 与 final RMSD/库名含回退标注）、磁盘上的阶段 JSON 产物、
  最近 3 条工具失败（state.errors + transcript 里 ok=false 的调用）、
  transcript 最后一条；`decisions.json` 记录每个判断与依据。
- `resume --project DIR`：重放 `agent/transcript.jsonl`，从断点继续
  （产物幂等 + 阶段级复用，R11/G8）。
- 报告：`projects/<项目>/reports/report.html`（Plotly 交互图 + 3Dmol 结构）
  与 `report.pdf`。

## 测试

```bash
env/bin/python -m pytest tests/ -m "not slow" -q     # 快速单测 (176 用例)
env/bin/python -m pytest tests/ -q                   # 含 slow（需 vina/GROMACS/RF 权重）
env/bin/python -m pytest tests/test_mdsim.py -m slow -q \    --basetemp=$PWD/data/fixtures/_ptmp             # MD e2e 建议本地盘 basetemp
```

> slow 构建 e2e 在本机 /tmp（tmpfs）上偶发失败：有后台清理进程会在
> gromacs 运行中删掉 pytest 的 tmp 目录。代码本身稳定，遇到偶发失败用
> 上面 `--basetemp` 指向本地盘重跑即可。

## 1HVI 全模块端到端（`--fast --auto`）

```bash
env/bin/python -m drugagent.cli run --target 1HVI --modules all --fast --auto
```

在 64 核 CPU 上约 3–5 h 完成全模块：靶点（1HVI 二聚体，99 残基/单体）→
小分子筛选（标准化 + Vina 对接 + GNINA 复打分 + agent 定命中）→
binder 设计（RFdiffusion 从头 + ESMFold 复合物打分）→ VHH 双轨 →
体系选择（agent）→ GROMACS EM + 5 ns×3 副本 MD →
RMSD/RMSF/Rg 分析 + HTML/PDF 报告。

过程中修复的关键问题（已并入代码）：

- **ESMFold 权重加载**：发布 checkpoint 只含 3923 个 trunk/adapter key，
  ESM2 由 `esm.pretrained` 单独加载；旧逻辑从随机初始化的模型 state_dict
  出发只回填 4 个 key，导致 ~99.9% 参数随机（pLDDT 恒为 ~50、CA 塌缩）。
  现改为以 checkpoint 权重为基底、仅重映射 4 个 IPA point 平铺 key。
- **RFdiffusion Fold 模板特征**：`InpaintSeq_Fold` checkpoint 需要
  d_t1d=28 / d_t2d=47；`SelfConditioning` 现在对无支架的 inpaint 补
  sec-struc (MASK) 与 block-adjacency (MASK) 特征。
- **GROMACS 副本启动**：steepest-descent EM 不总写 `em.cpt`；副本改从
  `em.gro` 起始并 `gen_vel=yes`；Berendsen 需显式 `compressibility`；
  grompp 加 `-maxwarn 1`（Berendsen 的“非严格系综”提示）。

## 已知注意事项

- **辅因子/血红素自动参数化（R8）**：`build_system` 构建前
  `_find_cofactors` 检出**被蛋白链包围的非蛋白残基**（重原子 ≥3、不在
  PROT/DNA/金属/水离子白名单、且该链含 ≥1 蛋白残基——纯配体单链不算），
  例如 HEM 类血红素、NAD/FAD（单残基小分子仍按配体走）。命中后：
  (1) `_reassign_cofactor_chains` 把辅因子残基改写到**新的自由链号**
  （pdb2gmx 对"同一链内蛋白+其他残基混编"报 "do not have a consistent
  type"，拆开即解）；(2) 蛋白侧 pdb2gmx、辅因子侧走 ACPYPE
  （`_run_acpype`：obabel PDB→mol2 后**必须统一 mol2 残基字段**——obabel
  把 `resname+resseq`（如 `HEM9`）写进原子残基列，antechamber 会检出
  多个残基 `{'HEM','MOL'}` 报错；`-b NAME` 原样保留为 itp 分子名）；
  (3) **嵌入金属**（`_embedded_metals`：辅因子残基内部、原子名属于
  METAL_PROPS 的原子，如 HEM 中央 FE）先**从 ACPYPE 输入剥离**——本环境
  无 sqm，antechamber 只能 Gasteiger 荷电，而 **Gasteiger 没有 Fe 参数**
  （"No Gasteiger parameter for atom (ID: 20, Name: FE)"），且孤立的中心
  金属会触发 antechamber "Atoms TOO scattered"（>3 Å 无邻居）；
  剥离后以**独立单原子离子**重新并入体系：`_append_embedded_ions` 写
  `embion_{el}.itp`（分子类型 + `[ atoms ]`），FF 未显式定义该元素时
  另写 `embion_{el}_types.itp`（纯 `[ atomtypes ]`）；
  (4) **GROMACS 2023 的 itp 指令顺序坑**：一旦解析到任何 `[ moleculetype ]`，
  后面再出现 `[ atomtypes ]` 就报 "Invalid order for directive atomtypes"
  ——所以 types include 插在第一个 `#include` 之后、分子 include 插在
  最后一个 `#include` 之后；(5) combined.gro 的离子行必须**插到盒向量行
  之前**（GRO 盒行必须是最后一行；直接 append 会顶掉盒行，而"只改
  最后一行"的修法则会覆盖刚写入的离子行——两个相邻的坑）。
  e2e（合成 HEM 复合物：8 ALA + 20 原子卟啉环 + 中央 FE）：检出 →
  拆链 → mol2 统一 → ACPYPE/Gasteiger → FE 独立离子 → EM 通过，
  全程无手工干预。已知局限：HEM 用 GAFF2/Gasteiger 荷电（无 sqm/AM1-BCC，
  血红素铁络合物电荷近似）；嵌入金属独立离子**无配位约束**（HEM 的 FE
  靠环平面约束间接定位，长 MD 可能漂出平面）；单原子孤立金属残基走
  `_build_ligand_system` 配体路径时会被丢弃（该路径对 METAL_RES 行既不
  归蛋白也不归配体——纯蛋白+金属体系不受影响，走 protein 路径）。
- **簇均值代表结构（R8）**：R7 的簇代表结构原来取簇**首帧**（`trjconv
  -dump`）——首帧只是簇边界的一个任意快照。R8 改为 `_cluster_mean_pdb`：
  MDAnalysis 读 TPR（拓扑）+ xtc，对簇内**全部帧**逐原子累加取均值，
  写出代表 PDB（MA 缺失回退首帧 dump）。本环境实测的三个 MA 2.10 +
  GROMACS 2023.1 xtc 坑：(a) MDA 把该 xtc 读成 **10× 尺度**（坐标与盒子
  一致地 ×10，Å 当 nm），而 PDB writer 原样写内部值当 Å——往返自洽，
  只要比较也留在 MDA 内部尺度即可（测试即按此对齐）；(b) **`AtomGroup
  .positions` 属性返回副本**（xtc 宇宙上连整个 universe 都是副本），
  项赋值 `positions[idx] = ...` 静默丢失——必须用 setter
  `u.atoms.positions = new_array`；(c) MA 2.10 没有
  `mda.lib.pdbwriter`，写 PDB 用 `AtomGroup.write('out.pdb')`。
  e2e（1BNA 0.1 ns 生产 + gromos 聚类）：代表 PDB 与 MDA 直算的簇均值
  逐原子偏差 < 0.5 Å。**逐原子 PBC 修正（R8 调试收尾）**：xtc 存的是
  未包裹的连续坐标，若某原子净扩散超过半个盒子，朴素平均会把均值抹到
  两个盒子图像之间——因此累加前对**簇首帧**做逐原子最小镜像
  （`d = p - ref; d -= box*round(d/box)`，等价 `gmx rms` 的内部做法；
  曾尝试的 `trjconv -pbc` 预处理步在 2023.1 与 `-fit` 互斥、输出组还会
  改变原子数，已整体移除，MDA 一步完成）。1HVI 单体 Cα 跨度 4.09 nm
  > 盒子 y 半边（3.68 nm）→ 基于 COM 的整体包裹有歧义，逐原子最小
  镜像才是正确的平均策略（只要单原子净位移 < 半盒）。
- **PBC 帧对齐（R8）**：MD 代表结构在**模拟盒子帧**、晶体 PDB 在
  **晶体学帧**，两者相差**逐原子整数盒子向量**（pdb2gmx 逐原子包裹
  进盒子，且 1HVI 的 CRYST1 晶胞 ≠ 生产盒子——NPT 平衡段把盒子从
  5.95×3.68×4.00 nm 膨胀到 6.21×3.84×4.18 nm）。全局 Kabsch 把这些
  整数向量平均成一个"糊掉"的变换（1HVI 二聚体上残差 190 Å，而真实
  构象差只有 ~3 Å）。修法：
  (1) `_mi_kabsch`——镜像算术全部放在**模型帧**（盒子晶格轴对齐）：
  迭代"逐原子选最近图像 → 刚体重拟合"，仅在接受分数改进时前进
  （盒子边界附近的原子会翻转图像选择，朴素循环会振荡），并多起点
  （整体裸拟合 + 首/中 10 残基核心，核心共享单一盒子向量故无歧义）；
  (2) `align_pdb_to_reference` 的**Cα 伙伴锚定**：单一刚体变换原则上
  消不掉逐原子盒子向量（1HVI 的 Cα 向量分布 4 种取值，写出的文件
  天然碎成多个盒子副本），改为每个原子写在"参考伙伴位置 + 模型帧
  局部偏移（残基内 < 3 Å，无盒子歧义）"，Cα 偏移带 PBC 修正并承载
  构象变化；无伙伴原子回退刚体变换。
  实测（1HVI rep vs 晶体）：190 Å → 对齐后 Cα 平均 3.7–4.0 Å /
  最大 5.5–6.2 Å；池化对接 MD 代表从 -4.5/-5.3（坏对齐）→
  重对接中（见 e2e）。合成 TDD：`test_mi_kabsch`（位置相干包裹 +
  核心单一向量 + 跨边界云）+ `test_cluster_representatives_smoke`
  改为按轨迹实际簇占比自适应。
- **生产段段错误与 barostat（R8）**：1HVI 全体系（8879 原子）生产 MD
  **两个副本都死在 step 10000 附近（20 ps）**——mdrun SIGSEGV，崩溃前
  能量表无 NaN，但连续出现 Berendsen PLOESS "pressure scaling more than
  1%"（mu 0.91–1.10）；同一 MDP 的 NPT 平衡段（带 posre）却完整跑完
  20 ps。判定：无约束体系弛豫到 ~20 ps 时 Berendsen 压控失稳，本 GROMACS
  2023.1 构建下直接段错误。R8 把平衡/生产 MDP 的压控统一改为
  **C-rescale（tau_p 2 ps）**（v-rescale 温控保持不变）。NVT 平衡段幂等
  补齐（`eq_nvt.gro` + log "Finished mdrun" 存在即跳过，NPT 从
  `eq_nvt.cpt` 续跑）——v5 重建时 NVT/NPT 幂等跳过直接进生产，省 ~30 分钟。
  另修**平衡段 nsteps 单位 bug**：模板曾按 `ps*1000/0.002` 计算
  （1 ps 当成 1000 ps 的步数），fast 的 10 ps NVT 实际跑了 10 ns
  （v2 重建 1h47m 的元凶）；正确为 `ps/0.002`（1 ps = 500 步），
  生产模板 `ns*1000/0.002` 本来就对。
  1HVI 全模块 smoke（2 ns × 2 副本，C-rescale）实测：最终 RMSD
  2.27 Å（修复前 5 ns×3 的旧轨迹 7.23 Å——旧流程带 10 ns 的错误
  NVT 加热 + Berendsen）、Rg 2.28 nm（旧 2.96 nm）、单一主导簇
  （99.8%）、二级结构 64.6%→64.9% 稳定、配体 RMSD 末段 2.1 Å
  （配体全程未脱靶）；R6 收敛判定 0 轮延长即通过（RMSD 平台
  0.2 Å + 主导簇 100%）。
- **对接侧防御（R8）**：`dock_one` 对刚体受体 PDBQT 做**图关键字剥离**
  ——若文件含 `ROOT/ENDROOT/BRANCH/ENDBRANCH/TORSDOF` 行，写
  `<stem>_nograph.pdbqt` 副本再对接（本 vina 构建对**刚体**受体的图
  关键字报 "Unknown or inappropriate tag found in rigid receptor"；
  柔性受体/柔性配体则必须保留图）。触发场景：`to_pdbqt(flex=True)`
  的 sanitizer 会给受体重新加回图关键字——调用方显式传 `flex=False`
  即可避免，`dock_one` 的剥离是最后一道保险。
- **磁盘**：环境+权重+库约 40 GB；MD 轨迹按 ns 增长（5 ns×3 副本约 2 GB）。
- **binder 序列来源**：RF 从头设计的 PDB 残基名为 GLY（不写序列）。
  现内置 ProteinMPNN（`data/tools/RFdiffusion/mpnn/`，vanilla v_48_010）：
  对设计骨架直接生成真实序列；RF inpaint 的多链输出会自动取全 GLY 的设计链。
  MPNN 缺失/失败时回退启发式序列（此时打分序列可能偏低，报告里会标注）。
  ESMFold 打分时用真实序列而非 PDB 里的 GLY 残基名。
- **ESMFold 打分口径**：binder/VHH 的 pLDDT 为 ESMFold 从序列重新折叠的
  置信度，不等于复现 RF 设计构象的置信度；de novo 设计 pLDDT 偏低属正常。
- **MD 参数**：agent 模式下 MDP 可由 agent 改写（模板 R8 起为 C-rescale/SPCE
  —— 1HVI 实测 Berendsen 在生产段 ~20 ps 双双 SIGSEGV，见下）；
  fast 模式短轨迹属稳定性验证，非生产级自由能/收敛结论。
- **MD 单位**：GROMACS 的 rms/gyrate 输出为 nm（rmsf 为 Å，cluster 截断 1.5 为 nm）；
  报告展示时 RMSD 已换算为 Å（×10），Rg 保持 nm。
- **VHH 对接**：完整 VHH 是"巨型配体"（~770 原子、200+ 可旋转键）。**R11/G10：
  默认刚性对接**（ESMFold 单模型只有一个构象，200+ 扭转的搜索空间只拖慢
  不增益；且该 vina 构建对大配体 `--cpu 64` 实测仍单核，见 ROUNDLOG R10/G10 证据），
  `options.vhh_dock_flex=true` 可开柔性（慢 1 个数量级以上）。本 vina 构建的
  柔性配体必须带 AD4 分子图关键字（ROOT/ENDROOT + TORSDOF），obabel 产物默认没有，
  `to_pdbqt(flex=True)` 已自动补齐（并修正 obabel 丢失的元素列）；PDBQT 缓存
  带 flex/rigid 一致性检查，模式不匹配自动重转。大配体自动降 exhaustiveness
  （<100 原子用 8，否则 1）。VHH 对接打分是"结合能力"的粗筛，不等于亲和力测定。
  **R11/G10-v2：CDR 片段对接**（`vhh_dock_cdr_only`，fast 默认开）——基准
  （`scripts/bench_vhh_dock.py`）显示对接成本 ~O(n^1.9) 于原子数（773 原子
  ~80-100 min，257 原子 6.25 min），且全长 VHH 分数由"撞墙"罚分主导
  （2.1e8 kcal/mol 量级）。故 fast 模式只对接 CDR/loop 片段：按残基 pLDDT<50
  的连续区切 1-3 个片段（pad 2 残基），综合分取最佳片段分（`fragment_scores`
  保留各片段分），实测 ~15× 提速。full 模式默认全长对接（更保守）。
- **VHH pLDDT 门槛（R11/G9）**：ESMFold 对 VHH 的 pLDDT 因 CDR3 loop 无序而
  集中在 30-35（100 条实测 p50=31.5、p90=34.5、>45 仅 1 条），旧 45/70 门槛
  让 fast 模式只剩 1 个对接样本。现 fast 35 / full 50（`vhh_plddt_min` 可覆盖），
  fast 屏 80 条。
- **小分子库回退**：DTP/PDBBind 镜像不稳（R11 实测：pdbbind.org.cn 的两个
  tar.gz URL 已 404，返回 138 字节 404 页；下载现需注册）。`resolve_library()`
  在 dtp/pdbbind 缺失或 <1MB 时自动回退 chembl35_small（50k），state 与报告
  标明 "fallback for dtp"；`setup --libraries pdbbind` 软失败不再打断整体
  setup，坏 tarball（<1MB）自动丢弃。
- **结构坑预检**：`analyze_pdb` 会扫描多 MODEL 叠加（NMR 系综/多构象）、备选位置
  altloc、缺失残基（含 SEQRES 对比）、金属离子、核酸、末端无序区等坑并给出建议
  action；`repair_structure(actions=[...])` 按 action 修好 PDB。`run_target_prep`
  会自动修 dedupe_models / keep_altloc_a 两个安全项，其余坑交由 agent 判断
  （decisions 留痕），报告第 1 节"结构预检"展示全部发现与修复。
- **MD 柔性诊断（R1）**：`analyze_replicas` 除整体 RMSD/RMSF/Rg/聚类外，还会：
  (1) 用 `gmx make_ndx splitch` 把 Protein 按拓扑分链（GROMACS TPR 会丢 PDB 链 ID，
  所以不能按链字符选组，`chain.ndx` 在 prepare 时自动建一次并复用），分别算
  **链间相对 RMSD**（fit 最大链）与**链内自拟合 RMSD**（fit 自身，不受参考链
  柔性/界面松弛放大）；(2) 用 MDAnalysis 逐帧做 DSSP 式氢键二级结构分类
  （O..N < 4.5 Å 且 C-O..N 角 > 100°，按 1HVI 晶体标定：69% 结构化 ≈ 已知
  β-三明治含量；全残基对、双向氢键，覆盖 β-三明治 30+ 的序列 offset），
  输出逐帧结构化占比 + 每残基稳定性；(3) 规则引擎把上面指标 + RMSF + 聚类
  汇总成中文"柔性解读"（整体稳定 / 链间域运动 vs 去折叠 / 高柔性区 /
  无主导构象态 / 二级结构丢失），写入 md summary 与报告第 5 节。
  `gmx_analyze` 新增 `kind=chain` / `kind=ss`。
  注意：MD 直接从 EM 起步，前 ~100-200 ps 属松弛段（Rg/链间 RMSD 会先跳变），
  解读时建议忽略早期段。
- **刚性漏斗局限**：默认对接/设计基于单一刚性构象，柔性只在 MD 阶段采样；
  对高度柔性靶点，高 RMSD 可能是构象采样/结构域运动而非失稳 —— 报告第 5 节的
  "柔性解读"会自动区分这两种情况（分链自拟合 RMSD + 二级结构 + RMSF + 聚类）。
- **MD 金属离子（R3）**：含金属的复合物（ZN/FE/MN/NI/CU/CO/MG/CA，单原子
  残基）进 MD 时会被误当成配体/蛋白残基，现自动处理：
  (1) `find_metals` 从复合物 PDB 检出金属（蛋白/配体拆分时剔除）；
  (2) 金属**并入所属蛋白链分子**（追加为该链 itp `[ atoms ]` 的最后一个
  原子，atomtype 由 `metal_{el}.itp`（仅 `[ atomtypes ]`，include 在力场
  include 之后）提供；电荷默认 +2 —— 由 genion 中和；sigma/epsilon=0，
  仅库仑作用，v1 局限）。之所以不做独立分子：GROMACS 的约束 section 只在
  **紧邻其前的 moleculetype 上下文**里解析，独立 ZNION 分子无法对蛋白原子
  下约束（grompp: "atom index out of bounds (1-1)"）；
  (3) 配位几何用 `[ distance_restraints ]`（写在所属链 itp 末尾，索引是
  **分子内相对**的；GROMACS 官方平坦底部势：r0–r1 = 晶体距离 ±0.15 Å 内
  零力，r1–r2 谐波，r2 以外线性；k = disre-fc(=2000) × fac；每条约束独立
  index 避免被当成 NOE 系综平均，type'=2 关闭时间平均）把金属与晶体内
  ≤2.5 Å 的蛋白 N/O/S 供体原子连起来。注意 GROMACS 2023 里没有
  `[ pair_restraints ]` 这个 section（grompp 报 Invalid directive）；
  MDP 需 `disre = Simple`（否则 grompp 静默丢弃约束）且 `nstdisreout = 0`
  （本 GROMACS 是 thread_mpi 构建，多 rank 下 disre 对输出不支持/EM 会
  段错误，故 EM mdrun 也带 `-ntmpi 1`，与生产 MD 一致）。
  拓扑索引映射 `_donor_index_map`
  用**坐标近邻+元素校验**匹配 PDB 供体 ↔ protein.gro 原子（pdb2gmx 会把 C
  端 O 改名 OC1/OC2、GRO 无链 ID、GRO 是 nm 单位 —— 名称/序号匹配都会踩坑，
  任一供体匹配不上就整体跳过约束并告警，不阻塞构建）。
  构建路径：`build_system` 两条分支（纯蛋白 / 蛋白+配体）都支持；combined.gro
  按拓扑分子顺序重排（各链 + 该链金属 + pdb2gmx 保留的晶体水 + 配体），
  与 `[ molecules ]` 顺序严格一致。无金属时行为与 R3 前完全一致
  （回归保护：1HVI 无金属构建 + 完整 ligand 路径 e2e）。
  已知局限：v1 只支持单原子金属残基（HEM 等配位环、多原子簇不处理，会按普通
  残基走）；+2 电荷对所有金属一刀切；约束是软性的，长时间 MD 下配位可能松弛。
  R2 起可用 `dock_conformer_set` 做**多构象+侧链柔性对接**缓解此局限（见下）。
- **多副本构象系综对接（R7）**：`pool_representatives` 跨所有副本
  汇集 GROMACS 聚类代表构象（每副本 top-3、占比 ≥ 5%），按占比降序
  贪心选取，**Cα Kabsch RMSD < 1 Å 判为冗余**（同一构象态在不同副本
  的重复出现）跳过，输出构象选择对接的受体系综。接入
  `dock_conformer_set`：crystal + 各副本簇代表（`md_r{rep}_c{cluster}`）
  逐一 Kabsch 对齐到晶体系 → PDBQT → 逐构象建 flex 文件 → Vina
  `--flex` 柔性对接 → `consensus_stats`（均值 = 共识值）。两个 R6 遗留
  坑的修复（实测固化）：
  (a) **代表结构提取必须用 `md_all.xtc`（合并轨迹）**——自动延长后
  聚类跑在合并时间轴上，簇首帧时刻可能超出原 md.xtc 的终点，从原
  轨迹 `-dump` 会报 "time outside trajectory"（静默回退到仅骨架，
  丢失侧链柔性）；
  (b) **trjconv -dump 的组名必须系统感知**——DNA-only tpr 没有
  "Protein" 组（GROMACS 按 moleculetype 名派生），硬编码 Protein
  提示会让每个 dump 失败回退；`_dump_group` 经 `make_ndx` 读 tpr 组
  列表后按 Backbone > DNA > RNA > Protein > System 选组，DNA 体系
  用 DNA 组（含碱基/糖/磷酸全原子，`full_atoms=True`）。
  e2e（1BNA 2 副本 × 0.5 ns）：池化代表结构全部 `full_atoms=True`
  （OP1/P 等磷酸碱基原子在场），占比降序，两两 Cα RMSD ≥ 1 Å
  （去重不变量）。R8 起代表结构默认取**簇均值结构**（MDAnalysis 对簇内
  全部帧逐原子取均值，见下），MDAnalysis 缺失时回退 `trjconv -dump` 首帧。
  已知局限：去重用 Cα 对 DNA 退化为骨架比较（可接受，去重只要求保守）。
  **R8 实测（1HVI smoke，3 构象池：晶体 + 2 个 MD 代表）**：flex 对接
  （flex，exhaustiveness 16，21 CPU，每构象 ~45 min）——晶体
  **-24.92**（flex 复现一致）；MD 代表：PBC 帧对齐修复**前**
  （全局 Kabsch 糊拟合）-8.29/-10.09，修复后（Cα 伙伴锚定）
  **-10.46/-13.58**——同一代表构象、同一参数，对齐修复净收益
  +2.2/+3.5 kcal/mol；consensus 均值 -14.44 → **-16.32**。MD 代表
  与晶体的差距（~11–14 kcal/mol）符合预期：代表构象口袋与晶体差
  ~4 Å，结合态偏好晶体几何；真正的构象选择价值要在"晶体失配"的
  配体上才能体现（见缺口）。池保留两个 MD 代表（原始 Cα RMSD 38 Å
  > 1 Å 去重阈值——两个不同盒子帧的代表天然"不同"，去重应对齐后
  再比，列入缺口）。
  （附）多副本并行踩坑：`run_replicas`/`extend_replicas` 的 worker
  函数必须是**模块级函数 + 纯数据 payload**——嵌套闭包经 joblib/
  cloudpickle 会按值序列化模块全局（loguru logger），而 pytest 的
  流捕获使 logger 持有不可 pickle 的捕获流（`EncodedFile`），单副本
  串行路径永远发现不了，`reps=2` 才炸。
- **MD 收敛判定与自动延长（R6）**：生产 MD 跑完先做收敛判定
  `_md_converged`（三个条件同时满足才算收敛）：
  (1) 轨迹长度 ≥ `md_converge_min_len_ps`（默认 50 ps）；
  (2) **RMSD 平台**——末 20% 时段的 RMSD 均值与前 80% 的偏差不超过
  `md_converge_rmsd_drift_nm`（默认 0.5 Å）；
  (3) **主导构象态**——最大簇占比 ≥ `md_converge_min_cluster`
  （默认 50%）。未收敛则 `extend_replicas` 给每个副本续跑
  `md_extend_ns`（默认 = `md_ns`）再重分析，最多 `md_max_extensions`
  轮（默认 2）。结果记录在 `md.json` 的 `extensions` 字段
  （rounds/converged/reason/total_ns）。
  续跑机制（2023.1 实测固化，坑很多）：
  (a) **延长 tpr 的 nsteps 是总步数**（原 + 延长）——`mdrun -cpi`
  把 tpr nsteps 当作目标终点步，只写延长步数会报 "checkpoint has
  already reached step N"；
  (b) mdrun 必须 `-cpi <末端 cpt> -noappend`：cpt 内部记录了原运行的
  输出文件名，deffnm 不同不加 `-noappend` 会拒绝启动；
  (c) `-cpi` 续跑的 xtc **时间轴连续**（如 500→1000 ps）且写成
  `<deffnm>.part000N.xtc` 分段文件；若不用 `-cpi` 而是新跑（grompp
  -t），xtc 从 t=0 重新开始；
  (d) **trjcat 按时间轴去重**——丢弃 ≤ 前一段末时刻的帧，所以新起
  t=0 的延长文件会被整体丢弃（0 帧保留）；用 `-cpi` 的连续时间轴
  正好只去重接点重复帧；
  (e) 延长 mdp 强制 `gen_vel = no` 继承 cpt 速度（避免接点处重新
  热化的速度跳变）。
  分析端自动优先读 `md_all.xtc`（合并轨迹）：RMSD/RMSF/Rg/配体
  RMSD/聚类/分链/SS 全部在合并后的连续时间轴上计算，烧入剔除仍只
  作用于生产段起始。e2e：1BNA 0.5 ns 生产 → 两轮 0.5 ns 延长 →
  合并轨迹 0–1500 ps（151 帧，接点帧正确去重），分析时间轴到
  990 ps 以上。已知局限：收敛判定是启发式（RMSD 平台 + 簇占比），
  不做自由能面分析；多副本各自独立延长（不共享随机数种子）。
- **MD 平衡段与烧入剔除（R5）**：EM 之后、生产 MD 之前插入两段
  平衡（每个构建一次，所有副本从平衡末端分叉）：
  (1) **NVT**（默认 50 ps，v-rescale 300 K）+ **NPT**（默认 100 ps，
  300 K / 1 bar，Berendsen，盒子弛豫），均带 pdb2gmx 生成的
  `posre_*.itp` 位置约束（蛋白/核酸重原子 1000 kJ/mol/nm²，力常数
  在 itp 段内逐原子定义）；(2) 生产 MD 从平衡末态
  （坐标+速度+盒子，`-t eq_npt.cpt`，`gen_vel=no`）开始。
  分析端 `burn_in_ps`（默认 100 ps）通过 gmx `-b` 从 RMSD/RMSF/Rg/
  配体 RMSD/聚类中剔除生产段起始的弛豫漂移。
  本轮踩坑固化：
  (a) **`posre-fc` 不是 MDP 参数**（grompp 报 "Unknown left-hand"）——
  位置约束力常数只在 posre itp 的 `[ position_restraints ]` 段；
  (b) **pdb2gmx 的 posre 文件没有 `[ moleculetype ]` 头**，其
  `[ position_restraints ]` 段按"紧跟的分子类型上下文"解析——必须
  把 `posre_X.itp` 的 include 紧跟在 `topol_X.itp` 之后（放最后会
  挂到 SOL 上 → "Atom index out of bounds"）；
  (c) **pdb2gmx 的 posre 侧文件写到 CWD** 而非 -p 同目录 →
  run_cmd 传 `cwd=build`；
  (d) 位置约束+压耦合需 `refcoord-scaling = all`，grompp 否则警告
  工件；
  (e) **`tc-grps = System` 而非 `Protein`**：GROMACS 默认 Protein
  组按 moleculetype 名派生，`DNA_chain_*` 不匹配 → DNA 体系
  "Group Protein not found"；
  (f) grompp 有 `-t` 文件且 `gen_vel=yes` 时速度仍会重新生成
  （-t 的速度被忽略）→ 副本分叉时强制 `gen_vel=no`。
  e2e：1BNA 构建 → NVT(1ps)→NPT(2ps) 平衡 → 1 ns 生产副本（自平衡
  态分叉）→ 0.5 ps 烧入剔除后分析全链通过。R8 起平衡/生产压控统一
  C-rescale（tau_p 2 ps）；自动延长/多构象系综策略见缺口清单。
- **MD 核酸（R4）**：`find_nucleic_acids` 检出含核酸残基的链（3 字母
  DA/DT/DC/DG 判 DNA；单字母 A/C/G/T 按 O2' 有无判 RNA，U 判 RNA）。
  构建规则：
  (1) **纯蛋白路径**（`is_ligand=False`）：DNA/RNA 链与蛋白链一起交给
  pdb2gmx 原生参数化（amber99sb-ildn 自带 b-DNA/RNA 参数；pdb2gmx 为
  每条核酸链生成 `topol_DNA_chain_X.itp` / `topol_RNA_chain_X.itp`，
  分子名 `DNA_chain_X`）；
  (2) **蛋白+配体路径**（`is_ligand=True`）：`_classify_chains` 把
  **≥10 个残基的核酸多聚体链**归到蛋白侧（原生参数化，比 GAFF2 物理
  上更正确），单残基核苷酸配体（NAD/ATP/FAD 等）仍走 ACPYPE/GAFF2；
  (3) R3 的金属集成对核酸链同样生效：`_chain_itp_stats` 按 itp 内
  `[ moleculetype ]` 的分子名解析（不再假定 `Protein_chain_*` 文件名），
  金属可以并入 DNA 链分子、约束供体集合扩到磷酸氧（OP1/OP2 或
  O1P/O2P，两种 PDB 命名都认）、糖环 O（O4'/O2'/O3'）与碱基 N/O
  （N1/N3/N7/N9、O6/O2）。
  原子类型去重：FF 的 `[ atomtypes ]`（ffnonbonded.itp）已显式定义的
  金属（CA/MG/CU 等）不再写自定义 `metal_*.itp`（grompp 会报
  "Atomtype ... defined previously"），直接用 FF 参数；仅 ions.itp
  隐式定义的（如 ZN）或完全没定义的（FE/MN/CO 等）才写自定义 itp。
  e2e 覆盖：1BNA 双螺旋单独构建、1HVI+1BNA（重标 C/D 链）混合构建、
  1BNA+合成 MG（约束落进 DNA 链 itp）。已知局限：反离子仍是
  genion -neutral 的单电荷离子（多聚阴离子 DNA 上 Mg2+ 更物理，未做）；
  RNA 2'-OH 的 O2' 判型依赖 PDB 含 O2' 原子。
- **柔性对接（R2）**：`make_flex_receptor` 从受体 PDBQT 里挑出"配体 cutoff 距离内
  残基的侧链"，按本 vina 构建（AD4 谱系）的 flex 格式写出 `BEGIN_RES/END_RES`
  块（每残基一个 AD4 扭曲张量树，从 3D 连通性推导，所有侧链二面角均可旋转）；
  `dock_conformer_set` 把同一配体对接到 **crystal + MD 聚类代表构象**（`--flex`
  侧链柔性），返回每构象分数 + **consensus 均分**（单构象异常不主导结论）。
  要点：flex 原子按**坐标**锚定到刚性受体，所以每个构象都要用自己的 PDBQT 重建
  flex 文件（工具已自动做）；本构建的 flex 解析器只认 `BEGIN_RES` 顶层标签
  （裸 ATOM 列表或 MODEL 都会报 "Unknown or inappropriate tag"）；受体 PDBQT
  必须用 obabel 的 AD4 列布局（不能直接过滤原始 PDB）。柔性对接比刚性慢数倍，
  建议只对少量命中/参考配体使用。
  - **MD 构象 = 全原子**：`cluster_representatives` 用 `gmx trjconv` 从轨迹按
    聚类首帧提取 Protein 全原子 PDB（`clusters_r{rep}.pdb` 的 MODEL 只是
    backbone-only 的**聚类质心**，与任何单帧都不重合，不能直接对接）。
  - **对齐**：MD 构象在模拟盒坐标系，对接盒子在晶体坐标系，`kabsch_transform`
    /`align_pdb_to_reference` 负责超叠：按 (链,残基名) 跨**所有残基编号偏移**
    收集 Cα 配对（吸收 ESMFold/GROMACS 重新编号/合并链，如 binder 拼到靶点链上
    任意起始号），再做**双射链分配**（同源二聚体 A→A 与 A→B 计数相同，按
    分配后的整体 RMSD 择优，防止对称拷贝污染拟合）+ 自适应离群剔除（拟合一致
    时才按 3 Å 修剪；构象差异大时保留全量，避免退化成 1-2 对）。
  - **已发现的坑（R2）**：`binder._ca_sequence` 曾用 `name[0]` 当单字母码，把
    ASP/GLU/LYS/ASN/ARG/GLN/TRP/TYR 全映射错 → ESMFold 打分的靶点序列被腐蚀
    （ARG→ALA 等 32 处）→ MD 体系实际跑在"ESMFold 重折叠的异序列靶点"上，
    与晶体折叠 RMSD ~13-27 Å，MD 构象对接分数因此偏低（见 R2 e2e 数字）。
    已修为完整三字母→单字母表；重跑 MD 后代表构象会回到正确序列/折叠。
- **GPU**：可选。当前 VM 无 `/dev/nvidia*` 时全部 CPU 运行（64 核足够）。
- **GNINA**：ELF 二进制依赖 `env/lib/python3.12/site-packages/nvidia/*/lib`（脚本已自动处理 LD_LIBRARY_PATH）。
- **力场**：GROMACS 2023.1 自带 amber99sb-ildn；配体用 ACPYPE GAFF2（AM1-BCC）；
  agent 可用 `gmx_env` 查看可用力场并选择。
- **DTP**：dtpbase.org 不稳定时自动回退 PDBBind / 本地 ChEMBL35。

- **R9 改进**：
  - **去重 PBC 感知**：`_ca_rmsd` 在 CRYST1 盒子存在时用 `_mi_kabsch`
    + 逐原子最小镜像残差（两个 MD 代表若处于不同盒子包裹，裸 Kabsch
    会把整数盒子向量平均成糊值，池化去重的 1 Å 阈值即被帧差污染）；
    无 CRYST1 时退化为裸 Kabsch（同帧比较不受影响）。
  - **平衡段幂等键含 MDP 指纹**：`run_equilibration` 的 NVT/NPT 复用
    判据加入 MDP 内容哈希（`eq_*.mdp.fp`，去注释行）——模板改动
    （如 barostat、时长）强制重跑对应平衡段；R8 前已完成的平衡段
    （无指纹文件）一次性重跑。
  - **options 真正生效**：`resolve_defaults(options)` 把 CLI/agent
    options 叠加到 Defaults（已知字段生效、未知键忽略、fast 先应用）；
    此前 `DEFAULTS.resolved(fast)` **静默丢弃**所有其他 option
    （如 `dock_exhaustiveness_final` 覆盖根本到不了工具）。
  - **池化对接成本**：MD 代表用 `dock_md_rep_exhaustiveness`
    （默认 12 / fast 8）对接——代表是 consensus 构象，不需要晶体级
    采样精度；flex 侧链 cutoff 改为 `flex_cutoff_ang`（默认 5 Å）
    可配。实测（1HVI，37 柔性残基）：exh 8 比 16 只省 ~10-15% 墙钟
    （flex 评估主导成本），但 r1 代表 best mode 从 -13.58 掉到
    -8.15 kcal/mol（r2 不变 -10.46）——默认取 12 折中；需要精度时
    直接覆盖该 option（R9 起 options 真正生效）。
  - **金属电荷（PDB 79-80 列）**：`find_metals` 读取 PDB "charge on
    atom" 字段（带符号/小数才采信，裸整数视为元素序号防误判）；
    未填时回退元素默认（METAL_PROPS，+2）。
  - **单原子金属"配体"不再丢失**：整条链都是单原子金属的"配体"
    （对接进来的 Zn2+/Mg2+）不再走 ACPYPE（单原子 mol2 脆弱），
    改走 standalone ion 通道（`_append_embedded_ions`）；同时修
    复多离子同元素时 [ molecules ] 按元素只写一行导致拓扑少数的
    问题（现在按离子逐个写）。
  - **死代码清理**：移除 `_core_kabsch` / `_pbc_compact_frame`
    （R8 末次重写后生产路径已由 `_mi_kabsch` + Cα 伙伴锚定取代，
    仅剩测试引用）。
  e2e（1HVI smoke）：池化对接 MD 代表按降档 exhaustiveness 8 重跑
  （vina 命令行路由正确；r2 -10.46 与 exh 16 一致，r1 -8.15，
  consensus -14.51——据此把默认从 8 调到 12）；smoke 全量复跑
  触发一次平衡段重跑（旧项目无 MDP 指纹）：NVT+NPT 以 **C-rescale**
  重平衡，盒子回到自然尺寸 **5.94×3.68×4.00 nm**（R8 的 6.21×
  3.84×4.18 nm 是 Berendsen 旧 MDP 的膨胀伪影）；生产段幂等复用
  （md.xtc 判据不变，生产盒子沿用旧值——列入缺口）。

## 缺口清单（R10 候选）
- **二价抗衡离子**：genion 只支持 Na/Cl；Mg2+/Ca2+ 作为抗衡离子（而非
  结构金属）时无法定量生成（`md_salt_m` 只算单价摩尔数）。
- **嵌入辅因子金属无协调约束**：HEM 的 Fe2+ 在 MD 中无 coordination
  restraints（GROMACS 无原生支持，需自定义 pair 势或 posre 到配位
  原子）；standalone ion 通道（R9）同样无约束。
- **HEM GAFF2 + Gasteiger（无 sqm）**：血红素电荷近似，长 MD 下 Fe 附近
  可能失稳（本 smoke 2 ns 未现）。
- **构象选择价值验证**：找一个晶体结合模式与 MD 代表口袋明显失配的配体，
  证明池化 docking 的 consensus 提升（1HVI 自身配体在晶体口袋，差距
  主要来自 MD 口袋构象差）。
- **多构象系综**：2 ns 单主导簇（99%）说明 1HVI 在此采样深度下构象
  单一；真正的多状态需要更长 MD 或副本系综（R6 自动延长已就绪但
  收敛判定即通过，未触发）。
- **vina 并行提示**：`--cpu` 超过物理核时 RSS 线性涨（37 柔性残基 +
  exh 16 时 ~23 GB），无上限提示。
- **flex 对接绝对成本**：37 柔性残基 + exh 16 ≈ 90–110 min/构象；
  晶体构象仍用 final exhaustiveness（16），MD 代表已降档（R9）。
- **池化去重的跨体系比较**：`_ca_rmsd` 现在 PBC 感知（R9），但跨
  不同盒子大小的体系（如重建后盒子变更）仍建议对齐到公共帧再比。
- **生产段幂等键不含 eq 状态**：R9 复跑中平衡段因 MDP 指纹重跑
  （C-rescale，新盒子 5.94×3.68×4.00 nm），但生产段沿用旧 md.xtc
  （旧 Berendsen 膨胀盒子 6.21×3.84×4.18 nm）——生产 tpr 的盒子
  与新 eq gro 不一致；修：生产幂等键加入 eq_npt.gro 的哈希。
