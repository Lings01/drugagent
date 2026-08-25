# DrugAgent — 药物发现一体化 Agent

**流水线是工具箱，LLM 是主程序。** DrugAgent 2.0 用一个 ReAct 主循环
（LLM + 原生 function calling）驱动整个药物发现流程：LLM 自己规划、调用
~41 个细粒度工具、自己写 MDP/选力场、失败时自己读日志排错重试；
4 个固定里程碑检查点 + LLM 随时发起的动态确认，`--auto` 时自动通过。

- 📖 **使用教程（零基础 → 独立使用）**：[TUTORIAL.md](TUTORIAL.md)
- 🏛 架构设计：[DESIGN.md](DESIGN.md)
- 📦 环境与坑的权威清单：[HANDOFF.md](HANDOFF.md)
- 📜 迭代历史（17 轮，每轮 计划/结果/反思）：[ROUNDLOG.md](ROUNDLOG.md)

## 技术栈

| 层 | 组件 | 版本/说明 |
|---|---|---|
| 语言/环境 | Python (conda) | 3.12，`env/` 独立环境 |
| LLM 主程序 | llama.cpp server（OpenAI 兼容端点） | 本地 `127.0.0.1:18080`，默认 `qwen3.8-27b-uncensored`（需 function calling；可用 `DRUGAGENT_LLM_BASE_URL`/`MODEL`/`API_KEY` 覆盖） |
| 靶点 | RCSB 下载、OpenBabel、ESMFold（vendor openfold + CPU 补丁） | 序列输入时自动建模 |
| 小分子 | RDKit、AutoDock Vina、GNINA 1.3.1 (ELF, CPU) | 库：DTP / ChEMBL35 / PDBBind（镜像不稳自动回退） |
| 蛋白设计 | RFdiffusion (hydra)、ProteinMPNN (vanilla v_48_010)、ESMFold + ESM2 | de novo binder / scaffold-guided VHH |
| MD | GROMACS 2023.1（自编译，amber99sb-ildn）、ACPYPE (GAFF2)、MDAnalysis 2.10 | 平衡段+多副本+自动延长+域级柔性诊断 |
| 报告 | Plotly（交互图）、3Dmol.js（3D 结构）、WeasyPrint（PDF） | `reports/report.html` + `.pdf` |
| 测试 | pytest | 218 快测 + 15 slow e2e |

硬件：64 核 CPU 即可（无 GPU 自动全 CPU）；磁盘 ~40 GB（环境+权重+库）。

## 模块（工具箱覆盖的功能）

| 模块 | 内容 | 关键工具 |
|---|---|---|
| A 靶点准备 | PDB 文件 / PDB ID / FASTA / 裸序列输入；完整性分析 + 结构坑预检（多 MODEL/altloc/缺失残基/金属/核酸/无序末端）；agent 判定；自动+手动修复 PDB；清洗；口袋检测；PDBQT 转换 | RCSB, obabel, ESMFold（仅序列输入） |
| B 小分子筛选 | 大库（DTP/ChEMBL/PDBBind 或自定义 SDF）→ 标准化 → 理化/ML 预过滤 → Vina 并行对接 → GNINA 复打分 → agent 定命中标准（参考配体重对接做阳性对照）；**柔性靶点工作流（R2/R5：MD 构象系综 → 多构象选择 + 侧链 `--flex`，consensus 平均分）** | RDKit, Vina, GNINA |
| C binder 设计 | RFdiffusion 从头设计 + ProteinMPNN 序列 + ESMFold 单体/复合物打分（界面 pLDDT）+ 几何接口度量（min 距离/接触对） | RFdiffusion, ProteinMPNN, ESMFold |
| D 纳米抗体 | 轨道A：VHH 库 → ESMFold 建模 → pLDDT 过滤（fast 35 / full 50）→ **刚性对接**（fast 默认 **CDR 片段对接**，~15× 提速，自适应盒子）→ 并行筛选；轨道B：RFdiffusion scaffold-guided 从头设计 + scaffold 保真度（`scaffold_rmsd_a`）+ 综合评分 | ESMFold, RFdiffusion, Vina |
| E MD 模拟 | 体系选择（配体/hit/binder/vhh/**apo**；修饰残基自动偏好 apo）→ pdb2gmx + ACPYPE → 溶胀/加离子/EM → **NVT→NPT 平衡（位置约束，C-rescale）+ 烧入剔除** → N ns × R 副本（自平衡分叉）→ **收敛判定 + 自动延长** → RMSD/RMSF/Rg/聚类 + **柔性诊断**（分链自拟合/二级结构/柔性区定位/**结构域 RMSD**/**域 vs-rest + 直径归一**/**刚性基线对照（R17：√t 折算 + 域 norm 基线）**/compact-unwrap 防紧凑盒 flapping）+ 金属离子协调 + 核酸原生参数化 + 辅因子/血红素自动参数化 | GROMACS 2023.1, ACPYPE, MDAnalysis |

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
env/bin/python -m drugagent.cli status --project projects/<项目>
env/bin/python -m drugagent.cli report --project projects/<项目>
```

零基础用户请从 [TUTORIAL.md](TUTORIAL.md) 开始（含逐参数解释、报告读法、
数字可信度指南、FAQ、术语表）。

常用 `run` 参数：`--library dtp|chembl35|pdbbind|<SDF路径>`，`--n-jobs 32`，
`--md-ns 100 --md-reps 3`，`--max-steps 300`（agent 步数预算），
`--no-llm`（确定性脚本模式），`--llm-base/--llm-model`（覆盖 LLM）。
MD 细调：`--md-salt`（离子浓度 M）、`--md-divalent MG --md-divalent-m 0.01`
（二价抗衡离子）、`--md-extend-ns`（自动延长步长 ns）、
`--md-max-extensions`（自动延长轮数）、`--md-burn-in-ps`（烧入段剔除）、
`--md-system`（强制 MD 体系）。VHH：`--vhh-plddt-min`、`--vhh-dock-flex`、
`--vhh-dock-cdr-only`。完整参数表见 TUTORIAL §9。

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
│   ├── calibration/                  # 标定用 PDB（1UBQ/1GFL/1CPS/1M17）
│   └── tools/
│       ├── gromacs/                  # GROMACS 2023.1 (自编译, amber99sb-ildn)
│       ├── vina/                     # autodock-vina
│       ├── gnina/                    # GNINA 1.3.1 ELF (CPU)
│       ├── RFdiffusion/              # RFdiffusion (hydra) + models/ 权重 + mpnn/
│       ├── vhh_scaffolds/            # 1EWN scaffold + secstruc 邻接 + NOTE.md
│       └── 3Dmol/                    # 3Dmol-min.js
├── drugagent/                        # 主包
│   ├── agent/                        # 2.0 核心: loop/提示词/41 个工具
│   ├── vendor/openfold/              # aqlaboratory/openfold + CPU 补丁
│   ├── modules/                      # A-E 五个模块（工具后端）
│   ├── graph.py                      # 1.0 LangGraph 状态机（兜底工具复用）
│   ├── cli.py                        # typer CLI (run/resume/rerun/status/report/setup)
│   └── report/                       # 交互式 HTML + PDF (WeasyPrint)
├── TUTORIAL.md                       # 使用教程（零基础）
├── DESIGN.md                         # 2.0 架构设计
├── HANDOFF.md                        # 环境与坑、当前状态、验证命令
├── ROUNDLOG.md                       # 迭代日志（每轮 计划/结果/反思）
├── projects/                         # 每次运行一个目录 (01_target…05_md, reports/, agent/)
├── tests/                            # pytest 套件 (218 快测 + 15 slow e2e)
└── logs/                             # 构建/运行日志
```

## 运行状态与断点

- `status`（R11 增强）：每个阶段的完成状态 + 关键数字（对接数/命中数/
  设计数/MD ns 与 final RMSD/库名含回退标注）、磁盘上的阶段 JSON 产物、
  最近 3 条工具失败（state.errors + transcript 里 ok=false 的调用）、
  transcript 最后一条；`decisions.json` 记录每个判断与依据。
- `resume --project DIR`：重放 `agent/transcript.jsonl`，从断点继续
  （产物幂等 + 阶段级复用，R11/G8）。
- `rerun --project DIR --stage X`：强制重跑单个阶段（G8；其余产物不动）。
- 报告：`projects/<项目>/reports/report.html`（Plotly 交互图 + 3Dmol 结构）
  与 `report.pdf`。

## 测试

```bash
env/bin/python -m pytest tests/ -m "not slow" -q --basetemp=$PWD/data/fixtures/_ptmp   # 快速单测 (218 用例)
env/bin/python -m pytest tests/ -q --basetemp=$PWD/data/fixtures/_ptmp                  # 含 slow（需 vina/GROMACS/RF 权重）
```

> slow 构建 e2e 在本机 /tmp（tmpfs）上偶发失败：有后台清理进程会在
> gromacs 运行中删掉 pytest 的 tmp 目录。代码本身稳定，遇到偶发失败用
> 上面 `--basetemp` 指向本地盘重跑即可。

## 1HVI 全模块端到端（参考）

```bash
env/bin/python -m drugagent.cli run --target 1HVI --modules all --fast --auto
```

在 64 核 CPU 上约 3–5 h 完成全模块：靶点（1HVI 二聚体，99 残基/单体）→
小分子筛选（标准化 + Vina 对接 + GNINA 复打分 + agent 定命中）→
binder 设计（RFdiffusion 从头 + ESMFold 复合物打分）→ VHH 双轨 →
体系选择（agent）→ GROMACS EM + 5 ns×3 副本 MD →
RMSD/RMSF/Rg 分析 + HTML/PDF 报告。完整实例见 `projects/r10_e2e/`
（含 5 个结构域的铰链信号分析）。

## 已知注意事项（摘要）

完整细节见 HANDOFF.md / ROUNDLOG.md 对应轮次；此处只列与使用相关的。

- **单位**：GROMACS 的 rms/gyrate 输出为 nm（rmsf 为 Å，cluster 截断 1.5 为 nm）；
  报告展示时 RMSD 已换算为 Å（×10），Rg 保持 nm。MD 分析链路已统一 nm 口径
  （R15 修复过一处 Å→nm 全链路 ×10 bug，有回归测试）。
- **紧凑盒子 PBC**：gmx 分析前自动 compact-unwrap（R17；蛋白接近盒尺寸时
  wrapped 坐标跨盒断裂会制造 RMSD 假突刺，GROMACS 论坛同款案例）。
- **辅因子/血红素**：MD 自动拆链 + ACPYPE/Gasteiger + 嵌入金属独立离子；
  已知近似：GAFF2/Gasteiger 电荷（无 sqm/AM1-BCC），嵌入金属无配位约束
  （长 MD 可能漂出平面）。
- **MD 金属离子**：单原子金属残基（ZN/FE/MN/NI/CU/CO/MG/CA）自动并入所属
  链 + 晶体配位距离约束（平坦底部势）；v1 局限：+2 电荷一刀切、软约束。
- **MD 核酸**：DNA/RNA 链原生参数化（amber99sb-ildn）；单残基核苷酸
  （NAD/ATP/FAD 等）走 ACPYPE 配体路径。
- **MD 平衡/收敛**：NVT 50 ps + NPT 100 ps（位置约束，C-rescale——本构建
  下 Berendsen 生产段 ~20 ps 会 SIGSEGV）；生产后收敛判定（RMSD 平台 +
  主导簇 ≥50%），未收敛自动延长（默认最多 2 轮）。
- **VHH 对接**：默认刚性 + CDR 片段（fast）——全长 VHH 是"巨型配体"
  （773 原子 ~100 min，成本 ~O(n^1.9)）；片段自适应盒子消除撞墙罚分
  （2.3e7 → 145 kcal/mol）。对接分是粗筛，非亲和力测定。
- **pLDDT 口径**：VHH 文库 pLDDT 因 CDR3 无序集中在 30-35（门槛 fast 35/
  full 50 因此标定）；binder/VHH 设计 pLDDT 是 ESMFold 重折叠置信度，
  de novo 偏低属正常；几何接口度量（min 距离/接触对）是位姿级判别。
- **scaffold 内容**：`vhh_scaffolds/` 的 1EWN 实为人 AAG 糖基化酶核心
  （目录名 `vhh_` 是历史前缀）；scaffold-guided 只条件 SS+邻接，设计骨架
  对 scaffold 漂移 15 Å 量级是预期（`scaffold_rmsd_a` 字段量化）。
  详见 `data/tools/vhh_scaffolds/NOTE.md`。
- **小分子库回退**：dtp/pdbbind 镜像不稳（实测 404/138 字节坏 tarball），
  缺失或 <1MB 自动回退 chembl35_small（50k），state 与报告标明 fallback。
- **结构坑预检**：`analyze_pdb` 扫多 MODEL/altloc/缺失残基/金属/核酸/无序区
  并给建议 action；`run_target_prep` 自动修两个安全项，其余交 agent 判断
  （decisions 留痕），报告第 1 节展示全部发现与修复。
- **修饰残基（R17）**：靶点"配体"全是修饰残基（CPM 交叉链接/金属/辅因子，
  64 名集合）时 MD 自动偏好 apo 通路（残基过滤后跑纯蛋白）；GFP 色原体
  CRO 类（骨架原子缺失）仍需残基级修复（R18 缺口）。
- **刚性漏斗局限**：默认对接/设计基于单一刚性构象，柔性只在 MD 阶段采样；
  报告第 5 节的"柔性解读"会自动区分"整体漂移=构象采样/域运动"与
  "局部去折叠"（分链自拟合 + 二级结构 + RMSF + 聚类 + 域 vs-rest）。
- **磁盘**：环境+权重+库约 40 GB；MD 轨迹按 ns 增长（5 ns×3 副本约 2 GB）。
- **GPU**：可选。无 `/dev/nvidia*` 时全部 CPU 运行（64 核足够）。
- **GNINA**：ELF 二进制依赖 `env/lib/python3.12/site-packages/nvidia/*/lib`
  （脚本已自动处理 LD_LIBRARY_PATH）。
- **binder 序列来源**：RF 设计 PDB 残基名为 GLY（不写序列），内置
  ProteinMPNN 生成真实序列（MPNN 缺失/失败回退启发式，报告标注）；
  ESMFold 打分用真实序列。

## 缺口清单（R18 候选，完整讨论见 ROUNDLOG 第 17 轮反思）

- **后台标定收口**：crambin（第二刚性点）/ CDK2（柔性参照）MD 完成后
  进标定表；r10_e2e 自动延长（二聚体相对运动 5 ns 未平台）完成后确认。
- **收敛判据与 unwrap 耦合**：用链自拟合（而非整体 RMSD）判"内部柔性
  漂移 vs 寡聚体相对运动"。
- **GFP/CRO 类修饰残基**：残基级修复（补缺骨架 O 或删受影响残基）。
- **二价抗衡离子**：genion 只支持 Na/Cl，Mg2+/Ca2+ 作抗衡离子无法定量
  （`md_salt_m` 只算单价摩尔数）。
- **嵌入辅因子金属无协调约束**（HEM 的 Fe2+，GROMACS 无原生支持）。
- **构象选择价值验证**：需晶体结合模式与 MD 代表口袋明显失配的配体。
- **生产段幂等键不含 eq 状态**（平衡重跑后生产 tpr 盒子可能与新 eq gro
  不一致）。
