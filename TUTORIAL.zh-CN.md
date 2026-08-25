# DrugAgent 2.0 使用教程（零基础 → 独立使用）

> 读者假设：会基本命令行操作，不知道分子对接是什么也没关系。
> 环境假设：已在 `/home/data/lrs/drug/drugagent` 完成部署（`env/`、`data/` 都在）。
> 🇬🇧 **English version**: [TUTORIAL.md](TUTORIAL.md)
> 配套文档：[README.zh-CN.md](README.zh-CN.md)（技术栈总览）、[DESIGN.md](DESIGN.md)（架构）、
> [HANDOFF.md](HANDOFF.md)（环境与坑的权威清单）、[ROUNDLOG.md](ROUNDLOG.md)（迭代历史）。

---

## 0. 这套系统能帮你做什么

给它一个蛋白（PDB ID / PDB 文件 / 甚至只是一条氨基酸序列），它替你完成：

| 你想做的事 | 模块 | 一句话解释 | 你会拿到 |
|---|---|---|---|
| **找药**：哪些小分子能结合我的蛋白 | `screen` | 从化合物库（几十万个小分子）里虚拟筛选 | 排名靠前的候选分子 + 结合分数 + 结合姿态（3D 结构文件） |
| **设计新蛋白**：设计一个小蛋白去结合/抑制我的蛋白 | `binder` | 从头设计 60-80 残基的新蛋白（RFdiffusion） | 设计结构 PDB + 氨基酸序列 + 置信度打分 |
| **设计抗体片段**：设计一个纳米抗体（VHH） | `vhh` | 双轨：从 VHH 库里筛 + 从头 scaffold 设计 | 候选 VHH + 打分 + 设计保真度指标 |
| **看蛋白动态**：我的蛋白稳不稳？哪个区是软的？ | `md` | 分子动力学模拟（GROMACS）+ 自动化柔性诊断 | RMSD/RMSF/结构域分析 + 自动中文解读（刚性/铰链/去折叠） |

一次运行可以只跑其中一项（`--modules screen`），也可以全跑（`--modules all`，
全链路：靶点 → 筛选 → binder → VHH → MD → 汇总报告）。

**交付物**：每次运行产出一个项目目录 `projects/<名字>/`，里面最重要的是
`reports/report.html`——一个浏览器打开就能看的交互式报告（3D 结构 + 图表 +
自动解读文字）。

---

## 1. 环境检查（30 秒）

```bash
cd /home/data/lrs/drug/drugagent

# ① python 环境（应该打印 3.12.x）
env/bin/python --version

# ①-b 安装 drugagent 命令（一次性；editable 安装——改代码即时生效；
#     软链进用户 PATH 后，任意目录直接 drugagent ...）
env/bin/pip install -e . --no-deps
ln -sf $PWD/env/bin/drugagent ~/.local/bin/drugagent

# ② （可选）本地 LLM 是否在线——只用 --no-llm 脚本模式则不需要
curl -s -m 5 http://127.0.0.1:18080/v1/models | head -c 200

# ③ 环境是否完整（幂等：缺什么补什么，已有跳过；首次部署必跑）
drugagent setup

# ④ 代码健康检查（~2 分钟，应显示 218 passed）
env/bin/python -m pytest tests/ -m "not slow" -q -p no:warnings --basetemp=$PWD/data/fixtures/_ptmp
```

说明：
- `--basetemp=$PWD/data/fixtures/_ptmp` 必须带：本机 `/tmp` 是 tmpfs 且有
  后台清理进程，测试临时文件放本地盘才稳。
- 磁盘需求：环境+权重+库约 40 GB；MD 轨迹按模拟时长增长（5 ns×3 副本约 2 GB）。
- 算力：64 核 CPU 即可（GPU 可选，无 GPU 时全 CPU 自动运行）。

---

## 2. 第一次运行：小分子对接（screen）

### 2.1 命令

```bash
drugagent run \
  --target 1HVI \
  --modules screen \
  --fast \
  --auto --no-llm \
  --name mydock
```

参数逐个解释：

| 参数 | 含义 | 你可以换成 |
|---|---|---|
| `--target 1HVI` | 靶点。PDB ID 会自动从 RCSB 下载 | 任何 PDB ID（如 `5P21`）；本地文件路径（`/path/protein.pdb`）；FASTA 文件；裸氨基酸序列（自动 ESMFold 建模） |
| `--modules screen` | 只跑小分子筛选 | `binder` / `vhh` / `md` / `screen,binder,md`（逗号组合）/ `all` |
| `--fast` | 验证规模（小库、短 MD，分钟~小时级） | 去掉 `--fast` = 生产规模（大库、100 ns MD，成本高很多） |
| `--auto` | 检查点自动通过（不等人在终端前） | 去掉 = 交互模式，agent 会在里程碑停下问你（见 §8） |
| `--no-llm` | 确定性脚本模式，不依赖 LLM 服务 | 去掉 = LLM 主程序模式（本地 llama.cpp 需在线） |
| `--name mydock` | 项目名（目录 `projects/mydock/`） | 任意名字；缺省用时间戳 |

> 零基础建议：前几次都带 `--auto --no-llm`，最稳、可复现。

### 2.2 后台在干什么（按顺序）

1. **靶点准备**：下载/读取 PDB → 完整性检查（多 MODEL 叠加？altloc？缺失残基？
   金属？核酸？无序末端？）→ 自动修复安全项 → 清洗 → 口袋检测 → 转 PDBQT。
   产物：`projects/mydock/01_target/`（`receptor.pdbqt`、`pocket.json`、`clean.pdb`）。
2. **库解析**：默认用本地主库 `nci_npatlas`（约 26.5 万化合物——NCI/DTP 开放化合物库 ∪ NPAtlas 天然产物，InChIKey 去重）。
   请求的库缺失/损坏时按 `nci_npatlas` → `chembl35_small`（约 5 万）自动回退，state/报告标注 "fallback for …"。
3. **标准化 + 预过滤**：RDKit 标准化、理化性质过滤（Lipinski 等）。
4. **对接**：Vina 并行对接（默认 32 进程，`--n-jobs` 可调）。
5. **复打分**：GNINA 对 top 结果复打分。
6. **命中判定**：取 top hits（agent 模式会参考配体重对接做阳性对照并留痕）。
7. **报告**：生成 `reports/report.html` + `report.pdf`。

### 2.3 看进度

```bash
drugagent status --project projects/mydock
```

输出每个阶段的 ✓/✗ + 关键数字（对接数/命中数/库名含回退标注）+ 产物清单 +
最近 3 条错误。

### 2.4 拿结果

```bash
xdg-open projects/mydock/reports/report.html     # 浏览器打开（推荐）
ls projects/mydock/02_screen/                     # 对接明细文件
```

报告第 2 节是筛选结果：候选分子排名表（SMILES、对接分数、复打分）+ 3D 结合
姿态（可旋转缩放）。

### 2.5 怎么读对接分数

- 单位 kcal/mol，**越低越好**（负得越多越"结合得紧"）。
- 它是**相对排名工具**：同一批分子之间比较可靠；跨体系/跨参数比较要小心。
- 刚性对接分数 ≠ 实验亲和力（没有自由能微扰），用于"从 5 万里挑 50 个"
  的富集，不是"测出 Kd"。

---

## 3. 项目目录解剖

```
projects/mydock/
├── state.json              # 阶段索引（status/resume 读它；别手改）
├── 01_target/
│   ├── clean.pdb           # 清洗后的蛋白
│   ├── receptor.pdbqt      # 对接用受体
│   └── pocket.json         # 口袋信息
├── 02_screen/              # 对接/打分/命中明细（JSON + PDBQT 姿态）
├── 03_binder/              # (若跑了 binder) 设计结构 + scored.json
├── 04_vhh/                 # (若跑了 vhh)
├── 05_md/                  # (若跑了 md)
│   ├── build/              # 构建产物（top/gro/tpr/平衡段）
│   ├── md_rep1/            # 副本 1 轨迹（md.xtc）
│   ├── md_rep2/  md_rep3/
│   └── analysis/           # RMSD/RMSF/聚类 xvg 与解读
├── reports/
│   ├── report.html         # ← 主交付物
│   └── report.pdf
└── agent/
    ├── transcript.jsonl    # LLM 完整工具调用流水（"它是怎么想的"）
    └── decisions.json      # 每个关键判断 + 理由
```

**幂等性**：同一条 `run` 命令重跑不会重算已完成的东西（工具产物级 + 阶段级
双层复用），重启后继续只补缺失部分。要强制重算某阶段用 `rerun`（§7）。

---

## 4. Binder 设计（从头设计一个小蛋白）

### 4.1 命令

```bash
drugagent run --target 1HVI --modules binder \
  --fast --auto --no-llm --name mybinder
```

### 4.2 后台流程

1. 靶点准备（同上）。
2. **RFdiffusion**：在靶点表面从头设计 60-80 残基的小蛋白骨架
   （`miniprotein` 模式）。
3. **ProteinMPNN**：给骨架生成真实氨基酸序列（RF 设计的 PDB 里残基名是
   GLY 占位，MPNN 序列才是真实序列）。
4. **ESMFold 打分**：用真实序列重新折叠，算 pLDDT / 界面置信度。
   设计条目还有**几何接口度量**（设计链-靶点最小 Cα 距离、<6 Å 接触对）——
   序列打分对"设计到底贴没贴到靶点"是盲的，几何接触是位姿级判别。
5. 报告：设计排名表 + 结构。

### 4.3 耗时与结果

- `--fast` 下约 **20-40 分钟**（RF 在 CPU 上跑，是全流程最慢的模块之一）。
- 结果看 `projects/mybinder/reports/report.html` 第 3 节；
  设计 PDB 在 `03_binder/`（`design_N.pdb` 单体、`design_N_complex.pdb`
  复合物）——拖进 PyMOL/ChimeraX 即可看。
- **pLDDT 怎么读**：de novo 设计的 pLDDT 偏低属正常（ESMFold 没见过这个
  序列）；重点看**界面 pLDDT**（设计链与靶点接触区）和几何接触指标
  （"目标距离 (Å)" 列，⚠未接触 标记表示设计没贴到靶点，要警惕）。
- 缓存：`03_binder/scored.json` 按设计文件 mtime 缓存，重跑只重验证
  被覆盖的设计（rerun 从 ~3 h 降到分钟级）。

---

## 5. VHH 纳米抗体（双轨）

```bash
drugagent run --target 1HVI --modules vhh \
  --fast --auto --no-llm --name myvhh
```

两条轨道并行：

- **轨道 A（文库筛选）**：VHH 库（`data/libraries/vhh_library.fasta` 或
  合成库）→ ESMFold 建模 → **pLDDT 过滤**（fast 门槛 35 / full 门槛 50，
  `--vhh-plddt-min` 可覆盖）→ **刚性对接**（fast 默认只对接 **CDR 片段**：
  按 pLDDT 低值区切 1-3 个片段 + 每片段自适应盒子，比全长快 ~15×）→ 并行筛选。
- **轨道 B（从头设计）**：RFdiffusion scaffold-guided（1EWN scaffold 做折叠
  条件）→ ProteinMPNN → 打分。

### 结果与解读

- 报告第 4 节：候选表（对接分/片段明细/pLDDT）+ **文库 pLDDT 分布直方图**
  （门槛虚线 + p50/p90 + 过门槛条数）+ 综合排名。
- **scaffold_rmsd_a**（轨道 B 设计）：设计骨架对 scaffold 骨架的 Kabsch
  RMSD（Å）。scaffold-guided 模式官方语义是"SS 模式 + 块邻接条件，允许
  序列与精细结构变化"，所以 **15 Å 量级的漂移是预期行为**（折叠拓扑保留、
  局部结构保留、长度精确匹配），该字段量化的是近似保真度而非失败。
  scaffold 内容说明见 `data/tools/vhh_scaffolds/NOTE.md`（1EWN 实为人 AAG
  糖基化酶核心，目录名 `vhh_` 是历史前缀）。
- VHH 对接分同样是粗筛分数（"结合能力"排序，非亲和力测定）。

---

## 6. MD 模拟（蛋白柔性诊断）

```bash
# 带配体的复合物（靶点自带配体或指定体系）
drugagent run --target 1HVI --modules md \
  --fast --auto --no-llm --name mymd

# 纯蛋白靶点（apo，蛋白单独 MD——柔性/铰链基线）
drugagent run --target data/calibration/1UBQ.pdb \
  --modules md --fast --auto --no-llm --name myapo

# 自定义时长/副本（去掉 --fast 后默认 100 ns×3；常用覆盖）
drugagent run --target 1HVI --modules md \
  --md-ns 5 --md-reps 3 --auto --no-llm --name mymd5
```

### 6.1 后台流程

1. **体系选择**：模拟哪个复合物（输入配体 / 筛选 hit / 设计的 binder / VHH /
   apo 蛋白；"配体"全是修饰残基时自动偏好 apo；`--md-system` 可强制）。
2. **构建**：pdb2gmx（amber99sb-ildn）→ 配体 ACPYPE/GAFF2 → 溶胀 → 加离子
   （`--md-salt`，默认 0.15 M）→ 能量最小化。
3. **平衡**：NVT（50 ps）→ NPT（100 ps），带位置约束，C-rescale 压控。
4. **生产 MD**：N ns × R 副本（自平衡末态分叉的独立副本），v-rescale 300 K。
   **收敛判定 + 自动延长**：RMSD 平台 + 主导簇占比不达标时自动续跑合并
   （`--md-extend-ns` 步长、`--md-max-extensions` 轮数）。
5. **分析**（全部自动）：
   - 核心 RMSD（backbone，相对起始构象）
   - 分链 RMSD：**链间相对**（fit 最大链后测其他链）+ **链内自拟合**
     （不受参考链柔性放大，判断"去折叠"用这个）
   - RMSF（每残基柔性）+ **柔性区定位**（连续高 RMSF 段 → 残基区间）
   - 回转半径 Rg、gmx 聚类（1.5 nm 截断）
   - **二级结构逐帧占比**（氢键规则，DSSP 式）
   - **结构域 RMSD**（SS 提取的结构域，逐域 Kabsch 自拟合 → 域内部形变）
   - **域-其余蛋白相对 RMSD**（铰链/变构刚体运动信号）+ **直径归一化
     final_norm**（跨域尺寸可比）
   - **自动中文解读**（刚性基线对照，见 6.2）

### 6.2 怎么读 MD 结果（重要）

报告第 5 节 + `analysis/` 里的 xvg 曲线。**不要只看一个总 RMSD 数字**，
现在有一套标定过的判读框架：

| 指标 | 含义 | 怎么判 |
|---|---|---|
| 核心 RMSD（backbone 自拟合） | 整体构象漂移 | 对照**刚性基线**：泛素（76 残基单体，公认稳定）2 ns 末端 ≤2.1 Å、域级 norm ≤0.08；按 √t 折算到你的时长（上限 = 2.1×√(ns/2) Å）。超了 → 进一步看下面两项区分原因 |
| 链内自拟合 RMSD | 单条链内部形变 | <0.6 nm（6 Å）→ 链没去折叠；>6 Å → 内部去折叠/大重排 |
| 链间相对 RMSD | 多链相对运动 | 远大于自拟合 → 二聚体/寡聚体相对运动（如 HIV-PR 5 ns 实测 1.15-3.88 nm，而链自拟合仅 0.05-0.14 nm） |
| 域 vs-rest final_norm | 结构域相对其余蛋白的刚性运动（铰链信号） | >0.08 超刚性基线 → 真实域运动/铰链；内部稳定(<3 Å) vs 大重排，报告分两档 |
| 柔性区（连续高 RMSF 段） | 局部柔性 loop | 报告列出残基区间；结合口袋位置判断是否影响结合 |
| 聚类占比 | 构象态数 | 单一主导簇（>90%）= 构象单一；多簇 = 构象系综 |
| SS 占比变化 | 二级结构保持 | 64%→64% 稳定；明显下降 → 去折叠迹象 |

**时间尺度提醒**：`--fast` 的 MD 是 2-5 ns 的**稳定性验证**，不是生产级
自由能/收敛结论；要研究铰链动力学请用 full 模式（默认 100 ns）或
`--md-ns 50` 起。MD 前 ~100-200 ps 是松弛段（曲线开头的跳变属正常，
分析端已默认剔除 100 ps 烧入段）。

### 6.3 耗时参考（64 核 CPU）

- 76 残基 monomer（ubiquitin）：2 ns×3 ≈ **18 分钟**
- 198 残基 dimer（HIV-PR）：5 ns×3 ≈ **1.5 小时**
- 290 残基 + 配体（CDK2）：2 ns×2 ≈ 30-40 分钟
- 自动延长触发时每副本再跑 `--md-extend-ns`（默认 = `--md-ns`）

---

## 7. 断点 / 重跑 / 状态

```bash
# 状态（随时可跑）
drugagent status --project projects/mydock

# 断点续跑（崩溃/中断后；重放 transcript + 产物复用）
drugagent resume --project projects/mydock

# 只强制重跑某个阶段（其余产物不动）
drugagent rerun --project projects/mydock --stage md
drugagent rerun --project projects/mydock --stage vhh --no-with-report

# 只重新生成报告
drugagent report --project projects/mydock
```

常用 stage 名：`target_prep` / `screening` / `binder` / `vhh` / `md` /
`report`。

**什么时候用 `rerun`**：改了代码/参数想重算某阶段；某个阶段产物可疑。
**什么时候用 `resume`**：进程被杀/机器重启，想接着跑。

---

## 8. 与 Agent 对话（交互模式）

默认（**不带** `--auto`）是交互模式：agent 在里程碑停下等你的输入。

```bash
drugagent run --target 1HVI --modules all --fast
# 前台运行，会在这些点停下：
#   [里程碑检查点: target]  靶点预检摘要 → 批准/修改/中止
#   [里程碑检查点: screening]
#   [里程碑检查点: design]
#   [里程碑检查点: md]      （附 MDP 计划：力场/盒子/盐/barostat）
#   [动态确认]              （LLM 随时发起：如选 binder 类型、选哪个复合物）
# 提示符：回答/选择 > 
```

- 输选项序号（`1`）= 批准；输 `x` / `中止` = 中止本次运行；
- **也可以直接打自由文本**（如"盐浓度改成 0.1，barostat 用 c-rescale"），
  会作为你的意见回给 LLM，它据此调整并 `record_decision` 留痕。
- LLM 的每个判断都写在 `agent/decisions.json`（stage / 选择 / 理由），
  完整工具调用流水在 `agent/transcript.jsonl`（报告第 6 节有摘要）。

### 模式选择速查

| 场景 | 模式 |
|---|---|
| 批量 / CI / 无人值守 | `--auto --no-llm` |
| 可复现优先 | `--no-llm`（±`--auto`） |
| 想让 agent 按你的偏好走（选参数、挑命中） | 交互模式（去 `--auto`，LLM 需在线） |
| LLM 服务没起 | `--no-llm`（确定性脚本模式） |

### 参数主权（LLM 模式）

力场、MDP（积分器/barostat/盐/盒子 margin）、模拟哪个复合物、取几个 hits、
binder 类型——都由 LLM 决定并留痕；`--no-llm` 模式用默认值。

---

## 9. 全部参数速查（run）

| 参数 | 默认 | 说明 |
|---|---|---|
| `--root` | 部署 projects/ | 项目输出根目录：项目建在 `<root>/<name>`——`--root .` 即在当前文件夹下跑（环境变量 `DRUGAGENT_PROJECTS_ROOT` 可固定默认） |
| `--target` | （必填） | PDB ID / PDB 路径 / FASTA / 裸序列 |
| `--modules` | `all` | `screen,binder,vhh,md` 逗号子集 |
| `--fast / --no-fast` | full | 验证规模 vs 生产规模 |
| `--auto / --no-auto` | 交互 | 检查点自动通过 |
| `--name` | 时间戳 | 项目名 |
| `--library` | `nci_npatlas` | `nci_npatlas`（主库：NCI/DTP ∪ NPAtlas，约 26.5 万）/ `chembl35` / `chembl35_small` / `pdbbind` / 自定义 SDF 路径；缺失/损坏自动经主库回退 |
| `--n-jobs` | 32 | 对接/打分并行度 |
| `--md-ns` | 100 (full) / 5 (fast) | MD 时长 ns |
| `--md-reps` | 3 | MD 副本数 |
| `--md-salt` | 0.15 | 离子浓度 M |
| `--md-divalent` / `--md-divalent-m` | 无 / 0.01 | 二价抗衡离子（MG/CA）及浓度 |
| `--md-extend-ns` | = md_ns | 自动延长步长 ns（0=用 md_ns） |
| `--md-max-extensions` | 2 | 自动延长最大轮数 |
| `--md-burn-in-ps` | 100 | 生产段烧入剔除 ps |
| `--md-system` | agent 选 | `input_ligand`/`screening_hit`/`binder`/`vhh`/`apo` 强制体系 |
| `--vhh-plddt-min` | fast 35 / full 50 | VHH 文库 pLDDT 门槛 |
| `--vhh-dock-flex` | 刚性 | VHH 柔性对接（慢一个数量级） |
| `--vhh-dock-cdr-only` | fast 开 | CDR 片段对接（~15× 提速） |
| `--max-steps` | 300 | LLM 步数预算 |
| `--llm-base` / `--llm-model` | 本地 18080/qwen3.8-27b | LLM 覆盖（也可用环境变量） |
| `--no-llm` | 关 | 确定性脚本模式 |

---

## 10. FAQ / 排错

**Q1：我的蛋白带配体/金属/DNA，能跑吗？**
- 配体：screen/binder 直接用；MD 走配体参数化路径（ACPYPE GAFF2）。
- 金属（ZN/FE/MG…单原子）：MD 自动并入所属链 + 配位距离约束。
- 核酸链：MD 原生参数化（amber99sb-ildn 自带 DNA/RNA）。
- 血红素等辅因子：MD 自动拆链 + ACPYPE + 嵌入金属独立离子（电荷近似、
  无配位约束，长 MD 注意）。
- **已知会失败的**：蛋白内嵌**修饰残基**（GFP 色原体 CRO 这类——pdb2gmx
  需要其全原子参数）：MD 构建会失败（R18 缺口）；对接/binder 不受影响。
  带交叉链接（如 crambin 的 CPM）的靶点现在会自动走 apo 通路（修饰残基
  归蛋白、过滤后跑纯蛋白）。

**Q2：报 "library dtp missing/corrupt (0 bytes)" 要紧吗？**
不要紧——`dtp` 是旧默认，指向的镜像已死（dtpbase.org，502）。现在默认
是本地主库 `nci_npatlas`；库缺失/损坏时按 `nci_npatlas` →
`chembl35_small`（5 万）自动回退，state/报告标注 `fallback for …`。
想补 PDBBind：`drugagent setup --libraries pdbbind`（软失败不打断 setup）。

**Q3：跑一半断了/机器重启了？**
`resume --project <目录>`。产物幂等（已完成的 mdrun/对接/PDBQT 自动复用），
不会从头烧。

**Q4：只想重算 MD（比如改了参数）？**
`rerun --project <目录> --stage md`（会按新参数重跑构建+MD+分析，其他阶段
产物不动）。

**Q5：对接分数 -12 和 -15，差多少？**
同一库内比较：差 3 kcal/mol 约等于 298 K 下结合常数差 ~7 倍（ΔG=RTlnK）。
跨库/跨参数不要硬比。

**Q6：binder pLDDT 才 40，设计失败了吗？**
不一定。de novo 设计 pLDDT 偏低属常态（ESMFold 从序列重折叠的置信度，
不等于复现 RF 设计构象的置信度）。重点看界面 pLDDT + 几何接触
（"目标距离" 列 + 接触对数）；报告第 3 节会综合排名。

**Q7：报告里 RMSD 5 Å 很大，蛋白散架了？**
先分三步判（§6.2 表）：链自拟合 <6 Å → 没散架；链间相对 >> 自拟合 →
是多聚体相对运动；域 norm >0.08 → 有真实铰链运动。`--fast` 短时长下
1-2 Å 漂移完全正常（刚性蛋白 2 ns 基线上限 ~2.1 Å）。

**Q8：没有 GPU 行不行？**
行。无 `/dev/nvidia*` 时全 CPU（RF/ESMFold/vina/GROMACS 都有 CPU 路径）。
64 核 CPU 全模块 fast e2e 约 3-5 h。

**Q9：只有序列没有结构？**
`--target` 直接给 FASTA 或裸序列（自动 ESMFold 建模；单序列用默认条件，
打分口径按"重折叠置信度"读）。

**Q10：用自己的小分子库？**
`--library /path/to/my.sdf`（SDF 格式；建议先 RDKit 标准化，管线也会做
二次标准化）。

**Q11：磁盘不够了清哪里？**
`projects/<旧项目>/05_md/md_rep*/`（xtc 轨迹，按 ns 增长）和
`projects/*/02_screen/` 的 PDBQT 缓存。`env/`、`data/` 别动（~40 GB）。
`logs/` 可清。

**Q12：外部命令失败的报错在哪看？**
异常信息自带日志尾部（最近 40 行），`status` 也会列最近 3 条工具失败。
原始日志在项目目录内（如 `05_md/md_rep1/md.log`）。

---

## 11. 数字可信度指南（什么能信、什么不能信）

| 数字 | 能信什么 | 不能信什么 |
|---|---|---|
| Vina 对接分 | 库内相对排名（富集） | 绝对亲和力；跨体系比较 |
| GNINA 复打分 | 同上（更物理，CASP 验证） | 同上 |
| ESMFold pLDDT（VHH 文库） | 文库内可折叠性排序 | 对 CDR3 loop 系统性偏低（30-35 集中是正常现象，不是 bug） |
| binder pLDDT / 界面分 | 相对排序 + 界面几何（接触对/距离） | 绝对"能结合"的断言（de novo 需实验验证） |
| scaffold_rmsd_a | scaffold-guided 保真度量化 | 小值才好的直觉（官方允许精细结构变化，15 Å 是预期量级） |
| MD RMSD（fast 2-5 ns） | 稳定性验证、柔性区定位、铰链信号 | 自由能/收敛结论（需 full 时长 + 收敛判定通过） |
| MD 柔性解读 | 刚性/铰链/去折叠的定性区分（有跨靶点刚性基线标定） | 定量构象动力学（时间太短） |

**总原则**：DrugAgent 的数字都是"排序 + 定性判断"级；要定量结论
（Kd、ΔG、速率常数）需要对应的实验/高级计算（FEP、umbrella sampling）。

---

## 12. 术语表

- **PDB ID**：蛋白质数据库编号（如 1HVI = HIV 蛋白酶），RCSB 可下载。
- **apo / 复合物（complex）**：不带配体的蛋白 / 带配体的蛋白。
- **口袋（pocket）**：蛋白表面适合结合小分子的凹槽，自动检测。
- **PDBQT**：Vina 对接专用格式（带原子类型/电荷/分子图）。
- **对接分数（kcal/mol）**：经验结合自由能估算，越低越好，相对量。
- **虚拟筛选（screening）**：把库里的分子逐个"试装"进口袋打分排序。
- **RFdiffusion**：扩散模型设计蛋白骨架（Rosalind 实验室）。
- **ProteinMPNN**：给定骨架生成序列（逆折叠）。
- **ESMFold**：序列 → 3D 结构 + pLDDT 置信度。
- **pLDDT**：局部结构预测置信度（0-100），高=可信。
- **VHH / 纳米抗体**：单域抗体（~120 残基），CDR 是互补决定区（结合位点）。
- **scaffold-guided**：用已知骨架的二级结构模式 + 块邻接做设计条件。
- **MD（分子动力学）**：牛顿力学积分原子运动，采样构象。
- **力场**：原子间作用的参数化函数（amber99sb-ildn 蛋白 / GAFF2 配体）。
- **RMSD**：结构均方根偏差，构象差异指标（nm 或 Å，1 nm = 10 Å）。
- **RMSF**：每残基的 RMSD（时间平均），柔性指标。
- **Kabsch**：两结构间的最优刚体对齐算法。
- **PBC（周期性边界）**：模拟盒子的"无缝"处理；轨迹分析需先 unwrap。
- **收敛判定**：RMSD 平台 + 主导簇占比，未收敛自动延长（R6）。
- **构象选择**：用 MD 采样到的多个构象逐一筛选（consensus 平均分）。
- **铰链（hinge）**：蛋白结构域间的柔性关节，变构/功能运动的关键。

---

## 13. 参考运行（可直接打开对照）

| 项目 | 内容 | 看点 |
|---|---|---|
| `projects/r10_e2e/` | HIV-PR 全模块（screen+binder+vhh+md，5 ns×3） | 完整报告 6 节；5 个结构域的 vs-rest 分析（0.8-2.4 Å 真实铰链信号）；VHH 设计 scaffold_rmsd_a=15.08 Å；vhh_30 frag1 −8.66 kcal/mol |
| `projects/rigid_ubq/` | 泛素 1UBQ apo 2 ns×3 | **刚性基线**（core RMSD 末端 0.11-0.21 nm，域 norm 0.04-0.08） |
| `projects/smoke_screen_demo/` | 1HVI 仅 screen（fast） | 最小对接运行范例 |

```bash
# 打开任一项目的报告
xdg-open projects/r10_e2e/reports/report.html
```

---

*教程基于第 17 轮代码（commit `4787a88`）。代码演进以 HANDOFF.md /
ROUNDLOG.md 为准；发现教程与代码不一致时以代码为准并提 issue 到
ROUNDLOG 反思段。（英文版：[TUTORIAL.md](TUTORIAL.md)）*
