# DrugAgent 2.0 — 架构设计

## 核心思想

**流水线 = 工具箱，Agent（LLM）= 主程序。**

旧版（1.0）：固定 LangGraph 状态机（target → screen → binder → vhh → md → report），LLM 只在
少数节点做"判断"（选口袋、定阈值、选体系）。流程顺序、参数、排错都是硬编码的。

新版（2.0）：LLM 通过 **function calling（ReAct 主循环）** 直接驱动整个流程：

- 每个原子操作（找口袋、对接、写 MDP、跑 mdrun、分析轨迹…）都是一个**工具**；
- LLM 自由决定调用顺序、参数、重试策略；
- 失败时 LLM 自己读日志、诊断、打补丁、重试（带预算，耗尽才回到人）；
- 力场、盒子、盐浓度、barostat、积分步长等参数由 LLM 自己写进 MDP（工具只校验/执行）；
- 人在 4 个固定里程碑 + LLM 主动发起的动态提问处介入（`--auto` 时自动通过）。

1.0 的确定性阶段函数**整体保留为工具**（`run_target_prep` / `run_screening` / …），
Agent 在需要"一段可靠的标准流程"时直接调用；需要定制时改用细粒度工具。

## 目录与模块

```
drugagent/
  agent/                 # 2.0 核心（新增）
    __init__.py          #   build_tools(), AgentLoop, system_prompt()
    loop.py              #   AgentLoop: ReAct 循环 + transcript + 预算 + 检查点
    prompts.py           #   中文系统提示词（目标模板 + 工作规约）
    tools_core.py        #   元工具: 文件/shell/决策记录/人工确认/finish
    tools_target.py      #   靶点准备工具
    tools_screen.py      #   小分子筛选工具
    tools_design.py      #   binder 设计工具
    tools_vhh.py         #   VHH 工具
    tools_md.py          #   GROMACS/MD 工具
    tools_report.py      #   报告工具
  modules/               # 1.0 阶段模块（工具的后端，基本不动）
    target_prep.py  screening.py  binder.py  vhh.py  mdmsim.py  esmfold_run.py
  report/report.py       # 报告生成（消费 state dict，接口不变）
  llm.py                 # AgentBrain（+新增 chat_tools: 原生 function calling）
  cli.py                 # run/resume/status/report（agent 语义）
```

## ReAct 主循环（agent/loop.py）

```
messages = [system(目标+规约), user(任务)]
loop:
  resp = LLM.chat(messages, tools=注册表)
  若 resp 无 tool_calls:
      补一条提示"请调用工具或 finish"（连续 2 次 → 视为 needs_human）
  对每个 tool_call:
      result = 执行工具(参数)      # 结果截断到 ~6KB，全量落盘 agent/out/
      messages.append(role=tool, result)
  transcript.jsonl 追加（每步落盘，支持 resume）
  步数 +1；超过 max_steps → 提示 LLM 收尾
终止条件:
  finish(status) 工具调用 / ask_human 阻塞（交互模式）/ 预算耗尽
```

- **transcript**：`agent/transcript.jsonl`，每行一个消息（assistant 的 tool_calls、tool 结果、
  id 全部保留）。`resume` = 重放 transcript 重建 messages，追加 "继续" 指令。
- **状态**：项目目录是唯一事实源。`state.json`（根）保持 1.0 的 state 形状
  （`target/options/target_prep/screening/binder/vhh/md/report/decisions/errors/status`），
  因此 `build_report(state)`、`status`、旧报告模板**零改动复用**。
- **决策留痕**：`record_decision(stage, choice, rationale)` 工具 → `decisions.json`
  （与 1.0 同格式，报告里照旧展示）。

## 工具清单（细粒度）

| 组 | 工具 |
|---|---|
| 元 | `list_dir` `read_file` `write_file` `edit_file` `shell` `record_decision` `ask_human` `checkpoint` `finish` `state_get` |
| 靶点 | `resolve_target` `analyze_pdb` `clean_pdb` `find_pockets` `pdb_to_pdbqt` `run_target_prep`(兜底) |
| 筛选 | `standardize_library` `prefilter_ligands` `dock_screen` `redock_ligand` `decide_hits` `run_screening`(兜底) |
| 设计 | `rf_design` `mpnn_sequence` `esmfold_score` `make_complex` `run_design`(兜底) |
| VHH | `vhh_screen` `vhh_design` `run_vhh`(兜底) |
| MD | `gmx_env` `md_prepare` `mdp_template` `grompp` `mdrun` `gmx_analyze` `md_summary` `run_md`(兜底) |
| 报告 | `build_report` |

工具约定：
- 所有工具返回 JSON 可序列化 dict；`ok: bool` + `error` 或结果字段；
- 长输出（日志、打分表）只回传摘要 + 文件路径，Agent 用 `read_file` 按需查看；
- 幂等：产物已存在且完整时直接复用（mdrun 检测 md.xtc + "Finished mdrun"）。

## 检查点（固定 + 动态）

- **固定 4 个**：`checkpoint(stage)`，stage ∈ target / screening / design / md（MD 跑前）。
  交互模式阻塞等待人（批准/修改/中止）；`--auto` 自动选第一个选项并记录。
- **动态**：`ask_human(question, options)` 任意时刻可调用；`--auto` 返回"自行决策"提示。

## CLI（agent 语义）

```
drugagent run --target 1HVI [--fast] [--auto] [--modules ...] [--md-ns 100] \
              [--md-reps 3] [--max-steps 300] [--library ...] [--no-llm]
drugagent resume --project DIR          # 重放 transcript 继续
drugagent status  [--project DIR]
drugagent report  --project DIR
drugagent setup   # 不变
```

`run` 默认 = agent 模式（取代 1.0 pipeline 入口）。`--no-llm` 时 Agent 退化为
"按固定顺序调用兜底阶段工具"的脚本模式（仍走同一循环，系统提示词换成确定性剧本）。

## 预算

- `max_steps`（LLM 轮数，默认 300）；每步可含多次工具调用；
- 工具超时由参数指定（`shell` 默认 600s，`mdrun` 默认 24h）；
- 预算耗尽 → 注入收尾提示 → 仍不收尾则 `finish(needs_human)`。
