# DrugAgent 2.0 — 交接说明（第 10 轮结束）

## 一句话状态
代码库处于"第 10 轮迭代完成"状态：用户优先级缺口 R1/R2/R3/R4/R5 的代码实现
全部到位（R1 区域级柔性、R2 策略引导为本轮补齐），快速测试 147/147 绿，
全模块 e2e（`projects/r10_e2e`，1HVI 全模块 --fast --auto --no-llm）收尾时
运行到 VHH 对接阶段（target/screening/binder 已 ok）。

## 关键路径
- 项目根：`/home/data/lrs/drug/drugagent`（唯一可写持久盘；**/tmp 是 tmpfs
  且有后台清理进程，pytest 慢测务必 `--basetemp=$PWD/data/fixtures/_ptmp`**）
- 环境：`env/bin/python`（py3.12 conda）；GROMACS 2023.1 在 `data/tools/gromacs`
- 测试：`env/bin/python -m pytest tests/ -m "not slow" -q --basetemp=$PWD/data/fixtures/_ptmp`
- e2e：`env/bin/python -m drugagent.cli run --target 1HVI --modules all --fast --auto --no-llm --name r10_e2e`
  （产物幂等复用：重跑同一 `--name` 会跳过已完成阶段产物）
- 轮次记录：`ROUNDLOG.md`（每轮计划/结果/反思，下轮计划从"反思/下轮缺口"开始）
- git 仓库已初始化（.gitignore 覆盖 env/data/projects/日志）

## 第 10 轮已完成（G1-G7，见 ROUNDLOG 详情）
1. G1 `run_cmd` 失败时异常带日志尾部（40 行/2KB）
2. G2 test_md_auto_extend 补 @pytest.mark.slow
3. G3 R1 收尾 `flexible_regions()`（RMSF 连续高值段→残基区间，进
   summary/解读/md_summary/报告）+ 报告 RMSF 单位修正（nm→Å）+
   gmx_analyze 支持 DNA_chain
4. G4 系统提示词新增"柔性靶点策略"段（短 MD 系综→make_flex_receptor→
   dock_conformer_set→consensus）
5. G5 CLI 新增 --md-salt/--md-divalent/--md-divalent-m/--md-extend-ns/
   --md-max-extensions/--md-burn-in-ps
6. G6 小分子库回退 `resolve_library()`（dtp 0 字节 → chembl35_small 50k，
   state/报告标明 fallback）
7. G7 VHH 对接并行 `dock_vhh_candidates()`（≤16 路 pmap；fast vhh_screen_n
   100→40）

## 环境坑（必读）
- **/tmp tmpfs**：前台命令 300s 超时 → 长任务 `setsid nohup ... > logs/x.log 2>&1 < /dev/null &`；
  pytest 用 `--basetemp` 指本地盘
- 持久 shell 会重置（300s 超时/OOM 触发），nohup 子进程一般幸存，但
  **detached 作业要 setsid 双保险**
- 本机 64 核 / 251GB 内存；LLM 本地 llama.cpp `http://127.0.0.1:18080/v1`
  model `qwen3.8-27b-uncensored`（支持 function calling，需先确认在跑：
  `curl -s -m 5 http://127.0.0.1:18080/v1/models`）
- dtpbase.org 镜像不稳（dtp.sdf.gz 0 字节、pdbbind.tar.gz 138 字节均损坏，
  G6 回退已覆盖 dtp；pdbbind 未覆盖）
- 编辑世界与 bash 世界可能不同 mount namespace（本轮未复现，若"文件已存在/
  不存在"矛盾时以 bash 为准，用 python 脚本改文件最稳）

## 下轮缺口优先级（ROUNDLOG 第 10 轮"反思"有完整版）
1. **G8 阶段级幂等**：run_screening/run_vhh/run_design 整段重跑（重启 e2e
   重跑了 35min 筛选）。stage json 完整 → 直接复用 + force 参数。
2. **G9 VHH 命中面太窄**：fast 下 40 条合成 VHH 仅 1 条过 pLDDT>45，
   对接样本=1。放宽阈值/加量/分层建模。
3. **G10 单全长 VHH 对接慢**（exh=1 仍 15-20min）：降 exhaustiveness/
   只柔 CDR/GNINA 复打分。
4. R1 真·结构域 RMSD（DSSP 域边界/NMDYN 式动态域，gmx_analyze kind="domain"）
5. R5 柔性工作流自动判据（无配体+高 RMSF → 自动建议短 MD 系综，进
   analyze_pdb issues）
6. status 命令增强（阶段产物+错误摘要）
7. pdbbind 库修复 / DTP 重试

## 收尾时 e2e 状态
`projects/r10_e2e`：01_target ✓（复用）、02_screening ✓（50k 回退库，
~3600 对接 + GNINA 复打分 + 命中判定）、03_binder ✓（RF 2 设计）、
04_vhh 进行中（1 条 VHH 过 pLDDT 过滤，单条对接 ~20min；之后 track B
de novo 2 设计 → MD 构建/平衡/3×5ns/分析/报告）。
验证命令：
```
tail -50 /home/data/lrs/drug/drugagent/logs/r10_e2e.log | strings
env/bin/python -m drugagent.cli status --project projects/r10_e2e
# 成功后:
ls projects/r10_e2e/reports/   # report.html / report.pdf
```
若进程已死：直接重跑同一条 run 命令（幂等复用），或 `drugagent resume --project projects/r10_e2e`。
