# DrugAgent 2.0 — 交接说明（第 11 轮结束）

## 一句话状态
代码库处于"第 11 轮迭代完成"状态：G8 阶段级幂等（run_* 完成即跳过 +
force）、G9 VHH 命中面加宽（pLDDT fast 35 / full 50，fast n=80）、G10 VHH
默认刚性对接（TORSDOF 0）、R5 柔性工作流自动判据（apo issue + apo×高 RMSF
建议）、status 命令增强、pdbbind 根因（镜像 404 注册制，软失败）。快速测试
175/175 绿。r10_e2e（1HVI 全模块 --fast --auto --no-llm）跑完收尾验证。

## 关键路径
- 项目根：`/home/data/lrs/drug/drugagent`（唯一可写持久盘；**/tmp 是 tmpfs
  且有后台清理进程，pytest 慢测务必 `--basetemp=$PWD/data/fixtures/_ptmp`**）
- 环境：`env/bin/python`（py3.12 conda）；GROMACS 2023.1 在 `data/tools/gromacs`
- 测试：`env/bin/python -m pytest tests/ -m "not slow" -q --basetemp=$PWD/data/fixtures/_ptmp`
- e2e：`env/bin/python -m drugagent.cli run --target 1HVI --modules all --fast --auto --no-llm --name r10_e2e`
  （产物幂等复用：重跑同一 `--name` 会跳过已完成阶段产物）
- 轮次记录：`ROUNDLOG.md`（每轮计划/结果/反思，下轮计划从"反思/下轮缺口"开始）
- git 仓库已初始化（.gitignore 覆盖 env/data/projects/日志）

## 第 11 轮已完成（见 ROUNDLOG 详情）
1. **G8 阶段级幂等**：`drugagent/agent/stages.py` — 6 阶段完成判据
   （target_prep: receptor_pdbqt+clean_pdb+pocket；screening: hit_decision；
   binder: designs+best；vhh: track_a+track_b；md: summary+replicas；
   report: report.html 存在）+ 缓存摘要。`run_target_prep/run_screening/
   run_design/run_vhh/run_md/build_report` 全加 `force` 参数：state 段完整
   → 整体跳过并返回缓存摘要（`reused: true`）。低层产物幂等不动。
2. **G9 VHH 命中面**：100 条已建模 VHH 的 pLDDT 分布 p50=31.5/p90=34.5/
   >45 仅 1 条（CDR3 loop 无序拉低整体 pLDDT）。`vhh_plddt_min` 进 Defaults
   （fast 35 / full 50，options/CLI `--vhh-plddt-min` 可覆盖）；fast
   `vhh_screen_n` 40→80；修 screen_vhh 写死 "pLDDT>70" 的日志。
3. **G10 VHH 刚性对接**：默认 `flex=False`（`vhh_dock_flex`/CLI
   `--vhh-dock-flex` 开柔性）。坑：此 vina 构建**所有** ligand 都要
   ROOT/ENDROOT 分子图（小分子也是），无图报 "Unknown or inappropriate
   tag" → 刚性 = 图保留 + `TORSDOF 0`（`make_rigid_pdbqt()`）；PDBQT 缓存
   按扭转数判 flex/rigid，模式不匹配自动重转。基准脚本
   `scripts/bench_vhh_dock.py`（data/fixtures/bench_vhh_dock/bench.json）。
4. **R5 自动判据**：`analyze_completeness` 对无配体/无金属结构加 info 级
   `apo_target` issue；`interpret_stability` 判据：apo（is_ligand=False，
   md_summary 从 state 注入）+ 平均 RMSF > 2.5 Å → 建议"短 MD 系综 + 柔性
   靶点工作流"。
5. **status 增强**：每阶段 ✓/✗ + 关键数字（对接数/命中/库名含 fallback 标注/
   设计数/MD ns+final RMSD）+ 磁盘阶段 JSON 产物清单 + 最近 3 条工具失败
   （state.errors + transcript ok=false）+ transcript 尾行。
6. **pdbbind**：根因 = pdbbind.org.cn 两个 tar.gz URL 已 404（138 字节文件
   就是 nginx 404 页，已删；下载改注册制 download.php）。`_setup_pdbbind`
   校验大小+tar 可读、丢 <1MB 坏件、软失败不中断 setup（resolve_library
   回退 chembl35_small 兜底）。

## 环境坑（必读）
- **/tmp tmpfs**：前台命令 300s 超时 → 长任务 `setsid nohup ... > logs/x.log 2>&1 < /dev/null &`；
  pytest 用 `--basetemp` 指本地盘
- 持久 shell 会重置（300s 超时/OOM 触发），nohup 子进程一般幸存，但
  **detached 作业要 setsid 双保险**
- 本机 64 核 / 251GB 内存；LLM 本地 llama.cpp `http://127.0.0.1:18080/v1`
  model `qwen3.8-27b-uncensored`（支持 function calling，需先确认在跑：
  `curl -s -m 5 http://127.0.0.1:18080/v1/models`）
- dtpbase.org 镜像不稳（dtp.sdf.gz 0 字节，G6 回退覆盖）；pdbbind.org.cn
  镜像 404（第 11 轮软失败处理）
- **vina（data/tools/vina，AD4 血缘 f458505-mod）**：ligand 必须带
  ROOT/ENDROOT+TORSDOF 分子图（`to_pdbqt`/`write_ligand_pdbqt` 已自动补齐）；
  大配体 `--cpu N` 实测仍单核（R10 证据）；exh=1 下全长 VHH 柔性 ~76min、
  刚性见 bench.json。vina 静默运行（无进度输出，完成后才写结果）
- 编辑世界与 bash 世界可能不同 mount namespace（若"文件已存在/不存在"矛盾
  时以 bash 为准，用 python 脚本改文件最稳）

## 下轮缺口优先级（ROUNDLOG 第 11 轮"反思"有完整版）
1. **R1 真·结构域 RMSD**：DSSP 域边界/NMDYN 式动态域，`gmx_analyze
   kind="domain"`（分链+柔性区已覆盖大部分需求，这是收尾项）。
2. G10 若刚性基准提速不显著：GNINA 复打分 / 截断框架只柔 CDR。
3. DTP 下载重试（dtpbase.org 偶活）；pdbbind 注册后手动补库。
4. 易用性：stage 级 `drugagent rerun --stage X` 子命令（force 的 CLI 化）。

## 验证命令
```
env/bin/python -m pytest tests/ -m "not slow" -q --basetemp=$PWD/data/fixtures/_ptmp
env/bin/python -m drugagent.cli status --project projects/r10_e2e
ls projects/r10_e2e/reports/   # report.html / report.pdf
```
