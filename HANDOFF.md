# DrugAgent 2.0 — 交接说明（第 13 轮结束）

## 一句话状态
代码库处于"第 13 轮迭代完成"状态：**track B 设计验证缓存**（scored.json
+ 顶层设计 glob 修复，rerun vhh 从 ~3 h 降到 7.5 min）、**R1-v2 域-其余
蛋白相对 RMSD**（铰链/变构检测 + 解读注 9 + 报告列/曲线）、**pLDDT 分布
直方图**（G9 门槛可视化）。快速测试 197/197 绿。r10_e2e 全状态最新：
vhh_30 综合分 1.0 居首（frag1 -8.66 kcal/mol）、MD 域表含相对其余列 +
注 9（dom3 构象重排 14.1 Å）。

## 关键路径
- 项目根：`/home/data/lrs/drug/drugagent`（唯一可写持久盘；**/tmp 是 tmpfs
  且有后台清理进程，pytest 慢测务必 `--basetemp=$PWD/data/fixtures/_ptmp`**）
- 环境：`env/bin/python`（py3.12 conda）；GROMACS 2023.1 在 `data/tools/gromacs`
- 测试：`env/bin/python -m pytest tests/ -m "not slow" -q --basetemp=$PWD/data/fixtures/_ptmp`
- e2e：`env/bin/python -m drugagent.cli run --target 1HVI --modules all --fast --auto --no-llm --name r10_e2e`
  （产物幂等复用：重跑同一 `--name` 会跳过已完成阶段产物）
- 轮次记录：`ROUNDLOG.md`（每轮计划/结果/反思，下轮计划从"反思/下轮缺口"开始）
- git 仓库已初始化（.gitignore 覆盖 env/data/projects/日志）

## 第 13 轮已完成（见 ROUNDLOG 详情）
1. **track B 设计验证缓存**：`vhh.score_designs()`（抽出）+
   `vhh_designs/scored.json` 按设计 PDB mtime 缓存（只缓存成功、错误
   重试、增量写盘）。**附带 bug 修复**：`vhh_design_*.pdb` glob 把打分
   副产物（`_binder`/`_complex`/`_binder_complex`…）当设计验证且自我繁殖
   （12 文件当 6 设计）——现只匹配 `vhh_design_N.pdb` 顶层设计。
   **实测 rerun vhh：~3 h → 7.5 min**（track A dock 缓存复用 + 设计
   验证缓存 + 1 个被 RF 覆盖的设计重验证）。
2. **R1-v2 域-其余蛋白相对 RMSD**：`mdsim._kabsch_transform()`（行向量
   约定 mobile @ M + t，M = U D Vt；修掉 R12 重构引入的转置 bug——
   准线性点云上转置 SVD 解非最优，brute-force 验证）；
   `domain_vs_rest_rmsd_series()`（每帧 fit 其余蛋白到 frame 0，同变换
   测该域 → 铰链/变构信号）。summary 新增 `domain_rmsd_vs_rest`
   （final/mean/series，多副本均值）；解读注 9 两分支；报告域表
   "末端相对其余 (Å)" 列 + 点线曲线。**实测 r10_e2e（3 副本）**：
   vs-rest 末端 8.4-14.1 Å（自拟合 4.5-5.9 Å）——HIV-PR apo 5 ns 半刚性
   域重排，dom3 (82-98) 最大 → 注 9 柔性分支。
3. **pLDDT 直方图**：screen_vhh 返回 `plddt_all`（全建模库，n=80）；
   报告 VHH 节直方图 + `vhh_plddt_min` 门槛虚线（随 fast/full 自适应）。

## 第 12 轮已完成（见 ROUNDLOG 详情）
1. **G10-v3 CDR 片段打磨**：(a) >20 残基 low-pLDDT run 自动切 ≤20 残基块
   （dock 成本 ~O(n^1.9)）；(b) **片段自适应盒子**（`vhh._frag_box`：中心=
   pocket 中心、边长=片段直径+8 Å、clip 12-30 Å）——实测 vhh_30 frag0
   6.25 min→1.3 min、2.3e7→145 kcal/mol（撞墙罚分消失）；(c) 无片段时
   fast 默认跳过（`vhh_dock_full_fallback=False`，score=NaN、pLDDT 仍参与
   排名；full 模式默认 True）；(d) `fragment_scores`/`n_fragments` 进 merge
   → 报告 VHH 表 "CDR片段" 列。
2. **R1 真结构域 RMSD**：`mdsim.find_ss_domains`（frame 0 SS 分类，≥8 连续
   结构化残基）+ `domain_rmsd_series`（每域 CA Kabsch **自拟合**到自身
   frame 0 → 逐帧域内部形变 RMSD）。`_ss_backbone_trajectory` 现在逐原子
   PBC unwrap + frame 0 两级 make-whole（残基内 min-image 一致 + 链走）——
   否则跨盒残基产生 50 Å+ 幻影距离。修了 Kabsch 约定 bug（H=AᵀB → R=U D
   Vᵀ）。summary 新增 `domains`/`domain_rmsd`（final/mean/series）；解读注 8
   （域构象漂移 >4 Å / 整体高 RMSD 但各域稳定→域间运动）；报告域表格+曲线。
   r10_e2e 实测 5 域、末端自拟合 RMSD 3.6-6.7 Å。
3. **`rerun --stage X`**：`drugagent rerun --project X --stage
   {target_prep,screening,binder,vhh,md,report} [--no-with-report]`——G8
   force 的 CLI 化；成功后默认重建报告。
4. **bench 收尾**：修 `scripts/bench_vhh_dock.py` 分数解析（VINA RESULT
   行取 score 字段而非行尾 rmsd_ub）；bench.json 最终 rigid 100.96 /
   flex 109.73 min（1.10×）、同条件分数差 <0.001 kcal/mol。
5. **e2e 重验**（`rerun --stage vhh`，一次验证 rerun CLI + G10 v3）：
   9 片段全 dock 成功（30 Å 自适应盒）；vhh_30 frag1 = **-8.66 kcal/mol**、
   vhh_30 综合分 1.0 居首；vhh_60/22/78 最佳片段 126-187。报告 CDR 片段
   列 ✓。
6. **DTP 重试**：dtpbase.org 仍 SSL EOF——fallback chembl35_small 兜底
   （软失败不阻断）。

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
3. **G10 VHH 刚性对接 + G10-v2 CDR 片段对接**：默认 `flex=False`
   （`vhh_dock_flex`/CLI `--vhh-dock-flex` 开柔性）。坑：此 vina 构建
   **所有** ligand 都要 ROOT/ENDROOT+TORSDOF 分子图（小分子也是），无图
   报 "Unknown or inappropriate tag" → 刚性 = 图保留 + `TORSDOF 0`
   （`make_rigid_pdbqt()`）；PDBQT 缓存按扭转数判 flex/rigid，模式不匹配
   自动重转。**基准**（`scripts/bench_vhh_dock.py`）：成本 ~O(n^1.9) 于
   原子数、扭转数无关（200 原子 TORSDOF 0 vs 20 = 3.4 vs 3.6 min）；全长
   773 原子 flex 76 min / rigid 101 min 且最优分几乎相同（2.1e8，撞墙
   主导）→ **G10-v2**：fast 默认 `vhh_dock_cdr_only`——`vhh_cdr_fragments()`
   按残基 pLDDT<50 连续区切 1-3 片段（pad 2，≥4 残基，整结构低则回退
   全长），分别对接，综合分=最佳片段分（`fragment_scores` 保留；实测
   257 原子片段 6.25 min vs 全长 ~80-100 min，~15×）。
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
