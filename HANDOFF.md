# DrugAgent 2.0 — 交接说明（第 16 轮结束）

## 一句话状态
代码库处于"第 16 轮迭代完成"状态：**纯蛋白 apo MD 通路**（target_prep
无配体崩溃修复 + select_system/build 的 apo 选项，ubq 1UBQ 2ns×3 首跑
成功）、**scaffold 保真度守卫**（scaffold_rmsd_a 字段；真实设计对 1EWN
骨架 15.08 Å 漂移——RF 官方文档确认 scaffold-guided 只条件 SS+邻接、
允许序列/精细结构变化）、**单位哨兵测试**（MDA vs gmx self-fit 比值带
[0.05,1.0]；过程中剥离三处 PBC/参考伪影）、**binder_rf_cautious**
（与 vhh 对称）、**跨靶点刚性基线**（ubq：self-fit ≤0.21 nm 干净分开
刚性/柔性；域级 final_norm 0.04-0.08 与 HIV-PR 0.041-0.109 重叠 →
域级 norm 需先扣刚性基线）。快速测试 212/212 绿。已知未修：gmx 链级
RMSD 用过期 tpr 作参考（rep1 "2.083" 应读 ~0.3）。

## 关键路径
- 项目根：`/home/data/lrs/drug/drugagent`（唯一可写持久盘；**/tmp 是 tmpfs
  且有后台清理进程，pytest 慢测务必 `--basetemp=$PWD/data/fixtures/_ptmp`**）
- 环境：`env/bin/python`（py3.12 conda）；GROMACS 2023.1 在 `data/tools/gromacs`
- 测试：`env/bin/python -m pytest tests/ -m "not slow" -q --basetemp=$PWD/data/fixtures/_ptmp`
- e2e：`env/bin/python -m drugagent.cli run --target 1HVI --modules all --fast --auto --no-llm --name r10_e2e`
  （产物幂等复用：重跑同一 `--name` 会跳过已完成阶段产物）
- 轮次记录：`ROUNDLOG.md`（每轮计划/结果/反思，下轮计划从"反思/下轮缺口"开始）
- git 仓库已初始化（.gitignore 覆盖 env/data/projects/日志）

## 第 16 轮已完成（见 ROUNDLOG 详情）
1. **P0a 纯蛋白 target_prep**：无配体靶点 `lig_pdb` 未初始化 →
   `UnboundLocalError`。修：`lig_pdb=None` 初始化；纯蛋白 e2e 走全程
   （ligand_* 为 None、pocket 表面回退）。
2. **P0b apo MD 通路**：`select_system` 无配体时追加 `apo` 选项；
   `build_complex_pdb` 加 apo 分支（仅靶点去水、保原生链号）。
3. **P1 scaffold 保真度**：`vhh.scaffold_fidelity()`（设计链 vs
   scaffold CA Kabsch，Å）+ `scaffold_rmsd_a` 字段。**发现**：真实
   r10_e2e 设计对 1EWN 骨架漂移 15.08 Å（noise_scale_ca=0 未刚性钉住）
   但局部保留（40-mer 滑窗 5.9 Å、200=200、SS 70% 链一致）；1EWN 实为
   人 AAG DNA 修复酶核心（非 VHH）；RF 官方文档：scaffold-guided 用
   SS+邻接条件、"允许序列与精细结构变化" → 漂移是预期，scaffold 序列
   作打分先验仍合理（近似）。
4. **P2 单位哨兵**：`tests/test_mdsim.py::test_md_unit_sentinel_mda_vs_gmx`。
   剥伪影三层：(a) rmsd_r1.xvg 参考是过期 tpr（与 xtc f0 差 ~21-23 Å）
   → 2.083 偏高（R17 修）；(b) 紧凑盒子 wrapped 边界 flapping（gmx
   0→1.36 nm vs MDA 0.13 nm）；(c) trjconv compact 成像平移 frame 0
   （flat 2.54 nm 伪影）。最终：trjconv compact 展开 + 自身 f0 参考 +
   gmx rms vs MDA CA self-fit，三副本基线比 0.13-0.40，带 [0.05,1.0]
   抓 10× 双向。
5. **P3 binder_rf_cautious**：`rfdesign(rf_cautious=True)` 在顶层
   design_*.pdb ≥ n_designs 时跳过 RF（默认 False 行为不变）；
   `run_design` 从 options 读 `binder_rf_cautious`。
6. **P4 跨靶点刚性基线**（ubq 1UBQ apo 2ns×3，~18 min）：整体 self-fit
   末端 ubq 0.110-0.208 vs HIV-PR 0.289-2.083（无重叠，干净分开）；
   域 vs-rest 0.076-0.183 vs 0.075-0.160（重叠）；final_norm 0.040-0.080
   vs 0.041-0.109（重叠）。→ 域级 norm 需扣刚性基线；self-fit 末端可
   作铰链判据；4 Å 绝对阈值远高于刚性基线（无误报）。GFP(1GFL) 因 CRO
   色原体（SER 65 缺 O）pdb2gmx 失败 → 修饰残基 R17 缺口。
7. **e2e**：212/212 快测绿（+7）。r10_e2e 状态未动。

## 第 15 轮已完成（见 ROUNDLOG 详情）
1. **P0 单位 bug（最重要）**：`mdsim._ss_backbone_trajectory` 返回
   MDA 原始 Å，下游按 nm → 域 RMSD ×10、`HBOND_DIST=4.5` 按 nm
   比较（ss_frac 0.99、域=整链）。修：出口 ÷10；`HBOND_DIST=0.45 nm`；
   合成 SS 夹具 ×0.1；回归测试 `test_ss_trajectory_units_are_nm`
   （真实轨迹：CA 相邻 0.25-0.50 nm、尺度 <8 nm）。**注意**：gmx 系
   指标（xvg）一直是 nm 正确，只有 MDA 域链路受影响。
2. **P1 scaffold 序列**：`vhh.design_vhh` 读 scaffold PDB 序列 →
   `score_designs(seqs=…)`；全 GLY 且等长时替换 seqB（`seq_used=
   "scaffold"`），缓存签名含 alt。假设：设计 200 残基=scaffold 200
   残基（长度匹配，未读 RF log 确认全序列保留）。
3. **P2 vhh_rf_cautious**：选项默认 true（旧行为）；false 时顶层
   vhh_design_*.pdb → `archive_<ts>/` 再跑 RF。
4. **P3 标定**：带配体参照（agent_smoke，2 ns×2）vs apo（5 ns×3）：
   final_norm 0.068-0.098 vs 0.041-0.109 重叠 → norm 不足以区分
   状态；绝对阈值 4 Å 在真值下无误报。跨靶点标定留 R16+。
5. **P4 binder 缓存 e2e**：`rerun --stage binder` 生成 03_binder/
   scored.json（签名含 MPNN 真实序列 + 权重标签）；新设计
   iface 52.5/39.6（真序列 → 分数分化，对照全 GLY 时代）。

## 第 14 轮已完成（见 ROUNDLOG 详情）
1. **P2 track B pLDDT 相同根因 + 几何接口**：两设计几何实质不同
   （目标相对帧位姿差 2.2 nm，min CA 距离 3.4 Å、接触对 30/31），
   相同 pLDDT 的根因 = RF 设计链全 GLY → ESMFold 只收序列 → 位点盲。
   新增 `vhh.design_interface_geom()`（min_dist_a / n_contacts /
   contact<8 Å）进设计条目 + ranked 候选（min_dist_a/contact）+
   报告候选表 "目标距离 (Å)" 列（⚠未接触 标记）。注意 PDB 坐标
   是 Å（勿 ×10——本轮手测脚本就踩过）。
2. **P1 binder 缓存**：`binder.score_designs` 同款 scored.json
   （03_binder/scored.json），签名 = [设计 mtime, 靶点 mtime,
   alt 序列, esmfold_version_tag]。
3. **P3 版本 + 直方图**：`esmfold_run.esmfold_version_tag()`
   （ckpt 名+大小）进 vhh/binder 缓存签名（权重更新 → 缓存整体
   失效，实测生效）；直方图标题 p50/p90 + 过门槛条数（plotly JSON
   把 `/` 转义 \u002f，标题避开斜杠）。
4. **P4 vs-rest 归一**：`mdsim.domain_diameters()`（frame 0 CA 直径）
   → 域条目 `diameter_nm`；summary vs-rest 附 `final_norm`/
   `mean_norm`（÷直径）；报告域表 "归一化 (÷直径)" 列。实测
   r10_e2e：直径 11.4-24.4 nm，final_norm 0.041-0.109。

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

## 下轮缺口优先级（R17，ROUNDLOG 第 16 轮"反思"有完整版）
1. **过期 TPR 参考（最大未修 bug）**：gmx 链级 RMSD（rmsd_rN.xvg +
   per-chain drift）以 md.tpr 为参考，tpr 是 2 ns 段产物、xtc 是延长
   段轨迹（起点漂移 ~21 Å）→ rep1 "2.083" 应读 ~0.3。修复：分析前从
   xtc 抽 frame 0 作参考，或延长前重新 grompp。
2. **铰链判据升级**：整体 self-fit 末端（第 16 轮证明可分刚性/柔性）
   并入 interpret 注 8/9，"整体+域级"两级判据；刚性基线（ubq）写进注释。
3. **修饰残基 + 第二刚性参照**：pdb2gmx 对 CRO 类修饰残基处理（或换
   无修饰残基刚性蛋白），补第二刚性点 + 已知铰链体系（无辅因子腺苷酸
   激酶），标定扩成小矩阵。
4. **scaffold 内容**：1EWN 是 AAG 非 VHH——换真 VHH scaffold 或改名
   目录并文档化"fold 条件而非序列条件"。
5. **展开策略一致性**：compact-com vs 逐原子 unwrap 在紧凑盒子下不等价
   （gmx rep2 0.984 vs MDA 0.130 nm）；修过期 TPR 时评估统一 unwrap。

## 验证命令
```
env/bin/python -m pytest tests/ -m "not slow" -q --basetemp=$PWD/data/fixtures/_ptmp
env/bin/python -m drugagent.cli status --project projects/r10_e2e
ls projects/r10_e2e/reports/   # report.html / report.pdf
```
