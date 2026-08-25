# DrugAgent — An Integrated Drug-Discovery Agent

**The pipeline is a toolbox; the LLM is the main program.** DrugAgent 2.0
drives the entire drug-discovery workflow with a single ReAct main loop
(LLM + native function calling): the LLM plans on its own, invokes ~41
fine-grained tools, writes its own MDP files / picks force fields, and
reads logs to self-diagnose and retry on failure; 4 fixed milestone
checkpoints + dynamic confirmations the LLM can raise at any time, all
auto-passed under `--auto`.

- 🇨🇳 **中文版 (Chinese)**: [README.zh-CN.md](README.zh-CN.md) · [TUTORIAL.zh-CN.md](TUTORIAL.zh-CN.md)
- 📖 **User tutorial (zero background → self-sufficient)**: [TUTORIAL.md](TUTORIAL.md)
- 🏛 Architecture: [DESIGN.md](DESIGN.md) *(in Chinese)*
- 📦 Environment & known pitfalls (authoritative list): [HANDOFF.md](HANDOFF.md) *(in Chinese)*
- 📜 Iteration log (17 rounds, each with plan / results / reflection): [ROUNDLOG.md](ROUNDLOG.md) *(in Chinese)*

## Tech Stack

| Layer | Component | Version / Notes |
|---|---|---|
| Language/env | Python (conda) | 3.12, standalone env in `env/` |
| LLM main program | llama.cpp server (OpenAI-compatible endpoint) | local `127.0.0.1:18080`, default `qwen3.8-27b-uncensored` (needs function calling; override via `DRUGAGENT_LLM_BASE_URL` / `MODEL` / `API_KEY`) |
| Target | RCSB download, OpenBabel, ESMFold (vendored openfold + CPU patches) | auto-modeling for sequence input |
| Small molecules | RDKit, AutoDock Vina, GNINA 1.3.1 (ELF, CPU) | libraries: **nci_npatlas master** (NCI/DTP open set ∪ NPAtlas natural products, ~265k, InChIKey-deduped) / ChEMBL35 / PDBBind (auto-fallback) |
| Protein design | RFdiffusion (hydra), ProteinMPNN (vanilla v_48_010), ESMFold + ESM2 | de novo binders / scaffold-guided VHH |
| MD | GROMACS 2023.1 (self-built, amber99sb-ildn), ACPYPE (GAFF2), MDAnalysis 2.10 | equilibration + replicas + auto-extension + domain-level flexibility diagnostics |
| Reports | Plotly (interactive charts), 3Dmol.js (3D structures), WeasyPrint (PDF) | `reports/report.html` + `.pdf` |
| Tests | pytest | 218 fast + 15 slow e2e (tests are local, see below) |

Hardware: 64 CPU cores is enough (no GPU → fully automatic CPU mode);
disk ~40 GB (env + weights + libraries).

## Modules (what the toolbox covers)

| Module | Contents | Key tools |
|---|---|---|
| A Target prep | PDB file / PDB ID / FASTA / raw sequence input; integrity analysis + structure-pitfall pre-check (multi-MODEL / altloc / missing residues / metals / nucleic acids / disordered termini); agent judgment; auto + manual PDB repair; cleaning; pocket detection; PDBQT conversion | RCSB, obabel, ESMFold (sequence input only) |
| B Small-molecule screening | Large library (nci_npatlas master / ChEMBL / PDBBind or custom SDF) → standardization → physchem/ML prefilter → parallel Vina docking → GNINA rescoring → agent-set hit criteria (co-crystallized ligand redock as positive control); **flexible-target workflow (R2/R5: MD conformational ensemble → multi-conformer selection + side-chain `--flex`, consensus averaging)** | RDKit, Vina, GNINA |
| C Binder design | RFdiffusion de novo design + ProteinMPNN sequences + ESMFold monomer/complex scoring (interface pLDDT) + geometric interface metrics (min distance / contact pairs) | RFdiffusion, ProteinMPNN, ESMFold |
| D Nanobody (VHH) | Track A: VHH library → ESMFold modeling → pLDDT filter (fast 35 / full 50) → **rigid docking** (fast default: **CDR-fragment docking**, ~15× speedup, per-fragment adaptive boxes) → parallel screening; Track B: RFdiffusion scaffold-guided de novo design + scaffold fidelity (`scaffold_rmsd_a`) + composite scoring | ESMFold, RFdiffusion, Vina |
| E MD simulation | System selection (ligand / hit / binder / VHH / **apo**; modified residues auto-prefer apo) → pdb2gmx + ACPYPE → solvation / ions / EM → **NVT→NPT equilibration (position restraints, C-rescale) + burn-in trimming** → N ns × R replicas (forked from equilibrated end-state) → **convergence check + auto-extension** → RMSD/RMSF/Rg/clustering + **flexibility diagnostics** (per-chain self-fit / secondary structure / flexible-region localization / **structural-domain RMSD** / **domain vs-rest + diameter normalization** / **rigid baseline comparison (R17: √t scaling + domain-norm baseline)** / compact-unwrap against tight-box flapping) + metal-ion coordination + native nucleic-acid parametrization + cofactor/heme auto-parametrization | GROMACS 2023.1, ACPYPE, MDAnalysis |

## Quick Start

```bash
cd /home/data/lrs/drug/drugagent

# 0) Install the `drugagent` command (one-time; editable — code edits
#    apply immediately; symlink puts it in your user PATH)
env/bin/pip install -e . --no-deps
ln -sf $PWD/env/bin/drugagent ~/.local/bin/drugagent

# 1) One-shot environment setup (idempotent, re-runnable)
drugagent setup

# 2) Fast end-to-end validation (1HVI target, all modules, small scale, unattended)
drugagent run --target 1HVI --modules all --fast --auto

# 3) Production run (interactive checkpoints: approve / modify / abort)
drugagent run --target 1HVI --modules screen,binder,vhh,md
drugagent resume --project projects/<project>

# 4) Status / reports
drugagent status --project projects/<project>
drugagent report --project projects/<project>
```

New users: start with [TUTORIAL.md](TUTORIAL.md) (per-parameter explanations,
how to read the report, number-reliability guide, FAQ, glossary).

**Run in the current folder**: add `--root .` (projects are created in
`./<name>/` instead of the deploy `projects/` dir; tools/weights stay put).
`status/resume/rerun/report` find such projects by bare name (checked
under `--root`, then the current directory, then the deploy dir). A
persistent default: `export DRUGAGENT_PROJECTS_ROOT=~/mywork`.

Common `run` options: `--library nci_npatlas|chembl35|pdbbind|<SDF path>` (default `nci_npatlas` — the NCI/DTP ∪ NPAtlas master library),
`--n-jobs 32`, `--md-ns 100 --md-reps 3`, `--max-steps 300` (agent step
budget), `--no-llm` (deterministic scripted mode), `--llm-base/--llm-model`
(LLM override). MD fine-tuning: `--md-salt` (ion concentration, M),
`--md-divalent MG --md-divalent-m 0.01` (divalent counterions),
`--md-extend-ns` (auto-extension step, ns), `--md-max-extensions`
(auto-extension rounds), `--md-burn-in-ps` (burn-in trimming),
`--md-system` (force the MD system). VHH: `--vhh-plddt-min`,
`--vhh-dock-flex`, `--vhh-dock-cdr-only`. Full option table: TUTORIAL §9.

Debugging hint: when an external command (gmx/vina/…) fails, the exception
already carries the tail of the log (last 40 lines) — just read the error,
no manual log digging.

## Agent Architecture (2.0)

- **Main loop** `drugagent/agent/loop.py`: each round the LLM
  (OpenAI-compatible endpoint) returns tool_calls → execute → feed results
  back → loop; terminates on the `finish` tool / budget exhaustion / human
  intervention.
- **Tools** `drugagent/agent/tools_*.py`: 41 tools = meta-tools
  (file/shell/decision/confirmation) + fine-grained tools for the five
  stages + five `run_*` whole-stage deterministic fallback tools (the 1.0
  pipeline kept intact).
- **Parameter sovereignty**: force fields listed by `gmx_env`, chosen by
  the agent; the MDP is written by the agent itself (`write_file`, the
  template is only a starting point); box / salt / barostat / dt all
  adjustable; every key decision is recorded via `record_decision`
  (`decisions.json`, shown in the report).
- **Self-debugging**: tool failure → agent reads the log (`read_file`) →
  diagnoses → patches with `edit_file` → retries; after 3 failed fixes for
  the same problem it calls `ask_human`.
- **Checkpoints**: 4 fixed (target / screening / design / md) + dynamic
  `ask_human`; `--auto` auto-passes and records.
- **State**: the project directory is the single source of truth.
  `state.json` stores the stage index; `agent/transcript.jsonl` records the
  full conversation + tool calls (resume = replay transcript and continue);
  two-layer idempotency: tool-artifact level (finished mdrun/docking/PDBQT
  are reused automatically) + **stage level (R11/G8: `run_*` tools skip the
  whole stage when its state.json section is complete; `force=true` forces
  a rerun — a crashed e2e no longer re-runs an already-finished 35-minute
  screening stage)**.

## LLM

Default `http://127.0.0.1:18080/v1` / `qwen3.8-27b-uncensored` (local
llama.cpp, must support function calling). Override with environment
variables: `DRUGAGENT_LLM_BASE_URL` / `DRUGAGENT_LLM_MODEL` /
`DRUGAGENT_LLM_API_KEY`.

## Directory Layout

```
/home/data/lrs/drug/drugagent/        # project root (persistent disk)
├── env/                              # python env (py3.12, conda) [local]
├── data/                             # [local]
│   ├── libraries/                    # small-molecule SDF / VHH library fasta
│   ├── weights/                      # ESMFold/ESM2 weights, torch/hf caches
│   ├── calibration/                  # calibration PDBs (1UBQ/1GFL/1CPS/1M17)
│   └── tools/
│       ├── gromacs/                  # GROMACS 2023.1 (self-built, amber99sb-ildn)
│       ├── vina/                     # autodock-vina
│       ├── gnina/                    # GNINA 1.3.1 ELF (CPU)
│       ├── RFdiffusion/              # RFdiffusion (hydra) + models/ weights + mpnn/
│       ├── vhh_scaffolds/            # 1EWN scaffold + secstruc adjacency + NOTE.md
│       └── 3Dmol/                    # 3Dmol-min.js
├── drugagent/                        # main package
│   ├── agent/                        # 2.0 core: loop/prompts/41 tools
│   ├── vendor/openfold/              # aqlaboratory/openfold + CPU patches
│   ├── modules/                      # modules A–E (tool backends)
│   ├── graph.py                      # 1.0 LangGraph state machine (fallback reuse)
│   ├── cli.py                        # typer CLI (run/resume/rerun/status/report/setup)
│   └── report/                       # interactive HTML + PDF (WeasyPrint)
├── TUTORIAL.md                       # user tutorial (zero background)
├── DESIGN.md                         # 2.0 architecture design (zh)
├── HANDOFF.md                        # environment & pitfalls, current state (zh)
├── ROUNDLOG.md                       # iteration log, per-round plan/results/reflection (zh)
├── projects/                         # one directory per run (01_target…05_md, reports/, agent/) [local]
├── tests/                            # pytest suite (218 fast + 15 slow e2e) [local, not in repo]
└── logs/                             # build/run logs [local]
```

## Status, Breakpoints & Reruns

- `status` (R11-enhanced): per-stage completion + key numbers
  (docking count / hit count / design count / MD ns & final RMSD / library
  name with fallback annotation), stage JSON artifacts on disk, the last 3
  tool failures (state.errors + ok=false calls in the transcript), last
  transcript entry; `decisions.json` records every judgment with rationale.
- `resume --project DIR`: replays `agent/transcript.jsonl` and continues
  from the breakpoint (artifact idempotency + stage-level reuse, R11/G8).
- `rerun --project DIR --stage X`: force-reruns a single stage (G8; other
  artifacts untouched).
- Reports: `projects/<project>/reports/report.html` (Plotly interactive
  charts + 3Dmol structures) and `report.pdf`.

## Tests

> The `tests/` directory is for local development (not published with the
> repo; `env/`, `data/`, `projects/` are likewise local deployment — see
> [TUTORIAL.md](TUTORIAL.md) §1 and HANDOFF.md for how to rebuild).

```bash
env/bin/python -m pytest tests/ -m "not slow" -q --basetemp=$PWD/data/fixtures/_ptmp   # fast unit tests (218 cases)
env/bin/python -m pytest tests/ -q --basetemp=$PWD/data/fixtures/_ptmp                  # incl. slow (needs vina/GROMACS/RF weights)
```

> Slow e2e builds intermittently fail on this machine's /tmp (tmpfs): a
> background cleaner removes pytest's tmp directory while gromacs is
> running. The code itself is stable — on a flaky failure, rerun with the
> `--basetemp` above pointing at a local disk.

## 1HVI Full-Module End-to-End (reference)

```bash
drugagent run --target 1HVI --modules all --fast --auto
```

On a 64-core CPU this completes all modules in ~3–5 h: target (1HVI dimer,
99 residues/monomer) → small-molecule screening (standardization + Vina
docking + GNINA rescoring + agent hit decision) → binder design (RFdiffusion
de novo + ESMFold complex scoring) → VHH dual track → system selection
(agent) → GROMACS EM + 5 ns × 3-replica MD → RMSD/RMSF/Rg analysis +
HTML/PDF reports. A complete instance: `projects/r10_e2e/` (includes
hinge-signal analysis of 5 structural domains).

## Known Caveats (summary)

Full details: HANDOFF.md / the corresponding round in ROUNDLOG.md (both in
Chinese); here only what matters to users.

- **Units**: GROMACS rms/gyrate output is nm (rmsf in Å, cluster cutoff 1.5
  in nm); the report displays RMSD converted to Å (×10), Rg in nm. The MD
  analysis chain is now unified on nm (R15 fixed an Å→nm bug that inflated
  every domain RMSD by ×10; regression-tested).
- **Tight-box PBC**: gmx analysis auto-runs compact-unwrap first (R17;
  when the protein nearly fills the box, wrapped coordinates break across
  box faces and create fake RMSD spikes — the classic case from the GROMACS
  forums).
- **Cofactors/heme**: MD auto-splits chains + ACPYPE/Gasteiger + embedded
  metals as standalone ions; known approximations: GAFF2/Gasteiger charges
  (no sqm/AM1-BCC), embedded metals have no coordination restraints (may
  drift out of plane in long MD).
- **MD metal ions**: single-atom metal residues (ZN/FE/MN/NI/CU/CO/MG/CA)
  are merged into the owning chain + crystal coordination distance
  restraints (flat-bottom potential); v1 limitations: +2 charge for all
  metals, soft restraints.
- **MD nucleic acids**: DNA/RNA chains are natively parametrized
  (amber99sb-ildn); single-residue nucleotides (NAD/ATP/FAD, etc.) go
  through the ACPYPE ligand path.
- **MD equilibration/convergence**: NVT 50 ps + NPT 100 ps (position
  restraints, C-rescale — with this build, Berendsen SIGSEGVs the
  production stage around 20 ps); after production, a convergence check
  (RMSD plateau + dominant cluster ≥ 50%), auto-extension if not converged
  (max 2 rounds by default).
- **VHH docking**: rigid + CDR fragments by default (fast) — a full-length
  VHH is a "giant ligand" (773 atoms ~100 min, cost ~O(n^1.9));
  per-fragment adaptive boxes eliminate wall-collision penalties
  (2.3e7 → 145 kcal/mol). Docking scores are coarse screening, not
  measured affinity.
- **pLDDT semantics**: VHH library pLDDT clusters at 30–35 because of
  CDR3 disorder (thresholds fast 35 / full 50 were calibrated to that);
  binder/VHH design pLDDT is ESMFold refolding confidence — low values are
  normal for de novo designs; geometric interface metrics (min distance /
  contact pairs) are the pose-level discriminator.
- **Scaffold content**: the 1EWN scaffold in `vhh_scaffolds/` is actually
  human AAG glycosylase core (the `vhh_` directory name is a historical
  prefix); scaffold-guided mode conditions only on SS + adjacency, so
  ~15 Å design-vs-scaffold drift is expected (quantified in the
  `scaffold_rmsd_a` field). See `data/tools/vhh_scaffolds/NOTE.md`.
- **Small-molecule libraries (R18)**: default `nci_npatlas` is a local
  master build — NCI/DTP open chemical repository (265,242, 3D, NSC
  ids) ∪ NPAtlas natural products (36,454, InChIKey), deduped by
  InChIKey (`scripts/merge_libraries.py`). The old `dtp` download mirror
  (dtpbase.org) is dead (502); legacy `--library dtp` auto-falls back to
  the master library, annotated in state/report as "fallback for dtp".
- **Structure-pitfall pre-check**: `analyze_pdb` scans multi-MODEL / altloc
  / missing residues / metals / nucleic acids / disordered regions and
  suggests actions; `run_target_prep` auto-fixes the two safe items, the
  rest go to the agent's judgment (recorded in decisions); report §1 shows
  all findings and repairs.
- **Modified residues (R17)**: when the target's "ligands" are all modified
  residues (CPM crosslinks / metals / cofactors, a 64-name set), MD
  auto-prefers the apo path (drop the residues, run the clean protein);
  GFP-chromophore-type CRO (missing backbone atoms) still needs
  residue-level repair (R18 gap).
- **Rigid-funnel limitation**: default docking/design uses a single rigid
  conformation; flexibility is sampled only in the MD stage. The report's
  §5 "flexibility interpretation" automatically distinguishes "global
  drift = conformational sampling / domain motion" from "local
  unfolding" (per-chain self-fit + secondary structure + RMSF + clustering
  + domain vs-rest).
- **Disk**: env + weights + libraries ~40 GB; MD trajectories grow with ns
  (5 ns × 3 replicas ≈ 2 GB).
- **GPU**: optional. Without `/dev/nvidia*`, everything runs on CPU
  (64 cores are enough).
- **GNINA**: the ELF binary depends on
  `env/lib/python3.12/site-packages/nvidia/*/lib` (the scripts handle
  LD_LIBRARY_PATH automatically).
- **Binder sequence source**: RF design PDBs carry GLY residue names
  (no sequence); built-in ProteinMPNN generates real sequences (fallback
  heuristic on MPNN failure, annotated in the report); ESMFold scoring
  uses the real sequences.
