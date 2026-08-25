# DrugAgent 2.0 User Tutorial (Zero Background → Self-Sufficient)

> Reader assumption: comfortable with basic command line; no prior
> molecular-docking knowledge required.
> Environment assumption: deployment already done at
> `/home/data/lrs/drug/drugagent` (`env/` and `data/` present).
> Companion docs: [README.md](README.md) (tech stack & overview) ·
> [DESIGN.md](DESIGN.md) (architecture, zh) · [HANDOFF.md](HANDOFF.md)
> (environment & pitfalls, authoritative, zh) · [ROUNDLOG.md](ROUNDLOG.md)
> (iteration history, zh).
> 🇨🇳 **中文版**: [TUTORIAL.zh-CN.md](TUTORIAL.zh-CN.md)

---

## 0. What Can This System Do for You

Give it a protein (PDB ID / PDB file / even just an amino-acid sequence)
and it will:

| You want to… | Module | In one sentence | You get |
|---|---|---|---|
| **Find drugs**: which small molecules bind my protein | `screen` | Virtual screen a compound library (hundreds of thousands of molecules) | Ranked candidate molecules + binding scores + binding poses (3D structure files) |
| **Design a new protein**: a small protein that binds/inhibits my protein | `binder` | De novo design of a 60–80 residue protein (RFdiffusion) | Design structure PDB + amino-acid sequence + confidence scores |
| **Design a nanobody**: a VHH | `vhh` | Dual track: library screening + scaffold-guided de novo design | Candidate VHHs + scores + scaffold-fidelity metrics |
| **See protein dynamics**: is my protein stable? which regions are soft? | `md` | Molecular dynamics (GROMACS) + automated flexibility diagnostics | RMSD/RMSF/domain analysis + automatic interpretation (rigid / hinge / unfolding) |

A single run can do just one of these (`--modules screen`) or all of them
(`--modules all`, full chain: target → screening → binder → VHH → MD →
summary report).

**Deliverable**: each run produces a project directory `projects/<name>/`;
the most important file is `reports/report.html` — an interactive report
(3D structures + charts + automatic interpretation text) you open in a
browser.

---

## 1. Environment Check (30 seconds)

```bash
cd /home/data/lrs/drug/drugagent

# ① python environment (should print 3.12.x)
env/bin/python --version

# ①-b install the `drugagent` command (one-time; editable install — code
#     edits apply immediately; the symlink puts it on your user PATH so you
#     can call `drugagent ...` from ANY directory)
env/bin/pip install -e . --no-deps
ln -sf $PWD/env/bin/drugagent ~/.local/bin/drugagent

# ② (optional) is the local LLM online — not needed if you only use --no-llm scripted mode
curl -s -m 5 http://127.0.0.1:18080/v1/models | head -c 200

# ③ environment completeness (idempotent: fills what's missing, skips what exists; must run on first deploy)
drugagent setup

# ④ code health check (~2 min, should show 218 passed)
env/bin/python -m pytest tests/ -m "not slow" -q -p no:warnings --basetemp=$PWD/data/fixtures/_ptmp
```

Notes:
- `--basetemp=$PWD/data/fixtures/_ptmp` is mandatory: this machine's
  `/tmp` is a tmpfs with a background cleaner that deletes test temp
  files; local disk keeps tests stable.
- Disk: env + weights + libraries ~40 GB; MD trajectories grow with
  simulated time (5 ns × 3 replicas ≈ 2 GB).
- Compute: 64 CPU cores is enough (GPU optional; without one, everything
  runs on CPU automatically).

---

## 2. First Run: Small-Molecule Docking (screen)

### 2.1 The command

```bash
drugagent run \
  --target 1HVI \
  --modules screen \
  --fast \
  --auto --no-llm \
  --name mydock
```

Every parameter explained:

| Parameter | Meaning | You can substitute |
|---|---|---|
| `--target 1HVI` | Target. PDB IDs are downloaded from RCSB automatically | Any PDB ID (e.g. `5P21`); a local file path (`/path/protein.pdb`); a FASTA file; a raw amino-acid sequence (auto-modeled with ESMFold) |
| `--modules screen` | Run only small-molecule screening | `binder` / `vhh` / `md` / `screen,binder,md` (comma-separated) / `all` |
| `--fast` | Validation scale (small library, short MD; minutes to hours) | drop `--fast` = production scale (big library, 100 ns MD; much more compute) |
| `--auto` | Checkpoints auto-pass (nobody needed at the terminal) | drop it = interactive mode; the agent stops at milestones and asks you (§8) |
| `--no-llm` | Deterministic scripted mode, no LLM service required | drop it = LLM-main-program mode (local llama.cpp must be online) |
| `--name mydock` | Project name (directory `projects/mydock/`) | any name; defaults to a timestamp |

> For your first runs, always use `--auto --no-llm`: most stable, fully
> reproducible.

### 2.2 What happens in the background (in order)

1. **Target prep**: download/read PDB → integrity check (multi-MODEL
   overlays? altloc? missing residues? metals? nucleic acids? disordered
   termini?) → auto-repair of the safe items → cleaning → pocket
   detection → PDBQT conversion. Artifacts: `projects/mydock/01_target/`
   (`receptor.pdbqt`, `pocket.json`, `clean.pdb`).
2. **Library resolution**: if the `dtp` library is missing/corrupt it
   auto-falls back to `chembl35_small` (~50k compounds; state and report
   are annotated "fallback for dtp").
3. **Standardization + prefilter**: RDKit standardization, physchem
   filters (Lipinski, etc.).
4. **Docking**: parallel Vina (32 jobs by default, `--n-jobs` to change).
5. **Rescoring**: GNINA rescoring of top results.
6. **Hit decision**: top hits are selected (in agent mode, the
   co-crystallized ligand is redocked as a positive control and the
   decision is recorded).
7. **Report**: `reports/report.html` + `report.pdf` generated.

### 2.3 Watch progress

```bash
drugagent status --project projects/mydock
```

Output: per-stage ✓/✗ + key numbers (docking count / hit count / library
name with fallback annotation) + artifact list + last 3 errors.

### 2.4 Get the results

```bash
xdg-open projects/mydock/reports/report.html     # open in a browser (recommended)
ls projects/mydock/02_screen/                     # docking detail files
```

Report §2 is the screening result: ranked candidate table (SMILES,
docking score, rescore) + 3D binding poses (rotatable/zoomable).

### 2.5 How to read docking scores

- Units kcal/mol, **lower is better** (more negative = tighter binding).
- It's a **relative ranking tool**: comparisons within the same batch of
  molecules are reliable; cross-system / cross-parameter comparisons need
  care.
- A rigid-docking score ≠ experimental affinity (no free-energy
  perturbation); it's for "enriching 50 out of 50,000", not "measuring a
  Kd".

---

## 3. Anatomy of a Project Directory

```
projects/mydock/
├── state.json              # stage index (read by status/resume; don't hand-edit)
├── 01_target/
│   ├── clean.pdb           # cleaned protein
│   ├── receptor.pdbqt      # docking receptor
│   └── pocket.json         # pocket info
├── 02_screen/              # docking/scoring/hit details (JSON + PDBQT poses)
├── 03_binder/              # (if binder ran) design structures + scored.json
├── 04_vhh/                 # (if vhh ran)
├── 05_md/                  # (if md ran)
│   ├── build/              # build artifacts (top/gro/tpr/equilibration)
│   ├── md_rep1/            # replica 1 trajectory (md.xtc)
│   ├── md_rep2/  md_rep3/
│   └── analysis/           # RMSD/RMSF/clustering xvg + interpretation
├── reports/
│   ├── report.html         # ← main deliverable
│   └── report.pdf
└── agent/
    ├── transcript.jsonl    # full LLM tool-call trace ("how it thought")
    └── decisions.json      # every key judgment + rationale
```

**Idempotency**: re-running the same `run` command does not recompute
finished work (two-layer reuse: tool-artifact level + stage level);
after a restart, only the missing parts are filled in. Force a stage
recompute with `rerun` (§7).

---

## 4. Binder Design (de novo small protein)

### 4.1 The command

```bash
drugagent run --target 1HVI --modules binder \
  --fast --auto --no-llm --name mybinder
```

### 4.2 What happens in the background

1. Target prep (as above).
2. **RFdiffusion**: designs a 60–80 residue small-protein backbone on the
   target surface (`miniprotein` mode).
3. **ProteinMPNN**: generates the real amino-acid sequence for the
   backbone (RF-design PDBs carry GLY placeholder residue names; the MPNN
   sequence is the real one).
4. **ESMFold scoring**: refolds from the real sequence and computes
   pLDDT / interface confidence. Design entries also carry **geometric
   interface metrics** (min Cα distance design-chain-to-target, contact
   pairs <6 Å) — sequence-level interface pLDDT is blind to "where the
   design actually sits"; geometric contacts are the pose-level
   discriminator.
5. Report: ranked design table + structures.

### 4.3 Cost and results

- With `--fast`, roughly **20–40 minutes** (RF runs on CPU; one of the
  slowest modules in the whole chain).
- Results: report §3 of `projects/mybinder/reports/report.html`; design
  PDBs in `03_binder/` (`design_N.pdb` monomer, `design_N_complex.pdb`
  complex) — drop into PyMOL/ChimeraX to view.
- **How to read pLDDT**: de novo designs with low pLDDT are normal
  (ESMFold has never seen this sequence). Focus on **interface pLDDT**
  (the design-chain/target contact region) and the geometric contact
  metrics (the "target distance (Å)" column; a ⚠no-contact flag means the
  design never touched the target — be wary).
- Caching: `03_binder/scored.json` caches by design-file mtime; reruns
  only revalidate designs overwritten by RF (rerun dropped from ~3 h to
  minutes).

---

## 5. VHH Nanobody (dual track)

```bash
drugagent run --target 1HVI --modules vhh \
  --fast --auto --no-llm --name myvhh
```

Two tracks in parallel:

- **Track A (library screening)**: VHH library
  (`data/libraries/vhh_library.fasta` or synthetic) → ESMFold modeling →
  **pLDDT filter** (fast threshold 35 / full 50, override with
  `--vhh-plddt-min`) → **rigid docking** (fast default: **CDR fragments
  only** — 1–3 fragments cut along low-pLDDT regions, each with its own
  adaptive box, ~15× faster than full-length) → parallel screening.
- **Track B (de novo design)**: RFdiffusion scaffold-guided (1EWN
  scaffold as fold condition) → ProteinMPNN → scoring.

### Results & interpretation

- Report §4: candidate table (docking score / fragment detail / pLDDT) +
  **library pLDDT distribution histogram** (threshold line + p50/p90 +
  count above threshold) + composite ranking.
- **scaffold_rmsd_a** (Track B designs): Kabsch RMSD (Å) of the designed
  backbone vs the scaffold backbone. Scaffold-guided mode's official
  semantics: "SS pattern + block adjacency as condition; sequence and
  fine structure may vary", so **~15 Å drift is expected behavior**
  (fold topology preserved, local structure preserved, length exact);
  this field quantifies approximation fidelity, not failure. Scaffold
  content: `data/tools/vhh_scaffolds/NOTE.md` (1EWN is actually human AAG
  glycosylase core; the `vhh_` directory name is a historical prefix).
- VHH docking scores are also coarse screening (ranking "binding
  capability", not measured affinity).

---

## 6. MD Simulation (protein flexibility diagnostics)

```bash
# Ligand-bound complex (target with its ligand, or a specified system)
drugagent run --target 1HVI --modules md \
  --fast --auto --no-llm --name mymd

# Pure protein target (apo — protein alone: flexibility/hinge baseline)
drugagent run --target data/calibration/1UBQ.pdb \
  --modules md --fast --auto --no-llm --name myapo

# Custom length/replicas (without --fast the default is 100 ns×3; common override)
drugagent run --target 1HVI --modules md \
  --md-ns 5 --md-reps 3 --auto --no-llm --name mymd5
```

### 6.1 What happens in the background

1. **System selection**: which complex to simulate (input ligand /
   screening hit / designed binder / VHH / apo protein; "ligands" that
   are all modified residues auto-prefer apo; force with `--md-system`).
2. **Build**: pdb2gmx (amber99sb-ildn) → ligand ACPYPE/GAFF2 → solvation →
   ions (`--md-salt`, default 0.15 M) → energy minimization.
3. **Equilibration**: NVT (50 ps) → NPT (100 ps), with position
   restraints, C-rescale pressure control.
4. **Production MD**: N ns × R replicas (independent replicas forked from
   the equilibrated end-state), v-rescale 300 K. **Convergence check +
   auto-extension**: if RMSD plateau + dominant-cluster fraction are not
   met, each replica is extended and merged automatically (`--md-extend-ns`
   step size, `--md-max-extensions` rounds).
5. **Analysis** (all automatic):
   - core RMSD (backbone, vs the initial structure)
   - per-chain RMSD: **inter-chain relative** (fit the largest chain,
     measure the others) + **intra-chain self-fit** (not inflated by
     reference-chain flexibility; use this to judge "unfolding")
   - RMSF (per-residue flexibility) + **flexible-region localization**
     (contiguous high-RMSF runs → residue intervals)
   - radius of gyration Rg, gmx clustering (1.5 nm cutoff)
   - per-frame secondary-structure fraction (hydrogen-bond rules, DSSP-like)
   - **structural-domain RMSD** (SS-extracted domains, per-domain Kabsch
     self-fit → intra-domain deformation)
   - **domain vs-rest-of-protein relative RMSD** (hinge/allosteric rigid
     motion signal) + **diameter-normalized final_norm** (comparable
     across domain sizes)
   - **automatic interpretation** (rigid-baseline comparison, §6.2)

### 6.2 How to read MD results (important)

Report §5 + the xvg curves in `analysis/`. **Don't judge from one
overall RMSD number** — there is now a calibrated reading framework:

| Metric | Meaning | How to judge |
|---|---|---|
| Core RMSD (backbone self-fit) | overall conformational drift | Against the **rigid baseline**: ubiquitin (76-residue monomer, canonically stable) ends ≤2.1 Å and domain norm ≤0.08 at 2 ns; scale by √t to your run length (upper bound = 2.1×√(ns/2) Å). Above it → use the next rows to find the cause |
| Intra-chain self-fit RMSD | internal deformation of a single chain | <0.6 nm (6 Å) → chain did not unfold; >6 Å → internal unfolding / major rearrangement |
| Inter-chain relative RMSD | relative motion of multiple chains | much larger than self-fit → dimer/oligomer relative motion (e.g. HIV-PR at 5 ns measured 1.15–3.88 nm while chain self-fit was only 0.05–0.14 nm) |
| Domain vs-rest final_norm | rigid motion of a domain vs the rest of the protein (hinge signal) | >0.08 exceeds the rigid baseline → genuine domain motion / hinge; report splits "internally stable (<3 Å)" from "major rearrangement" |
| Flexible regions (contiguous high-RMSF runs) | local flexible loops | the report lists residue intervals; judge against pocket location whether they affect binding |
| Cluster fractions | number of conformational states | single dominant cluster (>90%) = one state; multiple clusters = ensemble |
| SS fraction change | secondary-structure preservation | 64%→64% stable; a clear drop → unfolding signs |

**Time-scale caveat**: `--fast` MD is a 2–5 ns **stability check**, not a
production free-energy / convergence conclusion; for hinge dynamics use
full mode (default 100 ns) or `--md-ns 50`+. The first ~100–200 ps of MD
is a relaxation segment (early jumps in the curves are normal; the
analysis side trims a 100 ps burn-in by default).

### 6.3 Cost reference (64-core CPU)

- 76-residue monomer (ubiquitin): 2 ns × 3 ≈ **18 min**
- 198-residue dimer (HIV-PR): 5 ns × 3 ≈ **1.5 h**
- 290 residues + ligand (CDK2): 2 ns × 2 ≈ 30–40 min
- When auto-extension triggers, each replica runs another `--md-extend-ns`
  (default = `--md-ns`)

---

## 7. Breakpoints / Reruns / Status

```bash
# status (run anytime)
drugagent status --project projects/mydock

# resume from a breakpoint (after crash/restart; replays transcript + artifact reuse)
drugagent resume --project projects/mydock

# force-rerun one stage only (other artifacts untouched)
drugagent rerun --project projects/mydock --stage md
drugagent rerun --project projects/mydock --stage vhh --no-with-report

# regenerate only the report
drugagent report --project projects/mydock
```

Common stage names: `target_prep` / `screening` / `binder` / `vhh` /
`md` / `report`.

**When to use `rerun`**: you changed code/parameters and want one stage
recomputed; a stage's artifacts look suspect. **When to use `resume`**:
the process was killed or the machine rebooted and you want to continue.

---

## 8. Talking to the Agent (interactive mode)

Default (i.e. **without** `--auto`) is interactive: the agent stops at
milestones and waits for your input.

```bash
drugagent run --target 1HVI --modules all --fast
# runs in the foreground and stops at these points:
#   [milestone checkpoint: target]    target pre-check summary → approve/modify/abort
#   [milestone checkpoint: screening]
#   [milestone checkpoint: design]
#   [milestone checkpoint: md]        (includes the MDP plan: force field/box/salt/barostat)
#   [dynamic confirmation]            (raised by the LLM anytime: e.g. binder type, which complex)
# prompt: answer/choose >
```

- Enter an option number (`1`) = approve; enter `x` / `abort` = abort the
  run;
- **You can also type free text** (e.g. "set salt to 0.1, use c-rescale
  for the barostat") — it is passed back to the LLM as your opinion; the
  LLM adjusts accordingly and records a `record_decision`.
- Every LLM judgment is written to `agent/decisions.json` (stage / choice
  / rationale); the full tool-call trace is in
  `agent/transcript.jsonl` (report §6 has a summary).

### Mode selection cheat sheet

| Scenario | Mode |
|---|---|
| batch / CI / unattended | `--auto --no-llm` |
| reproducibility first | `--no-llm` (±`--auto`) |
| you want the agent to follow your preferences (parameter choices, hit selection) | interactive (drop `--auto`; LLM must be online) |
| LLM service not running | `--no-llm` (deterministic scripted mode) |

### Parameter sovereignty (LLM mode)

Force field, MDP (integrator / barostat / salt / box margin), which
complex to simulate, how many hits to keep, binder type — all decided by
the LLM and recorded. `--no-llm` mode uses defaults.

---

## 9. Full Option Reference (run)

| Option | Default | Description |
|---|---|---|
| `--target` | (required) | PDB ID / PDB path / FASTA / raw sequence |
| `--modules` | `all` | `screen,binder,vhh,md` comma-separated subset |
| `--fast / --no-fast` | full | validation scale vs production scale |
| `--auto / --no-auto` | interactive | auto-pass checkpoints |
| `--name` | timestamp | project name |
| `--library` | `dtp` | `dtp` / `chembl35` / `pdbbind` / custom SDF path (dtp auto-falls back to chembl35_small when corrupt) |
| `--n-jobs` | 32 | docking/scoring parallelism |
| `--md-ns` | 100 (full) / 5 (fast) | MD length, ns |
| `--md-reps` | 3 | MD replicas |
| `--md-salt` | 0.15 | ion concentration, M |
| `--md-divalent` / `--md-divalent-m` | none / 0.01 | divalent counterion (MG/CA) and concentration |
| `--md-extend-ns` | = md_ns | auto-extension step, ns (0 = use md_ns) |
| `--md-max-extensions` | 2 | max auto-extension rounds |
| `--md-burn-in-ps` | 100 | production-segment burn-in trimming, ps |
| `--md-system` | agent chooses | `input_ligand`/`screening_hit`/`binder`/`vhh`/`apo` to force the system |
| `--vhh-plddt-min` | fast 35 / full 50 | VHH library pLDDT threshold |
| `--vhh-dock-flex` | rigid | VHH flexible docking (an order of magnitude slower) |
| `--vhh-dock-cdr-only` | on in fast | CDR-fragment docking (~15× speedup) |
| `--max-steps` | 300 | LLM step budget |
| `--llm-base` / `--llm-model` | local 18080 / qwen3.8-27b | LLM override (env vars also work) |
| `--no-llm` | off | deterministic scripted mode |

---

## 10. FAQ / Troubleshooting

**Q1: My protein has a ligand / metals / DNA. Can it run?**
- Ligand: screen/binder use it directly; MD goes through the ligand
  parametrization path (ACPYPE GAFF2).
- Metals (single-atom ZN/FE/MG…): MD merges them into the owning chain +
  coordination distance restraints automatically.
- Nucleic-acid chains: natively parametrized in MD (amber99sb-ildn
  includes DNA/RNA).
- Heme-like cofactors: MD auto-splits chains + ACPYPE + embedded metals
  as standalone ions (charge approximation, no coordination restraints —
  watch long MDs).
- **Known to fail**: protein-internal **modified residues** (GFP
  chromophore CRO type — pdb2gmx needs full-atom parameters for them):
  the MD build fails (R18 gap); docking/binder are unaffected. Targets
  with crosslinks (e.g. crambin's CPM) now auto-route to the apo path
  (modified residues treated as protein, dropped, then run the clean
  protein).

**Q2: I see "library dtp missing/corrupt (0 bytes)". Important?**
No — it auto-falls back to chembl35_small (50k compounds), annotated as
`fallback for dtp` in state and report. To repair the library:
`drugagent setup --libraries dtp` (the dtpbase.org
mirror is flaky; `--libraries pdbbind` soft-fails without breaking setup).

**Q3: It died halfway / the machine rebooted?**
`resume --project <dir>`. Artifacts are idempotent (finished
mdrun/docking/PDBQT are reused automatically); no burning compute from
scratch.

**Q4: I only want to recompute MD (e.g. after changing parameters)?**
`rerun --project <dir> --stage md` (re-runs build + MD + analysis with the
new parameters; other stage artifacts untouched).

**Q5: Docking score −12 vs −15 — how big a difference?**
Within the same library: 3 kcal/mol ≈ ~7× difference in binding
constant at 298 K (ΔG=RTlnK). Don't force-compare across libraries or
parameters.

**Q6: Binder pLDDT is only 40 — is the design a failure?**
Not necessarily. Low pLDDT for de novo designs is normal (ESMFold
refolding confidence, not confidence of reproducing the RF-designed
conformation). Focus on interface pLDDT + geometric contacts (the
"target distance" column + contact-pair count); report §3 ranks
compositely.

**Q7: The report shows RMSD 5 Å — did the protein fall apart?**
Three-step check (table in §6.2): intra-chain self-fit <6 Å → not
unfolding; inter-chain relative >> self-fit → oligomer relative motion;
domain norm >0.08 → genuine hinge motion. At `--fast` short lengths,
1–2 Å drift is completely normal (rigid-protein 2 ns baseline upper bound
~2.1 Å).

**Q8: No GPU — does it work?**
Yes. Without `/dev/nvidia*` everything runs on CPU (RF/ESMFold/vina/
GROMACS all have CPU paths). Full-module fast e2e on 64 CPU cores takes
~3–5 h.

**Q9: I only have a sequence, no structure?**
`--target` accepts FASTA or a raw sequence directly (auto-modeled with
ESMFold; single sequence uses default conditions — read scoring as
"refolding confidence").

**Q10: I want to use my own small-molecule library?**
`--library /path/to/my.sdf` (SDF format; pre-standardizing with RDKit is
recommended, though the pipeline standardizes again).

**Q11: Running out of disk — what can I clean?**
`projects/<old>/05_md/md_rep*/` (xtc trajectories, grow with ns) and
PDBQT caches in `projects/*/02_screen/`. Don't touch `env/` or `data/`
(~40 GB). `logs/` is cleanable.

**Q12: Where are the errors when an external command fails?**
The exception already carries the log tail (last 40 lines); `status` also
lists the last 3 tool failures. Raw logs live inside the project
directory (e.g. `05_md/md_rep1/md.log`).

---

## 11. Number Reliability Guide (what to trust, what not to)

| Number | Trust it for | Don't trust it for |
|---|---|---|
| Vina docking score | relative ranking within a library (enrichment) | absolute affinity; cross-system comparison |
| GNINA rescore | same (more physical; CASP-validated) | same |
| ESMFold pLDDT (VHH library) | foldability ranking within the library | systematically low for CDR3 loops (a 30–35 cluster is normal, not a bug) |
| binder pLDDT / interface score | relative ranking + interface geometry (contacts / distances) | an absolute "it binds" claim (de novo designs need experimental validation) |
| scaffold_rmsd_a | quantifying scaffold-guided fidelity | the intuition that lower is better (the official mode allows fine-structure variation; ~15 Å is the expected magnitude) |
| MD RMSD (fast 2–5 ns) | stability check, flexible-region localization, hinge signals | free-energy / convergence conclusions (need full length + passed convergence check) |
| MD flexibility interpretation | qualitative rigid/hinge/unfolding distinction (calibrated against a cross-target rigid baseline) | quantitative conformational dynamics (time too short) |

**Rule of thumb**: DrugAgent's numbers are all "ranking + qualitative
judgment" grade; quantitative conclusions (Kd, ΔG, rate constants) need
corresponding experiments or advanced computations (FEP, umbrella
sampling).

---

## 12. Glossary

- **PDB ID**: Protein Data Bank identifier (e.g. 1HVI = HIV protease),
  downloadable from RCSB.
- **apo / complex**: protein without ligand / protein with ligand.
- **pocket**: a protein-surface recess suitable for binding small
  molecules; detected automatically.
- **PDBQT**: Vina's docking format (atom types / charges / molecular
  graph).
- **docking score (kcal/mol)**: empirical binding-free-energy estimate;
  lower is better; a relative quantity.
- **virtual screening**: "try-fitting" every library molecule into the
  pocket, score, and rank.
- **RFdiffusion**: diffusion-model protein-backbone design (Rosalind).
- **ProteinMPNN**: sequence generation for a given backbone
  (inverse folding).
- **ESMFold**: sequence → 3D structure + pLDDT confidence.
- **pLDDT**: local structure-prediction confidence (0–100); high =
  reliable.
- **VHH / nanobody**: single-domain antibody (~120 residues); CDRs are
  the complementarity-determining regions (binding site).
- **scaffold-guided**: design conditioned on a known backbone's SS
  pattern + block adjacency.
- **MD (molecular dynamics)**: Newtonian integration of atomic motion to
  sample conformations.
- **force field**: parametrized function of inter-atomic interactions
  (amber99sb-ildn for protein / GAFF2 for ligands).
- **RMSD**: root-mean-square deviation between structures; a
  conformational-difference metric (nm or Å, 1 nm = 10 Å).
- **RMSF**: per-residue time-averaged RMSD; a flexibility metric.
- **Kabsch**: optimal rigid-body alignment algorithm for two structures.
- **PBC (periodic boundary conditions)**: "seamless" simulation-box
  handling; trajectories must be unwrapped before analysis.
- **convergence check**: RMSD plateau + dominant-cluster fraction;
  auto-extension when not met (R6).
- **conformational selection**: screening against multiple MD-sampled
  conformers (consensus averaging).
- **hinge**: a flexible joint between protein domains; key for
  allosteric / functional motion.

---

## 13. Reference Runs (open and compare directly)

| Project | Content | What to look at |
|---|---|---|
| `projects/r10_e2e/` | HIV-PR full modules (screen+binder+vhh+md, 5 ns×3) | all 6 report sections; vs-rest analysis of 5 structural domains (0.8–2.4 Å genuine hinge signals); VHH design scaffold_rmsd_a = 15.08 Å; vhh_30 frag1 −8.66 kcal/mol |
| `projects/rigid_ubq/` | ubiquitin 1UBQ apo 2 ns×3 | the **rigid baseline** (core RMSD finals 0.11–0.21 nm, domain norm 0.04–0.08) |
| `projects/smoke_screen_demo/` | 1HVI screen only (fast) | minimal docking run example |

```bash
# open any project's report
xdg-open projects/r10_e2e/reports/report.html
```

---

*Tutorial reflects round-17 code (commit `4787a88`). Code evolution is
tracked in HANDOFF.md / ROUNDLOG.md; if the tutorial and code disagree,
the code wins — and note it in the ROUNDLOG reflection section.*
