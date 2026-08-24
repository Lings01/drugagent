"""Central configuration: paths, tools, LLM endpoint, simulation defaults."""
from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path


# --------------------------------------------------------------------------- #
# Root location
# --------------------------------------------------------------------------- #
def find_writable_root() -> Path:
    """Pick the first writable candidate for the project root."""
    candidates = [
        Path(os.environ.get("DRUGAGENT_ROOT", "")).expanduser()
        if os.environ.get("DRUGAGENT_ROOT")
        else None,
        Path("/home/data/lrs/drug/drugagent"),
        Path("/home/lrs/drugagent"),
        Path("/home/drugagent"),
        Path("/tmp/drugagent"),
        Path.home() / "drugagent",
    ]
    for c in candidates:
        if c is None:
            continue
        try:
            c.mkdir(parents=True, exist_ok=True)
            probe = c / f".probe_write.{os.getpid()}"
            probe.write_text("ok")
            probe.unlink(missing_ok=True)
            return c
        except OSError:
            continue
    raise RuntimeError("No writable location found for drugagent root")


ROOT = Path(os.environ.get("DRUGAGENT_ROOT", "/tmp/drugagent")).expanduser() \
    if os.environ.get("DRUGAGENT_ROOT") else find_writable_root()

DATA = ROOT / "data"
LIBRARIES = DATA / "libraries"
WEIGHTS = DATA / "weights"
MODELS = DATA / "models"
TOOLS = DATA / "tools"
PROJECTS = ROOT / "projects"
LOGS = ROOT / "logs"
ENV_DIR = ROOT / "env"

for _p in (DATA, LIBRARIES, WEIGHTS, MODELS, TOOLS, PROJECTS, LOGS):
    _p.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------- #
# LLM ("brain") endpoint - OpenAI-compatible (llama.cpp)
# --------------------------------------------------------------------------- #
@dataclass
class LLMConfig:
    base_url: str = os.environ.get("DRUGAGENT_LLM_BASE_URL", "http://127.0.0.1:18080/v1")
    api_key: str = os.environ.get("DRUGAGENT_LLM_API_KEY",
                                  "sk-d675373dcd907170f3a1f87d28e19a7c55d88cacd30ed177")
    model: str = os.environ.get("DRUGAGENT_LLM_MODEL", "qwen3.8-27b-uncensored")
    temperature: float = 0.2
    max_tokens: int = 3000
    timeout: float = 600.0

    @classmethod
    def from_env(cls) -> "LLMConfig":
        return cls()


# --------------------------------------------------------------------------- #
# External tool locations (auto-discovered)
# --------------------------------------------------------------------------- #
GMX_ENV_BIN = Path("/usr/local/anaconda3/envs/gmx/bin")


def which_tool(name: str, search_dirs: list[Path] | None = None) -> Path | None:
    dirs = [d for d in (search_dirs or []) if d.is_dir()]
    for d in dirs:
        p = d / name
        if p.is_file() and os.access(p, os.X_OK):
            return p
    return Path(shutil.which(name)) if shutil.which(name) else None


def gmx_bin(name: str) -> Path:
    """Locate a GROMACS binary (mdrun/pdb2gmx/...) in the known env or PATH."""
    p = GMX_ENV_BIN / name
    if p.is_file():
        return p
    found = shutil.which(name)
    if found:
        return Path(found)
    raise FileNotFoundError(f"GROMACS binary '{name}' not found")


# --------------------------------------------------------------------------- #
# Simulation / screening defaults
# --------------------------------------------------------------------------- #
@dataclass
class Defaults:
    # screening
    default_library: str = "dtp"
    screen_max_ligands: int = 50_000
    screen_min_ligands: int = 500
    screen_time_budget_h: float = 24.0
    dock_exhaustiveness_fast: int = 8
    dock_exhaustiveness_final: int = 16
    # R9: MD-rep conformers are screened for CONSENSUS, not precision —
    # a lower exhaustiveness keeps the pool cheap (flex docks are ~9x
    # slower per 16 units of exhaustiveness)
    dock_md_rep_exhaustiveness: int = 12
    # flex side-chain cutoff (A) around the ligand when building flex
    # receptor files per conformer
    flex_cutoff_ang: float = 5.0
    dock_rescore_topn: int = 200
    n_hits: int = 20
    # binder
    n_binder_designs: int = 8
    binder_len: int = 70
    # vhh
    vhh_lib_size: int = 20_000
    vhh_screen_n: int = 1000
    vhh_de_novo_designs: int = 8
    # R11/G9: mean-structure pLDDT gate for track A. ESMFold pLDDT for
    # VHHs (real or synthetic) concentrates at 30-35 because CDR3 loops
    # are disordered in the model — a 70 gate would screen out almost
    # everything. 50 keeps the full-mode bar meaningful while fast mode
    # relaxes further (see fast dict) because docking, not modeling, is
    # the fast-mode cost driver.
    vhh_plddt_min: float = 50.0
    # R11/G10: dock full VHHs as RIGID bodies by default. A folded domain
    # from a single ESMFold model has one conformation; letting hundreds of
    # torsions roam turns a triage screen into a 30-min-per-candidate job
    # (this vina build also runs single-core on big ligands). Set
    # vhh_dock_flex=true via options for a (much slower) flexible screen.
    vhh_dock_flex: bool = False
    # md
    md_ns: float = 100.0
    md_reps: int = 3
    md_timestep_ps: float = 0.002
    md_save_ps: float = 10.0
    md_box_margin_nm: float = 1.0
    md_salt_m: float = 0.15
    # R10: divalent counterions (MG/CA); None disables. genion adds the
    # cation + 2x Cl- to stay neutral; the ion atomtype is injected
    # (GROMACS ion.itp ships only NA/CL/K)
    md_divalent: str | None = None
    md_divalent_m: float = 0.01
    # R5: equilibration stages (once per build, not per replica) +
    # production burn-in trimmed from the analysis
    md_eq_nvt_ps: float = 50.0
    md_eq_npt_ps: float = 100.0
    md_burn_in_ps: float = 100.0
    # R6: convergence-gated auto-extension (each round extends every
    # replica by md_extend_ns, or md_ns when 0)
    md_max_extensions: int = 2
    md_extend_ns: float = 0.0
    md_converge_min_len_ps: float = 50.0
    md_converge_rmsd_drift_nm: float = 0.05
    md_converge_min_cluster: float = 0.5
    md_forcefield: str = "amber19sb"
    md_ligand_forcefield: str = "gaff2"
    # fast mode (CI-like)
    fast: dict = field(default_factory=lambda: {
        "screen_max_ligands": 2000,
        "screen_time_budget_h": 2.0,
        "n_binder_designs": 2,
        "vhh_lib_size": 2000,
        # R11/G9: widen the hit surface — measured pLDDT of 100 modeled
        # VHHs: p50=31.5, p90=34.5, only 1/100 > 45. A 45 gate left a
        # single docking candidate (Module D had no statistical power).
        # 35 passes ~6-10% of candidates; R11/G10 rigid docking makes the
        # extra candidates cheap. 80 keeps --fast a validation scale while
        # giving the dock step a real sample.
        "vhh_screen_n": 80,
        "vhh_plddt_min": 35.0,
        "vhh_de_novo_designs": 2,
        "md_ns": 5.0,
        "md_reps": 3,
        # short equilibration in fast mode (big systems: 50/100 ps is
        # hours of wall time)
        "md_eq_nvt_ps": 10.0,
        "md_eq_npt_ps": 20.0,
        "md_max_extensions": 1,
        "dock_exhaustiveness_fast": 4,
        "dock_exhaustiveness_final": 8,
        "dock_md_rep_exhaustiveness": 8,
        "dock_rescore_topn": 50,
    })

    def resolved(self, fast: bool) -> "Defaults":
        d = Defaults()
        if fast:
            for k, v in self.fast.items():
                setattr(d, k, v)
        return d


def resolve_defaults(options: dict | None = None,
                     *, fast: bool | None = None) -> "Defaults":
    """R9: Defaults with agent/CLI options layered on top.

    Previously `DEFAULTS.resolved(options.get("fast"))` silently DROPPED
    every other option (e.g. a user override of dock_exhaustiveness_final
    never reached the tools). Known dataclass fields are honored (None
    values ignored), unknown keys (name, max_steps, ...) are skipped."""
    opts = options or {}
    if fast is None:
        fast = bool(opts.get("fast"))
    d = Defaults().resolved(fast)
    for k, v in opts.items():
        if v is None or k == "fast" or not hasattr(Defaults, k):
            continue
        try:
            object.__setattr__(d, k, v)
        except Exception:
            setattr(d, k, v)
    return d


DEFAULTS = Defaults()


# --------------------------------------------------------------------------- #
# re-exports (historically imported from .config by some modules)
# --------------------------------------------------------------------------- #
from .utils import download_file, fetch_pdb, n_cores  # noqa: E402,F401
