"""One-click environment setup (idempotent)."""
from __future__ import annotations

import os
import shutil
import subprocess
import tarfile
import zipfile
from pathlib import Path

from loguru import logger

from .config import (ENV_DIR, LIBRARIES, ROOT, TOOLS, WEIGHTS)
from .utils import download_file, fetch_pdb

RF_REPO_SRC = Path("/home/data/lrs/RFdiffusion-main")
RF_WEIGHTS = {
    "Base_ckpt.pt": "http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt",
    "Complex_base_ckpt.pt": "http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt",
    "InpaintSeq_Fold_ckpt.pt": "http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt",
}
VHH_SCAFFOLD_PDB = "1EWN"


def _done_marker(name: str) -> Path:
    return ROOT / ".setup_done" / name


def _mark(name: str) -> None:
    m = _done_marker(name)
    m.parent.mkdir(parents=True, exist_ok=True)
    m.write_text("ok")


def _done(name: str) -> bool:
    return _done_marker(name).exists()


# --------------------------------------------------------------------------- #
# python env
# --------------------------------------------------------------------------- #
def setup_python_env() -> None:
    if (ENV_DIR / "bin" / "python").is_file():
        logger.info("python env already present")
        return
    conda = "/usr/local/anaconda3/envs/rna/bin/conda"
    env = dict(os.environ)
    env.update({"CONDA_PKGS_DIRS": str(ROOT / "conda_pkgs"),
                "PIP_USER": "0", "PYTHONNOUSERSITE": "1",
                "PIP_CACHE_DIR": str(ROOT / "pip_cache")})
    Path(ROOT / "conda_pkgs").mkdir(parents=True, exist_ok=True)
    subprocess.run([conda, "create", "-p", str(ENV_DIR), "python=3.12", "pip",
                    "-y", "-c", "conda-forge", "--override-channels"],
                   env=env, check=True)
    pip = ENV_DIR / "bin" / "pip"
    pkgs = ["langgraph", "langgraph-checkpoint-sqlite", "langchain-openai",
            "typer", "loguru", "requests", "joblib", "numpy", "pandas",
            "scipy", "scikit-learn", "rdkit", "biopython", "MDAnalysis",
            "plotly", "kaleido", "weasyprint", "mordred", "tqdm", "pytest",
            "httpx", "rich", "omegaconf", "fair-esm", "acpype"]
    subprocess.run([str(pip), "install", "-q", *pkgs], env=env, check=True)
    subprocess.run([str(pip), "install", "-q", "torch",
                    "--index-url", "https://download.pytorch.org/whl/cpu"],
                   env=env, check=True)
    subprocess.run([conda, "install", "-p", str(ENV_DIR), "-c", "conda-forge",
                    "-y", "--override-channels", "openbabel", "sqm"],
                   env=env, check=True)
    _mark("python_env")


# --------------------------------------------------------------------------- #
# tools: vina / gnina / 3dmol
# --------------------------------------------------------------------------- #
def setup_tools() -> None:
    if _done("tools"):
        logger.info("tools already installed")
        return
    # vina: conda-forge
    vina_bin = TOOLS / "vina" / "bin" / "vina"
    if not vina_bin.is_file():
        conda = "/usr/local/anaconda3/envs/rna/bin/conda"
        env = dict(os.environ, CONDA_PKGS_DIRS=str(ROOT / "conda_pkgs"))
        (TOOLS / "vina").mkdir(parents=True, exist_ok=True)
        subprocess.run([conda, "create", "-p", str(TOOLS / "vina"),
                        "-c", "conda-forge", "--override-channels",
                        "autodock-vina", "-y"], env=env,
                       stdout=subprocess.DEVNULL)
        # conda-forge 'vina' package provides `vina`
        for cand in (TOOLS / "vina" / "bin" / "vina",
                     TOOLS / "vina" / "bin" / "autodock-vina"):
            if cand.is_file():
                (TOOLS / "vina" / "bin" / "vina").symlink_to(cand)
    # gnina: CPU build from GitHub releases
    gnina = TOOLS / "gnina" / "gnina"
    if not gnina.is_file():
        url = ("https://github.com/gnina/gnina/releases/download/"
               "v1.3.2/gnina.1.3.2")
        raw = TOOLS / "gnina_raw"
        try:
            download_file(url, raw)
            if raw.suffix == ".zip" or zipfile.is_zipfile(raw):
                with zipfile.ZipFile(raw) as z:
                    z.extractall(TOOLS / "gnina")
            else:
                try:
                    with tarfile.open(raw) as t:
                        t.extractall(TOOLS / "gnina")
                except tarfile.ReadError:
                    (TOOLS / "gnina").mkdir(parents=True, exist_ok=True)
                    shutil.copy(raw, gnina)
            for f in (TOOLS / "gnina").rglob("gnina*"):
                if f.is_file() and "gnina" in f.name and not f.name.endswith(
                        (".txt", ".md")):
                    f.chmod(0o755)
                    if f.name != "gnina":
                        gnina.symlink_to(f) if not gnina.exists() else None
        except Exception as e:  # noqa: BLE001
            logger.warning(f"gnina download failed (optional): {e}")
    # 3Dmol.js
    js = TOOLS / "3Dmol" / "3Dmol-min.js"
    if not js.is_file():
        js.parent.mkdir(parents=True, exist_ok=True)
        try:
            download_file("https://3dmol.org/build/3Dmol-min.js", js,
                          expected_size_min=50_000)
        except Exception as e:  # noqa: BLE001
            logger.warning(f"3Dmol download failed (optional): {e}")
    _mark("tools")


# --------------------------------------------------------------------------- #
# GROMACS
# --------------------------------------------------------------------------- #
def setup_gromacs() -> None:
    gmx = TOOLS / "gromacs" / "bin" / "gmx"
    spc = TOOLS / "gromacs" / "share" / "gromacs" / "top" / "spc216.gro"
    if gmx.is_file() and spc.is_file():
        logger.info("GROMACS 2023.1 already built")
        return
    src = Path("/home/data/lrs/AI_Design_Antibody/gromacs-2023.1")
    if not src.is_dir():
        logger.warning("GROMACS source not found; will fall back to old env")
        return
    bdir = TOOLS / "gromacs-build"
    bdir.mkdir(parents=True, exist_ok=True)
    cmd = ["cmake", "-S", str(src), "-B", str(bdir),
           f"-DCMAKE_INSTALL_PREFIX={TOOLS / 'gromacs'}",
           "-DCMAKE_BUILD_TYPE=Release", "-DGMX_HDF5=OFF", "-DGMX_API=OFF",
           "-DGMX_X11=OFF", "-DGMX_BUILD_MPIRUN=OFF"]
    subprocess.run(cmd, check=True,
                   stdout=open(bdir / "cmake.log", "w"), stderr=subprocess.STDOUT)
    subprocess.run(["cmake", "--build", str(bdir), "-j", "64"], check=True,
                   stdout=open(bdir / "build.log", "a"), stderr=subprocess.STDOUT)
    subprocess.run(["cmake", "--install", str(bdir)], check=True,
                   stdout=open(bdir / "install.log", "a"), stderr=subprocess.STDOUT)
    _mark("gromacs")


# --------------------------------------------------------------------------- #
# RFdiffusion
# --------------------------------------------------------------------------- #
def setup_rfdiffusion() -> None:
    repo = TOOLS / "RFdiffusion"
    if not repo.is_dir():
        if RF_REPO_SRC.is_dir():
            shutil.copytree(RF_REPO_SRC, repo,
                            ignore=shutil.ignore_patterns("*.pyc", "__pycache__",
                                                         "example_outputs"))
        else:
            raise FileNotFoundError(f"RFdiffusion source not found: {RF_REPO_SRC}")
    models = repo / "models"
    models.mkdir(exist_ok=True)
    for name, url in RF_WEIGHTS.items():
        dst = models / name
        if not dst.exists():
            download_file(url, dst, retries=2)
    # VHH scaffold (1EWN) preprocessing
    sdir = TOOLS / "vhh_scaffolds"
    if not (sdir / "vhh1ewn_ss.pt").exists():
        sdir.mkdir(parents=True, exist_ok=True)
        pdb = fetch_pdb(VHH_SCAFFOLD_PDB, sdir)
        (sdir / "vhh1ewn.pdb").write_text(pdb.read_text())
        # single chain A only
        txt = "".join(l for l in pdb.read_text().splitlines(True)
                      if not (l.startswith(("ATOM", "HETATM")) and l[21] != "A"))
        (sdir / "vhh1ewn.pdb").write_text(txt + "END\n")
        py = ENV_DIR / "bin" / "python"
        subprocess.run(
            [str(py), str(repo / "helper_scripts" / "make_secstruc_adj.py"),
             "--input_pdb", str(sdir / "vhh1ewn.pdb"),
             "--out_dir", str(sdir)],
            check=True, stdout=subprocess.DEVNULL, cwd=sdir)
        # rename outputs to vhh1ewn_*
        for f in sdir.glob("1EWN_*.pt"):
            f.rename(sdir / f.name.replace("1EWN_", "vhh1ewn_"))
        (sdir / "scaffold_list.txt").write_text("vhh1ewn\n")
    _mark("rfdiffusion")


# --------------------------------------------------------------------------- #
# libraries
# --------------------------------------------------------------------------- #
def setup_libraries(which: str = "all") -> None:
    names = (["dtp", "chembl", "pdbbind", "vhh"] if which == "all"
             else [which])
    for n in names:
        if n == "dtp":
            _setup_dtp()
        elif n == "vhh":
            _setup_vhh_library()
        elif n == "pdbbind":
            _setup_pdbbind()


def _setup_dtp() -> None:
    dst = LIBRARIES / "dtp.sdf"
    if dst.exists() and dst.stat().st_size > 100_000_000:
        return
    urls = [
        "http://www.dtpbase.org/download/All_Druglike_Compounds.sdf.gz",
        "https://www.dtpbase.org/download/All_Druglike_Compounds.sdf.gz",
    ]
    last = None
    for u in urls:
        try:
            gz = LIBRARIES / "dtp.sdf.gz"
            download_file(u, gz, retries=1, expected_size_min=100_000_000)
            with tarfile.open(gz, "r:gz") as t, open(dst, "wb") as fo:
                for m in t.getmembers():
                    if m.isfile():
                        with t.extractfile(m) as src:
                            shutil.copyfileobj(src, fo)
            return
        except Exception as e:  # noqa: BLE001
            last = e
            logger.warning(f"DTP download failed: {e}")
    raise RuntimeError(f"DTP library download failed: {last} "
                       "(fallback: `drugagent setup --libraries pdbbind`)")


def _setup_pdbbind() -> None:
    dst = LIBRARIES / "pdbbind.sdf"
    if dst.exists():
        return
    # PDBBind provides tar.gz of SDFs; try known URL, else build small library
    urls = [
        "https://www.pdbbind.org.cn/pdbbind/DownloadData/PDBBind_v20201216.tar.gz",
        "http://www.pdbbind.org.cn/pdbbind/DownloadData/PDBBind_v20201216.tar.gz",
    ]
    for u in urls:
        try:
            raw = LIBRARIES / "pdbbind.tar.gz"
            download_file(u, raw, retries=1)
            with tarfile.open(raw) as t:
                members = [m for m in t.getmembers() if m.name.endswith(".sdf")]
                with open(dst, "wb") as fo:
                    for m in members[:20]:
                        with t.extractfile(m) as src:
                            shutil.copyfileobj(src, fo)
            return
        except Exception as e:  # noqa: BLE001
            logger.warning(f"PDBBind download failed: {e}")
    raise RuntimeError("PDBBind download failed")


def _setup_vhh_library() -> None:
    dst = LIBRARIES / "vhh_library.fasta"
    if dst.exists():
        return
    from .modules.vhh import generate_vhh_library, save_library
    seqs = generate_vhh_library(20_000)
    save_library(seqs, dst)
    logger.info(f"VHH library: {len(seqs)} sequences -> {dst}")


# --------------------------------------------------------------------------- #
# weight prefetch (optional, big)
# --------------------------------------------------------------------------- #
def setup_weights() -> None:
    """Prefetch ESMFold weights (ESM2-3B + esmfold_3B_v1)."""
    import urllib.request
    esmfold_pt = WEIGHTS / "esmfold_3B_v1.pt"
    if not esmfold_pt.exists():
        url = "https://dl.fbaipublicfiles.com/fair-esm/models/esmfold_3B_v1.pt"
        download_file(url, esmfold_pt, retries=2)
    _mark("weights")


def run_setup(libraries: str = "all", gromacs: bool = True,
              tools: bool = True, rfdiffusion: bool = True,
              weights: bool = False) -> None:
    logger.info(f"=== DrugAgent setup (root={ROOT}) ===")
    setup_python_env()
    if tools:
        setup_tools()
    if gromacs:
        setup_gromacs()
    if rfdiffusion:
        setup_rfdiffusion()
    if libraries:
        setup_libraries(libraries)
    if weights:
        setup_weights()
    logger.info("=== setup complete ===")
