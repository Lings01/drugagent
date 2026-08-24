"""Shared utilities: subprocess, IO, downloads, PDB helpers, parallelism."""
from __future__ import annotations

import json
import os
import re
import subprocess
import time
from pathlib import Path
from typing import Any, Callable, Iterable

import requests
from joblib import Parallel, delayed
from loguru import logger

# silence RDKit warnings
os.environ.setdefault("RDKIT_WARNINGS", "0")


# --------------------------------------------------------------------------- #
# subprocess
# --------------------------------------------------------------------------- #
def run_cmd(
    cmd: list[str] | str,
    *,
    cwd: str | Path | None = None,
    env: dict[str, str] | None = None,
    log_file: str | Path | None = None,
    check: bool = True,
    timeout: float | None = None,
    capture: bool = True,
    stdin: str | None = None,
) -> subprocess.CompletedProcess:
    """Run a command with logging. Returns CompletedProcess."""
    cmd_s = cmd if isinstance(cmd, str) else " ".join(str(c) for c in cmd)
    logger.info(f"$ {cmd_s}")
    kwargs: dict[str, Any] = dict(cwd=str(cwd) if cwd else None, env=env, timeout=timeout)
    if stdin is not None:
        kwargs["input"] = stdin.encode()
    fh = None
    if log_file is not None:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        fh = open(log_file, "a")
        fh.write(f"\n$ {cmd_s}\n")
        fh.flush()
        kwargs["stdout"] = fh
        kwargs["stderr"] = subprocess.STDOUT
        capture = False
    else:
        kwargs["stdout"] = subprocess.PIPE
        kwargs["stderr"] = subprocess.STDOUT
    try:
        if isinstance(cmd, str):
            proc = subprocess.run(cmd, shell=True, **kwargs)
        else:
            proc = subprocess.run([str(c) for c in cmd], **kwargs)
    finally:
        if fh is not None:
            fh.close()
    if capture and proc.stdout is not None:
        out = proc.stdout.decode(errors="replace")
        if len(out) > 4000:
            out = out[:2000] + "\n...[truncated]...\n" + out[-2000:]
        if out.strip():
            logger.debug(f"stdout: {out}")
    if check and proc.returncode != 0:
        tail = ""
        if log_file is not None:
            # R10: a failed *logged* command must surface its own output —
            # otherwise the agent sees "Command failed (1): gmx solvate ..."
            # with no idea WHY and burns a read_file step (or the user
            # gives up). Include the log tail (capped) in the error.
            try:
                lines = Path(log_file).read_text(errors="replace").splitlines()
                tail = "\n".join(lines[-40:])
                if len(tail) > 2000:
                    tail = ("…[截断, 完整日志见 " + str(log_file) + "]\n"
                            + tail[-2000:])
                tail = f"\n--- {Path(log_file).name} (tail) ---\n{tail}"
            except OSError:
                pass
        elif capture and proc.stdout is not None:
            raw = proc.stdout.decode(errors="replace")
            tail = ("\n" + raw[-2000:]) if raw.strip() else ""
        raise RuntimeError(
            f"Command failed ({proc.returncode}): {cmd_s}{tail}")
    return proc


# --------------------------------------------------------------------------- #
# IO
# --------------------------------------------------------------------------- #
def jsave(path: str | Path, obj: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=str)


def jload(path: str | Path) -> Any:
    with open(path) as fh:
        return json.load(fh)


def download_file(
    url: str,
    dest: str | Path,
    *,
    retries: int = 3,
    timeout: int = 60,
    expected_size_min: int | None = None,
) -> Path:
    """Download with resume + retries."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(1, retries + 1):
        try:
            headers = {"Range": f"bytes={tmp.stat().st_size}-"} if tmp.exists() else {}
            with requests.get(url, stream=True, headers=headers, timeout=timeout,
                              allow_redirects=True) as r:
                if r.status_code == 416:  # already complete
                    break
                r.raise_for_status()
                mode = "ab" if (tmp.exists() and tmp.stat().st_size > 0 and r.status_code == 206) else "wb"
                if mode == "wb":
                    tmp.write_bytes(b"")
                with open(tmp, mode) as fh:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        if chunk:
                            fh.write(chunk)
            if expected_size_min and tmp.stat().st_size < expected_size_min:
                raise RuntimeError(f"downloaded file too small: {tmp.stat().st_size}")
            tmp.rename(dest)
            logger.info(f"downloaded {url} -> {dest} ({dest.stat().st_size} bytes)")
            return dest
        except Exception as e:  # noqa: BLE001
            logger.warning(f"download attempt {attempt} failed: {e}")
            if attempt == retries:
                raise
            time.sleep(3 * attempt)
    return dest


# --------------------------------------------------------------------------- #
# parallelism
# --------------------------------------------------------------------------- #
def pmap(func: Callable, items: Iterable, *, n_jobs: int = 1, **kwargs):
    """Parallel map with logging-friendly fallback."""
    items = list(items)
    if n_jobs <= 1:
        return [func(x) for x in items]
    return Parallel(n_jobs=n_jobs, backend="loky", verbose=0)(
        delayed(func)(x) for x in items
    )


def n_cores() -> int:
    return os.cpu_count() or 8


# --------------------------------------------------------------------------- #
# PDB helpers
# --------------------------------------------------------------------------- #
def fetch_pdb(pdb_id: str, dest_dir: str | Path | None = None, fmt: str = "pdb") -> Path:
    """Fetch a structure from RCSB (pdb or mmCIF)."""
    dest_dir = Path(dest_dir) if dest_dir else Path(".")
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{pdb_id}.{fmt}"
    if not dest.exists():
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.{fmt}"
        try:
            download_file(url, dest)
        except Exception:  # noqa: BLE001
            # fall back to the local fixture cache (offline mode)
            from .config import ROOT
            cached = ROOT / "data" / "fixtures" / f"{pdb_id.upper()}.{fmt}"
            if cached.is_file():
                import shutil
                shutil.copyfile(cached, dest)
            else:
                raise
    return dest


def pdb_ligands(pdb_path: str | Path) -> list[str]:
    """Return unique HETATM residue names (excluding common waters/ions)."""
    skip = {"HOH", "WAT", "H2O", "K", "NA", "CL", "MG", "ZN", "CA", "MN", "FE",
            "NI", "CU", "SO4", "PO4", "GOL", "EDO", "ACT", "PEG", "DOD", "CRO", "UNX"}
    names = set()
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("HETATM"):
                name = line[17:20].strip()
                if name and name.upper() not in skip:
                    names.add(name)
    return sorted(names)


def pdb_chains(pdb_path: str | Path) -> dict[str, int]:
    """chain -> number of ATOM records."""
    counts: dict[str, int] = {}
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and not line.startswith("HETATM"):
                chain = line[21]
                counts[chain] = counts.get(chain, 0) + 1
    return counts


def centroid_from_pdb(pdb_path: str | Path, resnames: list[str] | None = None,
                      chains: list[str] | None = None) -> tuple[float, float, float]:
    """Geometric centroid of atoms (optionally restricted to residues/chains)."""
    xs = ys = zs = n = 0.0
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = line[17:20].strip()
            chain = line[21]
            if resnames and resname not in resnames:
                continue
            if chains and chain not in chains:
                continue
            xs += float(line[30:38]); ys += float(line[38:46]); zs += float(line[46:54]); n += 1
    if n == 0:
        raise ValueError("no atoms matched centroid query")
    return xs / n, ys / n, zs / n


def pdb_extent(pdb_path: str | Path, center: tuple[float, float, float],
               resnames: list[str] | None = None) -> float:
    """Max distance (A) of matching atoms from center."""
    cx, cy, cz = center
    max_d = 0.0
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if resnames and line[17:20].strip() not in resnames:
                continue
            d = ((float(line[30:38]) - cx) ** 2 + (float(line[38:46]) - cy) ** 2
                 + (float(line[46:54]) - cz) ** 2) ** 0.5
            max_d = max(max_d, d)
    return max_d


def is_pdb_id(s: str) -> bool:
    return bool(re.fullmatch(r"[A-Za-z0-9]{4}", s.strip()))


def ensure_parent(path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p
