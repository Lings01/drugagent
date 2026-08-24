"""ESMFold structure-prediction driver (fair-esm 2.0, ESM2-3B).

Uses the vendored `openfold` package; weights cached under data/weights.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
from loguru import logger

_model = None
_model_device = None


def _vendor_dir() -> Path:
    cands = [Path(__file__).resolve().parent.parent / "vendor",
             Path(__file__).resolve().parent.parent.parent / "drugagent" / "vendor"]
    for c in cands:
        if c.is_dir():
            return c
    return cands[0]


def _ensure_path():
    vd = _vendor_dir()
    if str(vd) not in sys.path:
        sys.path.insert(0, str(vd))


def _ckpt_path(name: str) -> Path:
    from ..config import WEIGHTS
    return WEIGHTS / "torch_cache" / "hub" / "checkpoints" / name


def get_model(device: str = "cpu"):
    """Load ESMFold; uses local weights when present, else downloads."""
    global _model, _model_device
    if _model is not None and _model_device == device:
        return _model
    _ensure_path()
    from ..config import WEIGHTS
    os.environ.setdefault("TORCH_HOME", str(WEIGHTS / "torch_cache"))
    os.environ.setdefault("HF_HOME", str(WEIGHTS / "hf_cache"))
    logger.info("Loading ESMFold (first call downloads ~4 GB of weights)...")
    import torch
    from esm.esmfold.v1.esmfold import ESMFold

    pt = _ckpt_path("esmfold_3B_v1.pt")
    if not pt.exists():
        from esm.esmfold.v1.pretrained import esmfold_v1
        m = esmfold_v1()
    else:
        data = torch.load(str(pt), map_location="cpu", weights_only=False)
        m = ESMFold(esmfold_config=data["cfg"]["model"])
        # Start from the checkpoint weights (NOT the freshly initialised
        # model, which would leave ~99.9% of parameters random). The
        # released checkpoint stores IPA point weights flat
        # (linear_q_points.weight) while the vendored openfold nests them
        # under a PointProjection (.linear); remap those 4 keys.
        # ESM2 itself (esm.*) is loaded separately by esm.pretrained in
        # ESMFold.__init__, so it is intentionally absent from the ckpt.
        sd = dict(data["model"])
        model_keys = set(m.state_dict().keys())
        remapped = 0
        for model_key, ckpt_key in (
            ("trunk.structure_module.ipa.linear_q_points.linear.weight",
             "trunk.structure_module.ipa.linear_q_points.weight"),
            ("trunk.structure_module.ipa.linear_q_points.linear.bias",
             "trunk.structure_module.ipa.linear_q_points.bias"),
            ("trunk.structure_module.ipa.linear_kv_points.linear.weight",
             "trunk.structure_module.ipa.linear_kv_points.weight"),
            ("trunk.structure_module.ipa.linear_kv_points.linear.bias",
             "trunk.structure_module.ipa.linear_kv_points.bias"),
        ):
            if model_key in model_keys and ckpt_key in sd:
                sd[model_key] = sd.pop(ckpt_key)
                remapped += 1
        missing, unexpected = m.load_state_dict(sd, strict=False)
        logger.info(f"Loaded {len(sd)} checkpoint keys "
                    f"({remapped} flat->nested remapped); "
                    f"missing={len(missing)} unexpected={len(unexpected)}")
    if device == "cpu":
        # ESM is cast to half precision in __init__; fp16 is slow on CPU
        m.esm.float()
    m.to(device)
    m.eval()
    _model, _model_device = m, device
    return m


def predict(
    sequences: str | list[str],
    *,
    num_recycles: int = 3,
    device: str = "cpu",
) -> dict:
    """Predict structure(s).

    sequences: one sequence, or a list of chain sequences (multimer, chains are
    joined by a 25-G linker internally).

    Returns dict with:
      pdb (str), plddt (np [L]), pae (np [L,L] | None), ptm (float | None),
      chain_index (np [L]), aatype (np [L])
    """
    import torch

    _ensure_path()
    from esm.esmfold.v1.misc import output_to_pdb

    if isinstance(sequences, str):
        sequences = [sequences]
    model = get_model(device)
    # ESMFold v1: a list is a BATCH of monomers; a multimer must be a single
    # string with chains separated by ":" (25-G linker inserted internally).
    multimer = ":".join(sequences)
    logger.info(f"ESMFold predicting {len(sequences)} chain(s), recycles={num_recycles}")
    with torch.no_grad():
        out = model.infer(multimer, num_recycles=num_recycles)

    pdb_strings = output_to_pdb(out)
    pdb = pdb_strings[0] if pdb_strings else ""
    plddt = out["plddt"][0].detach().cpu().numpy()  # [L, n_atoms]
    plddt_res = plddt.mean(axis=1)
    # multimer inputs contain 25-residue poly-G linkers with no atoms;
    # res_present marks the residues that actually appear in the PDB
    res_present = (out["atom37_atom_exists"][0].sum(dim=1) > 0)
    res_present = res_present.detach().cpu().numpy()
    aatype = out["aatype"][0].detach().cpu().numpy()
    chain_index = (out["chain_index"].detach().cpu().numpy()
                   if "chain_index" in out else np.zeros(len(aatype), dtype=int))
    pae = (out["predicted_aligned_error"][0].detach().cpu().numpy()
           if "predicted_aligned_error" in out else None)
    ptm = float(out["ptm"][0]) if "ptm" in out else None
    return {
        "pdb": pdb,
        "plddt": plddt_res,
        "mean_plddt": float(np.mean(plddt_res)),
        "min_plddt": float(np.min(plddt_res)),
        "pae": pae,
        "ptm": ptm,
        "aatype": aatype,
        "chain_index": chain_index,
        "res_present": res_present,
    }


def interface_metrics(pdb_path: str | Path, res_plddt: np.ndarray,
                      cutoff_ang: float = 8.0) -> dict:
    """Interface pLDDT/PAE for a 2-chain (or N-chain) PDB.

    Interface residue = has an atom of another chain within cutoff_ang.
    res_plddt must align 1:1 with the residues in the PDB (in file order).
    """
    from drugagent.utils import pdb_chains

    import numpy as np

    atoms: list[tuple[str, str, np.ndarray]] = []  # (chain, resi, xyz)
    res_order: list[tuple[str, int]] = []
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                resi = line[22:26].strip()
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                atoms.append((chain, resi, xyz))
                if not res_order or res_order[-1] != (chain, resi):
                    res_order.append((chain, resi))
    if len(res_order) != len(res_plddt):
        logger.warning(f"interface_metrics: plddt len {len(res_plddt)} != residues {len(res_order)}")
    n = min(len(res_order), len(res_plddt))
    # map each residue to its mean atom position
    res_pos: dict[tuple[str, int], list[np.ndarray]] = {}
    for (chain, resi, xyz) in atoms:
        key = (chain, int(float(resi)))
        res_pos.setdefault(key, []).append(xyz)
    res_cent: dict[tuple[str, int], np.ndarray] = {
        k: np.mean(v, axis=0) for k, v in res_pos.items()
    }
    keys = list(res_cent.keys())[:n]
    centers = np.array([res_cent[k] for k in keys])
    chains_of = np.array([k[0] for k in keys])
    pld = np.array(res_plddt[:n])

    interf: list[bool] = []
    for i in range(len(keys)):
        is_if = False
        for j in range(len(keys)):
            if keys[i][0] != keys[j][0]:
                if np.linalg.norm(centers[i] - centers[j]) <= cutoff_ang:
                    is_if = True
                    break
        interf.append(is_if)
    interf = np.array(interf)
    out = {
        "n_interface_residues": int(interf.sum()),
        "interface_plddt_mean": float(plddt_mean(pld[interf])) if interf.sum() else None,
        "interface_plddt_min": float(pld[interf].min()) if interf.sum() else None,
    }
    return out


def plddt_mean(x: np.ndarray) -> float:
    return float(np.mean(x)) if len(x) else 0.0


def write_pdb(out: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(out["pdb"])
    return path
