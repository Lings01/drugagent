import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from drugagent.utils import download_file  # noqa: E402


@pytest.fixture(scope="session")
def hivp_pdb(tmp_path_factory) -> Path:
    """HIV-1 protease + saquinavir (A77), PDB 1HVI (local cache preferred)."""
    from drugagent.config import ROOT
    cached = ROOT / "data" / "fixtures" / "1HVI.pdb"
    if cached.is_file():
        return cached
    d = tmp_path_factory.mktemp("pdb")
    p = d / "1HVI.pdb"
    if not p.exists():
        download_file("https://files.rcsb.org/download/1HVI.pdb", p)
    return p


@pytest.fixture(scope="session")
def dna_pdb() -> Path:
    """B-DNA dodecamer duplex, PDB 1BNA."""
    from drugagent.config import ROOT
    cached = ROOT / "data" / "fixtures" / "1BNA.pdb"
    assert cached.is_file(), "1BNA fixture missing: data/fixtures/1BNA.pdb"
    return cached


@pytest.fixture(scope="session")
def has_net() -> bool:
    import requests
    try:
        r = requests.get("https://files.rcsb.org", timeout=8)
        return r.status_code in (200, 301, 302)
    except Exception:  # noqa: BLE001
        return False


@pytest.fixture
def tiny_pdb(tmp_path: Path) -> Path:
    """A synthetic 40-residue alpha helix (CA + N/C/O) for offline tests."""
    import math
    lines = []
    pqr = 3.6
    rise = 1.5
    for i in range(40):
        x = 2.3 * math.cos(2 * math.pi * i / 3.6)
        y = 2.3 * math.sin(2 * math.pi * i / 3.6)
        z = i * rise
        resi = i + 1
        lines.append(
            f"ATOM  {4*i+1:5d}  CA  ALA A{resi:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  ")
        lines.append(
            f"ATOM  {4*i+2:5d}  N   ALA A{resi:4d}    "
            f"{x+0.5:8.3f}{y:8.3f}{z+0.4:8.3f}  1.00  0.00           N  ")
        lines.append(
            f"ATOM  {4*i+3:5d}  C   ALA A{resi:4d}    "
            f"{x-0.5:8.3f}{y:8.3f}{z-0.4:8.3f}  1.00  0.00           C  ")
    lines.append("END\n")
    p = tmp_path / "tiny.pdb"
    p.write_text("".join(lines))
    return p


@pytest.fixture
def small_sdf(tmp_path: Path) -> Path:
    """Tiny SDF library (5 small molecules)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    smis = ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O",  # benzene, aspirin
            "CN1C=NC2=C1C(=O)NC(=O)N2C",  # caffeine
            "c1ccc2c(c1)cc(cc2)-c3ccccc3",  # biphenyl-ish
            "CCN(CC)CC"]  # TFE
    blocks = []
    for i, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m, randomSeed=42)
        mb = Chem.MolToMolBlock(m)
        blocks.append(mb + f"$$$$\n")
    p = tmp_path / "lib.sdf"
    p.write_text("".join(blocks))
    return p
