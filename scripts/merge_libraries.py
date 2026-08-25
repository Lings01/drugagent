#!/usr/bin/env python3
"""Build the master screening library `nci_npatlas.sdf` as the UNION of
NCI-Open (NCI/DTP open chemical repository) and NPAtlas (natural products),
deduplicated by canonical InChIKey.

Provenance
----------
- data/libraries/nci_open_2012-05-01.sdf
    NCI/DTP DIS export (Release: December 2010, exported 2012-05-01),
    265,242 records, 3D coordinates, NSC ids. Source file:
    https://dctd.cancer.gov/data-tools-biospecimens/data
- data/libraries/npatlas_2024_09.sdf
    NPAtlas (Natural Product Atlas) 2024-09 download, 36,454 records,
    2D coordinates, names + InChI/InChIKey. Source: https://www.npatlas.com

Semantics
---------
- Union: a molecule present in both files appears ONCE (first-seen wins;
  NCI is processed first). In-file duplicates are removed the same way.
- Original record text is preserved verbatim (properties, coordinates);
  CRLF is normalized to LF.
- A record that RDKit cannot parse is kept (conservative: it may still
  parse through the pipeline's tolerant block parser) but is NOT counted
  as a dedup key (its InChIKey is None).

Usage:  env/bin/python scripts/merge_libraries.py
Output: data/libraries/nci_npatlas.sdf
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
LIBS = ROOT / "data" / "libraries"
SRC_NCI = LIBS / "nci_open_2012-05-01.sdf"
SRC_NPA = LIBS / "npatlas_2024_09.sdf"
OUT = LIBS / "nci_npatlas.sdf"


def iter_blocks(path: Path):
    """Yield raw record text (without trailing $$$$), CRLF normalized."""
    buf: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.replace("\r\n", "\n").replace("\r", "\n")
            if line == "$$$$\n" or line.rstrip("\n") == "$$$$":
                if buf:
                    yield "".join(buf)
                    buf = []
                continue
            buf.append(line)
    if buf:
        yield "".join(buf)


def inchikey_of(block: str) -> str | None:
    """InChIKey of the molblock inside an SDF record. RDKit's
    MolFromMolBlock wants the canonical 3-line pre-header (title,
    comment, blank) before the V2000/V3000 counts line, so try the
    counts line and the few lines before it."""
    lines = block.split("\n")
    for i, line in enumerate(lines):
        if "V2000" in line or "V3000" in line:
            for start in (i - 3, i - 2, i - 1, i):
                if start < 0:
                    continue
                try:
                    m = Chem.MolFromMolBlock(
                        "\n".join(lines[start:]), removeHs=True, sanitize=True)
                except Exception:
                    m = None
                if m is not None:
                    try:
                        return Chem.MolToInchiKey(m)
                    except Exception:
                        return None
            return None
    return None


def main() -> int:
    t0 = time.time()
    for p in (SRC_NCI, SRC_NPA):
        if not p.exists():
            print(f"missing source: {p}", file=sys.stderr)
            return 1

    # Pass A: NPAtlas — small; load its InChIKey set + blocks in memory.
    npa_blocks: list[str] = []
    npa_keys: list[str | None] = []
    t = time.time()
    for blk in iter_blocks(SRC_NPA):
        npa_keys.append(inchikey_of(blk))
        npa_blocks.append(blk)
    npa_set = {k for k in npa_keys if k}
    print(f"[npatlas] {len(npa_blocks)} records, {len(npa_set)} unique keys "
          f"({time.time() - t:.0f}s)", flush=True)

    # Pass B: NCI — stream; write uniques directly to OUT.
    seen: set[str] = set(npa_set)  # NCI record colliding with NPAtlas -> skip
    n_total = n_uniq = n_dup = n_nopa = 0
    t = time.time()
    with open(OUT, "w", encoding="utf-8") as out:
        for blk in iter_blocks(SRC_NCI):
            n_total += 1
            k = inchikey_of(blk)
            if k is None:
                n_nopa += 1
            elif k in seen:
                n_dup += 1
                continue
            else:
                seen.add(k)
            n_uniq += 1
            out.write(blk if blk.endswith("\n") else blk + "\n")
            out.write("$$$$\n")
            if n_total % 20000 == 0:
                print(f"[nci] {n_total} processed, {n_uniq} kept, "
                      f"{n_dup} dup ({time.time() - t:.0f}s)", flush=True)
        # Then NPAtlas uniques (its own in-file dups removed too).
        seen_npa: set[str] = set()
        n_npa_written = 0
        for blk, k in zip(npa_blocks, npa_keys):
            if k is None:
                pass
            elif k in seen_npa:
                continue
            else:
                seen_npa.add(k)
            n_npa_written += 1
            out.write(blk if blk.endswith("\n") else blk + "\n")
            out.write("$$$$\n")

    print(f"[nci]     total={n_total} kept={n_uniq} dup={n_dup} "
          f"unparsed={n_nopa} ({time.time() - t:.0f}s)")
    print(f"[npatlas] kept={n_npa_written} (in-file dups removed)")
    total_out = n_uniq + n_npa_written
    print(f"[merge]  output={OUT} records={total_out} "
          f"({OUT.stat().st_size / 1e9:.2f} GB) "
          f"total={time.time() - t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
