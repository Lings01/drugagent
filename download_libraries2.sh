#!/bin/bash
# Library v2: local ChEMBL 35 sqlite -> SDF (primary, offline), DTP retry (bonus)
set -uo pipefail
ROOT=/home/data/lrs/drug/drugagent
LIB=$ROOT/data/libraries
WORK=$LIB/chembl35
mkdir -p $WORK
SRC=/home/data/lrs/chembl_35_sqlite.tar.gz
DB=$WORK/chembl_35.db

# 1) extract sqlite if needed
if [ ! -f "$DB" ]; then
  tar -xzf "$SRC" -C "$WORK" 2>/dev/null
  found=$(find "$WORK" -maxdepth 2 -name '*.db' -o -maxdepth 2 -name '*.sqlite' 2>/dev/null | head -1)
  [ -n "$found" ] && mv -f "$found" "$DB"
fi
[ -f "$DB" ] || { echo "SQLITE_NOT_FOUND"; exit 1; }
echo SQLITE_READY

# 2) SDF export: molfile + canonical_smiles + id  (ChEMBL 35 ~2M mols)
python3.12 - <<'EOF'
import sqlite3, sys
from pathlib import Path
root = Path("/home/data/lrs/drug/drugagent/data/libraries")
db = root / "chembl35/chembl_35.db"
out = root / "chembl35.sdf"
con = sqlite3.connect(str(db))
cur = con.execute("SELECT molecule_id, canonical_smiles, molfile FROM molecule WHERE molfile IS NOT NULL")
n = 0
with open(out, "w") as fo:
    for mid, smi, molfile in cur:
        if not molfile or "M  END" not in molfile:
            continue
        header = f"{mid}\n  DrugAgent\n\n"
        # keep molfile but ensure SMILES prop
        fo.write(header + molfile.rstrip() + "\n  >  <canonical_smiles>\n" + (smi or "") + "\n\n$$$$\n")
        n += 1
print("SDF molecules:", n, "->", out)
EOF
echo CHEMBL_SDF_DONE
ls -la $LIB/chembl35.sdf

# 3) bonus: DTP retry in background of this job (long timeout)
for u in "http://www.dtpbase.org/download/All_Druglike_Compounds.sdf.gz" \
         "https://www.dtpbase.org/download/All_Druglike_Compounds.sdf.gz"; do
  if [ ! -f $LIB/dtp.sdf ] || [ $(stat -c%s $LIB/dtp.sdf 2>/dev/null || echo 0) -lt 100000000 ]; then
    curl -sL --retry 2 --retry-delay 10 --max-time 2400 -o $LIB/dtp.sdf.gz "$u" && \
      gzip -t $LIB/dtp.sdf.gz 2>/dev/null && gunzip -f $LIB/dtp.sdf.gz && echo DTP_DONE
  fi
done
echo LIB_V2_FINISHED
