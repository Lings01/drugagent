#!/bin/bash
# small-molecule library: DTP primary, PDBBind fallback
set -uo pipefail
ROOT=/home/data/lrs/drug/drugagent
LIB=$ROOT/data/libraries
mkdir -p $LIB
echo "trying DTP..."
if ! { [ -f $LIB/dtp.sdf ] && [ $(stat -c%s $LIB/dtp.sdf) -gt 100000000 ]; }; then
  for u in "http://www.dtpbase.org/download/All_Druglike_Compounds.sdf.gz" \
           "https://www.dtpbase.org/download/All_Druglike_Compounds.sdf.gz"; do
    curl -sL --retry 1 --max-time 1800 -o $LIB/dtp.sdf.gz "$u" && \
      gzip -t $LIB/dtp.sdf.gz 2>/dev/null && \
      gunzip -f $LIB/dtp.sdf.gz && echo DTP_DONE && break
  done
fi
if [ -f $LIB/dtp.sdf ] && [ $(stat -c%s $LIB/dtp.sdf) -gt 100000000 ]; then
  echo LIB_DONE
  exit 0
fi
echo "DTP failed, trying PDBBind..."
rm -f $LIB/dtp.sdf.gz
for u in "https://www.pdbbind.org.cn/pdbbind/DownloadData/PDBBind_v20201216.tar.gz" \
         "http://www.pdbbind.org.cn/pdbbind/DownloadData/PDBBind_v20201216.tar.gz"; do
  if curl -sL --retry 1 --max-time 1800 -o $LIB/pdbbind.tar.gz "$u" && tar -tf $LIB/pdbbind.tar.gz >/dev/null 2>&1; then
    python3.12 - <<'EOF'
import tarfile
from pathlib import Path
lib = Path("/home/data/lrs/drug/drugagent/data/libraries")
with tarfile.open(lib / "pdbbind.tar.gz") as t, open(lib / "pdbbind.sdf", "wb") as fo:
    members = [m for m in t.getmembers() if m.name.endswith(".sdf")]
    for m in members:
        with t.extractfile(m) as src:
            fo.write(src.read())
print("PDBBind sdf built from", "members")
EOF
    echo PDBBIND_DONE
    break
  fi
done
[ -f $LIB/pdbbind.sdf ] || echo LIB_FAILED
echo LIB_STAGE_FINISHED
