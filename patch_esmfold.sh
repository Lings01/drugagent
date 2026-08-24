#!/bin/bash
# fair-esm site-packages py3.12 dataclass patches (run after env is ready)
set -euo pipefail
ROOT=/home/data/lrs/drug/drugagent
SP=$ROOT/env/lib/python3.12/site-packages/esm/esmfold/v1
python3.12 - <<'EOF'
from pathlib import Path
sp = Path("/home/data/lrs/drug/drugagent/env/lib/python3.12/site-packages/esm/esmfold/v1")
# patch 1: esmfold.py — mutable dataclass default
p = sp / "esmfold.py"
t = p.read_text()
if "trunk: T.Any = FoldingTrunkConfig()" in t:
    t = t.replace("trunk: T.Any = FoldingTrunkConfig()",
                  "trunk: T.Any = field(default_factory=FoldingTrunkConfig)")
    if "from dataclasses import dataclass" in t and ", field" not in \
            [l for l in t.splitlines() if l.startswith("from dataclasses import")][0]:
        t = t.replace("from dataclasses import dataclass",
                      "from dataclasses import dataclass, field", 1)
    p.write_text(t)
    print("patched esmfold.py")
elif "field(default_factory=FoldingTrunkConfig)" in t:
    print("esmfold.py already patched")
# patch 2: trunk.py — mutable dataclass default
p = sp / "trunk.py"
t = p.read_text()
if "structure_module: StructureModuleConfig = StructureModuleConfig()" in t:
    t = t.replace(
        "structure_module: StructureModuleConfig = StructureModuleConfig()",
        "structure_module: StructureModuleConfig = field(default_factory=StructureModuleConfig)")
    lines = [l for l in t.splitlines() if l.startswith("from dataclasses import")]
    if lines and "field" not in lines[0]:
        t = t.replace(lines[0], lines[0] + ", field", 1)
    elif not lines:
        t = t.replace("import torch", "import torch\nfrom dataclasses import field", 1)
    p.write_text(t)
    print("patched trunk.py")
elif "field(default_factory=StructureModuleConfig)" in t:
    print("trunk.py already patched")
EOF
echo ESMFOLD_PATCHED
