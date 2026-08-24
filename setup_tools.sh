#!/bin/bash
# Tools: conda-pkgs seed, vina, gnina, 3Dmol, RFdiffusion+weights, VHH scaffold
set -x
ROOT=/home/data/lrs/drug/drugagent
CONDA=/usr/local/anaconda3/envs/rna/bin/conda
export CONDA_PKGS_DIRS=$ROOT/conda_pkgs
mkdir -p $CONDA_PKGS_DIRS $ROOT/logs
# 1) seed conda package cache (saves ~17G of network)
cp -rn /home/data/lrs/.conda/pkgs/. $CONDA_PKGS_DIRS/ 2>/dev/null || true
echo PKGS_SEEDED
# 2) vina (separate conda env)
rm -rf $ROOT/data/tools/vina
$CONDA create -p $ROOT/data/tools/vina -c conda-forge --override-channels vina -y -q > $ROOT/logs/vina.log 2>&1 && echo VINA_DONE
# 3) gnina v1.3.1 ELF binary
mkdir -p $ROOT/data/tools/gnina
curl -sL --retry 2 --max-time 900 "https://github.com/gnina/gnina/releases/download/v1.3.1/gnina1.3.1" \
  -o $ROOT/data/tools/gnina/gnina && chmod +x $ROOT/data/tools/gnina/gnina && echo GNINA_DONE
# 4) 3Dmol.js
mkdir -p $ROOT/data/tools/3Dmol
curl -sL --retry 2 "https://3dmol.org/build/3Dmol-min.js" -o $ROOT/data/tools/3Dmol/3Dmol-min.js && echo 3DMOL_DONE
# 5) RFdiffusion repo + weights
cp -rn /home/data/lrs/RFdiffusion-main $ROOT/data/tools/RFdiffusion 2>/dev/null || cp -r /home/data/lrs/RFdiffusion-main $ROOT/data/tools/RFdiffusion
echo RF_COPY_DONE
mkdir -p $ROOT/data/tools/RFdiffusion/models
cd $ROOT/data/tools/RFdiffusion/models
curl -sL --retry 2 -o Base_ckpt.pt "http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt" && echo RF_W1_DONE
curl -sL --retry 2 -o Complex_base_ckpt.pt "http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt" && echo RF_W2_DONE
curl -sL --retry 2 -o InpaintSeq_Fold_ckpt.pt "http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt" && echo RF_W3_DONE
cd $ROOT
# 6) VHH scaffold 1EWN (chain A only) + secstruc adjacency
mkdir -p $ROOT/data/tools/vhh_scaffolds
curl -sL --retry 2 "https://files.rcsb.org/download/1EWN.pdb" -o $ROOT/data/tools/vhh_scaffolds/1EWN.pdb && echo 1EWN_DONE
python3.12 - <<'EOF'
from pathlib import Path
p = Path('/home/data/lrs/drug/drugagent/data/tools/vhh_scaffolds/1EWN.pdb')
txt = "".join(l for l in p.read_text().splitlines(True)
              if not (l.startswith(("ATOM", "HETATM")) and l[21] != "A"))
Path('/home/data/lrs/drug/drugagent/data/tools/vhh_scaffolds/vhh1ewn.pdb').write_text(txt + "END\n")
print("vhh1ewn.pdb written")
EOF
# RF scripts/sitecustomize shim (np.long + nvtx noop)
cat > $ROOT/data/tools/RFdiffusion/scripts/sitecustomize.py <<'EOF'
# DrugAgent shims for running RFdiffusion in the shared venv:
# - numpy 1.24+ removed np.long/np.ulong (old scipy lazy-import paths)
# - torch CPU build lacks the _nvtx C binding; se3_transformer calls
#   torch.cuda.nvtx ranges -> replace with no-ops
import sys

try:
    import numpy as _np
    if not hasattr(_np, "long"):
        _np.long = int
        _np.ulong = int
except Exception:
    pass

try:
    import torch.cuda.nvtx as _nv

    class _NoopNVTX:
        @staticmethod
        def rangePushA(*a, **k):
            return 0

        @staticmethod
        def rangePop(*a, **k):
            return None

        @staticmethod
        def markA(*a, **k):
            return None

    _nv._nvtx = _NoopNVTX()
except Exception:
    pass
EOF
# se3_transformer nvtx import shim
cd $ROOT/data/tools/RFdiffusion/env/SE3Transformer && python3.12 - <<'EOF'
import pathlib
shim = '''try:
    from torch.cuda.nvtx import range as nvtx_range
    with nvtx_range("_probe_"):
        pass
except Exception:  # torch CPU build lacks _nvtx C binding
    from contextlib import contextmanager

    @contextmanager
    def nvtx_range(*args, **kwargs):
        yield
'''
n = 0
for p in pathlib.Path("se3_transformer").rglob("*.py"):
    t = p.read_text()
    if "from torch.cuda.nvtx import range as nvtx_range" in t:
        t = t.replace("from torch.cuda.nvtx import range as nvtx_range", shim.rstrip())
        p.write_text(t)
        n += 1
print("patched se3 files:", n)
EOF
find . -name "*.pyc" -delete
echo TOOLS_DONE
