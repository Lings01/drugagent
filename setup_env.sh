#!/bin/bash
# DrugAgent python env (final recipe, adapted to persistent workspace root)
set -euo pipefail
ROOT=/home/data/lrs/drug/drugagent
ENV=$ROOT/env
CONDA=/usr/local/anaconda3/envs/rna/bin/conda
export PIP_USER=0 PIP_CACHE_DIR=$ROOT/pip_cache PYTHONNOUSERSITE=1
export CONDA_PKGS_DIRS=$ROOT/conda_pkgs
mkdir -p "$PIP_CACHE_DIR" "$CONDA_PKGS_DIRS"
if [ ! -x "$ENV/bin/pip" ]; then
  $CONDA create -p "$ENV" python=3.12 pip -y -c conda-forge --override-channels \
    > $ROOT/logs/conda_env.log 2>&1 || { tail -5 $ROOT/logs/conda_env.log; exit 1; }
fi
P=$ENV/bin/pip
$P install -q --upgrade pip
$P install -q langgraph langgraph-checkpoint-sqlite langchain-openai typer loguru requests joblib \
  numpy pandas scikit-learn rdkit biopython MDAnalysis plotly kaleido weasyprint mordred \
  tqdm pytest httpx rich omegaconf fair-esm acpype ml-collections dm-tree modelcif einops \
  "numpy<2" "scipy==1.11.4"
echo "PIP_STAGE_DONE"
$P install -q torch==2.2.1 --index-url https://download.pytorch.org/whl/cpu
echo "TORCH_DONE"
$P install -q dgl -f https://data.dgl.ai/wheels/repo.html torchdata==0.7.1
$P install -q "e3nn==0.3.3" pynvml wandb decorator pyrsistent socksio dllogger
echo "EXTRAS_DONE"
$P install -q nvidia-cublas-cu12 nvidia-cuda-runtime-cu12 nvidia-cufft-cu12 nvidia-nvtx-cu12 \
  nvidia-nvtx-cu11 nvidia-cusolver-cu11 nvidia-cudnn-cu11
echo "NVIDIALIBS_DONE"
$P install -q openbabel || $CONDA install -p "$ENV" -c conda-forge -y -q --override-channels openbabel sqm
echo "CONDA_EXTRAS_DONE"
# SE3Transformer (RFdiffusion) editable, no deps
(cd $ROOT/data/tools/RFdiffusion/env/SE3Transformer && $P install -q --no-deps -e .)
echo "SE3_INSTALL_DONE"
$ENV/bin/python -c "import torch; print('torch', torch.__version__)"
$ENV/bin/python -c "import dgl; print('dgl', dgl.__version__)"
$ENV/bin/python -c "from se3_transformer.model import SE3Transformer; print('se3 ok')"
echo "ENV_SETUP_DONE"
