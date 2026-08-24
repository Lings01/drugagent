#!/bin/bash
set -euo pipefail
ROOT=/home/data/lrs/drug/drugagent
mkdir -p $ROOT/logs $ROOT/data/tools/gromacs-build
export CC=/usr/bin/gcc CXX=/usr/bin/g++
cmake -S /home/data/lrs/AI_Design_Antibody/gromacs-2023.1 -B $ROOT/data/tools/gromacs-build \
  -DCMAKE_INSTALL_PREFIX=$ROOT/data/tools/gromacs \
  -DCMAKE_BUILD_TYPE=Release -DGMX_HDF5=OFF -DGMX_API=OFF -DGMX_X11=OFF \
  -DGMX_BUILD_MPIRUN=OFF -DGMX_DEFAULT_MANBUILD=OFF > $ROOT/logs/gmx_build.log 2>&1
cmake --build $ROOT/data/tools/gromacs-build -j 64 >> $ROOT/logs/gmx_build.log 2>&1
cmake --install $ROOT/data/tools/gromacs-build >> $ROOT/logs/gmx_build.log 2>&1
echo GMX_DONE
