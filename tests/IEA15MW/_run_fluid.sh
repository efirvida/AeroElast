#!/bin/bash
# BEM fluid participant — parameterized run script for yaw sweep
# Usage: _run_fluid.sh <yaml_path> [workdir]
#   workdir  explicit working directory for preCICE sockets and output paths.
#            Must match the solid participant's CWD so both find the same socket.
#            Default: directory of yaml_path (aeroelast default behaviour).

set -euo pipefail

YAML="$1"
WORKDIR_ARG=()
[[ -n "${2:-}" ]] && WORKDIR_ARG=(--workdir "$2")

module purge
module load glu gcc/14.2.0_sequana

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_btl_openib_allow_ib=0
ulimit -n 65535 >/dev/null 2>&1 || true

export VENV_DIR=$SCRATCH/venv
export PATH=${VENV_DIR}/bin/:$PATH

GCC14_LIB=$(gcc -print-file-name=libstdc++.so.6 | xargs dirname)
export LD_LIBRARY_PATH=${GCC14_LIB}:${VENV_DIR}/lib/:${VENV_DIR}/lib64/:/usr/lib64/psm2-compat:/scratch/app/openmpi/4.1.4_gnu/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:${VENV_DIR}/lib/pkgconfig:${VENV_DIR}/lib64/pkgconfig

# BEM is not parallelised — no OMP/PETSC tuning needed
aeroelast "$YAML" "${WORKDIR_ARG[@]}"
