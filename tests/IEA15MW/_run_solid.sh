#!/bin/bash
# FEM solid participant — parameterized run script for yaw sweep
# Usage: _run_solid.sh <yaml_path>
# Must be called from the desired working directory (preCICE sockets and
# output.folder paths will land relative to CWD).
#
# Environment variables (optional overrides):
#   OMP_NUM_THREADS   default: 4  (reduced from 8 to allow 5 parallel cases)
#   RUST_LOG          default: info,aeroelast_solvers::petsc::fsi::profiling=debug

set -euo pipefail

YAML="$1"

module purge
module load glu gcc/14.2.0_sequana

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_btl_openib_allow_ib=0
ulimit -n 65535 >/dev/null 2>&1 || true

export VENV_DIR=$SCRATCH/venv
export PATH=${VENV_DIR}/bin/:$PATH

# GCC 14.2 libstdc++ must come FIRST so GLIBCXX_3.4.32 resolves correctly
# (libprecice.so.3 was compiled against GCC 14; the system /lib64 only has 3.4.25)
GCC14_LIB=$(gcc -print-file-name=libstdc++.so.6 | xargs dirname)
export LD_LIBRARY_PATH=${GCC14_LIB}:${VENV_DIR}/lib/:${VENV_DIR}/lib64/:/usr/lib64/psm2-compat:/scratch/app/openmpi/4.1.4_gnu/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:${VENV_DIR}/lib/pkgconfig:${VENV_DIR}/lib64/pkgconfig

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OPENBLAS_NUM_THREADS=${OMP_NUM_THREADS:-4}
export RUST_LOG=${RUST_LOG:-info,aeroelast_solvers::petsc::fsi::profiling=debug}
export PETSC_OPTIONS="-pc_factor_mat_solver_type mumps"

aeroelast "$YAML"
