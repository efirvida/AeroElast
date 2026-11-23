#!/usr/bin/env bash
# Build CalculiX 2.20 with 8-byte integers (ccx_i8), without PaStiX.
#
# Why: the 32-bit ccx overflows in mastructure (u_realloc, negative size) on the
# es 0.12 m IEA-15 blade mesh (303 711 quadratic nodes), so the mesh-convergence
# test cannot reach its finest point (odd/tasks/gmsh-official-mesher.md L-S77).
#
# Recipe (CCX README.INSTALL, "8-byte integers"): Fortran with -fdefault-integer-8,
# C with -DLONGLONG, ARPACK in integer*8 (so its BLAS/LAPACK too), SPOOLES stays
# integer*4. The existing venv/bin/ccx is left untouched; the result is
# $VENV/bin/ccx_i8, selected per run with CCX_BIN.
#
# Usage: scripts/build_ccx_i8.sh            (needs network; ~15-30 min)
set -euo pipefail

VENV=${VENV:-/scratch/leahk/eduardo.donestevez/venv}
ROOT=${ROOT:-/scratch/leahk/eduardo.donestevez/build/ccx-2.20-i8}
PREFIX=$ROOT/prefix
JOBS=${JOBS:-8}
LAPACK_VER=3.12.0
ARPACK_TAG=3.9.1

here=$(cd "$(dirname "$0")" && pwd)
# shellcheck disable=SC1091
source "$here/aeroenv.sh" >/dev/null 2>&1 || true
export CC=gcc FC=gfortran
mkdir -p "$ROOT" "$PREFIX/lib"
cd "$ROOT"

FFLAGS_I8="-O2 -fPIC -fdefault-integer-8 -fallow-argument-mismatch"

# 1. Reference BLAS + LAPACK, integer*8.
if [ ! -f "$PREFIX/lib/liblapack.a" ]; then
  [ -f v$LAPACK_VER.tar.gz ] || wget -q https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v$LAPACK_VER.tar.gz
  rm -rf lapack-$LAPACK_VER build-lapack && tar xzf v$LAPACK_VER.tar.gz
  cmake -S lapack-$LAPACK_VER -B build-lapack -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF \
    -DCBLAS=OFF -DLAPACKE=OFF -DCMAKE_Fortran_FLAGS="$FFLAGS_I8" >/dev/null
  cmake --build build-lapack -j "$JOBS" >/dev/null
  cmake --install build-lapack >/dev/null
fi

# 2. arpack-ng, integer*8 (INTERFACE64), static, serial.
if [ ! -f "$PREFIX/lib/libarpack.a" ]; then
  rm -rf arpack-ng build-arpack
  git clone -q --depth 1 --branch $ARPACK_TAG https://github.com/opencollab/arpack-ng.git
  cmake -S arpack-ng -B build-arpack -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS=OFF -DINTERFACE64=ON \
    -DMPI=OFF -DICB=OFF -DEXAMPLES=OFF -DTESTS=OFF \
    -DBLAS_LIBRARIES="$PREFIX/lib/libblas.a" -DLAPACK_LIBRARIES="$PREFIX/lib/liblapack.a" \
    -DCMAKE_Fortran_FLAGS="$FFLAGS_I8" >/dev/null
  cmake --build build-arpack -j "$JOBS" >/dev/null
  cmake --install build-arpack >/dev/null
fi

# 3. CalculiX 2.20, -DLONGLONG, no PaStiX. SPOOLES (int4) from the venv.
[ -f ccx_2.20.src.tar.bz2 ] || wget -q http://www.dhondt.de/ccx_2.20.src.tar.bz2
[ -d CalculiX/ccx_2.20/src ] || tar xjf ccx_2.20.src.tar.bz2
cd CalculiX/ccx_2.20/src
cat > Makefile_i8_nopastix <<EOF
CFLAGS = -O2 -fopenmp -fpermissive -w -I$VENV/include/spooles -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DLONGLONG
FFLAGS = -O2 -fopenmp -w -fdefault-integer-8 -fallow-argument-mismatch
CC=gcc
FC=gfortran
.c.o :
	\$(CC) \$(CFLAGS) -c \$<
.f.o :
	\$(FC) \$(FFLAGS) -c \$<
include Makefile.inc
SCCXMAIN = ccx_2.20.c
OCCXF = \$(SCCXF:.f=.o)
OCCXC = \$(SCCXC:.c=.o)
OCCXMAIN = \$(SCCXMAIN:.c=.o)
LIBS = $VENV/lib/spooles.a $PREFIX/lib/libarpack.a $PREFIX/lib/liblapack.a $PREFIX/lib/libblas.a -lpthread -lm -lc
ccx_2.20_i8: \$(OCCXMAIN) ccx_2.20.a
	\$(FC) -O2 -fopenmp -static-libgfortran -static-libquadmath -static-libgcc -o \$@ \$(OCCXMAIN) ccx_2.20.a \$(LIBS) -Wl,-Bstatic -lgomp -Wl,-Bdynamic -ldl
ccx_2.20.a: \$(OCCXF) \$(OCCXC)
	ar vr \$@ \$?
EOF
rm -f ccx_2.20_i8
make -j "$JOBS" -f Makefile_i8_nopastix ccx_2.20_i8 >"$ROOT/make-ccx.log" 2>&1
install -m 755 ccx_2.20_i8 "$VENV/bin/ccx_i8"
echo "installed $VENV/bin/ccx_i8"

# 4. Same build with PARDISO (MKL ILP64, static) as the default solver.
# SPOOLES stays integer*4 and overflows while factoring the es 0.125 m blade
# (SIGSEGV in I2Ohash_insert <- FrontMtx_splitUpperMatrices, SLURM 11757677, L-S78).
# Built with -DPARDISO, CalculiX picks PARDISO by default (ini_cal.c), so decks
# need no SOLVER= change. Skipped when MKL is absent.
MKL=${MKL:-/opt/intel/oneapi/mkl/2022.0.2/lib/intel64}
if [ -f "$MKL/libmkl_gf_ilp64.a" ]; then
  cd "$ROOT/CalculiX/ccx_2.20"
  rm -rf src_pardiso && mkdir src_pardiso
  (cd src && cp *.c *.f *.h Makefile.inc date.pl src_pardiso_dummy 2>/dev/null ../src_pardiso/ || true)
  cd src_pardiso
  sed -e "s#-DLONGLONG\$#-DLONGLONG -DPARDISO -DMKL_ILP64 -I$MKL/../../include#" \
      -e "s#^LIBS = \(.*\)\$#LIBS = \1 -Wl,--start-group $MKL/libmkl_gf_ilp64.a $MKL/libmkl_gnu_thread.a $MKL/libmkl_core.a -Wl,--end-group#" \
      ../src/Makefile_i8_nopastix > Makefile_i8_pardiso
  grep -q -- '-DPARDISO' Makefile_i8_pardiso && grep -q mkl_core Makefile_i8_pardiso
  make -j "$JOBS" -f Makefile_i8_pardiso ccx_2.20_i8 >"$ROOT/make-ccx-pardiso.log" 2>&1
  install -m 755 ccx_2.20_i8 "$VENV/bin/ccx_i8_pardiso"
  echo "installed $VENV/bin/ccx_i8_pardiso"
fi
