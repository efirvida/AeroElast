---
name: build-aeroelast
description: 'Trigger: compilar aeroelast, build aeroelast, instalar aeroelast, maturin, rebuild Rust, pip install aeroelast, GLIBCXX error, assemble_m_lumped, AttributeError PyMeshAssembler. Compile and install the aeroelast Python+Rust package on the SDumont HPC cluster.'
argument-hint: 'Describe the action: build, verify, or diagnose'
---

# Build Aeroelast

Compile and install the `aeroelast` Python+Rust hybrid package on SDumont.

## When to Use

- After pulling Rust changes (new methods, refactors in `crates/`)
- `AttributeError: 'PyMeshAssembler' has no attribute '...'` — stale `.so`
- `GLIBCXX_3.4.3x not found` — wrong `libstdc++` in `LD_LIBRARY_PATH`
- Fresh venv setup on a new node

## Hard Rules

- ALWAYS load `gcc/14.2.0_sequana` before compiling — libprecice and the `.so` need `GLIBCXX_3.4.32+`
- ALWAYS put `GCC14_LIB` **first** in `LD_LIBRARY_PATH` at runtime — system `/lib64` only has 3.4.25
- NEVER use `pip install -e .` alone — it does NOT recompile Rust; use `maturin develop --release`
- `maturin` requires `VIRTUAL_ENV` to be set; activate the venv or export it explicitly
- Compile from `crates/aeroelast-py/` — that's where `Cargo.toml` lives for the PyO3 crate

## Build Procedure

```bash
module purge
module load glu gcc/14.2.0_sequana

export VENV_DIR=$SCRATCH/venv
export VIRTUAL_ENV=$VENV_DIR
export PATH=${VENV_DIR}/bin/:$PATH
GCC14_LIB=$(gcc -print-file-name=libstdc++.so.6 | xargs dirname)
export LD_LIBRARY_PATH=${GCC14_LIB}:${VENV_DIR}/lib/:${VENV_DIR}/lib64/:/usr/lib64/psm2-compat

cd /scratch/$USER/fem-shell/crates/aeroelast-py
maturin develop --release
```

Expected finish: `🛠 Installed aeroelast-py-0.1.0` in ~40–90s.

## Verify After Build

```bash
python3 -c "
import _aeroelast
methods = [m for m in dir(_aeroelast.PyMeshAssembler) if 'assemble' in m]
print('OK — methods:', methods)
"
```

Must print `assemble_m_lumped` in the list. If it doesn't, the old `.so` is loading.

## Diagnose Stale .so

Two `.so` files can conflict:

| File | Description |
|------|-------------|
| `$VENV/site-packages/_aeroelast.so` | Legacy flat module — may be outdated |
| `$VENV/site-packages/_aeroelast/_aeroelast.cpython-312-*.so` | Current package `.so` |

If the flat `_aeroelast.so` exists and is older than the package `.so`, remove it:

```bash
ls -lh $SCRATCH/venv/lib/python3.12/site-packages/_aeroelast*.so
# If the flat .so is older, delete it:
rm $SCRATCH/venv/lib/python3.12/site-packages/_aeroelast.so
```

## Runtime LD_LIBRARY_PATH Pattern (for run.sh / runAll.srm)

```bash
module load glu gcc/14.2.0_sequana
GCC14_LIB=$(gcc -print-file-name=libstdc++.so.6 | xargs dirname)
export LD_LIBRARY_PATH=${GCC14_LIB}:${VENV_DIR}/lib/:${VENV_DIR}/lib64/:/usr/lib64/psm2-compat
```

This resolves dynamically — no hardcoded paths needed.
