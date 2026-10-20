# Validation environment

The environment used by the validation harness. Historical rows of the validation
store (`docs/validation/`) retain their recorded measurement versions; the current
OpenFAST installation does not reproduce the old 4.2.1 environment.
Reproduce it before quoting a number, because several references are version-sensitive
(CalculiX element output, OpenFAST AeroDyn tables, the ccblade BEM solution).

## Pinned versions

| component | version | where |
| --- | --- | --- |
| Python | 3.12.14 | conda env `aeroelast-dev` |
| Ruff | 0.16.0 | lint, also in CI |
| CalculiX (ccx) | 2.23 | `~/miniconda3/envs/aeroelast-dev/bin/ccx` |
| OpenFAST | 5.0.0 | `/scratch/leahk/eduardo.donestevez/conda-envs/openfast/bin/openfast`; set `OPENFAST_BIN` explicitly |
| ccblade | 1.3.1 | editable install; see the fragility note below |
| neuralfoil | 0.3.3 | polar generation |
| preCICE | 3.4.0 | `libprecice.pc` in the conda env |
| PETSc / petsc4py | 3.25.5 | conda env |
| SLEPc / slepc4py | 3.25.1 | conda env |
| NumPy / SciPy | 2.5.3 / 1.18.1 | conda env |
| maturin | 1.15.0 | builds `_aeroelast` |
| Rust | stable | `cargo test -p aeroelast-core` |

The CalculiX version matters: the modal mesh-convergence test of §4.8b measured a
high-mode worst gap of 2.12% under CCX **2.23**, while issue #8 measured 4.03% under
CCX **2.20**. A missing external tool makes the affected rows **skip**, not fail; a tool that is present and
refuses its input is a different case and shows up as an error, which is the signal it should be.

**The cluster has only OpenFAST v5.0.0, not the historical 4.2.1 binary.** Use
`OPENFAST_BIN=/scratch/leahk/eduardo.donestevez/conda-envs/openfast/bin/openfast`;
`aerodyn_driver` sits next to it. `find_openfast_bin` resolves `OPENFAST_BIN` first,
then sibling/home `openfast` environments, then PATH. Its fallback locations do
not cover this cluster installation, so set the variable explicitly.

The parity harness writes the v5 driver layout, including the mandatory SeaState
section (`CompSeaSt=0`, `SeaStFile="unused"`, with `MHK=0`). On the temporary
AeroDyn primary copy it removes the disabled `Buoyancy` switch and adds zero
`TwrCp`/`TwrCa` tower coefficients required by v5. The vendored deck is unchanged.
The format reference is the OpenFAST `v5.0.0` source:
[`AeroDyn_Driver_Subs.f90`](https://github.com/OpenFAST/openfast/blob/v5.0.0/modules/aerodyn/src/AeroDyn_Driver_Subs.f90)
and [`AeroDyn_IO.f90`](https://github.com/OpenFAST/openfast/blob/v5.0.0/modules/aerodyn/src/AeroDyn_IO.f90).
`AbortLevel` and `ModCoupling` are not inputs to this v5 AeroDyn driver.

The comparison is quasi-steady and tower-free (`UA_Mod=0`, `TwrPotent=0`,
`TwrShadow=0`, `TwrAero=False`, `Wake_Mod=1`, `BEM_Mod=1`, `DBEMT_Mod=0`).
**F6 correction:** the old harness retained the deck's `DBEMT_Mod=2`; the 4.2.1
store numbers therefore included dynamic inflow. Its time scale (~R/U, ~13 s)
exceeds the rotor period (~8 s), filtering yaw/shear induction variation that
quasi-steady per-azimuth CCBlade resolves. Explicitly disabling it implements
K1's quasi-steady intent; old margins are history, not a like-for-like baseline.
No skew, tip/hub-loss or other induction settings are tuned for the new run.

From the worktree root, use the bootstrap and local source path (to avoid importing
another worktree's installed code):

```bash
OPENFAST_BIN=/scratch/leahk/eduardo.donestevez/conda-envs/openfast/bin/openfast PYTHONPATH=$PWD/src scripts/aeroenv.sh python -m pytest -q -p no:cacheprovider -rA -s tests/validation/bem/test_bem_openfast_parity.py
```

Compare the binary version before quoting a number. A missing binary skips;
a binary that refuses its input remains a visible setup error.

## Run the suite

```bash
source ~/miniconda3/etc/profile.d/conda.sh && conda activate aeroelast-dev
export CCX_BIN=~/miniconda3/envs/aeroelast-dev/bin/ccx
# pyproject.toml puts -v in addopts; override it for a quiet, machine-readable run.
python -m pytest -o addopts="" -q -rxX
```

That is the run that produces "429 tests: 416 passed, 13 xfailed, 0 failed, 0 errors"
on the current tree. The 13 `xfail` are the documented validity limits, not failures;
their bounds and drivers are in the **Validity Envelope** section of the matrix.

Rust core:

```bash
cargo test --manifest-path crates/Cargo.toml -p aeroelast-core
```

## Rebuild `_aeroelast`

The PyO3 extension needs HDF5 and preCICE discoverable at build time. On this machine
the HDF5 development files live in the conda env but there is no `hdf5.pc` in its
`pkgconfig`, so point the build script at the env and add the preCICE pkg-config path:

```bash
export HDF5_DIR=$CONDA_PREFIX
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
maturin develop --release
```

Without `HDF5_DIR` the `hdf5-metno-sys` build script panics; without the
`PKG_CONFIG_PATH` the `precice` build script cannot find `libprecice`.

## Environment fragilities (keep out of the cited path)

- **ccblade editable into `/tmp`.** The installed ccblade was an editable install whose
  finder mapped `ccblade` to `/tmp/CCBlade/ccblade`; a `/tmp` cleanup removed it and the
  BEM tests failed with `ModuleNotFoundError`. Rebuild it into a stable directory:

  ```bash
  git clone --depth=1 https://github.com/WISDEM/CCBlade.git ~/CCBlade
  cd ~/CCBlade && pip install --no-build-isolation --no-deps -e .
  ```

- **OpenFAST lives in a second env.** The tests resolve the binary themselves; if the
  `openfast` env is absent the A1/A2/A3 parity rows skip.
- **`pdftoppm`, not `pdftotext`, for PDF numbers.** The recovered scans have no text
  layer over their mathematics; every PDF-sourced cell is a lead until re-read with
  vision.

## CI

`.github/workflows/rust-core.yml` runs `cargo test -p aeroelast-core` and `ruff` on every
push/PR. It does **not** run the external-solver rows, which need the heavy stack of
`.github/workflows/makefile.yml`; a full local run is the command above.
