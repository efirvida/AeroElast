# Validation environment

The exact environment in which the rows of the validation store (`docs/validation/`)
were measured.
Reproduce it before quoting a number, because several references are version-sensitive
(CalculiX element output, OpenFAST AeroDyn tables, the ccblade BEM solution).

## Pinned versions

| component | version | where |
| --- | --- | --- |
| Python | 3.12.14 | conda env `aeroelast-dev` |
| Ruff | 0.16.0 | lint, also in CI |
| CalculiX (ccx) | 2.23 | `~/miniconda3/envs/aeroelast-dev/bin/ccx` |
| OpenFAST | 4.2.1 | conda env `openfast`; resolved by the tests via `OPENFAST_BIN`/`PATH` |
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

**OpenFAST has drifted from the table above.** It says 4.2.1 and the binary that answers here
reports OpenFAST-v5.0.0 (GCC 15.3.0, single precision, built 2026-09-19). That matters because
write_aerodyn_dvr in tests/support/openfast_bem.py writes an **OpenFAST 4.x** driver input, and v5
makes mandatory fields that 4.x defaulted: it refuses the deck at AbortLevel first, and at
ModCoupling once AbortLevel is supplied, both with FATAL ERROR and no simulation run. Every test in
tests/validation/bem/test_bem_openfast_parity.py that runs the driver therefore **errors in setup**
-- twelve of the fourteen -- while the two that do not run it pass.

They are not skipped, and they should not be. A skip says the dependency is absent; here it is
present and refusing the input, so the honest signal is an error until the writer is ported to v5.
That is a porting task, not a missing file, and the rows of group 11 carry no measurement because
of it.

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
