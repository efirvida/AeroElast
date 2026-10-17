# Feature: self-consistent dynamic arbiter for the IEA-15 v1.1 blade (OpenFAST v5 + BeamDyn), runs only
Branch: main. Parent: `odd/tasks/validation-l3-b3-audit.md` ("B3-2 RESULT"). Analysis is a SEPARATE session (maintainer, 2026-10-17).

## Specs
- S1: "si claro tienes mi autorizacion para hacer todas las corridas que entiendas necesarias"
- S2: "pero el analisis lo hacemos en otra sesion de trabajo" (this unit produces runs and raw numbers only; no store edits, no verdicts, no conclusions about the shell).
- S3: from the previous report, accepted by S1: "OpenFAST v5 con BeamDyn sobre el deck v1.1" (blade structure from the BeamDyn 6x6 deck instead of ElastoDyn's mode-shape representation), "S-4 1F/2F y reabrir los veredictos de V-02".

## Tasks
- [ ] T1 locate or obtain an OpenFAST v5.0.0 binary (and a BeamDyn-enabled IEA-15 v1.1 deck); record paths.
- [ ] T2 blade-only modal runs, BeamDyn blade: parked (0 rpm) and 7.56 rpm, collective (MBC3) frequencies for 1F, 1E, 2F, 2E (+ torsion if present).
- [ ] T3 same two cases with the ElastoDyn blade of the same deck as a control (must reproduce the stored S-4 targets 0.544/0.739/1.610 parked, 0.5666/0.7461/1.6372 rated within ~1 %; if not, report before anything else).
- [ ] T4 raw result table + handoff block in this file (paths, job ids, commands), no analysis.

## Log
- L1: "continuamos con la validacion ya terminaros las simulaciones de S4" (session origin); L2: "si claro tienes mi autorizacion para hacer todas las corridas que entiendas necesarias"; L3: "pero el analisis lo hacemos en otra sesion de trabajo"

## Results (raw, unanalysed)

2026-10-17; repository baseline `d3420a4`. **Partial campaign: parked BeamDyn and both ElastoDyn controls completed; rotating BeamDyn failed numerically.** No shell comparison, validation-store update, tolerance change, or verdict reopening was performed. S1/S2 covered; S3 partially covered (BeamDyn parked numbers delivered; rotating numbers remain unavailable). T1/T3/T4 delivered; T2 incomplete.

### Executable and input provenance

- Binary: `/scratch/leahk/eduardo.donestevez/tools/openfast-env/bin/openfast`.
- `openfast -v`: **OpenFAST-v5.0.0**, GCC 14.3.0, 64 bit, **single precision**, OpenMP No, compiled Apr 11 2026 22:27:10. Captured in `/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/of_bd/openfast-version.log`. No build was necessary.
- Converted v5 primary deck: `/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/IEA-15-240-RWT-Monopile/IEA-15-240-RWT-Monopile.fst`.
- Converted ElastoDyn input: `/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/IEA-15-240-RWT-Monopile/IEA-15-240-RWT-Monopile_ElastoDyn.dat`.
- BeamDyn geometry input: `/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/IEA-15-240-RWT/IEA-15-240-RWT_BeamDyn.dat`. Its blade input needed two v5 format updates: remove the obsolete pitch-actuator section; add the modal-damping fields that v5 reads even with damping disabled.
- Structural 6x6 input: `/scratch/leahk/eduardo.donestevez/fem-shell/tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat`. SHA256 `bbc76b78ba2a09fbe863663dbb807e9cdbe000af423441580874e694bce28608`; identical to the converted deck's and bundled repository's blade property files. The distributed stiffness/mass text in the working copy is byte-identical; only the damping flag and v5 modal-damping metadata differ.
- Bundled repository found at `/scratch/leahk/eduardo.donestevez/IEA-15-240-RWT`; `git describe --tags --always` reports `v1.1.18-4-ge4993d6`. This is a v1.1-family checkout, **not a verified exact v1.1.0 tag**. The explicitly requested structural file is the authority for these runs.
- Machine-readable provenance and job map: `/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/of_bd/provenance.json`.

### Configuration and extraction

- Three blades; `CompElast=2` for BeamDyn, `1` for ElastoDyn; inflow, aero, servo, sea state, hydro, substructure, mooring, ice and soil disabled. Tower, platform, yaw, drivetrain, generator and pitch DOFs disabled. Rotor speed prescribed through `RotSpeed` with `GenDOF=False`; all pitches zero. Original cone, shaft tilt, hub radius and reference-axis geometry retained.
- **Gravity = 0 m/s²**. BeamDyn physical damping disabled (`damp_type=0`). ElastoDyn control retains its original blade damping; tabulated frequencies are natural frequencies `abs(lambda)/(2*pi)`, not damped frequencies.
- BeamDyn: quartic interpolation (`order_elem=4`), 15 or 30 members/elements (61 or 121 unique nodes per blade), Gaussian quadrature (`quadrature=1`), `RotStates=True`, `QuasiStaticInit=True`. Three geometry keypoints per member, generated with SciPy `CubicSpline` through the original 50 geometry points, uniform element boundaries in reference z.
- Original cases: `DT=0.005 s`, `ModCoupling=3`, numerical `rhoinf=0`. Parked: one linearization at 1 s, `TMax=1.01 s`. Rotating: 12 times starting at 30 s, separated by `60/7.56/12 s` (30-degree spacing), `TMax=37.285 s`. `Linearize=True`, `CalcSteady=False`, `LinInputs=LinOutputs=0`, `LinOutJac=LinOutMod=False`, `OutFmt="ES24.16E3"`.
- Separate rotating numerical attempts retain the same physical configuration: `_dt001` uses `DT=0.001 s`, numerical `rhoinf=0.8`, 12 times starting at 0.1 s, `TMax=7.385 s`; `_loose` additionally sets `ModCoupling=1`. Neither reached its first linearization. No tolerance was widened.
- Postprocessing uses the installed `/scratch/leahk/eduardo.donestevez/venv` `openfast_toolbox.linearization.mbc.fx_mbc3`. It averages the transformed A matrices and extracts the collective state indices from `RotTripletIndicesStates2[:,0]` and their velocities. Eigenvalues of that collective block give the table below; complete MBC spectra and dominant-state descriptions are saved in each successful case's `modes.json`.
- Collective-to-other-state off-block norm divided by collective-block norm: ED parked `1.20456e-8`; ED rotating `7.26626e-8`; BD15 parked `1.90640e-12`; BD30 parked `1.81961e-13`.
- Primary inputs are mirrored under `tests/IEA15MW/of_bd_modal/`. OpenFAST runs the identical `case.fst` under `logs/check/of_bd/<case>/` because it writes beside the primary file. Outputs are `.lin` files plus logs/JSON and tiny **header-only** `.out` files (no time-history rows: empty module OutLists, `TStart=TMax`, `DT_Out=100`). No generated files were written into the working input directories.

### Raw modal frequencies (Hz)

Bending columns are 1F, 1E, 2F, 2E in that order. BeamDyn labels follow the dominant X/Y displacement states; the torsion entry is dominated by rotation in Z. The ElastoDyn blade has only one edgewise modal DOF, so it cannot supply a 2E control.

| Case | rpm | FE elements | 1F / 1E / 2F / 2E (Hz) | First torsion-dominated (Hz) | OpenFAST job | Post job |
|---|---:|---:|---|---:|---|---|
| `ed_0` | 0 | modal representation | 0.543187284 / 0.740340749 / 1.611795996 / unavailable | unavailable | 11617427, COMPLETED | local, completed |
| `ed_7.56` | 7.56 | modal representation | 0.566646499 / 0.746140815 / 1.637324651 / unavailable | unavailable | 11617429, COMPLETED | local, completed |
| `bd15_0` | 0 | 15 | 0.505298359 / 0.694587606 / 1.479336188 / 2.139833302 | 4.071353457 | 11617436, COMPLETED | 11617442, COMPLETED |
| `bd30_0` | 0 | 30 | 0.502270103 / 0.699466965 / 1.479818506 / 2.145863008 | 4.032181038 | 11617440, COMPLETED | 11617443, COMPLETED |
| `bd15_7.56` | 7.56 | 15 | unavailable: run failed before linearization | unavailable | 11617439, FAILED | not run |
| `bd30_7.56` | 7.56 | 30 | unavailable: run failed before linearization | unavailable | 11617441, FAILED | not run |
| `bd15_7.56_dt001` | 7.56 | 15 | unavailable: run failed before linearization | unavailable | 11617444, FAILED | not run |
| `bd15_7.56_loose` | 7.56 | 15 | unavailable: run failed before linearization | unavailable | 11617446, FAILED | not run |

Parked element-count differences, `100*(f_30/f_15 - 1)` (%), without an acceptance verdict: **1F -0.599301; 1E +0.702483; 2F +0.032604; 2E +0.281784; first torsion -0.962147**. Rotating element-count convergence is unavailable.

Control differences from the stored targets (%), ordered 1F / 1E / 2F: parked **-0.149396 / +0.181427 / +0.111553**; 7.56 rpm **+0.008207 / +0.005470 / +0.007614**. All six are within the requested approximately 1% control criterion. This statement concerns only the OpenFAST control, not the shell.

### Absolute case and log paths

All rows use the binary recorded above. Each log directory also contains `case.fst` (the actual primary input passed to the executable), `job.json` (exact sbatch command), and successful `.lin` / `modes.json` outputs where available.

| Case | Working primary deck | Log/output directory |
|---|---|---|
| `ed_0` | `/scratch/leahk/eduardo.donestevez/fem-shell/tests/IEA15MW/of_bd_modal/ed_0/case.fst` | `/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/of_bd/ed_0/` |
| `ed_7.56` | `/scratch/leahk/eduardo.donestevez/fem-shell/tests/IEA15MW/of_bd_modal/ed_7.56/case.fst` | `/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/of_bd/ed_7.56/` |
| `bd15_0` | `/scratch/leahk/eduardo.donestevez/fem-shell/tests/IEA15MW/of_bd_modal/bd15_0/case.fst` | `/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/of_bd/bd15_0/` |
| `bd30_0` | `/scratch/leahk/eduardo.donestevez/fem-shell/tests/IEA15MW/of_bd_modal/bd30_0/case.fst` | `/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/of_bd/bd30_0/` |
| `bd15_7.56` | `/scratch/leahk/eduardo.donestevez/fem-shell/tests/IEA15MW/of_bd_modal/bd15_7.56/case.fst` | `/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/of_bd/bd15_7.56/` |
| `bd30_7.56` | `/scratch/leahk/eduardo.donestevez/fem-shell/tests/IEA15MW/of_bd_modal/bd30_7.56/case.fst` | `/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/of_bd/bd30_7.56/` |
| `bd15_7.56_dt001` | `/scratch/leahk/eduardo.donestevez/fem-shell/tests/IEA15MW/of_bd_modal/bd15_7.56_dt001/case.fst` | `/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/of_bd/bd15_7.56_dt001/` |
| `bd15_7.56_loose` | `/scratch/leahk/eduardo.donestevez/fem-shell/tests/IEA15MW/of_bd_modal/bd15_7.56_loose/case.fst` | `/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/of_bd/bd15_7.56_loose/` |

Failure evidence retained:

- Initial BD15 parked input-format attempts: jobs **11617430** (`openfast.log`: invalid numerical input reading `NNodeOuts`) and **11617434** (`openfast-2.log`: invalid numerical input reading `n_modes`). The working inputs were converted as described above. Successful attempt **11617436** is in `openfast-3.log`, with `job-3.json` recording its command.
- Original rotating BD15 **11617439** and BD30 **11617441**: `BD_CheckRotMat:DBDSQR did not converge` at **0.01 s**, before any linearization.
- BD15 `_dt001` **11617444**: tight-coupling nonconvergence warnings at steps 2–4, then `BD_CheckRotMat:DBDSQR did not converge` at **0.005 s**.
- BD15 `_loose` **11617446**: `BD_DynamicSolutionGA2:Solution does not converge after the maximum number of iterations` at **0.027 s**. Correction attempts stopped here.
- SLURM wrapper outputs: `/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/of_bd/<case>-<jobid>.out`; post outputs use `<case>-post-<jobid>.out`. No jobs remain outstanding.

## Handoff

**Deliverable state:** raw parked BeamDyn modal runs at 15/30 elements and two successful ElastoDyn controls; **no rotating BeamDyn frequencies**. S1/S2 respected, S3 incomplete. Do not treat failed rotating attempts as modal measurements. No validation store, existing tests, tolerances, `src/`, or git index/commits were changed.

**Script:** `tools/of_bd_blade_modal.py` builds isolated inputs, submits serial OpenFAST jobs through `sbatch -p sequana_cpu --ntasks=1 --time=120`, and postprocesses `.lin` files with MBC3. It does not import aeroelast. Use the existing venv; no aeroenv bootstrap is needed. Completed linearizations are protected from submission overwrite. `--attempt N` preserves separate job records/logs for initialization retries. A fresh checkout is needed for a clean full rerun without overwriting completed artifacts.

Rerun recipe from a fresh checkout (substitute one of the exact case identifiers above):

```bash
PY=/scratch/leahk/eduardo.donestevez/venv/bin/python
OPENBLAS_NUM_THREADS=1 "$PY" tools/of_bd_blade_modal.py build --case bd15_0
"$PY" tools/of_bd_blade_modal.py submit --case bd15_0
# Poll the returned job id; inspect logs/check/of_bd/bd15_0/openfast.log.
# Only after normal completion, submit heavy postprocessing on a compute node:
sbatch -p sequana_cpu --ntasks=1 --mem=48G --time=120 \
  -o logs/check/of_bd/bd15_0-post-%j.out \
  --wrap "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 $PY tools/of_bd_blade_modal.py post --case bd15_0 > logs/check/of_bd/bd15_0/post.log 2>&1"
```

The completed current BD15 parked case used attempt 3, so current-artifact postprocessing requires `post --case bd15_0 --attempt 3`; BD30 parked uses the default attempt 1. Reprocessing overwrites only the derived `modes.json`, not `.lin`. Exact original OpenFAST submission commands are in the per-case job JSON files. Post jobs **11617442** and **11617443** used `--mem=48G`, `--ntasks=1`, `--time=120`, partition `sequana_cpu`, one BLAS thread, and `afterok` dependencies on their corresponding run jobs.

Observed verification: Ruff (`--no-cache`) passed; all eight generated cases passed `openfast_toolbox.io.FASTInputFile` parameter checks (including DOF isolation, gravity, rotor speeds, linearization counts and FE orders/counts); the 6x6 distributed-property text was checked byte-identical; the public submit CLI rejected an attempted overwrite of an existing `.lin` file with nonzero exit and unchanged file hash/job records. No repository tests were edited. This is a runs-only campaign rather than a package behavior change; no formal RED/GREEN test suite or package/physics suite was run. The failed OpenFAST runs remain reported above, not hidden as passing verification.

**Open:** obtain a numerically successful blade-only BeamDyn 7.56 rpm periodic linearization at 15 and 30 elements, then its collective MBC3 frequencies and rotating element-count differences. Determine the cause of the recorded integration failures in a separately scoped continuation; do not silently adjust physics or tolerances. Analysis of these raw results, shell comparisons, validation-store edits and any reopening of S-4/V-02 verdicts remain deferred to the maintainer's separate analysis session. Logs/linearizations are under ignored `logs/check/of_bd/` (about 560 MiB of `.lin` data alone); preserve that directory with this handoff. There are no outstanding jobs or scratch directories.
