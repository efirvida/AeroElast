# Reduced-Order VLM vs BEM Validation for IEA-15 MW

This note records the current validation status of the quasi-steady VLM backend
against the existing BEM reference path on the IEA-15 MW rotor.

## Why this exists

The current Python VLM backend is cheap enough to validate only on a reduced
full-rotor aerodynamic lattice. Running it directly on the structural
`RotorMesh` path is not a practical validation route today.

The command below builds a reduced rotor from `BladeAero`, runs the reduced VLM
and the current VLM-plus-polar correction prototype, and compares both against
the existing BEM implementation:

```bash
module load gcc glu
/scratch/leahk/eduardo.donestevez/venv/bin/python3 tools/validate_vlm_vs_bem_sweep.py
```

To run all available validation branches side by side (baseline, hybrid, and
the non-default experimental lift closure), use:

```bash
module load gcc glu
/scratch/leahk/eduardo.donestevez/venv/bin/python3 tools/validate_vlm_vs_bem_sweep.py \
  --model all --experimental-lift-gain 0.15
```

The sweep CLI now reports compact objective metrics directly in the summary:

- `MAE(|1-T_ratio|)`, `MAE(|1-Q_ratio|)`, `MAE(|1-P_ratio|)` on reference points
- zero-pitch `MAE(|1-Q_ratio|)`
- zero-pitch TSR `CP_RMSE`
- `Cp`-peak mismatch (`Δrpm`, `ΔCP`)

For lift-closure debugging, use the sectional diagnostic helper:

```bash
module load gcc glu
/scratch/leahk/eduardo.donestevez/venv/bin/python3 tools/diagnose_vlm_sectional.py \
  --wind-speed 10.659 --omega-rpm 7.559 --pitch-deg 0.0
```

This prints stripwise incident-vs-total angle of attack, polar coefficients,
VLM trailing-panel circulation-based lift indicator, and baseline strip force
contributions. Use `--csv <path>` for reproducible post-processing.

The reusable comparison helpers live in
`src/aeroelast/solvers/aero/validation.py`.

## Validation setup

- Rotor: IEA-15 MW WindIO blade definition
- VLM mesh: reduced full rotor, default `12 x 2` panels per blade
- BEM reference: `aeroelast.solvers.bem.engine.BEMSolver`
- Reference operating points: the same wind-speed / RPM / pitch schedule used in
  `tests/test_iea15mw_v01_rotor_performance.py`
- Extra sweep: fixed `10.659 m/s`, zero pitch, varying RPM to probe TSR behavior

## Current findings

At the representative IEA-15 MW operating points, the baseline reduced VLM
still looks like this:

- thrust ratios stayed within `0.618 - 1.124`
- torque ratios stayed within `0.058 - 0.636`
- the backend reproduces the order of magnitude of axial loading better than the
  order of magnitude of power extraction
- `Ct` trends are qualitatively reasonable
- `Cp` trends are not acceptable yet in the zero-pitch operating region

At the fixed-wind zero-pitch TSR sweep:

- the BEM reference reaches its `Cp` peak near `7 rpm`
- the current VLM backend keeps increasing through the top of the tested RPM
  range instead of reproducing that peak
- the zero-pitch torque ratios collapse to roughly `0.06`, which is the clearest
  signal that the current backend is not ready for power prediction

An experimental `hybrid` validation mode now also exists in
`tools/validate_vlm_vs_bem_sweep.py`. It reconstructs a panel-local chord/span
frame from the reduced mesh, evaluates a self-induction-free incident flow at
each panel, and adds sectional profile drag from `BladeAero` polars on top of
the baseline VLM panel force.

That prototype is useful, but only in a narrow sense today:

- it improves torque ratios for the pitch-controlled reference points
- it also improves the zero-pitch reference torque ratios slightly
- it breaks the strictly monotone `Cp` ramp of the baseline TSR sweep
- it still underpredicts power badly and shifts the reduced-order `Cp` peak too far down in RPM

So it should be treated as a validation scaffold for future force-model work,
not as a production aerodynamic closure.

A second non-default branch now exists for lift-closure research:

- model key: `experimental`
- closure: trailing-panel delta-`Cl` correction relative to VLM trailing-panel
  circulation-derived `Cl`
- control: `--experimental-lift-gain`

This branch is for R&D diagnostics only and is intentionally separated from the
stable `hybrid` default so baseline behavior remains reproducible.

A third experimental branch is available for targeted zero-pitch mismatch work:

- model key: `experimental-gated`
- closure: trailing-panel delta-`Cl`, but only positive (`Cl_polar > Cl_vlm`)
  and gated by incident angle magnitude
- controls: `--experimental-lift-gain` and `--experimental-alpha-gate-deg`

This branch is designed to avoid high-TSR over-correction while recovering
missing lift where the reduced VLM underpredicts sectional loading.

A fourth experimental branch extends the gated variant with explicit TSR
damping so the lift correction weakens above a chosen TSR:

- model key: `experimental-gated-tsr`
- controls:
  - `--experimental-lift-gain`
  - `--experimental-alpha-gate-deg`
  - `--experimental-tsr-cutoff`
  - `--experimental-tsr-falloff`
  - `--experimental-radial-tip-scale` (optional)
  - `--experimental-radial-power` (optional)

Use this branch when the gated model improves torque error but still pushes the
`Cp` peak too far toward high RPM.

## Important caveat

The structural-shell route is not a valid validation path today. A probe using a
coarse structural `RotorMesh` (`1815` panels) eventually finished only after
about `1241 s` and returned absurd integrated loads, on the order of `1e12 N`
and `1e14 N.m`.

That means the reduced-order aerodynamic lattice is currently the only practical
way to validate this backend.

## What this means technically

The current backend is useful as an intermediate rotor-loading model only if the
target is qualitative axial behavior.

It is not yet validated for:

- rotor torque
- rotor power
- `Cp` peak location
- above-rated and zero-pitch power behavior

## Recommended next step

Do not spend more time tuning post-processing around the present inviscid VLM
response.

The next useful improvement is still a better lift reconstruction for the
reduced model. The incident-flow `alpha` reconstruction is now in place and the
drag term is no longer blind, but a chordwise force reconstruction based on the
current multi-panel VLM circulation still overpredicts torque. The remaining
work is therefore in the lift closure, not in polar access or incident-angle
recovery.