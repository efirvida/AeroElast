# Feature: attribute the wall-flow coupled twist movement (#30, item P2 of #18)

Status: T1-T2 done (instrument `dccf9e2`, record in this note). T3, T3b and T4 pending.
Owner: this session (2026-10-10); T3b continues in a clean session - see the handoff at the end.
Related: issue **#30** (roadmap item `P2` of **#18**), issue **#16** (blocked by #30),
`odd/tasks/bem-wall-flow-activation.md` (T5/T6 hold the three measurements),
`docs/validation/gaps.yaml` id `moment_realization_over_delivers`, group 31
(`tests/validation/parity/test_thin_walled_tube_moment_realization.py`), gap
`reference_projection_geometry_approximation`, issue #17 (per-process node ordering).

## The measured state (do not re-measure; it is already on the record)

`odd/tasks/bem-wall-flow-activation.md` §T5/T6 and `docs/validation_closures.md` §`#16` hold
the numbers. In one line: with the same per-strip resultant (`|dF| = 1.403e-16`,
`|dM| = 5.984e-16` in the activation guard) the coupled tip section rotation moves from
`+2.583` to `-27.404` deg (`ROTZ`) / `+1.208` to `-9.270` deg (best fit), flapwise tip
displacement moves 0.179 %, and the aero state moves 0.23 % of thrust. Run C (element-bearing
mesh, no property map) reproduces A to 8-9 digits, so the mechanism is the **realisation
itself**. The section probe (`$SCRATCH/bfs16/probe_section_gj.py`) tied (5.002 vs 5.454 deg
against 8.582 deg of Bredt) and cannot explain the coupled factor.

## T1 - the transfer audit: nearest-neighbor is the identity here (DONE 2026-10-10)

Hypothesis raised by the user: the node-based `nearest-neighbor` exchange loses the
distributed pattern (a known complaint in the OpenFOAM coupled cases), and the ~30 deg is
that. Measured on the three existing run directories, from preCICE's own mapping statistics
(`$SCRATCH/bfs16/{A,B,C}/logs/fluid.log:30-36`, identical in all three):

```text
Computing "nearest-neighbor" mapping from mesh "Fluid-Mesh" to mesh "Solid-Mesh" in "write" direction.
Mapping distance min:0 max:0 avg: 0 var: 0 cnt: 27609
Computing "nearest-neighbor" mapping from mesh "Solid-Mesh" to mesh "Fluid-Mesh" in "read" direction.
Mapping distance min:0 max:0 avg: 0 var: 0 cnt: 27609
```

Both meshes register **27609** vertices and **every** correspondence is at distance **0** in
**both** directions, so the mapping is the exact identity permutation, not an interpolation:

| fact | evidence |
| --- | --- |
| fluid registers node-only, 27609 nodes, 0 elements | `A/logs/fluid.log:6-8` (`Filtered to node set 'allOuterShellNods': 27609 nodes`, `Mesh ready: 27609 nodes, 0 elements`) |
| both sides are the same node set of the same generator | `fluid_yaw_0.yaml` `coupling_node_set: allOuterShellNods`; solid YAML `boundaries: [allOuterShellNods]` |
| the two written meshes are bit-identical as point sets | `meshio` nearest-distance check on `A/fluid_mesh.vtu` vs `A/solid_mesh.vtu`: 32325/32325 points, `max = 0.0` both directions |
| no connectivity is declared anywhere | no `set_mesh_triangles` / `set_mesh_edges` in `src/`; `src/aeroelast/solvers/fsi/base.py:111` and `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs:723` call `set_mesh_vertices` only |

**Verdict.** For the BEM→rotor rotor path the exchange loses **nothing**: a zero-distance
bijection is the identity matrix, so a higher-order mapping (`nearest-projection` needs
connectivity, RBF needs neither) would return the same vector up to round-off. The ~30 deg is
**not** a mapping loss. The user's prior stays valid where it came from - the OpenFOAM cases,
whose fluid mesh is genuinely non-matching - and that is a separate observation, not #30.

Corollary that narrows #30 usefully: the *same* force vector reaches both the min-norm and the
wall-flow solid solve. The difference is therefore entirely in the applied **distribution** and
in how the structure answers it, with no transfer step in between.

## The readings that remain open

1. **The wall flow is physical.** The strip's torsion enters the closed skin cell as a couple,
   the blade twists by its `GJ`, and the minimum-norm field was leaking the moment into
   non-torsional skin deformation. Against this: on the closed tube (group 31) the minimum-norm
   field over-delivers by `rate/Bredt = 31.69` while the wall flow sits at `1.0031`, so the
   over-delivery is real - but the coupled sign is *reversed* on the blade.
2. **The wall flow is a drilling-DOF artefact.** The tangential wall flow excites the
   zero-thickness shell's rotational (drilling) DOFs, which carry no physical stiffness, and the
   ~30 deg is that. First thing to audit, per the issue.
3. **The load frame no longer follows the deformation.** The realisation is built from the
   *undeformed* strip offsets and `span_direction: [0,0,1]` while the solved blade has 16.3 m of
   flapwise deflection and a rotated section; the wall flow is a strongly anisotropic
   (wall-tangential dipole) pattern, the minimum-norm field is nearly isotropic, so the same
   frame error costs the wall flow far more. This is the same class as #19/#26 and the residual
   gap `reference_projection_geometry_approximation`.

Note the arithmetic that all three have to explain: both fields deliver the **same per-strip
force and moment**, so their difference is a self-equilibrated pattern per strip; yet the
consequence is a *distributed* twist growing to the tip. A linear structure answers a
self-equilibrated pattern locally, so either the response is not linear-in-effect (drilling,
corotational frame) or the moment is not actually the same in the *deformed* frame.

## Tasks

**T1 - the transfer audit.** Done above; a null result with the preCICE evidence. Closes when the
numbers and the "not the mechanism" verdict are in this note (this section). Evidence: the three
`logs/fluid.log` mapping blocks, the `meshio` point-set check.

**T2 - isolate the mechanism on the actual load field.** Recompute the two applied nodal fields
(min-norm and wall-flow) on the production mesh at the campaign's own setting, then solve the
**same** structural model twice, once per field, and compare the section twist. Variants: static
linear on the production mesh (skin + webs) and on the skin-only coupling mesh; report the twist
profile, the drilling rotation content and the moment of the difference field about the
*undeformed* and the *deformed* span axis. Closes when the table says whether a linear structure
reproduces the sign flip and its magnitude. If it does not, the mechanism is in the nonlinear /
corotational path and T2 says so explicitly.

**Hypothesis T2 tests (section elongation, code-verified).** `_distribute` is a true constrained
minimum-norm solve, `f = Aᵀ (A Aᵀ)⁻¹ b` (`force_projection.py:788-869`). For a pure span couple the
closed form is `f_j = omega x d_j`, so its **magnitude is proportional to `|d_j|`**, the distance
from the strip centroid, and its direction is perpendicular to the radius - which on a *circular*
ring always coincides with the wall tangent and carries the same magnitude everywhere, i.e. exactly
Bredt for uniform node spacing. On an **elongated** airfoil ring (IEA-15MW: chord ~5 m, thickness
~0.5 m) `|d_j|` runs from ~0.25 m at mid-chord to ~2.5 m at the LE/TE, so the same couple is
delivered by a shear flow that is uniform-*tangential* but **~10x non-uniform in magnitude**, while
Saint-Venant requires the uniform `q = T/(2A)` the wall flow applies. That predicts: the two
realisations nearly agree on a tube (the shape where the min-norm pattern *is* Bredt) and diverge
badly on the blade (the shape where it is not). It also predicts that the blade's difference is a
*local section-distortion* mechanism, not a torsion-magnitude one - which is what T2 has to separate.

*Evidence that the structural side already sees this*: group 31
(`tests/validation/parity/test_thin_walled_tube_moment_realization.py`, re-run 2026-10-10, 1.98 s,
PASSES) measures the production wall flow on the validated closed tube against Bredt's `T/GJ` inside
the validated interior window: `2.75983e-05` against `2.75132e-05`, **0.3091 %** against a 5 % bound
(`rows/31-tube_moment_realization.yaml`). The minimum-norm field on that same tube is the `31.69`
of the T6 note - a number this note does **not** yet reproduce, and part of T3.

### T2 - result (2026-10-10): the wall flow IS Bredt on the deck's own section, and the elongation guess is refuted

Probe: `tools/diagnose_section_moment_realization.py` (standalone, reuses the tube fixture's helpers;
command `scripts/aeroenv.sh python tools/diagnose_section_moment_realization.py`, exit 0). Control =
the group-31 rectangle; treatment = the deck's own outer-skin ring at `z ~= 40 m` (51 nodes, from
`BladeMesh` on `tests/IEA-15-240-RWT.yaml`), extruded to `L = 30.268 m` (`L/dim = 6`, the rectangle's
ratio), `nz = 40`, window `(0.4L, 0.9L)`, `T = 1.0e4 N.m`.

| shape | realisation | BC | rate [rad/m] | rate/Bredt | distortion rms/dim |
| --- | --- | --- | --- | --- | --- |
| rectangle B=1, H=0.6 | wall flow, hand `_ring_shear_flow` | self-eq | 2.759826e-5 | **1.00309** | 1.270e-7 |
| rectangle | wall flow, prod. `realise_section_load` | self-eq | 2.759826e-5 | **1.00309** | 1.270e-7 |
| rectangle | **min-norm `_distribute`** | self-eq | 8.719852e-4 | **31.69331** | 1.0649e-3 |
| rectangle | wall flow / min-norm | clamped | 2.728111e-5 / 2.863174e-4 | 0.99156 / 10.40654 | 5.204e-6 / 6.0596e-4 |
| **airfoil 51 nodes** | wall flow, prod. | self-eq | 2.075047e-5 | **0.95943** | 4.1267e-6 |
| **airfoil 51 nodes** | **min-norm** | self-eq | 3.582830e-5 | **1.65658** | 7.4327e-5 |
| airfoil 51 nodes | wall flow / min-norm | clamped | 2.078209e-5 / 3.392130e-5 | 0.96089 / 1.56841 | 4.0818e-6 / 7.7050e-5 |

Airfoil reference: `A = 5.0655 m2`, `GJ = 4.623668e8 N.m2`, `theta' = 2.162785e-5 rad/m`, with each
wall's own deck `A66` (which varies 16.8x around the perimeter); `A66` enters as a harmonic sum.
Convergence: rectangle `n_seg` 4/8/12 gives wall `1.00309/0.99987/0.99890` and min-norm
`31.69331/31.39697/30.52327`; airfoil 29/51 nodes gives wall `0.96639/0.95943` and min-norm
`1.39723/1.65658`, **still rising**, so the airfoil min-norm ratio is a lower bound. Caveat stated by
the probe: the FEM wall is an **isotropic substitute** with `G*t = A66`, because the deck's layup
frame (`xDir` = span) is not the hand-built element's local frame - feeding the raw composite `cm`
rotated `A66` and gave a spurious wall/Bredt of ~1.74. `E`/`nu` are substitutes.

**Verdicts.** (1) The wall flow sits at Bredt on *both* shapes (`0.96-1.00`, distortion 4e-6), so on
the torsion path it is the physical realisation, on the deck's own section, against a closed-form
reference - reading (b) (a drilling artefact) is **refuted on this path**, and the coupling between
torsion and the wall pattern is not the problem. (2) The minimum-norm field over-delivers on every
non-circular section measured (`31.7x` rectangle, `1.57-1.66x` and rising airfoil), so the **parked
production default is the wrong one** on closed cells - the store's `moment_realization_over_delivers`
is now quantified on the deck's own geometry. (3) The elongation guess in this note is **refuted**:
the airfoil (chord/thickness `3.05`) is more elongated than the rectangle (`1.67`) yet diverges far
*less*. The driver is the **radius-versus-wall-tangent mismatch at the discretisation nodes** (sharp
corners), and section distortion tracks it (`~8400x` amplification on the rectangle against `~18x` on
the airfoil).

**What this does to the issue.** The coupled movement can no longer be read as "the torsion response
to the realisation": statically the two realisations differ by a factor of about `1.6` with the *same
sign*, while the coupled runs differ by a factor of about `13` with a **sign flip**. That is a
quantified contradiction, and it says the coupled anomaly lives **outside** the pure-torsion path.

**Metric caveat found while writing this up.** `ROTZ` is the 6th nodal DOF (`checkpoint.py:34-35`),
written from the full DOF vector; over the mesh it spans `-83.09..+80.31` deg (A) and
`-164.41..+75.87` deg (B), with the extremes at the tip, and it is 2-3x larger than the
displacement-based best-fit rotation at the same station. It is a nodal field, **not** a section
twist, so `ROTZ` and the best fit must not be used interchangeably in any row this issue writes.

**T3 - the arbiter, properly posed.** The validated closed-tube fixture is the arbiter, not a new
construction: `test_thin_walled_tube_torsion.py` already has the interior window
`WINDOW = (0.4L, 0.9L)`, the Bredt `T/GJ` reference and the validated bound. The addition needed is
the **minimum-norm** realisation through the same fixture, which turns the T6 note's `31.69` into a
pinned measurement (or refutes it) - and the probe above has already produced both numbers in the
self-equilibrated and clamped configurations, so the promoted test has its expected values.

**T3b - the candidate the torque experiment cannot see.** In the wall-flow branch the residual
`_distribute` call receives only the **non-span** part of `M_strip` (`residual_moment = M_strip -
M_shear`), while in the min-norm branch `_distribute` receives the **whole** `M_strip`. So the same
`F_strip` is spread over the section by two *different* optimisations, i.e. the two realisations also
differ in the **flapwise force pattern**, not only in the torsion. With 16.3 m of flapwise deflection
that difference acts through the bending-torsion coupling, which a pure-torque experiment by
construction cannot see. Test: apply the actual strip loads (`F` and `M`) through both realisations on
the production mesh (skin + webs) and compare the flapwise nodal pattern, the resulting twist in a
corotational static solve at the run's own deflection, and the same with the aero feedback disabled.
Closes when the sign flip is reproduced statically or the mechanism is shown to need the coupled loop.

**T4 - the verdict and the closure.** Whichever reading is demonstrated, on the case that has the
reference. Already settled by T2: the wall flow is the physical realisation on a case with a
reference. Deliverable: the row(s) in `docs/validation/rows/`, the gap narrowed in `gaps.yaml`, the
section in `docs/validation_closures.md`, and the comment that closes #30 (which unblocks #16's
activation decision). Closes when register and issue agree - and, if T3b reproduces the sign flip, the
activation decision is the wall flow *and* a declared residual for the coupled rotation field.

*Not absorbed:* making the OpenFOAM coupling declare mesh connectivity and use
`nearest-projection` (the user's hint, valid where the meshes do *not* coincide). If it earns its
own issue, it goes to #18's table, not here.

## Traps

- **Environment.** Any command importing `aeroelast` needs `scripts/aeroenv.sh <cmd>` (GCC 14
  `CXXABI_1.3.15` otherwise). The tree is shared with another session
  (`CLAUDE.md -> AGENTS.md` staged, `scripts/`, `odd/tasks/token-efficiency.md` untracked):
  never `git add -A`, commit with explicit pathspecs.
- **Never run `coherence`** on the validation store (#24).
- **A campaign is not a measurement.** The 5 s A/B cost a 5 h queue slot for three runs. Every
  hypothesis here first gets a static or single-participant probe; a coupled re-run only decides
  what the probes cannot.
- **The coupling mesh is skin-only** (`allOuterShellNods`); the solved mesh has the webs. A probe
  that solves the coupling mesh alone is not the production structure - say which one it solved.
- **`element_size: 1.0` closes 1-3 cells per ring**, `0.25` closes exactly one. A multi-cell
  route on the coarse mesh is not the campaign's own state.

## Commits (work units)

| commit | what |
| --- | --- |
| `dccf9e2` | `test(tools): probe the section-moment realisation on the deck's own airfoil ring (#30)` - the instrument, 566 lines, ruff clean, run exit 0 |
| this commit | `docs(odd): record T1-T2 of #30 - the transfer is exact and the wall flow is Bredt` |

Both were committed with **explicit pathspecs**; this tree is shared with another session (staged
`CLAUDE.md -> AGENTS.md`, untracked `scripts/`, `odd/tasks/token-efficiency.md`), so never
`git add -A`.

## Handoff - resume recipe for T3b, in a clean session

**Goal.** Attribute the coupled movement of #16's A/B runs: tip section rotation `+2.583 ->
-27.404` deg (`ROTZ`) / `+1.208 -> -9.270` deg (best fit) at an identical per-strip resultant. T1
and T2 are closed; T3b is the open measurement. T4 is the closure.

**Read first, in this order, and nothing else.**
1. this note's `T1` and `T2 - result` sections - every number and both refutations are there;
2. `gh issue view 30` for the closes-when;
3. `docs/validation/gaps.yaml` id `moment_realization_over_delivers`.
Memory: `mem_search` with `wall flow Bredt airfoil ring` (scope project), then the topic key
`odd/wall-flow-twist-attribution/tasks`.

**The measurement T3b has to make.** In the wall-flow branch `_distribute` receives only the
**residual** non-span moment (`residual_moment = M_strip - M_shear`); in the min-norm branch it
receives the **whole** `M_strip`. So the same `F_strip` is spread over the section by two different
optimisations, and the realisations also differ in the **flapwise** force pattern - which a
pure-torque experiment cannot see. Test: apply the actual strip loads (`F` and `M`) through both
realisations on the production mesh (skin + webs) and compare (a) the flapwise nodal pattern per
strip, (b) the resulting twist in a corotational static solve at the run's own deflection, (c) the
same with the aero feedback disabled. It reproduces the sign flip, or it shows the mechanism needs
the coupled loop - either outcome is a result.

**Where the input data is.** Runs `$SCRATCH/bfs16/{A,B,C}`; probes
`$SCRATCH/bfs16/{probe_section_gj,probe_cells,probe_resultant,probe_ring_counts}.py`,
`compare_ab.py`, `field_ab.py`, `run_ab.srm`; the 5 s case `$SCRATCH/bfs16/case5s`. The 5 s A/B
cost a 5 h queue slot - do not re-run it before T3b says it is needed.

**Traps.** `scripts/aeroenv.sh` for every import (GCC 14, `CXXABI_1.3.15`); never `git add -A`;
never run `coherence` (#24); `ROTZ` is the 6th **nodal** DOF and not a section twist; the coupling
mesh is skin-only while the solved mesh carries the webs - always state which one a probe solved.
