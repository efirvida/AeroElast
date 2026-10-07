# Feature: the radii feedback's measuring stick (#13)

Status: in progress
Owner: this session (2026-10-07)
Related: issue #13 (roadmap item 1 of #18), store gap `force_projection_axial_extension`,
`odd/tasks/composite-bend-twist-verdict.md` §22.6/§22.7, issue #19 (item 0).
Deliberately off the critical path of #19: everything here runs on the production participant
offline, with no coupled run.

## Why this exists

The premise of #13: the one-way radii feedback re-loads the rotor (`+1.24%` thrust / `+0.95%`
power) and part of that re-loading is an artefact of **measuring a path length** on a linearly
deflected curved axis, not a physical axial extension. Its own "what would close it" lists two
routes (accept-and-document, or a second-order geometric solve) that both need a production
measurement, and both are blocked because no healthy coupled run exists.

The roadmap comment on #18 (2026-10-07) proposed a third route: refute the premise, since

1. `r_def` in `BEMFSIParticipant._compute_deformed_geometry` is a **span projection**
   (`r_def[k] = mean(x_def . s) + offset`, `s` a fixed unit vector) and never measures a path;
2. production computes no path length anywhere (`grep` for `cumsum | arclen | path_len |
   norm(diff)` across `src/aeroelast/solvers/` and `crates/aeroelast-solvers/src/` returns
   nothing);
3. the converged campaign's own `bem_sectional.csv` shows `r_def` contracting, the opposite
   direction from the `+2.112%` the diagnostic measured **on its own metric**.

(1) and (2) hold. (3) does not carry the argument: the sign is a property of the load case, not
of the metric, and on production's own fixture it is the other way around (E2). So the
refutation has to rest on the metric identity plus a **measured counterfactual**, which is what
this feature produces.

## Evidence as of 2026-10-07

### E1 — production's radius is the span projection (code)

- `src/aeroelast/solvers/bem/fsi_participant.py:632`:
  `r_def[k] = float(np.mean(deformed_coords[idx] @ s)) + self._mesh_datum_offset`, with
  `deformed_coords = self._ref_coords + displacements` and `s = self._span_dir` fixed at init.
- `_ref_r = blade_aero.r` (hub-referenced AeroDyn stations); the module docstring states the
  datum and the expected physics: for a pure flapwise deflection `delta` at radius `r`, the
  projected span shortens by order `delta^2 / (2r)`.
- The Rust FSI rotor loop receives no radius field at all; `rotor.py` passes only the scalar
  `rotor_radius` from configuration.
- `grep -rnE "arclen|path_len|path length|curve length|cumsum|norm\(np\.diff"` over
  `src/aeroelast/solvers/` and `crates/aeroelast-solvers/src/`: no hits.

=> No path term can enter the feedback. The `+2.112%` path growth is real and lives in the
diagnostic's measuring stick.

### E2 — the roadmap's sign argument does not hold (measured)

Fixture `tests/validation/blade/test_blade_deloading_vs_reference.py` at HEAD (production
participant, real IEA-15MW mesh, real AeroDyn `BladeAero`, one shell solve under the
participant's own projected loads): **4 passed**. Printed table, verbatim rows:

```text
  deformed radii only, reference twist   the geometric re-loading (axial-stretch artefact)   +1.24%    +0.95%    +0.5093          +0.0000
  tip nodal displacement: axial (span) +1.1352 m, radial (in-plane) +17.3065 m (u_tip = [ 1.6387 17.2287  1.1352])
```

| quantity | value | reading |
| --- | ---: | --- |
| tip `r_def - _ref_r` (the table's `tip dr`) | **+0.5093 m** | the radius **grows** |
| tip nodal spanwise displacement | **+1.1352 m** | the linear solve **extends** the span in this load case |
| radii-only row | **+1.24% / +0.95%** | coherent with a growing radius, not with a contraction |

Campaign `$SCRATCH/frontiersin_results_corotational_100s/yaw_0` (2026-09-30, pre-fix),
`fluid/<t>/bem_sectional.csv` column `r[m]`, 53 strips, t = 0.01 -> t = 100:

| statistic | value |
| --- | ---: |
| tip | **-0.605 m** |
| sum over strips | **-10.363 m** |
| positive strips | **0 of 53** (contracts monotonically) |

=> The sign of the radius increment flips between a one-way static solve and a converged
coupled one. It is set by the load case, not by the metric, so "production contracts while the
path grows" cannot be the refutation.

### E3 — what T1 measures

1. `r_def(0) - _ref_r`: the **definition bias** between the mesh's span projection and the
   AeroDyn station radii (present with zero displacement).
2. `r_def(u) - r_def(0) = mean(u . s)` per strip: the **deformation increment** production
   actually feeds.
3. `dr_path[k]`: the same increment as the path metric would measure (cumulative deformed
   centroid polyline minus the reference one), i.e. the second-order transverse term.
4. BEM rows against one baseline (`bem(_ref_r, _ref_twist)`): the existing production row
   (`bem(r_def, ref_twist)`, `+1.24%`), the increment-only row with the projection metric
   (`_ref_r + dr_proj`), and the increment-only row with the path metric (`_ref_r + dr_path`).
   The two increment-only rows split the `+1.24%` into definition bias and deformation, and
   the path row says what the artefact would add if production measured a path.
5. Invariance to pin: `r_def(u) - r_def(0) == mean(u . s)` per strip, to float precision.

### E4 — the measurement (T1, 2026-10-07)

`tests/validation/blade/test_blade_deloading_vs_reference.py`, production participant, real
IEA-15MW mesh + real AeroDyn `BladeAero` at the rated point, one shell solve under the
participant's own projected loads. **5 passed.** All rows against the same
construction baseline `_rebuild_bem_solver(_ref_r, _ref_twist)` = 2.541662 MN / 16.389372 MW:

| row fed to the BEM | d thrust | d power | tip dr [m] |
| --- | ---: | ---: | ---: |
| radii, definition bias only (zero displacement) | **-0.38%** | **-0.17%** | **-0.5969** |
| radii, deformation increment only (projection metric) | **+1.89%** | **+1.86%** | **+1.1062** |
| radii, deformation increment only (path metric) | **+3.63%** | **+3.56%** | **+2.1121** |
| radii feedback as production feeds it (bias + increment) | +1.24% | +0.95% | +0.5093 |

Unchanged: twist only `-26.01% / -15.61%`, production `-25.31% / -14.77%`.

- **Metric identity (pinned assert):** `dr_proj[k] == mean(displacements[strip k] . span_dir)`
to `1e-9`. The feedback carries the strip's own spanwise displacement, to first order.
- **Path metric**: `dr_path` tip = `+2.112061 m` against `dr_proj` = `+1.106156 m` (**1.909x**);
the second-order transverse term is **47.6%** of the path increment. `L_ref` = 116.181290 m,
`L_def` = 118.293351 m, `+1.8179%` (a different construction from the store's old
`117.2256 -> 119.7010 m`, which came from the deleted diagnostic).
- **Definition bias (new defect, measured):** `r_def_rigid - _ref_r` spans `+0.397959 m`
(root strip) to `-0.596870 m` (tip strip). It alone moves the BEM `-0.38% / -0.17%`, so the
deformation-driven `+1.89% / +1.86%` is presented to the rotor as `+1.24% / +0.95%`. The
deltas are three separate BEM evaluations, not additive.

Mechanism: the rigid path short-circuits before any deformed geometry is measured
(`fsi_participant.py:910-912`, `disp_max < 1e-12` -> pre-built reference solver on
`blade_aero.r`), while the deformed path rebuilds the stations from the mesh projection
(`:632`) with the datum anchored at one point (`:402`,
`_mesh_datum_offset = _ref_r[0] - _mesh_span.min()`). The rigid-path guard therefore compares
the short-circuit path against itself and cannot see the difference.

=> The premise of #13 is refuted **with a magnitude**: a path metric would overstate the
increment 1.91x and would feed `+3.63% / +3.56%` instead of the measured `+1.24% / +0.95%`,
so the artefact's share of the production row is **zero**. What the row's wording hid is not
an artefact but the new definition bias.

## Tasks

- [x] T1 Measure, on the production participant: the metric split (bias / increment /
  path increment), `L_ref` vs `L_def`, and the BEM rows above. **Done 2026-10-07**;
  numbers in E4. Surface: `tests/validation/blade/test_blade_deloading_vs_reference.py`.
- [x] T2 Rewrite the store entry from the measured numbers. **Done 2026-10-07**:
  `force_projection_axial_extension` leaves `docs/validation/gaps.yaml`; the new defect
  `radii_datum_definition_bias` enters it with the measured numbers; the closure section and
  its re-run trigger land in `docs/validation_closures.md` (header date 2026-10-07); the
  `out_of_scope` reason for the deloading module in `docs/validation/groups.yaml` follows the
  four radii rows; and the dangling reference to the removed id in
  `tests/test_blade_ccx_parity.py` now points at the closure section. Surfaces:
  `docs/validation/gaps.yaml`, `docs/validation_closures.md`,
  `tests/validation/blade/test_blade_deloading_vs_reference.py`, `docs/validation/groups.yaml`,
  `tests/test_blade_ccx_parity.py`.
- [x] T3 Verify. **Done 2026-10-07** (`gentle-ai-verify`, read-only): `5 passed`; 16 of 17
  numeric claims in the two store texts matched the live output, and the one mismatch was the
  run's wall time (24.67 s against the 27.6 s the prose cited), which was removed from all
  three files instead of kept; `validation_matrix.py check` 0 errors / 0 warnings; `ruff check`
  and `ruff format --check` clean; the code facts of E1 re-confirmed at `:397`, `:402`, `:632`,
  `:910-912`; the removed gap id is absent from `gaps.yaml`; both YAML files parse.
- [ ] T4 Close #13 with the measured evidence and annotate #19 with the corrected sign table.
  The #19 note also closes its T4: job `11610151` (`twistfix_cosine_5s`, the campaign's own
  case at HEAD) converged 34 of its 493 windows (6.9%) with the rest on the 30-iteration
  ceiling, thrust ≈ −0.7 MN and not settling — the case is ruled out and the non-contraction
  is a code regression, as the roadmap comment predicted.

## Follow-ups, out of scope here

- If T1 shows the definition bias dominates the `+1.24%`, that is a **new** defect (mesh span
  projection vs AeroDyn station radii, a datum/definition question) and gets its own issue and
  its own store entry. It is not a path-measure artefact and it is not #13.
- The campaign's `r[m]` trend was read as evidence of the code path only; the campaign is
  pre-fix and no number from it is citable.
- `tools/validation_matrix.py coherence` cannot run in this venv (`test_ko2017_performance.py`
  collects no nodes here); unrelated to this feature.
- S-7 WIP (`tests/validation/blade/test_iea15mw_s7_torsion.py` uncommitted delta, store reports
  the module as neither grouped nor declared) is a separate work unit; untouched here.
