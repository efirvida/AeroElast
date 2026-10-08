# Feature: validate the rated load's application on a case with an exact answer (issue #14)

Status: planned 2026-10-07 (T1-T4 open; machinery mapped, nothing written yet)
Owner: next session
Related: issue #14 (roadmap item 2 of #18), `odd/tasks/rated-twist-sign-convention.md` (the sign
fixes this builds on), `odd/tasks/shear-centre-arbitration.md` (which eliminated the shear-centre
hypothesis and left the application as the mechanism).

## Why this exists

The withdrawn rated-twist magnitude has one prerequisite left, stated by the test itself: the load
**application** must be validated on a case with an exact answer. Everything else about #14 is now
settled by measurement:

- The stiffness side is bounded by #20 (shell `+6.4%` globally, `+23.9%` sectionally against the
  BeamDyn deck), so stiffness cannot produce a `1.61x`.
- The shear-centre hypothesis is **eliminated**: the deck puts its shear centre within `0.018 c` of
  its own reference axis, so the eccentricity the rated twist rides is the **application point**
  (`x_AC`, up to `0.2545 c` at the root), not the section's shear centre.
- Two sign defects in the application were found and fixed, and the applications now agree with the
  production projector in sign.

What remains is that the same load set moves the tip twist by a factor of **6.17** (section rotation
across the four `_rated_load_cases` applications) purely by where the force is placed, and that the
shell's rated twist is substantially a **sectional deformation** (`distortion / |omega| = 0.534` at
the tip ring). A magnitude that depends that much on an unvalidated application is not a measurement
of the model. So the application is either validated — on a geometry whose answer is exact — or the
magnitude stays withdrawn with the mechanism proven.

## What is mapped (2026-10-07, read-only)

### The exact case exists

- `tests/validation/parity/test_thin_walled_tube_torsion.py` (moved there from
  `tests/test_thin_walled_tube_torsion.py`; the twist test's docstring still cites the old path) is
  the "exact thin-walled-tube case". A closed rectangular tube (b = 1.0 m, h = 0.6 m, L = 6-24 m,
  t = 0.01 m or 0.5 mm) under `T = 1e4 N.m`, compared against the Bredt-Batho closed form
  `theta' = T / GJ` with `GJ = 4 A^2 A66 / perimeter`. Its helpers: `_tube_mesh()`, `_iso_prop()`,
  `_bredt_isotropic()`, `_wall_traction()`, `_ring_shear_flow()`.
- `tests/validation/parity/test_thin_walled_tube_moment_realization.py` (group 31) validates the
  production `realise_section_load` against the same Bredt rate; `docs/validation/rows/31-...yaml`
  records ~0.31% error for both the hand-written shear flow and the production path.

So the **moment** realisation is validated. The **distributed force plus its in-plane moment** — the
`at_ac` pattern, which is what dominates the rated twist — is not.

### What is reusable and what is missing

Geometry-generic helpers in `tests/validation/blade/test_blade_rated_twist.py`: `_physical_stations`,
`_in_plane_edges`, `_unit_shear_flow` (returns a unit perimeter shear flow and the moment it
realises), `_ring_kinematics` (the affine fit that yields `omega`). Blade-specific: `_section_couples`
(needs the `shell` dict) and `_rated_load_cases`.

**Gap:** nothing applies a transverse force at a chosen chord fraction plus the corresponding
in-plane moment to a tube. `_wall_traction` gives the force part and `_ring_shear_flow` the
out-of-plane torsion, but the `at_ac` pattern itself exists only inside the blade's station
machinery. That helper is the unit's first deliverable.

## The acceptance criterion (decide it before writing the test)

A closed rectangular tube has an exact answer and a **known shear centre: the geometric centre of the
section**. A uniform transverse load `f` per unit length carried through that centre by a perimeter
traction must therefore twist nothing, and the same load with a transfer moment `e * f` per unit
length must reproduce the exact clamped-free twist of a linearly varying internal torque
`T(z) = e f (L - z)`:

    theta(L) = e f L^2 / (2 GJ)          with GJ the Bredt-Batho value the module already pins

So the validation is: **the pattern's tip rotation reproduces the closed form for the offset it was
given.** The bound is the suite's 5% rule, fixed before the run because a bound fitted afterwards is
the defect this project keeps finding. Three outcomes, and all three are results:

1. It reproduces it ⇒ the application is validated, and the blade's magnitude becomes promotable.
2. It does not, and the discrepancy scales with the offset ⇒ the application is the defect, proven
   on a geometry with an exact answer rather than argued from the blade.
3. It does not, and the discrepancy appears only for a thin wall and a single-wall load path ⇒ that
   is the shear-lag/distortion mode the module already documents at 125.7x section shear.

## Verification, and what it narrowed (2026-10-07)

An independent read-only verifier re-checked everything at `7bfe577`: the module's five tests pass,
the printed table matches every claimed number, and the closed form was **re-derived independently**
from the module's own constants (`G = 2.1e11 / 2.6 = 8.0769231e10`, `GJ = 0.45 G t`, `theta = e f L^2 /
(2 GJ)` reproducing `9.904762e-4 / 1.980952e-3` and `1.980952e-2 / 3.961905e-2`), `check` 0 errors,
`status` 0 undeclared, `ruff` clean, and every number in this note's task section matching live output.

Three things it narrowed, all of which change how the result should be stated:

1. **The test validates the pattern, not the production path.** `_ring_traction` and
   `_distributed_pattern_load` are constructed *in the test*; the blade's `at_ac` case and the
   production `ForceProjector` are **not executed here**. Group 31 covers the production *moment*
   realisation; the production *distributed force* side is not exercised by this test. So "the
   application is validated" means "this pattern, on a geometry with an exact answer": carrying it to
   the blade's code path is an inference, and the promotion step (T3) is where the production path
   would get its own evidence.
2. **The control's gate is 2%, not round-off.** The *measurement* is round-off (`1.76e-15`,
   `2.67e-12` rad) but the *assertion* tolerates 2% of the closed form, which for `e = 0.2 m` is a
   resultant up to ~4 mm off-centre. The static torque ruler (1e-9 of `e f L`) is the sharp part: for a
   +y-only load it pins the resultant's chordwise position exactly, so the physics claim rests on the
   ruler rather than on the 2% rotation gate.
3. **Group 6's rows understate the file in three ways, not one**: the new test is reported unclaimed,
   `test_wall_traction_excites_a_section_distortion` is *also* unclaimed (pre-existing), and the
   clamped-root row is unmapped with `measured.status: not_measured`. `regression --group 6` reports
   `4 same, 2 unclaimed, 1 unmapped; source digest: moved`. The re-extract that would fix it resets the
   group's eleven hand-written prose fields (the extractor writes its own placeholder text), so it
   needs a prose-preserving flow or a hand merge; not done blind.

Also confirmed: the linearity check would catch a *fixed* misplacement (ratios `(e + delta) / e` at the
further offset), while a *proportional* one is caught by the magnitude check — the pair covers both.
And `theta_fit`'s 2.1% / 12.5% deviation matches the clamped-root table's own `0.9916` / `0.8825`
pattern, i.e. the documented contamination.

## Tasks

**T1 and T2 are done, and the answer is that the application is validated.**
`tests/validation/parity/test_thin_walled_tube_torsion.py` gained the rated pattern on the tube: a
uniform transverse load carried through the section's shear centre by a perimeter traction, plus the
transfer moment `e * f` per unit length, against the closed form `theta(L) = e f L^2 / (2 GJ)` with
`GJ` from Bredt. Measured (clamped-free tube, `f = 1e5 N/m`, `L = 6 m`):

| wall | offset `e` | closed form | `theta_z` | ratio | `theta_fit` | control: traction alone |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 mm | 0.20 m | 9.9048e-4 | 9.8816e-4 | **0.9977** | 0.9786 | 1.76e-15 rad |
| 10 mm | 0.40 m | 1.9810e-3 | 1.9763e-3 | **0.9977** | 0.9786 | ~0 |
| 0.5 mm | 0.20 m | 1.9810e-2 | 1.9788e-2 | **0.9989** | 0.8750 | 2.67e-12 rad |
| 0.5 mm | 0.40 m | 3.9619e-2 | 3.9577e-2 | **0.9989** | 0.8750 | ~0 |

Three things this settles:

- **The pattern reproduces an exact answer to 0.23% (10 mm wall) and 0.11% (thin wall)**, by the
  metric the module itself converged (`theta_z`). The bound was the suite's 5% rule, fixed before the
  run.
- **The control is exact**: a traction whose resultant passes through the section's centre twists the
  tube by round-off (`1e-15`, `1e-12` rad). The pattern's geometry is right, not merely close.
- **The response is linear in the offset**, identically at `e = 0.2` and `e = 0.4`, which is what the
  closed form requires.

`theta_fit` comes out 2.1% and 12.5% low, worst on the thin wall, which is the section in-plane shear
contamination the module already documents and why it asserts only `theta_z`. Not a new defect.

**So the blade's magnitude is no longer blocked by the application.** The spread across the four
`_rated_load_cases` applications (6.17x) is what *inconsistent* force distributions cost: the tube
shows that the consistent traction plus its transfer moment is exact, and the blade's three other
cases distribute the same force differently. The validated case is `at_ac`.

- [x] T1 The tube-side pattern helper. **Done**: `_ring_traction` (edge-length weighted, resultant at
  the perimeter centroid) and `_distributed_pattern_load` (tributary weights, traction plus transfer
  moment), reusing `_ring_shear_flow`, `_theta_z`, `_theta_fit`, `_solve_clamped`, `_bredt_isotropic`.
- [x] T2 Validate against the closed form at two offsets and two thicknesses. **Done**: the table
  above, `1 passed`; the torque ruler is asserted (the traction alone is torsion-free to 1e-9
  relative, the total equals `e f L`) and the control is a physics assertion, not a bound.
- [ ] T3 Act on the outcome, which is now a **decision** rather than a blocker: the magnitude is
  promotable for the validated application. Promoting it means a store row with its own review (and
  the twist test's docstring, which withdraws the magnitude "until the application is validated on a
  case with an exact answer", stops being accurate as written). The maintainer decides whether to
  promote now or leave the residual reported with the validation on the record. **Store mechanics
  owed either way**: the tube file is already a group (`docs/validation/groups.yaml`), so the two new
  analytical comparisons need their rows (`extract --group <id> --write`, then
  `regression --group <id> --write` to record the measured values); until that runs, `regression` for
  that group reports the new printed residuals as unmapped, which is why it is listed here rather
  than left implicit. **Done 2026-10-07, without `extract`**: the row was written by hand (the
  extractor's id is the test name truncated to 40 characters, its label is `rtol at line N`) because a
  re-extract resets that group's eleven hand-written prose fields. One structural finding came out of
  it: **an assertion inside a loop over four cases is one comparison to the store and four printed
  residuals to the run**, so the row could not map (`4 printed residual(s) for 1 asserted
  comparison(s)`) - the four assertions are now written out at their own call sites, and a helper does
  not substitute because the extractor only follows same-module helpers the group declares. Result:
  `regression --group 6 --write` -> `4 new_baseline, 4 same, 1 unclaimed, 1 unmapped`, `check` 197
  rows / 255 comparisons / 0 errors, module `5 passed`. **Two stalenesses remain and both are
  pre-existing, not this unit's**: the clamped-root row is unmapped for the same loop-assertion reason
  (2 printed for 1 asserted) and `test_wall_traction_excites_a_section_distortion` is unclaimed (the
  `regression` command not honouring `non_validation_tests`, issue #21's class). Both deserve their
  own unit rather than a fix folded into a validation feature.
- [x] T4 (store half) Row and measurement done; verification done (the three findings above). The
  remaining `unmapped`/`unclaimed` entries in group 6 are pre-existing stalenesses with their own
  follow-ups.

## Risks, stated up front

- **The tolerance is the whole test.** The existing tube test's own accuracy (~0.31% in group 31)
  is the natural anchor for the bound; a bound fitted to whatever the new pattern happens to produce
  would be the defect this project keeps finding. State the bound before the run.
- **A tube's shear centre is at its geometric centre, a blade's is not.** A pass on the tube
  validates the *pattern* (force plus moment applied consistently on a ring), not the blade's choice
  of where to put the force — that choice is the aerodynamic centre, and WU-C established it is the
  term that matters.
- **Do not confuse the two applications.** Group 31 already covers the moment realisation; this unit
  is about the force-plus-moment pattern. If the new helper ends up re-testing the shear flow, stop:
  that part is already pinned.
- **Blast radius.** The magnitude is withdrawn and no store row depends on it, so a failed validation
  changes no number; a successful one promotes a residual, which is a decision with its own review.

## Measurements kept outside the repo

`$SCRATCH/s7_diag/` holds the sign work's and WU-C's logs; this unit's go there too.
