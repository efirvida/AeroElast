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
section**. Applying a transverse force `F` at an offset `e` from that centre plus the corresponding
in-plane moment must therefore produce the exact twist rate of a pure torque `e * F`:

    theta' = (e * F) / GJ          with GJ the Bredt-Batho value the existing tube test already pins

So the validation is: **the pattern's tip twist reproduces the closed form for the offset it was
given.** Three outcomes, and all three are results:

1. It reproduces it (to a stated tolerance) ⇒ the application is validated. The blade's magnitude
   becomes promotable, and promoting it is its own decision with the store row to match.
2. It does not, and the discrepancy scales with the offset ⇒ the application is the defect, now
   proven on a geometry with an exact answer rather than argued from the blade.
3. It does not, and the discrepancy appears only for a thin wall and a single-wall load path ⇒ that
   is the shear-lag/distortion mode the existing tube test already documents at 125.7x section shear
   and 2050% metric disagreement, and the honest conclusion is that the pattern cannot be validated
   by a shear-flow-equivalent criterion.

## Tasks

- [ ] T1 Write the tube-side application helper: a transverse force at a chosen chord fraction of a
  tube ring, plus the in-plane moment that goes with the pattern, in the same wall-traction style
  `_wall_traction` uses. Reuse the generic helpers; do not re-derive the perimeter machinery.
- [ ] T2 Validate against the closed form: the tip twist for a force at offset `e` equals
  `(e * F) / GJ` from `_bredt_isotropic`, at two offsets and two wall thicknesses (the thin one is
  the interesting one), with the numbers printed and the tolerance justified from the existing tube
  test's own accuracy rather than fitted.
- [ ] T3 Act on the outcome: if validated, say so in the store (the blade row's promotion is then a
  decision, not a blocker) and update the twist test's docstring, which currently withdraws the
  magnitude "until the application is validated on a case with an exact answer". If it fails, record
  the failure with its numbers in the same place, and the magnitude stays withdrawn for a reason
  that is now measured.
- [ ] T4 Independent read-only verification, then comment on #14 with the outcome, which is what
  finally lets it choose between closing as documented non-transferability and promoting the
  magnitude.

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
