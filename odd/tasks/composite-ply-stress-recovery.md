# Feature: composite outer-fibre stress recovery is not delivered yet (#27, item 4 of #18)

Status: planned 2026-10-09; T1-T2 open (first write not made)
Owner: this session (2026-10-09)
Related: issue #27 (roadmap item 4, `P1`, successor of #3), issue #18 (roadmap and
order of record), `docs/validation/gaps.yaml` id `composite_stress_recovery`,
`docs/validation/rows/18-composite_layup.yaml` and `groups.yaml` group 18 (the
stiffness half, closed), `odd/tasks/composite-bend-twist-verdict.md` (the #3 record).
New issue required by the roadmap rule: the smeared recovery found here is an unknown
found mid-item, so it gets its own issue and its own row in #18 instead of being
absorbed into #27.

## Why this exists

Issue #27 was written on the assumption that the missing thing is the **reference**:
"`*SHELL SECTION, COMPOSITE` ignores `OUTPUT=3D`, so CalculiX cannot act as the
independent stress judge it is for isotropic shells". Measured on 2026-10-09, both
halves of that assumption are wrong, in opposite directions:

* CalculiX **can** judge the outer fibre (E2), it just reports one value per node on
  the midsurface, and `OUTPUT=3D` is inert for composites (E3);
* the recovery in the tree **does not produce a ply stress at all** (E1). It produces
  the homogenised, thickness-averaged equivalent stress of the section.

So the item cannot close on a row today: there is nothing on our side to compare
against the reference. The user's decision (2026-10-09) is to build the capability
test-first and then close #27 with the CalculiX + closed-form row, instead of closing
item 4 on a weaker claim.

## Evidence as of 2026-10-09

Reproduce with `module load glu gcc/14.2.0_sequana` plus
`LD_LIBRARY_PATH=/scratch/app_sequana/gcc/14.2.0/lib64` and the probes in
`/tmp/ccx_composite_probe/` (`probe.py`, `probe2.py`, `probe3.py`); see the instrument
notes at the bottom before re-deriving any of this.

### E1 - the delivered recovery is smeared, not ply-level

Membrane tension (`L=1.0`, `B=0.1`, four plies of `0.00125 m`, `E1=120 GPa`,
`E2=10 GPa`, `G12=5 GPa`, `G23=3 GPa`, `nu12=0.3`, `F=1000 N`), centre element, free
field:

| quantity | `[0/90]s` | `[90/0]s` |
| --- | --- | --- |
| CLT closed form, outer ply | **3.696 MPa** | **0.304 MPa** |
| CLT closed form, second ply | 0.304 MPa | 3.696 MPa |
| CLT thickness mean `A/h . eps0` | 2.000 MPa | 2.000 MPa |
| `StressRecovery`, TOP / MIDDLE / BOTTOM | **2.0000 / 2.0000 / 2.0000** | **2.0000 / 2.0000 / 2.0000** |

Stable across `8x2` and `16x4`, and independent of the stack: no `z` dependence in
membrane, no dependence on the ply whose fibre is outermost. Error against the real
outer fibre: **-45.9 %** for `[0/90]s`, **+558 %** for `[90/0]s`. In bending the error is
smaller by accident: our 616.10 MPa against CalculiX 595.59 MPa (-3.33 %) on
`[0/90/90/0]`.

Cause, not a scale slip: `src/aeroelast/postprocess/stress_recovery.py:24-32` defines
`sigma_m = (C / h) . (B_m . u_e)` with `C = Cm() / h`, and for the composite element
`Cm()` returns the laminate's integrated `A`, so `C / h` is the **homogenised section
modulus**. The documented `h` trap is the mechanism of the smear.

### E2 - CalculiX's composite-shell FRD value is the outer-fibre ply stress

Same deck through `write_ccx_mesh(..., quadratic=True)` with
`*SHELL SECTION, COMPOSITE` (S8R), membrane tension, **free-field** node at the
mid-length:

| laminate | CCX 8x2 | CCX 16x4 | CLT outer ply |
| --- | --- | --- | --- |
| `[0/90]s` | 3.6984 MPa | 3.7005 MPa | 3.696 MPa (+0.07 %) |
| `[90/0]s` | 0.2991 MPa | 0.2994 MPa | 0.304 MPa (-1.6 %) |

The value follows the closure ply, so it is the outer fibre and not a thickness mean
(which would be 2.000 MPa for both stacks).

### E3 - `OUTPUT=3D` is inert for composite sections

Bending-dominated `[0/90/90/0]` cantilever, same deck, `OUTPUT=2D` against `OUTPUT=3D`:
max `|SZZ|` 595.588 MPa both, max von Mises 574.521 MPa both, **relative difference
0.000e+00**. This confirms and sharpens the note already committed in
`tests/validation/parity/test_shell_stress_ccx_parity.py` (899.67 MPa unchanged);
that note was measured on a membrane-dominated case, this one on pure bending.

### E4 - the clamped edge is a restraint singularity, not the answer

On the membrane case the maximum interior `SXX` **grows with refinement** (4.4413 MPa
at 8x2, 5.0905 MPa at 16x4) while the free field stays at 3.700. The row must read the
centre of the patch, never the global maximum, or it will validate a mesh-dependent
singularity.

## Contract note for #27

The issue's `closes-when` ("a converged-vs-converged comparison against CalculiX S8R
with `*SHELL SECTION, COMPOSITE` ... plus a row with its reference and its bound") stays
as written; only its premise changes: the CalculiX route exists and the missing half is
ours. The premise in the issue body and in the `composite_stress_recovery` gap reason
must be corrected in the same push as the closure.

## Tasks

* **T1 - RED: ply-level outer-fibre stress against an independent closed form.**
  New `tests/validation/parity/test_composite_ply_stress_parity.py`. The reference must
  be re-implemented in the test (the `_hand_clt_abd` pattern of
  `test_composite_layup_parity.py:375-438`: own `Q`, own `Qbar(theta)`, own `A`), not
  imported from `aeroelast.core.laminate`. Cases: `[0/90]s` and `[90/0]s` membrane,
  asserted against 3.696 MPa and 0.304 MPa. RED today because the smeared value misses
  the bound by 46 % / 558 %. Solve through the **production** path (the standing rule:
  measure through the production path, never through a test-only construction), i.e.
  the Rust assembly plus the production static solver as
  `test_shell_stress_ccx_parity.py::test_production_solver_matches_the_scipy_replica`
  does, feeding `StressRecovery`.
* **T2 - GREEN: ply-resolved location, not a parallel API.**
  **Design decision (2026-10-09): fix the existing API instead of adding a second one.**
  `compute_element_stresses(location=...)` and `compute_nodal_stresses(location=...)`
  promise the stress at a through-thickness location; for a composite the stress at
  `TOP` **is** the outer ply's stress in its own axes, so the smeared answer is a bug in
  that promise, not a missing feature. Resolve the ply that contains the requested `z`
  and evaluate `sigma = Qbar(theta_ply) . eps(z)` with `eps(z)` from the existing strain
  recovery (check the Voigt shear convention against the Rust `compute_stress_field`
  before trusting it). The ply stack comes from the element's
  `CompositeShellProperty.laminate`, because the Rust element only holds ABD. For an
  isotropic section the new path must reproduce the old numbers to machine precision:
  that is the no-regression guard (T7).
* **T3 - coverage beyond the discriminator, plus the additive ply-wise API.** Pure
  bending of a symmetric laminate (`sigma(z)` linear per ply, TOP and BOTTOM of opposite
  sign), an angle-ply case (`[+-45]s`) where the flight to ply axes matters, and a new
  additive API that returns every ply (index, angle, `z_bottom`/`z_top`, stress in global
  and ply axes) evaluated at the ply mid-depth. Additive means it must not change the
  existing entry points.
* **T4 - CalculiX row.** Converged-vs-converged outer-fibre comparison against S8R with
  `*SHELL SECTION, COMPOSITE` on each code's own mesh sequence, read in the free field
  (E4), bound in the style of the isotropic sibling (`TOL_CCX = 0.05`). The row's notes
  carry the measured E2/E3 characterisation of what that CalculiX value is.
* **T5 - store.** Row under group 18 (or a new group with its own `groups.yaml` entry),
  `gaps.yaml` `composite_stress_recovery` removed or sharpened to whatever remains
  unvalidated, gates re-run with `tools/validation_matrix.py` (`check`, `regression`).
  Do not run the `coherence` write path: it rewrites the whole row file and wipes the
  hand-written prose and `measured` blocks (issue #24, still open).
* **T6 - closure and publication.** New issue for the smeared-recovery defect with the
  E1 table, comment on #27, the #18 body rows, the `docs/validation_closures.md`
  section. Push, issue and comment need the user's explicit go-ahead.
* **T7 - guards.** The isotropic outer-fibre test
  (`test_shell_stress_ccx_parity.py`) stays green and unchanged in its numbers.

## Instrument notes (they cost time, do not re-derive)

* `import aeroelast` fails on SDumont without
  `module load glu gcc/14.2.0_sequana` and its `libstdc++` on `LD_LIBRARY_PATH`; a
  verifier subagent burned 30 min on a hanging bash call for exactly this.
* The workable composite CalculiX deck recipe is
  `tests/validation/parity/test_orthotropic_shell_parity.py` (`_write_ccx_inp` +
  `write_ccx_mesh(..., quadratic=True)`); S8R is required, S4 does not accept
  `COMPOSITE`.
* `tests/support/ccx_io.py` gives `parse_frd_stress`, `von_mises_from_voigt`;
  `ccx_bin_or_skip` lives in `tests/conftest.py` (`CCX_BIN` env var, then `PATH`).
* The FRD coordinates of a shell deck are midsurface nodes; the through-thickness
  location of the reported value is not in the file, so E2 had to be pinned against the
  closed form rather than read off.
