# Feature: the rated twist's sign convention, flipped between the two load applications (#14)

Status: delivered 2026-10-07 (T1-T4 done; T5 verification and T6 pending)
Owner: this session (2026-10-07)
Related: issue #14 (roadmap item 2 of #18), issue #20 (item 8, closed: the torsional
stiffness direction), `tests/validation/blade/test_blade_rated_twist.py`,
`tests/validation/blade/test_blade_twist_anchor_beam.py`,
`tools/diagnose_zhou_loads_reverse.py`.

## Why this exists

#14 asks whether our rated tip twist against Zhou et al. 2025 is a structural or a
definitional difference. The route chosen was to arbitrate the sign convention that
`tools/diagnose_zhou_loads_reverse.py` uses, because that tool prints our shell and our beam
under the same load and the two numbers come out with opposite signs (`+5.7822` against
`-2.0790` deg), which makes the pair unreadable and blocks any magnitude against Zhou.

Arbitrating it found a real defect, not a printing nit.

## Evidence as of 2026-10-07

### E1 - `omega` and `theta_z` are the same signed quantity

`_ring_kinematics` fits the ring's in-plane field to `u_x = a11 x + a12 y + tx`,
`u_y = a21 x + a22 y + ty` and takes `omega = 0.5 * (a21 - a12)`
(`tests/validation/blade/test_blade_rated_twist.py:537`). For a rigid rotation `alpha` about
`+z`, `a12 = -alpha` and `a21 = +alpha`, so `omega = alpha`: the same signed rotation the
`theta_z` DOF (index 5) carries. Measured confirmation in the same run (`mp_only` case):
`theta_z = -4.8337` with `omega = -4.7851`, same sign. They cannot legitimately disagree in
sign for a rotate-dominated section, and at the tip the file itself measures
`distortion / |omega| = 0.359`, so distortion is not what separates them.

### E2 - the frame derivation fixes the sense, and it is `omega > 0`

`docs/validation_closures.md` states the convention once, from the deck and the load frame:
the leading edge sits at `+x` and the load frame's downwind thrust at `+y`, so a rigid `+z`
rotation moves the leading edge downwind and **nose-down is `omega > 0`**. Corroborated by the
de-loading sign (a nose-down twist lowers the angle of attack and the thrust, and the measured
de-loading is negative: `-25.31% / -14.77%` at the rated point), by
`tools/diagnose_sign_chain.py` (which arbitrated the chain against the deck's own aerofoils),
and by both rated realisations agreeing in sign: `+8.1048` deg (minimum norm) and `+9.6669` deg
(multi-cell with properties).

### E3 - the two load applications of the same nose-down moment disagree in sign

Measured at HEAD (`pytest tests/validation/blade/test_blade_rated_twist.py -k
"production_loads_is_nose_down or couple_sets_the_physical_sign"`, `3 passed` in 21 s):

| path | what it applies | tip result | calls it |
| --- | --- | ---: | --- |
| the production projector | the rated load, multi-cell realisation | `omega = +9.6669` | nose-down |
| `_section_couples` (the test-local application) | the same nose-down `Mp` as force couples | `-7.3487` | nose-down |

Both cannot be right. The same `nodes` file asserts all three of these:

- `:255` a nose-down couple gives `theta_z < 0` => nose-down is `omega < 0`;
- `:772` the four `_rated_load_cases` applications give `omega < 0` => nose-down is `omega < 0`;
- `:1005` the production path gives `omega > 0` => nose-down is `omega > 0`.

E2 decides it: the production side is the arbitrated one, and `_section_couples` carries the
applied `Mp` with the opposite sign to the production projector's chain.

### E4 - this is why the tool's two numbers look opposite

`tools/diagnose_zhou_loads_reverse.py` drives the shell through the **production projector**
(`:163`) and the beam through `m_b = Mp_b + ((pitch - 0.25) * chord_b - x_shear) * Np_b - ...`
(`:168`), i.e. the same sign family as `_section_couples` and as
`test_blade_twist_anchor_beam.py:174` (`m = Mp + (x_ac - xS) * Np - (0 - yS) * Tp`). So the tool
prints one path in the arbitrated convention (`omega = +5.7822`) and the other in the flipped
one (`phi = -2.0790`) without saying so. The tool's own docstring ("our `omega` is positive
nose-down") is right; its beam line is the odd one, and the pair is presented as comparable.

The tree already half-knows this: the docstring of
`test_section_moment_as_a_couple_sets_the_physical_sign` records that "the `+0.98` deg
attributed to the BeamDyn anchor is the opposite sense to both the literature and the
aerodynamics". That observation is the same defect seen from the anchor side; what was missing
is that the flipped side is the *shared test-local application*, not the anchor alone.

### E5 - blast radius

No stored row moves. Group 27's rows are the rigid BEM table (BEM-only), our applied loads
against our own BEM resultants, and the production twist as a self-reference; the Zhou twist
magnitude is **withdrawn**, not a row. The four-application test's spread and distortion-ratio
assertions are sign-independent in magnitude; only its sign assertions and the anchor's
`tip_artefact_deg < 0` encode the flipped sense. So this is a defect with no citable number
riding on it, which is why it can be fixed in place.

**E5 was too optimistic about one row.** T2c below moves four stored measured values: the
`kind: self` row `blade_rated_twist.test_rated_aero_loads_reproduce_the_bem_resultants.single`
compares the applied moment against `sum(Mp * dr)`, which is the convention the applications
had before T2, so the expectation travels with the application.

### E6 - the second defect: the aerodynamic centre is on the wrong chord side

`tools/diagnose_leading_edge.py` measures the leading edge at the **high-x end** of every ring
(nine stations, two independent methods agreeing; only the degenerate root ring differs), which
is the closure log's convention: leading edge at `+x`, load frame's downwind thrust at `+y`, so
nose-down is `omega > 0`.

`_rated_load_cases` places the aerodynamic centre at `xs.min() + 0.25 * chord`, i.e. 0.25 c
from the **low-x** end. Measured on the rated mesh (element size 1.0, physical rings):

| k | z [m] | x_min | x_max | x_ac(code) | x_ac(LE) | x_c | (ac_code - x_c) | (ac_LE - x_c) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.0 | -2.577 | 2.605 | -1.281 | +1.309 | +0.021 | **-1.302** | **+1.288** |
| 12 | 9.6 | -2.957 | 2.318 | -1.638 | +0.999 | -0.252 | -1.386 | +1.251 |
| 24 | 19.1 | -3.553 | 2.054 | -2.151 | +0.652 | -0.384 | -1.767 | +1.036 |
| 42 | 33.4 | -3.735 | 1.713 | -2.373 | +0.351 | -0.661 | -1.712 | +1.011 |
| 185 | 117.0 | -0.316 | 0.184 | -0.191 | +0.059 | -0.061 | -0.130 | +0.120 |

The transfer a downwind normal force makes when it is moved from the ring centroid to the
centre of pressure is `(x_ac - x_c) * Np`, and `M_z > 0` is nose-down here. The code's column is
negative at every station and the leading-edge column positive, so the **force-driven twist is
inverted** by the same class of error as E3. The helper's docstring records that the *term's*
sign was fixed (":333: `(x_ac - x_c) * F`, not `(x_c - x_ac) * F`") while the *position* was left
on the wrong chord end; a half-fix, and the reason the four applications still measured
nose-up after T2 (`at_ac` = `-0.193` deg against the production path's `+8.1`..`+9.7`).

With E3 fixed and E6 unfixed, the two contributions nearly cancel (the `Mp` part `+10.28` deg
and the force part `-10.46` deg at the tip), which is why the residual showed up as a small
negative number rather than as an obvious inversion.

## Tasks

**State: the fix is complete and green.** Both defects are fixed (T2, T2b), the resultant
invariant's reference travels with the application (T2c), the assertions and docstrings that
named the flipped sense are corrected (T3), and the tool prints one convention and runs as the
issue documents it (T4). `12 passed` across both blade modules, `ruff check` and
`ruff format --check` clean, `validation_matrix.py check` 0 errors / 0 warnings and `status`
0 files neither grouped nor declared. Group 27 was re-measured: **2 changed, 5 same**, both
movements far inside their bound. T5 (an independent read-only verifier) and T6 (the #14
comment) are what remain.

- [x] T1 Run the decisive experiment before changing anything. **Done 2026-10-07**: a pure
  uniform nose-down `Mp` (integral `-1.0029e6 N.m`) through the production projector gives
  `omega = +29.0724` deg and through `_section_couples` gives `-38.6420`, and flipping `Mp`
  flips both exactly. The flipped term is localized; `omega` and `theta_z` agree in sign
  inside each path, as E1 predicts.
- [x] T2 (done) Fix the flipped `Mp`: `_section_couples` and the `Mp` term of
  `_rated_load_cases` now realise `M_z = -Mp`; the in-code comment states why. The anchor table
  and the tool convert once at their comparison boundary, keeping the beam's textbook
  convention. **Still owed**: E6's second defect and the resultant expectation below.
- [x] T2b (done) Fix E6: the aerodynamic centre is now `xs.max() - 0.25 c`, and the four
  applications were re-measured: `mp_only` `+4.8337`, `at_ac` `+15.4921`, `spread_omega`
  `6.170`, `distortion[at_ac]/distortion[mp_only]` `18.367`, `ratio_omega` `2.7327`,
  `ratio_theta_z` `4.3034`. `_rated_load_cases` places the aerodynamic centre at `xs.min() + 0.25 c`,
  which with the leading edge at high-x is the trailing-edge side, so the force transfer's sign
  is inverted. It must be `xs.max() - 0.25 c`. Verify by re-measuring the four applications and
  checking `at_ac` lands nose-down and in the production path's order of magnitude.
- [x] T2c (done) Resolve the resultant expectation: `I_Mz = -I_Mp` is now the reference, and
  the row was re-measured with the maintainer's authorization (`regression --group 27 --write`:
  **2 changed, 5 same**; `at_ac` normal 0.0767% -> 0.1676% and tangential 0.3623% -> 0.1933%,
  both against a 0.5% bound; `mp_only resultant moment` did not move because its reference
  travelled with the application).
  `test_rated_aero_loads_reproduce_the_bem_resultants` compares the applied moment against
  `sum(Mp * dr)`, the convention the applications had before T2. With the frame's `M_z = -Mp`
  the applied resultant is `-sum(Mp * dr)`, so the check's expectation moves with the
  application or the two invariants contradict. Row `blade_rated_twist.test_rated_aero_loads_
  reproduce_the_bem_resultants.single` is stored as `kind: self`, so this **moves four stored
  measured values** and needs `regression --group 27 --write` plus the maintainer's go-ahead.
- [x] T3 (done) Correct what asserts the flipped sense. The couple test's `+7.3487`, the
  `at_ac` sense, the moment-carrying applications' sign assert, the anchor table's ratios
  (`ring/beam` now `4.5622`, was negative) and the docstrings are all corrected. One
  over-broad assert had to be narrowed rather than flipped: "every application gives
  nose-down" is false, because the two applications that place the normal force at the
  perimeter centroid give the opposite sense - which is the line-of-action difference the
  `at_ac` case exists to measure, not a defect.
- [x] T4 (done) Tool: one printed convention, the beam negated once at the print, and the
  invocation fixed. Run as `python tools/diagnose_zhou_loads_reverse.py` (the way the issue
  documents it) it now prints shell `+5.7822`, beam `+2.0790` and Zhou `+3.6000` deg in one
  nose-down-positive convention: `1.61x` and `0.58x`. Was: shell `+5.7822` against beam
  `-2.0790`, unreadable.
- [ ] T5 Verify with an independent read-only verifier: both blade modules green, the new
  cross-path check non-vacuous (it fails on the pre-fix code), group 27's rows re-measured and
  consistent, the store's `check` 0 errors, `ruff` clean, and every quoted number re-read from
  live output.
- [ ] T6 Comment on #14 with the arbitrated convention, the corrected tool output and the
  consequences for its three routes.

## Delivered so far

Committed as work units: the two sign defects with the cross-path consistency check, the
anchor's comparison conversion, the tool, and the re-measured group 27.

## Follow-ups, out of scope here

- **`tools/diagnose_pitching_moment_sense.py`'s docstring is stale in the way the S-7 table was.**
  It states "Today they disagree: -0.8696 deg against +0.0773 deg" and hypothesises that
  flipping `Mp` is the fix; the moment-realisation pair now reads `+8.1048` and `+9.6669` (both
  nose-down, per the closure log), so the disagreement it documents was closed by the
  realisation fix and its premise ("the projector applies `Mp` without converting that sense")
  is contradicted by the production path's measured nose-down result. Its falsifiable test is
  still worth keeping; its prose is not. Not fixed here because the wording question (which side
the polars' `Cm` convention really is) deserves its own measurement rather than another edit.
- The load-application spread (6.170x on the section rotation across the four cases) is #14's
  own finding and is not resolved by having the signs right.
- Whether the withdrawn magnitude can be promoted still needs the application validated on a
  case with an exact answer.

## Measurements kept outside the repo

`$SCRATCH/s7_diag/`: `rated_sign_live.log` (E3), `twist_zhou_live.log` (the four applications),
`zhou_reverse_live.log` (the tool at 0.5 m), `zhou_reverse_es10.log` (the same at 1.0 m), and
the T1 cross-path log once it exists.

## Follow-ups, out of scope here

- The load-application spread (4.3x on the tip twist, 5.5x on the section rotation) is a
  separate finding of #14 and is not resolved by fixing the sign.
- Whether the withdrawn magnitude can be promoted needs the application validated on a case
  with an exact answer; the file names that as its own unit and this feature does not do it.
