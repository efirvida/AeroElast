# Feature: the shell blade's per-station shear centre against the BeamDyn deck's (WU-C of #14)

Status: planned 2026-10-07 (T1-T5 open; inputs mapped, nothing written yet)
Owner: next session
Related: issue #14 (roadmap item 2 of #18), `odd/tasks/rated-twist-sign-convention.md` (the unit
that named WU-C and fixed the two sign defects this builds on), issue #20 (item 8: the structural
bound).

## Why this exists

`odd/tasks/rated-twist-sign-convention.md` fixed the rated twist's two sign defects and left one
number readable that was not before: under **their** own reconstructed load our **shell** twists
`1.61x` Zhou's `-3.60 deg` and our **beam** `0.58x`, in one convention and with the same sign. So
our two structural models still disagree by **2.8x** under one load, and the rated magnitude our
shell reports is dominated by an eccentricity: the aerodynamic normal force acts at the 25 %-chord
aerodynamic centre while the section's shear centre sits elsewhere, and the resulting transfer
torque dominates the pitching moment (`mp_only` `+4.7851` deg against `at_ac` `+9.8376` deg at the
tip, and `distortion/|omega| = 0.534` at that ring).

Whether that eccentricity is a **geometry difference** (the mesh's shear centre genuinely sits
where the deck's does not) or a **load-application artefact** is what WU-C decides. It is the unit
the twist test's own docstring names as the next one, and it is the only remaining lever on #14's
"is the 1.61x structural or definitional" question that the tree can pull without the paper's
unpublished inputs.

## What is mapped (2026-10-07, read-only)

### The deck side exists

`K66toPropsDecoupled(K66[i], convention="BeamDyn")` returns `xS` at index **8** and `yS` at **9**
(`tests/validation/blade/test_blade_twist_anchor_beam.py:154,165`), in **metres**, relative to the
**cross-section origin O** (`openfast_toolbox/converters/beam.py:206`). It is the same `xS` the
anchor beam's `m` uses (`:174`) and the tool's beam uses.

### The mesh side does not

No code in `src/`, `tools/` or `tests/` computes a shear centre, elastic axis or tension centre of
the shell mesh. What exists is:

- `SectionalExtractor.section_stiffness(st, membership)` (`src/aeroelast/postprocess/sectional.py`,
  returning at `:276`): a Lagrange-multiplier homogenization that returns **four scalars** —
  `EA`, `GJ`, `EI_flap`, `EI_edge`.
- `SectionalExtractor.ei_from_static_response(properties, load=1e6)`: `(z, EI_flap, EI_edge)`.

Neither returns the coupling terms a shear centre needs. `ComputeShearCenter(K)`
(`openfast_toolbox/converters/beam.py:265-268`) takes a **6x6**, splits it into `K1 = K[0:3, 0:3]`
and `K3 = K[0:3, 3:6]`, solves `K1 Y = -K3` and returns `[x_S, y_S] = [-Y[1,2], Y[0,2]]`, i.e. it
reads the axial-torsion and bending-shear couplings the extractor never computes.

### The datum and the frame, which is where the two prior defects lived

The extractor measures its own frame from the mesh by PCA: chord = the major transverse spread,
flap = the minor, with the comment "so no twist sign convention is assumed"
(`src/aeroelast/postprocess/sectional.py:43`, `:198`). The deck's `xS` is relative to the section
origin in the BeamDyn frame. **The mapping between those two frames is not defined anywhere in one
function**, so WU-C has to state and pin it — the same class of work as
`_resolve_twist_sense` and the `M_z = -Mp` conversion, and the same class of defect if it is left
implicit.

### A number in prose with nothing behind it

`tests/validation/blade/test_blade_rated_twist.py:735` states "the mesh's measured shear centre sits
at 0.477 of the chord". No code computes it, no store row carries it, and it appears nowhere else.
It is a docstring number nobody measures — the same defect class as the S-7 table (#20) and
`tools/diagnose_pitching_moment_sense.py`'s stale prose. WU-C either measures it or strikes it.

### What already tried this

- `tools/diagnose_load_point_sign.py` compares the shell's rotation under a downwind lift at the
  aerodynamic centre (`Mp = 0`) against the beam's with only the anchor's arm term, i.e. it checks
  the arm's **sign**, not a mesh shear centre.
- `tools/diagnose_moment_reference.py` compares three balance points — the arithmetic-mean
  centroid, the aerodynamic centre, and a **shear centre it never implemented**. WU-C is partly
  finishing that.

## Tasks

**T1 is done and it is a negative result, which is why it was a prototype.** `tools/diagnose_shear_centre.py`
assembles the 6x6 by prescribing six unit generalized strains on the upper face of a clamped slice,
and validates it two ways. Measured (element_size 1.0, z = 5.97 / 44.16 / 77.60 / 115.21 m):

| check | result |
| --- | --- |
| `EA`, `GJ`, `EI_flap`, `EI_edge` against `section_stiffness` | **rel.diff 0.00e+00 .. 2.19e-16** — machine precision, all four, all four stations |
| `max\|K - K^T\| / max\|K\|` (Maxwell-Betti) | **0.50 / 0.92 / 0.99 / 0.78** — fails |

So the frame, the datum, the DOF mapping and the scaling are right, and the **couplings are not**: a
fully clamped lower face absorbs work outside the section's generalized coordinates, so the mutual
terms are not the duals of the diagonal ones. **This scheme cannot produce the 6x6**, and no shear
centre is computed from it. The cross-check earned its keep twice on the way: it caught a double `dz`
scaling on the rotation states (the diagonals were off by exactly `dz`), and it is what lets the
negative result be stated precisely instead of "nothing matched".

- [x] T1 Prototype in the repo as a diagnostic, with the cross-check as the acceptance gate.
  **Done 2026-10-07**: `tools/diagnose_shear_centre.py`, diagonals at machine precision, couplings
  rejected by reciprocity, full matrix printed so it can be re-mapped without re-running.
- [ ] T2 Pick a formulation whose mutual terms are reciprocal **by construction** and rebuild the
  6x6 with it: a proper influence-coefficient / Saint-Venant treatment, or the energy route, or the
  BECAS-style approach `openfast_toolbox` was built around. Then promote it into
  `SectionalExtractor` with a test that keeps the two checks T1 established (diagonals against
  `section_stiffness`, and symmetry) as the gate.
- [ ] T3 State and pin the **PCA chord/frame to BeamDyn `x`** mapping, and the ordering
  `ComputeShearCenter` assumes (it documents neither), then compare `xS/c` per station on both sides
  in chord fractions. Report the figure that replaces the unbacked `0.477`.
- [ ] T4 Decide: geometry or application? Then assert the mesh-side shear centre as its own
  reference, or record the gap if the two sides genuinely disagree.
- [ ] T5 Verify independently and comment on #14 with the verdict.

## Risks, stated up front

- **The frame transform is the defect-prone part.** Two of this project's last four defects were a
  convention that was implicit (`_twist_mesh_to_bem`'s sense, `Mp`'s sign) and a half-fix (the
  term's sign right, the point's position wrong). The PCA-to-BeamDyn mapping gets pinned with a
  case whose answer is known before it is used.
- **A 6x6 from a homogenization is not free.** The Lagrange-multiplier route is what the extractor
  already trusts for its four scalars, but the couplings must be validated against those scalars
  (T1's cross-check) rather than assumed.
- **The blast radius is contained.** No stored row depends on the mesh shear centre today; the
  magnitude is withdrawn. So this unit adds a measurement and can be abandoned without unpicking
  anything.

## Measurements kept outside the repo

`$SCRATCH/s7_diag/` holds the sign work's logs; T1's are `shear_centre_proto.log` (the first
revision, whose diagonals showed the exact `dz` factor) and `shear_centre_proto2.log` (after the
fix: diagonals at machine precision, symmetry failing).
