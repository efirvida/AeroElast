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

- [ ] T1 Prototype the ring 6x6 in `$SCRATCH` before touching `src/`: reuse the Lagrange-multiplier
  approach `section_stiffness` uses, assemble the full symmetric 6x6 for one ring, and get
  `(x_S, y_S)` out of `ComputeStiffnessProps().ComputeShearCenter`. Cross-check the prototype's
  `EA`/`GJ`/`EI_*` against the existing `section_stiffness` on the same ring, which is a free
  correctness check the prototype must pass before anyone believes its couplings.
- [ ] T2 If T1 holds, add the 6x6 extraction where it belongs (a new `SectionalExtractor` method,
  `src/aeroelast/postprocess/sectional.py`), with the coupling terms named and the frame stated
  once. Production code, so it carries its own test: the diagram terms against the scalars T1
  cross-checked, and a rotation-invariance or a known-section case if one is cheap.
- [ ] T3 State and pin the **PCA chord/frame to BeamDyn `x`** mapping, then compare `xS/c` per
  station on both sides in chord fractions (the deck's is in metres against the section origin, the
  extractor's against the slice centroid, so the normalisation is part of the unit). Report the
  per-station difference and the chord-fraction figure that replaces the unbacked `0.477`.
- [ ] T4 Decide: is the eccentricity geometry or application? Then either assert the mesh-side
  shear centre in the store as its own reference (making the `at_ac` placement defensible), or
  record the gap if the two sides genuinely disagree and the deck's `xS` is the one to trust.
- [ ] T5 Verify with an independent read-only verifier (the two sign defects this builds on were
  caught by prose and by a cross-path check, so assume nothing), and comment on #14 with the
  verdict, which is what finally decides between closing it as documented non-transferability and
  promoting the magnitude.

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

`$SCRATCH/s7_diag/` already holds the sign work's logs; T1's prototype and its log go there too.
