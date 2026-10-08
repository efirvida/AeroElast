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

**Done, and the answer is that the shear centre was never the eccentricity.** The cheap route first:
the deck's own `xS` against the point the rated load is applied at, both in the deck's datum (the
pitch axis, the frame the anchor beam's `(x_ac - xS) * Np` term already uses and whose sign was
arbitrated on a measurement). `tools/diagnose_shear_centre.py --deck-only`:

| frac | chord [m] | pitch_axis | `xS` [m] | `xS/c` | `x_AC/c` | arm `/c` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.000 | 5.200 | 0.505 | 0.006 | **0.0011** | 0.2545 | 0.2534 |
| 0.150 | 5.647 | 0.376 | 0.040 | **0.0071** | 0.1258 | 0.1186 |
| 0.300 | 5.367 | 0.312 | 0.077 | **0.0142** | 0.0624 | 0.0482 |
| 0.550 | 3.904 | 0.288 | 0.040 | **0.0102** | 0.0379 | 0.0276 |
| 0.850 | 2.525 | 0.323 | 0.044 | **0.0175** | 0.0731 | 0.0555 |
| 1.000 | 0.500 | 0.368 | 0.006 | **0.0115** | 0.1182 | 0.1067 |

`|xS/c|` never exceeds **0.018**, so the deck's shear centre sits essentially on its own reference
axis, and the eccentricity the rated twist rides is the **application point** (`x_AC`, up to 0.2545 c
at the root), not the shear centre.

**That decides WU-C: the eccentricity is not geometry of the shear centre.** A mesh-side shear centre
could differ from the deck's by a percent or two of chord at most, and that cannot move a twist which
differs between our two models by a factor 2.8. So the remaining explanation is the **load
application** (and the shell's sectional deformation, `distortion/|omega| = 0.534` at the tip ring),
which is exactly what the four-application test already reports — now with the shear-centre
hypothesis eliminated rather than merely unexamined.

**The 6x6 is therefore de-scoped.** It would only quantify a term already known to be second order,
and it costs a reciprocal-by-construction formulation (T2b below). The prototype's negative result
stays as the evidence that the naive route is closed.

- [x] T1 Prototype the 6x6 with the cross-check as the acceptance gate. **Done 2026-10-07**:
  `tools/diagnose_shear_centre.py`; the four diagonals reproduce `section_stiffness` to machine
  precision (rel.diff `0.00e+00 .. 2.19e-16`), the couplings fail Maxwell-Betti (`0.50 / 0.92 /
  0.99 / 0.78`), so this scheme cannot produce a section stiffness matrix. Negative result, kept.
- [x] T2 (cheap route) The deck's shear centre against the application point. **Done 2026-10-07**:
  `--deck-only`, the table above; verdict in the same section.
- [ ] T2b (de-scoped unless needed) The reciprocal-by-construction 6x6: a proper
  influence-coefficient / Saint-Venant treatment, the energy route, or the BECAS-style approach.
  Only worth it if a mesh-side shear centre is needed for something else, because WU-C no longer
  depends on it. If it is built, T1's two checks are the permanent gate.
- [ ] T3 The unbacked `0.477`. The deck puts the pitch axis between 0.288 and 0.505 of the chord and
  the shear centre within 0.018 c of it, and the shell's sectional stiffness tracks the deck (S-1),
  so `tests/validation/blade/test_blade_rated_twist.py:735`'s "the mesh's measured shear centre sits
  at 0.477 of the chord" is far from anything measured and has no code behind it. Either measure it
  or strike it; it is the same defect class as the S-7 table (#20) and the stale
  `diagnose_pitching_moment_sense.py` docstring.
- [ ] T4 Land the verdict: comment on #14 with the elimination (the eccentricity is application, not
  shear-centre geometry), which is what lets the issue choose between closing as documented
  non-transferability and promoting the magnitude; and record in the store what WU-C established.
- [ ] T5 Independent read-only verification of the deck profile and the negative result.

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
