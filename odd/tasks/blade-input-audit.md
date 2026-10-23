# Blade input audit: mesh vs yaml, field by field

Does the 0.25 m shell mesh carry what `tests/IEA-15-240-RWT.yaml` says? This audits the yaml →
NuMAD → section → element translation, the one cause category that neither L1/L2 nor the CalculiX
parity can exclude (CCX shares the geometry writer and the CLT). Diagnostic only: no fix and no
store row in this feature. Context: `l3-physical-validation.md`, `twist-discrepancy-closure.md` (P1 J1).

## Specs

- S1 (user, verbatim): "en lo que esperamos vamos a hacer la auditoria campo por campo de la malla contra el yaml"
- S2 Independent reference. The yaml is read by a standalone reader (PyYAML + numpy) that imports
  nothing from `aeroelast.models.blade` or its NuMAD layer. Span interpolation is linear on the
  yaml's own grids. The `nd_arc` convention is taken from the WindIO ontology and cited, not
  assumed.
- S3 Stations: r/R = 0.30, 0.60, 0.90 (the P3a stations) and 0.97 (where K66 collapses, P1 J1).
  Mesh: official yaml, `element_size = 0.25`, the production `Blade` path.
- S4 Fields and pre-registered tolerances (gated unless marked):

  | id | field | mesh side | yaml side | tolerance |
  |---|---|---|---|---|
  | F1 | chord, twist, pitch-axis position, prebend, sweep | the ring of outer nodes at the station | `outer_shape_bem` | chord 0.5 %; twist 0.1°; positions 1 cm or 0.5 % c |
  | F2 | airfoil shape | ring polygon | the yaml airfoils interpolated as the yaml defines (`airfoil_position`) | max normal distance 0.5 % c. If the yaml leaves the blending undefined, the scheme is reported and the field is not gated |
  | F3 | laminate stack per ring element: material, ply thickness, fiber angle, order | the element's section properties | the layers active at (r, `nd_arc` of the element centroid) | total and per-ply thickness 1 %; angle 0.1°; material and order exact |
  | F4 | webs: attachment `nd_arc` on both sides, and web stack | web elements at the station | `webs`, `web*` layers | position within es/2 = 0.125 m; stack as F3 |
  | F5 | integrals per unit length: mass Σρt·ds, axial Σ(E₁-projected)t·ds | ring elements | the yaml stack on the yaml geometry | 1 %. Deck M11/K11 printed as context, not gated |
  | F6 | material constants passed to the element (E1, E2, G12, ν12, ρ) | element material | `materials` | 1e-9 relative |
  | F7 | reference surface: where the nodes sit (OML or mid), and the offset used | mesh nodes and element offset | yaml OML plus the stack thickness | reported, not gated |

  An element whose centroid lies within es/2 of a layer's `start_nd_arc`/`end_nd_arc` is a
  "boundary" element. It is reported apart and not gated, because the boundary discretisation is
  ±es/2 by construction.
- S5 Decision, fixed now. Any gated field outside its tolerance on a non-boundary element is an
  **input-translation defect candidate**. It is reported with its location, then estimated for its
  effect on S-1 (EA, EI, GJ at that station) in a later, separately pre-registered step. It is not
  fixed here. If all gated fields pass, the input translation is cleared at these four stations, and
  the S-1/S-2 gaps point to model class or reference quality (P3, C1).
- S6 Outputs: `tools/audit_blade_input.py` and its unit tests; `logs/check/input-audit/report.json`
  and `stations.csv`; the summary table committed here.

## Tasks

- T1 | S1-S4 | parent | done | pre-registration (this file), committed before any tool edit
- T2 | S2-S4, S6 | worker | done | `tools/audit_blade_input.py`, 24 unit tests (RED: collection error; GREEN 24/24), ruff clean
- T3 | S4-S6 | parent | done | login-node run (~25 s), table below; S5: defect candidates A-E, none fixed
- T4 | S5 | - | pending, new session | pre-register the effect estimate of A and B on S-1 (EA, EI_edge, GJ at 0.30/0.60/0.90), then measure it

## Result (2026-10-23, revision of the tool commit; `logs/check/input-audit/report.json`, `stations.csv`, 1568 checks)

Counts are pass/fail/boundary/not-gated.

| r/R | F1 shape | F2 airfoil | F3 stack | F4 webs | F5 integrals | F6 materials | F7 surface |
|---|---|---|---|---|---|---|---|
| 0.30 | 5/0/0/0 | 0/0/0/1 | 294/72/106/52 | 20/1/0/0 | 2/0/0/6 | 25/0/0/25 | 0/0/0/3 |
| 0.60 | 5/0/0/0 | 0/0/0/1 | 199/17/104/36 | 21/0/0/0 | 0/2/0/6 | 25/0/0/25 | 0/0/0/3 |
| 0.90 | 5/0/0/0 | 1/0/0/0 | 80/44/68/24 | 17/4/0/0 | 0/2/0/6 | 15/0/0/15 | 0/0/0/3 |
| 0.97 | 4/1/0/0 | 1/0/0/0 | 96/8/68/24 | 1/0/0/0 | 0/2/0/6 | 15/0/0/15 | 0/0/0/3 |

F5, mesh vs the yaml point value:

| r/R | mass | axial |
|---|---|---|
| 0.30 | -0.10 % | -0.48 % |
| 0.60 | -5.15 % | -2.28 % |
| 0.90 | -6.86 % | -10.2 % |
| 0.97 | +2.73 % | +2.03 % |

Fields that pass:
- Outer shape: twist within 0.035°, positions within tolerance.
- Airfoil: where gated, within 0.29 % c.
- In-plane materials: exact.
- Fibre angles: all 0 in the yaml and in the mesh.
- Material and ply order: they differ only on boundary elements.
- Nodes: they sit on the OML with the laminate centred on them (F7, reported only).

**S5 verdict: input-translation defect candidates exist.** None was fixed. They are listed with the
worker's numbers, and only D was spot-checked by the parent:
- **A. Partly covered layers are smeared over the whole NuMAD region** (`stackdb._map_continuous_layers`):
  - `Shell_skin` in `TE_REINF` is -15.5 % at 0.30 and -13.3 % at 0.60;
  - `Shell_skin_inner` in `TE_PANEL` is -9 % at 0.30, -23 % at 0.90 and -30 % at 0.97;
  - TE-reinforcement glass mass at 0.60 is 47.1 vs 61.0 kg/m (-23 %), which is the whole F5 mass miss there.

  This bears on EI_edge, which S-1 measures as -7.5 % signed (softer).
- **B. Thickness is held per NuMAD bay**, not evaluated at the element. At 0.90 the spar caps are 10.18 vs
  11.43 mm (-11.0 %). Bay-averaging the yaml brings F5 at 0.90 to -0.9 % mass and -0.6 % axial, so
  B accounts for the 0.90 miss. Mesh refinement cannot remove it: the bay comes from the NuMAD
  stations, not from `element_size`. Rule 9 convergence studies therefore do not see it.
- **C. Chord at 0.97 is +0.77 %** (1.8897 vs 1.8752 m). Pchip versus linear explains part of it, and
  0.5 % is unexplained.
- **D. Out-of-plane constants** (`build_rust_properties`, model.py:721):
  - the Rust order is (g12, g23, g13; nu12, nu23, nu31), and windIO (xy, xz, yz) is passed in;
  - glass_uni g13/g23 come out swapped (0.57 %);
  - nu23 is 0.35 vs 0.264, and 0.0 vs 0.17 for glass_triax.

  How the shell uses these constants (transverse shear) has not been checked.
- **E. The web1 suction-side attachment at 0.30 is 0.157 m off** (tolerance 0.125 m).

Literal-application notes, each recorded in the report meta:
- The ring is the planar cut at z.
- The TE flat-back elements are excluded from F3/F5 (0.001-1.6 % of the ring mass).
- r/R is the blade-length fraction.
- The `nd_arc` convention is inferred from the airfoil order and the `side` tags. The windIO schema
  text is not vendored, so it is not cited, which departs from S2.
- The position tolerance is the larger of the two values.

**Reading.** "Implementation is clean" no longer holds for the input chain. A and B reach
the stiffness of the TE and the spar caps, the members that set EI_edge and EA. Whether they explain
S-1/S-2 is a quantitative question for T4. The answer is not assumed here.

## Log

- L1 (2026-10-23) User: "en lo que esperamos vamos a hacer la auditoria campo por campo de la malla contra el yaml". This follows the four-way split agreed in the same session: implementation, input translation, model class, reference quality. Input translation goes first because it is the cheapest discriminator.
