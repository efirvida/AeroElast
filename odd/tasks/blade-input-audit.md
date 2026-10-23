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
- T2 | S2-S4, S6 | worker | pending | the standalone reader and the audit tool, with unit tests
- T3 | S4-S6 | parent | pending | run at es 0.25, record the table, apply the S5 decision

## Log

- L1 (2026-10-23) User: "en lo que esperamos vamos a hacer la auditoria campo por campo de la malla contra el yaml". This follows the four-way split agreed in the same session: implementation, input translation, model class, reference quality. Input translation goes first because it is the cheapest discriminator.
