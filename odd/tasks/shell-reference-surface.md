# [modelling] shell reference surface: OML nodes and a laminate centred on them

*Opened 2026-10-11, and parked by the maintainer's decision: the blade keeps its skin nodes on the
outer mould line (OML), which is the position FSI needs. This file is where to start when the
structural model is generalised to other problems, or when L3's stiffness excess is attacked.
Measurements: `odd/tasks/s7-torsion-ratio.md` section E8.*

## The state of the code

- **Nodes.** The NuMAD-derived blade generator (`src/aeroelast/core/mesh/generators.py` ->
  `models/blade/numad`) places the skin nodes on the OML. At eta 0.3 / 0.5 / 0.7 the skin ring's
  largest node-to-node extent equals the yaml chord to within `3e-4`.
- **Laminate.** `Laminate::new` (`crates/aeroelast-core/src/materials/laminate.rs`) puts the plies
  from `z = -h/2` to `+h/2` about the node surface, with ply 0 at the bottom (`-h/2`).
- **Element.** MITC4+/D consumes the mid-surface ABD. ADR-6
  (`openspec/changes/archive/2026-09-29-mitc4plusd-faithful/design.md`, spec gap G7) forbids
  `to_shell_constitutive_with_offset` on the element path. The paper's kinematics
  `x = x_m + t a/2 V_n` puts the reference surface at the mid-surface, and shifting the ABD
  without shifting the kinematics has no derivation in either paper.
- **Consequence.** Half of every skin laminate lies outside the aerofoil. The torsion cell
  encloses the OML area, and the spar caps sit t/2 further from the neutral axis.

## What it costs (es = 0.25, official blade)

| quantity | as built (OML) | skin moved to mid-surface | reference |
| --- | ---: | ---: | --- |
| global twist ratio, 0.3-0.9 R | 0.8670 (GJ +15 %) | 0.9732 | BeamDyn K66, integral to z |
| flap compliance at 0.9 L | 3.79e-5 (17-19 % stiff) | 4.41e-5 (4-6 % stiff) | 6x6 rotated by ElastoDyn twist |
| edge compliance at 0.9 L | 2.14e-5 (7.5-10 % soft) | 2.22e-5 (+1.6 to 3.4 %) | same |

First-order Bredt predicts a GJ factor of `(A_mid/A_OML)^2 = 0.864`; the counterfactual measures
0.891-0.904. Not yet verified: that the deck's sectional tool stacks plies inward from the OML.

## The three ways to model it

| option | nodes | laminate | stiffness | FSI coupling surface | cost |
| --- | --- | --- | --- | --- | --- |
| `oml` (default today) | OML | centred | GJ +12-15 %, flap ~18 % stiff | exact (wet surface) | none |
| `mid` | inset by t/2 along the in-section normal | centred | matches the deck to the residual above | 2-5 cm inside the wet surface | low: a generator step, already prototyped |
| `oml_offset` | OML | stacked inward (reference surface at the laminate's outer face) | correct | exact | high: reopens ADR-6 |

- **`mid`** is prototyped in `tools/diagnose_oml_reference_surface.py`. Each skin node moves in its
  section plane by the mean half-thickness of its skin elements, along their averaged normal.
  - The skin normals point outward: the vote is unanimous at es = 1.0. The 175 dissenters at
    es = 0.5 sit on the concave aft pressure side, and a point-in-polygon test confirms the shift
    lands inside.
  - `TE_FLAT` nodes must stay on the OML (`--keep-te-flat`). Otherwise the HP and LP shifts cross
    in the thin trailing edge and 90 elements invert.
  - The shear webs stay connected because they share the moved nodes.
  - Before it becomes a generator flag, check that the L5 guards which read geometry from
    `allOuterShellNods` hold: the AC datum, the wall-flow activation (group 37) and the
    deformed-geometry guards.
- **`oml_offset`** is the target for FSI work, because the wet surface stays exact. It is a
  formulation item on L1, not a mesh option:
  - Derive the offset term in `u_b`/`x_b`, or show that integrating the through-thickness moments
    over the shifted interval is consistent with the element's kinematics.
  - Verify against two arbiters: a closed-form Bredt for a thin tube modelled at its outer radius
    with plies inward (`GJ = 2 pi G t R_mid^3`; group 6 already has a tube), and CalculiX
    `*SHELL SECTION, OFFSET`.
  - Settle the ply order first. With outward normals and ply 0 at `-h/2`, a layup listed
    outer-ply-first is stacked inside out. That flips B for every asymmetric skin and corrupts
    `D - 2 z0 B`. It also affects ply-stress recovery (group 36) under any option.

## Where it matters beyond the blade

Any thin-walled section whose laminate is thick relative to its depth: other blades, spars, box
beams, fuselage or nacelle skins. It also matters whenever a geometry is meshed on its wetted or
outer surface for coupling. `(A_mid/A_OML)^2` is the quick estimate; the error grows as
`t/depth` grows.

## Open items

- [ ] Verify the deck's layer convention (IEA 15 MW report, or the provenance of the yaml's
      `elastic_properties_mb`).
- [ ] Verify how BeamDyn interpolates the 6x6 between key points. This moves the S-7 reference
      by 2-3 %.
- [ ] The ply order against the element normal (see `oml_offset` above).
- [ ] The dynamic flap is in tension already at HEAD: V-02's 1st flap is -3.8 % while the static
      flap is about 18 % too stiff, so mass or the modal set-up is part of the question.
- [ ] Composite MITC3 triangles get no span-angle offset (`node_ids.len() == 4` in
      `crates/aeroelast-py/src/assembler.rs`): 43 at es = 0.5, 88 at es = 0.25.
