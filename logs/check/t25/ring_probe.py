"""T25: why 50 of 694 precomputed ring sections are None on the branch (test_wall_flow_activation).
Rebuilds the test's production fixture and classifies every None ring: group < 3 nodes, ValueError from
ring_section (message), or unusable cell area. Prints strip, ring, z range and group size.
usage (cwd = worktree, PYTHONPATH = <worktree>.pp:<worktree>): python ring_probe.py"""
import collections, numpy as np
import tests.validation.bem.test_wall_flow_activation as T
import aeroelast.solvers.bem.force_projection as FP
from aeroelast.cli.run_bem_fsi import _build_mesh
from aeroelast.models.blade.model import Blade
from aeroelast.solvers.bem.fsi_participant import build_from_config

mesh, viz = _build_mesh(T._production_cfg(), T.YAML)[:2]
model = Blade(str(T.YAML), element_size=T.ELEMENT_SIZE); model.generate_mesh()
props = model.get_element_properties()
ids = {e.id for e in mesh.elements}
print(f"coupling mesh {len(mesh.nodes)} nodes, {len(mesh.elements)} elements; props for {len(props)} ids, "
      f"coupling elements without props {len(ids - set(props))}")
cfg = {"blade_file": str(T.YAML), "bem": {"wind_speed": 10.59, "omega": 7.56 * 2 * np.pi / 60, "pitch": 0.0,
       "air_density": 1.225, "dynamic_viscosity": 1.81206e-5, "hub_height": 150.0, "shear_exp": 0.0, "default_re": 1e7}}
p = build_from_config(mesh, cfg, viz_mesh=viz, element_properties=props)._projector
coords = np.array([n.coords for n in mesh.nodes]); span = np.array(T.SPAN_DIR)
why = collections.Counter(); k_flat = 0
for k, strip in enumerate(p._strips):
    for g, (group, sec) in enumerate(zip(p._strip_ring_groups[k], p._strip_ring_sections[k])):
        if sec is None:
            idx = strip.node_indices[np.asarray(group, dtype=np.intp)]
            z = coords[idx] @ span
            if len(group) < 3:
                r = "group<3"
            else:
                try:
                    cells, *_ = FP.ring_section(mesh, idx, span, props, coords=coords)
                    r = "unusable-area" if not FP._cells_have_usable_area(cells, strip.offsets[np.asarray(group)], span) else "?"
                except ValueError as e:
                    r = "ValueError: " + str(e)[:120]
            why[r.split(":")[0] if r.startswith("ValueError") else r] += 1
            print(f"flat {k_flat:4d} strip {k:3d} ring {g:2d} n {len(group):4d} z {z.min():8.3f}..{z.max():8.3f}  {r}")
        k_flat += 1
print("total rings", k_flat, dict(why))
