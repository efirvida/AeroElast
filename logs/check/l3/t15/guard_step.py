"""T18: is the pure-torque twist step at z ~ 9.551 (frame-tear / torsion-stiffness guards) a ring-fit artefact?

Same realisation as the guard fixture (es 0.25, wall shear flow +T tip / -T root, rigid modes removed); the
section rotation per station is read three ways: the guards' affine-fit rotation (all nodes), the rigid fit on
all nodes, and the rigid fit on the box nodes (caps + webs, T15). Prints the steps around the blend zone.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from scipy.sparse import coo_matrix  # noqa: E402
from _aeroelast import PyMeshAssembler  # noqa: E402

from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.postprocess.sectional import box_node_mask, rigid_rotation_fit  # noqa: E402
import tests.validation.blade.test_blade_iea15mw_validation as bv  # noqa: E402
import tests.validation.blade.test_blade_section_frame_tear as ft  # noqa: E402
import tests.validation.blade.test_blade_section_torsion_stiffness as torsion  # noqa: E402
import tests.validation.parity.test_thin_walled_tube_torsion as tube  # noqa: E402

Node._id_counter = 0
MeshElement._id_counter = 0
model = Blade(str(ft.YAML), element_size=0.25)
model.generate_mesh()
mesh, props = model.mesh, model.get_element_properties()
coords = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes], dtype=float)
stations = np.asarray(ft._physical_stations(coords), dtype=float)
rings = {zz: np.where(np.abs(coords[:, 2] - zz) < ft.STATION_GAP_TOLERANCE)[0] for zz in stations}
asm = PyMeshAssembler.from_model(bv._to_rust_mesh(mesh, props), props, list(bv.SPAN_DIRECTION), None)
r, c, v = asm.assemble_k()
K = coo_matrix((np.asarray(v), (np.asarray(r), np.asarray(c))), shape=(asm.dofs_count,) * 2).tocsr()
f = torsion._contour_shear_flow(coords, torsion._outer_contour(mesh, coords, rings[stations[-1]]), +torsion.TORQUE)
f += torsion._contour_shear_flow(coords, torsion._outer_contour(mesh, coords, rings[stations[0]]), -torsion.TORQUE)
u = tube._solve_rigid_removed(K, f, tube._rigid_modes(coords))
ut = u.reshape(-1, 6)[:, :3]
box = box_node_mask(mesh)
phi_aff, phi_all, phi_box, nbox = [], [], [], []
for zz in stations:
    nodes = rings[zz]
    cen = coords[nodes].mean(axis=0)
    phi_aff.append(torsion._ring_section_rotation(coords, u, nodes))
    phi_all.append(rigid_rotation_fit(coords[nodes], cen, ut[nodes])[2])
    nb = nodes[box[nodes]]
    nbox.append(len(nb))
    phi_box.append(rigid_rotation_fit(coords[nb], coords[nb].mean(axis=0), ut[nb])[2] if len(nb) >= 6 else np.nan)
phi_aff, phi_all, phi_box = map(np.asarray, (phi_aff, phi_all, phi_box))
for name, p in (("affine(all)", phi_aff), ("rigid(all)", phi_all), ("rigid(box)", phi_box)):
    d = np.diff(p)
    bad = np.where(d <= 0)[0]
    print(f"{name:>12}: non-increasing steps {len(bad)} at z {[round(stations[k + 1], 3) for k in bad]}; "
          f"tip {p[-1]:.6e}")
print("z, ring_n, box_n, dphi affine / rigid-all / rigid-box (station k -> k+1)")
for k in range(len(stations) - 1):
    if 8.5 <= stations[k + 1] <= 13.5:
        print(f"{stations[k + 1]:8.3f} {len(rings[stations[k + 1]]):4d} {nbox[k + 1]:4d} "
              f"{phi_aff[k + 1] - phi_aff[k]:+.3e} {phi_all[k + 1] - phi_all[k]:+.3e} {phi_box[k + 1] - phi_box[k]:+.3e}")
