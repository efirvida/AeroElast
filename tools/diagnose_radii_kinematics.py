"""Does the radii feedback carry the wrong-sign second-order term on the production state?

T1 of ``odd/tasks/radii-linear-kinematics.md`` (store gap ``radii_linear_kinematics_axial_term``).
The participant feeds ``r_def[k] = _ref_r[k] + mean(u[strip k] . s)``. Along the strip-centroid
line (trapezoid rule) the fed ``u_z`` splits into the linear prebend/sweep kinematics
``lin = -int (y0' u_y' + x0' u_x') dz`` plus an axial-strain residual, while a geometrically
linear solve omits the second-order shortening ``nl = -1/2 int (u_y'^2 + u_x'^2) dz``.
``|nl| / lin > 1`` means the fed term has the wrong physical sign.

States solved (official IEA-15MW mesh, element_size 1.0, rigid-blade rated aero load):
  control : K u = F_aero  (the de-loading fixture; reproduces lin +1.0973, nl -2.1534 m)
  prod    : (K + K_G + K_SP) u = F_aero + F_cf(X0), 7.56 rpm, axis +y, centre 0 (as every
            case YAML), K_G from the static centrifugal solve with free_dofs
  plus the aero-only and centrifugal-only parts on the stiffened operator.
Omitted: gravity (edgewise; only raises u_x'^2) and the twist feedback's de-loading.

    scripts/aeroenv.sh python -W ignore tools/diagnose_radii_kinematics.py
"""
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from aeroelast.core.assembler import MeshAssembler  # noqa: E402
from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.elements import ElementFamily  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant  # noqa: E402
import tests.validation.blade.test_blade_deloading_vs_reference as fx  # noqa: E402
from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402

RPM = 7.56
OMEGA = RPM * 2 * np.pi / 60
AXIS = np.array([0.0, 1.0, 0.0])

Node._id_counter = 0
MeshElement._id_counter = 0
model = Blade(str(fx.YAML), element_size=1.0)
model.generate_mesh()
mesh = model.mesh
props = model.get_element_properties()
asm = MeshAssembler(mesh, {"elements": {"element_family": ElementFamily.SHELL,
                                        "span_direction": (0.0, 0.0, 1.0),
                                        "properties": props}})
n, dpn = asm.dofs_count, asm.dofs_per_node
r, c, v = asm._rust.assemble_k()
K = coo_matrix((np.asarray(v), (np.asarray(r), np.asarray(c))), shape=(n, n)).tocsr()
root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
free = np.array(sorted(set(range(n)) - {dpn * i + d for i in root for d in range(dpn)}), dtype=np.int64)
coords = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes])
print(f"mesh z range {coords[:,2].min():.3f} .. {coords[:,2].max():.3f} m, nodes {len(coords)}")

kg = asm.assemble_geometric_stiffness(omega=OMEGA, rotation_axis=AXIS,
                                      rotation_center=np.zeros(3), free_dofs=free)
ai, aj, av = kg.getValuesCSR()
Kg = csr_matrix((av, aj, ai), shape=(n, n))
mr, mc, mv = asm._rust.assemble_m()
mr, mv = np.asarray(mr), np.asarray(mv)
nodal_m = np.bincount(mr[mr % dpn == 0] // dpn, weights=mv[mr % dpn == 0], minlength=len(coords))
P = -OMEGA**2 * (np.eye(3) - np.outer(AXIS, AXIS))
rows = np.repeat(np.arange(len(coords)) * dpn, 9) + np.tile(np.repeat(np.arange(3), 3), len(coords))
cols = np.repeat(np.arange(len(coords)) * dpn, 9) + np.tile(np.tile(np.arange(3), 3), len(coords))
Ksp = coo_matrix((np.outer(nodal_m, P.ravel()).ravel(), (rows, cols)), shape=(n, n)).tocsr()
f_cf = np.asarray(asm._rust.centrifugal_load(OMEGA, list(AXIS), [0.0, 0.0, 0.0], asm._rho_per_elem))
print(f"blade mass (lumped) {nodal_m.sum():.0f} kg; sum F_cf = {f_cf[0::6].sum():.4g}, "
      f"{f_cf[1::6].sum():.4g}, {f_cf[2::6].sum():.4g} N")

blade_aero = build_blade_aero_from_aerodyn(fx.AD_PRIMARY)
part = BEMFSIParticipant(mesh, blade_aero, fx.BEM_CONFIG)
f_aero_n, _ = part._compute_forces(np.zeros((len(coords), 3)))
f_aero = np.zeros(n)
for d in range(3):
    f_aero[d::6] = f_aero_n[:, d]


def solve(A, f):
    u = np.zeros(n)
    u[free] = spsolve(A[np.ix_(free, free)].tocsc(), f[free])
    return np.column_stack([u[0::6], u[1::6], u[2::6]])


def cumtrapz(y, x):
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])


ref = np.stack([part._ref_coords[i].mean(axis=0) for i in part._strip_node_indices])
z0 = ref[:, 2]
x0p, y0p = np.gradient(ref[:, 0], z0), np.gradient(ref[:, 1], z0)


def report(name, U):
    um = np.stack([U[i].mean(axis=0) for i in part._strip_node_indices])
    uxp, uyp = np.gradient(um[:, 0], z0), np.gradient(um[:, 1], z0)
    lin = -cumtrapz(y0p * uyp + x0p * uxp, z0)
    nl = -0.5 * cumtrapz(uyp**2 + uxp**2, z0)
    r_def, _ = part._compute_deformed_geometry(U)
    k = int(np.argmax(z0))
    print(f"\n[{name}] tip strip z0 = {z0[k]:.3f} m")
    print(f"  u_x {um[k,0]:+.4f}  u_y {um[k,1]:+.4f}  u_z(fed) {um[k,2]:+.4f}  "
          f"dr_fed {r_def[k]-part._ref_r[k]:+.4f}")
    print(f"  lin {lin[k]:+.4f}  nl {nl[k]:+.4f}  axial residual {um[k,2]-lin[k]:+.4f}  "
          f"|nl|/lin {abs(nl[k])/lin[k]:.3f}")
    print(f"  true 2nd-order position u_z_true ~ {um[k,2]+nl[k]:+.4f}")
    for frac in (0.5, 0.75, 0.9):
        j = int(np.argmin(abs(z0 - frac * z0[k])))
        print(f"   z0 {z0[j]:6.2f}: uz {um[j,2]:+.4f} lin {lin[j]:+.4f} nl {nl[j]:+.4f} "
              f"ratio {abs(nl[j])/lin[j] if lin[j] else float('nan'):.3f}")
    return um[k]


report("control: K, aero only", solve(K, f_aero))
Ks = K + Kg + Ksp
report("prod: K+K_G+K_SP, aero + F_cf", solve(Ks, f_aero + f_cf))
report("prod aero part: K+K_G+K_SP, aero only", solve(Ks, f_aero))
report("prod cf part: K+K_G+K_SP, F_cf only", solve(Ks, f_cf))
