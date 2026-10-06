"""Is the shear-centre arm's sign in the beam formula the physical one?

Issue #12's structural anchor drives a 1D beam with the sectional moment

    m(z) = Mp + (x_ac - xS) Np - (y_ac - yS) Tp,   x_ac = (pitch - 0.25) chord

and that beam gives -1.7765 deg where Zhou's coupled GEBT gives -3.60, which is
how "our torque is about 2x short" was measured.  Its arm term integrates to
+1.288109e5 N.m against the polars' Mp of -2.736575e5 N.m: the arm *cancels* 47%
of the pitching moment and carries the opposite sign.

The standard formulation has `M_ea = M_ac + (x_ac - x_ea) L` with both offsets
from the leading edge and L the lift, so with the aerodynamic centre ahead of the
elastic axis its transfer term is **negative** - it *adds* to the nose-down
pitching moment.  If the anchor's is positive, the sign is inverted and the
"deficit" is a sign error, not an aerodynamic one.

The arbiter is convention-free because the shell's own statics is verified
(`tools/diagnose_moment_reference.py`: balancing about the aerodynamic centre
reproduces production exactly, and the realised moment is checked against an
independent ruler).  So:

  * the shell, driven by the same loads with `Mp` zeroed, applies a downwind lift
    at the aerodynamic centre and answers with the physical sign;
  * the beam, driven by `Mp = 0` and its arm term alone, answers with the anchor's
    sign.

If those two disagree the anchor's arm sign is wrong.
"""

import numpy as np
from _aeroelast import PyMeshAssembler
from openfast_toolbox.converters.beam import K66toPropsDecoupled
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

import tests.validation.blade.test_blade_iea15mw_validation as bv
import tests.validation.blade.test_blade_rated_twist as t
from aeroelast.core.mesh.entities import MeshElement, Node
from aeroelast.models.blade.model import Blade
from aeroelast.solvers.bem.engine import BEMSolver
from aeroelast.solvers.bem.force_projection import ForceProjector
from tests.support.openfast_bem import build_blade_aero_from_aerodyn
from tests.validation.blade.test_blade_twist_anchor_beam import (
    BEAMDYN_BLADE,
    ELASTODYN_BLADE,
    GKT_MEDIAN_FACTOR,
    _numeric_rows,
)

SPAN = np.array([0.0, 0.0, 1.0])

rows = _numeric_rows(BEAMDYN_BLADE, "DISTRIBUTED PROPERTIES")
n_st = len(rows) // 13
frac = np.array([rows[i * 13][0] for i in range(n_st)])
K = np.array([np.asarray(rows[i * 13 + 1 : i * 13 + 7]) for i in range(n_st)])
deck = [K66toPropsDecoupled(K[i], convention="BeamDyn") for i in range(n_st)]
x_shear = np.array([p[8] for p in deck])
y_shear = np.array([p[9] for p in deck])
gkt = np.array([p[5] for p in deck])
ed = np.asarray(_numeric_rows(ELASTODYN_BLADE, "DISTRIBUTED BLADE PROPERTIES"))
pitch_frac, pitch_axis = ed[:, 0], ed[:, 1]

Node._id_counter = 0
MeshElement._id_counter = 0
model = Blade(str(t.YAML), element_size=0.5)
model.generate_mesh()
mesh = model.mesh
if mesh is None:
    raise SystemExit("no mesh")
props = model.get_element_properties()
coords = mesh.coords_array
assembler = PyMeshAssembler.from_model(
    bv._to_rust_mesh(mesh, props), props, list(bv.SPAN_DIRECTION), None
)
n_dof = assembler.dofs_count
rows_k, cols_k, vals_k = assembler.assemble_k()
Kmat = coo_matrix(
    (np.asarray(vals_k), (np.asarray(rows_k), np.asarray(cols_k))), shape=(n_dof, n_dof)
).tocsr()
root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
free = np.array([i for i in range(n_dof) if i not in {6 * r + d for r in root for d in range(6)}])
Kff = Kmat[np.ix_(free, free)]
phys = t._physical_stations(coords)
tip = np.where(np.abs(coords[:, 2] - phys[-1]) < t.STATION_GAP_TOLERANCE)[0]

blade_aero = build_blade_aero_from_aerodyn(t.AD_PRIMARY)
bem = BEMSolver(blade_aero, rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0).compute(
    t.V_RATED, t.RPM_RATED, t.PITCH_RATED
)
if bem.Mp is None:
    raise SystemExit("no pitching moment")
proj = ForceProjector(mesh, blade_aero, span_direction=SPAN, element_properties=props)

import dataclasses  # noqa: E402  (kept next to its only use)


def shell_omega(result) -> float:
    forces = proj.project(result)
    f = np.zeros(n_dof)
    f[0::6] = forces[:, 0]
    f[1::6] = forces[:, 1]
    f[2::6] = forces[:, 2]
    u = np.zeros(n_dof)
    u[free] = np.asarray(spsolve(Kff, f[free]), dtype=float).ravel()
    return float(np.rad2deg(t._ring_kinematics(coords, u, tip)["omega"]))


def beam_twist(m_of_z, gkt_of_z) -> float:
    r = float(blade_aero.hub_radius) + frac * float(blade_aero.blade_length)
    floor = float(np.median(gkt_of_z)) / GKT_MEDIAN_FACTOR
    ok = np.minimum(gkt_of_z[:-1], gkt_of_z[1:]) >= floor
    dr = np.diff(r)
    torque = np.zeros_like(m_of_z)
    for i in range(len(m_of_z) - 2, -1, -1):
        torque[i] = torque[i + 1] + (0.5 * (m_of_z[i] + m_of_z[i + 1]) * dr[i] if ok[i] else 0.0)
    theta = np.zeros_like(m_of_z)
    for i in range(1, len(m_of_z)):
        theta[i] = theta[i - 1] + (
            0.5 * (torque[i - 1] / gkt_of_z[i - 1] + torque[i] / gkt_of_z[i]) * dr[i - 1]
            if ok[i - 1]
            else 0.0
        )
    return float(np.rad2deg(theta[-1]))


r_hub = np.asarray(blade_aero.r, dtype=float)
r = float(blade_aero.hub_radius) + frac * float(blade_aero.blade_length)
Np = np.interp(r, r_hub, bem.Np)
Tp = np.interp(r, r_hub, bem.Tp)
Mp = np.interp(r, r_hub, bem.Mp)
chord = np.interp(r, r_hub, blade_aero.chord)
pitch = np.interp(frac, pitch_frac, pitch_axis)
x_ac = (pitch - 0.25) * chord
arm = (x_ac - x_shear) * Np - (0.0 - y_shear) * Tp

zero_mp = dataclasses.replace(bem, Mp=np.zeros_like(np.asarray(bem.Mp)))

print("case B - a downwind lift at the aerodynamic centre, no pitching moment:")
print(f"  shell (physical statics)   omega = {shell_omega(zero_mp):+9.4f} deg")
print(f"  beam, anchor arm term      phi   = {beam_twist(arm * 0 + arm, gkt):+9.4f} deg")
print()
print("case C - the polars' pitching moment alone:")
print(f"  shell                      omega = {shell_omega(dataclasses.replace(bem, Np=np.zeros_like(np.asarray(bem.Np)), Tp=np.zeros_like(np.asarray(bem.Tp)))):+9.4f} deg")
print(f"  beam, Mp alone             phi   = {beam_twist(Mp, gkt):+9.4f} deg")
print()
print(f"full load: shell omega {shell_omega(bem):+9.4f} deg   beam phi {beam_twist(Mp + arm, gkt):+9.4f} deg")
print(
    "physical expectation: with the aerodynamic centre ahead of the shear centre, a downwind "
    "lift must pitch the section NOSE-DOWN, and so must the polars' (negative) Cm."
)
