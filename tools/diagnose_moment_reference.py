"""Where must the section moment be balanced for the twist to match the reference?

The mechanism module (`test_blade_twist_mechanism.py`, issue #9 WU7) established that
the blade carries no laminate bend-twist coupling at all - all 696 sections are
balanced, D16 = D26 = 0 - and that the twist is a *load-path* response: for one and
the same root moment, moving the line of action of the resultant along the chord
swings the tip twist by tens of degrees.

The projector balances each strip's load about its **arithmetic-mean node centroid**,
which on this mesh sits 0.47 c behind the leading edge, and then adds the transfer
`-cross(centroid - ac, F)`. That transfer is +7.607e5 N.m against the polars' -2.766e5,
so it dominates the applied torque by 2.75x. Balancing the same force and the same
aerodynamic pitching moment about the **aerodynamic centre** instead is the classical
formulation, and the two differ only in where the resultant line of action sits.

This measures the tip section rotation under three realisations of the SAME strip
forces and the SAME aerodynamic pitching moment:

  * about the arithmetic-mean centroid, with the transfer term (production today);
  * about the aerodynamic centre, no transfer (the classical aero formulation);
  * about the shear centre, if the deck's anchor stations are available.

The reference the issue quotes: the deck's decoupled 1D beam gives -1.7765 deg and
Zhou's coupled GEBT gives -3.60 deg, both nose-down.
"""

import numpy as np
from _aeroelast import PyMeshAssembler
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

import tests.validation.blade.test_blade_iea15mw_validation as bv
import tests.validation.blade.test_blade_rated_twist as t
from aeroelast.core.mesh.entities import MeshElement, Node
from aeroelast.models.blade.model import Blade
from aeroelast.solvers.bem.engine import BEMSolver
from aeroelast.solvers.bem.force_projection import ForceProjector, _Strip
from tests.support.openfast_bem import build_blade_aero_from_aerodyn

SPAN = np.array([0.0, 0.0, 1.0])

Node._id_counter = 0
MeshElement._id_counter = 0
model = Blade(str(t.YAML), element_size=1.0)
model.generate_mesh()
mesh = model.mesh
if mesh is None:
    raise SystemExit("no mesh")
props = model.get_element_properties()
coords = mesh.coords_array
asm = PyMeshAssembler.from_model(
    bv._to_rust_mesh(mesh, props), props, list(bv.SPAN_DIRECTION), None
)
n = asm.dofs_count
rows, cols, vals = asm.assemble_k()
K = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))), shape=(n, n)).tocsr()
root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
fixed = {6 * i + d for i in root for d in range(6)}
free = np.array([i for i in range(n) if i not in fixed], dtype=np.int64)
Kff = K[np.ix_(free, free)]
phys = t._physical_stations(coords)
tip = np.where(np.abs(coords[:, 2] - phys[-1]) < t.STATION_GAP_TOLERANCE)[0]

aero = build_blade_aero_from_aerodyn(t.AD_PRIMARY)
bem = BEMSolver(aero, rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0).compute(
    t.V_RATED, t.RPM_RATED, t.PITCH_RATED
)
if bem.Mp is None:
    raise SystemExit("no pitching moment")
mp = np.asarray(bem.Mp, dtype=float)
proj = ForceProjector(mesh, aero, span_direction=SPAN, element_properties=props)


def project_about(reference: str) -> np.ndarray:
    """Realise the same F and the same Mp, balanced about `reference`."""
    forces = np.zeros((coords.shape[0], 3))
    for k, strip in enumerate(proj._strips):
        if len(strip.node_indices) == 0:
            continue
        f_strip = (
            float(bem.Np[k]) * strip.dr * proj._strip_normal_dirs[k]
            + float(bem.Tp[k]) * strip.dr * proj._strip_chord_dirs[k]
        )
        m_ac = float(mp[k]) * strip.dr * proj._strip_moment_axis_sign[k] * SPAN
        if reference == "centroid":
            ref = strip.centroid
        elif reference == "ac":
            ref = strip.centroid - proj._strip_ac_offsets[k]
        else:
            raise ValueError(reference)
        st = _Strip(
            node_indices=strip.node_indices,
            r_center=strip.r_center,
            dr=strip.dr,
            centroid=np.asarray(ref, dtype=float),
            offsets=coords[strip.node_indices] - ref,
        )
        forces[strip.node_indices] = ForceProjector._distribute(st, f_strip, m_ac)
    return forces


def omega_of(forces: np.ndarray) -> float:
    f = np.zeros(n)
    f[0::6] = forces[:, 0]
    f[1::6] = forces[:, 1]
    f[2::6] = forces[:, 2]
    u = np.zeros(n)
    u[free] = np.asarray(spsolve(Kff, f[free]), dtype=float).ravel()
    return float(np.rad2deg(t._ring_kinematics(coords, u, tip)["omega"]))


for label, ref in (("about the mean centroid (production)", "centroid"), ("about the AC", "ac")):
    forces = project_about(ref)
    torque = float(np.cross(coords - coords[tip].mean(axis=0), forces).sum(axis=0) @ SPAN)
    print(f"{label:36s} tip omega {omega_of(forces):+8.4f} deg   realised torque {torque:+.6e} N.m")
print("reference: deck decoupled beam -1.7765 deg, Zhou coupled GEBT -3.60 deg (nose-down)")
