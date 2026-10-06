"""Is the aerodynamic pitching moment applied with the right sense?

`_section_ends` picks the leading edge at the END WITH THE LARGEST projection on
`_strip_chord_dirs[k]`, so the resolved chord axis points from the trailing edge
towards the leading edge. The pitching moment in `bem_result.Mp` comes from the
polars, whose convention is the aerofoil one (chord from leading to trailing edge,
positive nose-up), and the projector applies it as `Mp * dr * span_hat` without
converting that sense. If the two conventions disagree the term is inverted.

The falsifiable test: two *different* realisations of the same requested section
moment - the minimum-norm field and the verified per-cell wall flow - must give tip
rotations of the same sign, because both are faithful to that moment. Today they
disagree: -0.8696 deg against +0.0773 deg. Flipping the sign of Mp should make
them agree if the sense is the defect.
"""

import dataclasses

import numpy as np
from _aeroelast import PyMeshAssembler
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

import tests.validation.blade.test_blade_iea15mw_validation as bv
import tests.validation.blade.test_blade_rated_twist as t
from aeroelast.core.mesh.entities import MeshElement, Node
from aeroelast.models.blade.model import Blade
from aeroelast.solvers.bem.engine import BEMSolver, BEMResult
from aeroelast.solvers.bem.force_projection import ForceProjector
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
solver = BEMSolver(aero, rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0)
bem = solver.compute(t.V_RATED, t.RPM_RATED, t.PITCH_RATED)
if bem.Mp is None:
    raise SystemExit("no pitching moment in the BEM result")
mp_reference = np.asarray(bem.Mp, dtype=float)

plain = ForceProjector(mesh, aero, span_direction=SPAN)
with_props = ForceProjector(mesh, aero, span_direction=SPAN, element_properties=props)


def total_applied_torque(result: BEMResult, projector: ForceProjector) -> float:
    """Sum over strips of Mp*dr and the aerodynamic-centre lever arm."""
    mp = result.Mp
    total = 0.0
    for k, strip in enumerate(projector._strips):
        f_strip = (
            float(result.Np[k]) * strip.dr * projector._strip_normal_dirs[k]
            + float(result.Tp[k]) * strip.dr * projector._strip_chord_dirs[k]
        )
        m_ac = float(mp[k]) * strip.dr if mp is not None else 0.0
        total += m_ac - float(np.cross(projector._strip_ac_offsets[k], f_strip) @ SPAN)
    return total


def omega_of(projector: ForceProjector, result: BEMResult) -> float:
    forces = projector.project(result)
    force = np.zeros(n)
    force[0::6] = forces[:, 0]
    force[1::6] = forces[:, 1]
    force[2::6] = forces[:, 2]
    u = np.zeros(n)
    u[free] = np.asarray(spsolve(Kff, force[free]), dtype=float).ravel()
    return float(np.rad2deg(t._ring_kinematics(coords, u, tip)["omega"]))


for label, result in (
    ("Mp as delivered", bem),
    ("Mp sign flipped", dataclasses.replace(bem, Mp=-mp_reference)),
):
    min_norm = omega_of(plain, result)
    multi = omega_of(with_props, result)
    agree = "AGREE" if np.sign(min_norm) == np.sign(multi) else "DISAGREE"
    torque = total_applied_torque(result, plain)
    print(
        f"{label:16s} applied torque {torque:+.6e} N.m | "
        f"min-norm {min_norm:+8.4f} deg | multi-cell {multi:+8.4f} deg | {agree}"
    )
