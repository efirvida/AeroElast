"""Which way must the BEM tangential force point?  Static FSI fixed point for both senses (T7b).

A blade moves toward its leading edge, and a positive Tp (positive torque, positive power)
acts along that motion.  On the official IEA-15MW mesh the projector measures the leading
edge at +x and the solid's omega x r is +x, so the physical ``tangential_direction`` is +x.
Every production fluid YAML sets -x.  This script iterates BEM on the deformed geometry
against K + K_G + K_SP (rated, rigid-blade start, relaxation 0.5) for both senses and
prints the flexible thrust, flap tip and the elastic tip twist the BEM sees.

Measured 2026-10-12 (odd/tasks/coupled-divergence.md T7b): -x gives +5.55 % thrust,
16.6 m flap, -0.86 deg (nose-up); +x gives -21.6 %, 10.0 m, +6.28 deg (nose-down).
"""
import sys
import time
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


Ks = K + Kg + Ksp
cfg_prod = dict(fx.BEM_CONFIG, span_direction=[0, 0, 1], normal_direction=[0, 1, 0],
                tangential_direction=[-1, 0, 0])
part_p = BEMFSIParticipant(mesh, blade_aero, cfg_prod)
fn_p, _ = part_p._compute_forces(np.zeros((len(coords), 3)))
f_p = np.zeros(n)
for d in range(3):
    f_p[d::6] = fn_p[:, d]
dr = np.array([s.dr for s in part._projector._strips])
cfg_px = dict(cfg_prod, tangential_direction=[1, 0, 0])
for name, cfg in (("tangential -x (production)", cfg_prod), ("tangential +x (motion/LE)", cfg_px)):
    pt = BEMFSIParticipant(mesh, blade_aero, cfg)
    U = np.zeros((len(coords), 3))
    t0 = time.time()
    hist = []
    for it in range(14):
        fn, res = pt._compute_forces(U)
        T = 3 * float(np.sum(np.asarray(res.Np)[:len(dr)] * dr))
        hist.append(T / 1e6)
        if it == 0:
            T0 = T
        f = np.zeros(n)
        for d in range(3):
            f[d::6] = fn[:, d]
        Un = solve(Ks, f + f_cf)
        U = Un if it == 0 else 0.5 * U + 0.5 * Un
    _, tw = pt._compute_deformed_geometry(U)
    k = int(np.argmax(coords[:, 2]))
    print(f"[{name}] rigid T {T0/1e6:.4f} MN  flexible T {T/1e6:.4f} MN ({100*(T/T0-1):+.2f}%)  tip uy {U[k,1]:+.2f} m  tip elastic twist(BEM) {np.degrees(tw[-1]-pt._ref_twist[-1]):+.2f} deg  ({time.time()-t0:.0f}s)")
    print("  last iterations [MN]:", " ".join(f"{x:.4f}" for x in hist[-3:]))
