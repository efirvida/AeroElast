"""Arbitrate the blade's twist against CalculiX under the identical nodal load.

The shell gives +8.1048 deg of tip section rotation under the production
aerodynamic load while the deck's decoupled 1D beam gives -1.7765 deg for the same
sectional moment, a 4.56x amplification. That number decides whether the extra
twist is real 3D physics or our shell's load path.

The arbiter is an independent FE code fed **the exact same nodal force vector** and
clamped the same way: if CalculiX lands near +8.1048 deg the shell is validated
with that load and the residual gap is the aero side; if it does not, the load path
is the defect. A few percent is expected and not a failure - CalculiX runs S8R
against our linear MITC4 - while a sign flip or a factor of several would be.

The per-node field travels through `write_ccx_mesh(load_field=...)`, which emits one
`*CLOAD` per non-zero component on the mesh node's own deck label; the FRD's node
numbers are the mesh index plus one, because that is how the writer labels `*NODE`.
"""

import argparse
import shutil
import tempfile
from pathlib import Path

import numpy as np
from _aeroelast import PyMeshAssembler
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

import tests.validation.blade.test_blade_iea15mw_validation as bv
import tests.validation.blade.test_blade_rated_twist as t
from aeroelast.core.mesh.entities import MeshElement, Node
from aeroelast.core.mesh.io.writers import write_ccx_mesh
from aeroelast.models.blade.model import Blade
from aeroelast.solvers.bem.engine import BEMSolver
from aeroelast.solvers.bem.force_projection import ForceProjector
from tests.support.ccx_io import fail_ccx, parse_frd_disp, run_ccx
from tests.support.openfast_bem import build_blade_aero_from_aerodyn

SPAN = np.array([0.0, 0.0, 1.0])

_parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
_parser.add_argument(
    "--element-size",
    type=float,
    default=1.0,
    help="blade mesh size in metres (default 1.0; the arbitration gets sharper at 0.5)",
)
_args = _parser.parse_args()

Node._id_counter = 0
MeshElement._id_counter = 0
model = Blade(str(t.YAML), element_size=_args.element_size)
model.generate_mesh()
mesh = model.mesh
if mesh is None:
    raise SystemExit("no mesh")
props = model.get_element_properties()
coords = mesh.coords_array
n_nodes = coords.shape[0]

assembler = PyMeshAssembler.from_model(
    bv._to_rust_mesh(mesh, props), props, list(bv.SPAN_DIRECTION), None
)
n_dof = assembler.dofs_count
rows, cols, vals = assembler.assemble_k()
K = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))), shape=(n_dof, n_dof)).tocsr()
root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
free = np.array([i for i in range(n_dof) if i not in {6 * r + d for r in root for d in range(6)}])
Kff = K[np.ix_(free, free)]

phys = t._physical_stations(coords)
tip = np.where(np.abs(coords[:, 2] - phys[-1]) < t.STATION_GAP_TOLERANCE)[0]

blade_aero = build_blade_aero_from_aerodyn(t.AD_PRIMARY)
bem = BEMSolver(blade_aero, rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0).compute(
    t.V_RATED, t.RPM_RATED, t.PITCH_RATED
)
forces = ForceProjector(
    mesh, blade_aero, span_direction=SPAN, element_properties=props
).project(bem)

# ---- AeroElast: the same vector through its own clamped solve -----------------
f = np.zeros(n_dof)
f[0::6] = forces[:, 0]
f[1::6] = forces[:, 1]
f[2::6] = forces[:, 2]
u = np.zeros(n_dof)
u[free] = np.asarray(spsolve(Kff, f[free]), dtype=float).ravel()
omega_ae = float(np.rad2deg(t._ring_kinematics(coords, u, tip)["omega"]))

# ---- CalculiX: the same vector as one *CLOAD per node -------------------------
node_ids = [nd.id for nd in mesh.nodes]
field = {
    node_ids[i]: [float(forces[i, 0]), float(forces[i, 1]), float(forces[i, 2])]
    for i in range(n_nodes)
    if np.any(forces[i])
}
ccx_bin = shutil.which("ccx")
if ccx_bin is None:
    raise SystemExit("ccx not found in PATH")

with tempfile.TemporaryDirectory(prefix="ccx_twist_arbitration_") as tmp:
    deck = Path(tmp) / "blade.inp"
    write_ccx_mesh(
        mesh,
        str(deck),
        properties=props,
        boundary_nodeset="RootNodes",
        solver_type="LinearStatic",
        quadratic=True,
        span_direction=(0.0, 0.0, 1.0),
        load_field=field,
    )
    proc = run_ccx(deck, ccx_bin)
    if proc.returncode != 0:
        fail_ccx(proc, deck)
    frd = deck.with_suffix(".frd")
    if not frd.exists():
        raise SystemExit("CCX produced no FRD")
    disp = parse_frd_disp(frd, list(range(1, n_nodes + 1)))

print(f"nodes {n_nodes}  dofs {n_dof}  loaded nodes {len(field)}")
print(f"AeroElast (MITC4, same vector)  tip section rotation {omega_ae:+9.4f} deg")
u_ccx = np.zeros(n_dof)
matched = 0
for k in range(n_nodes):
    d = disp.get(k + 1)
    if d is None:
        continue
    u_ccx[6 * k : 6 * k + 3] = d
    matched += 1
if matched == 0:
    raise SystemExit("no FRD displacement matched the mesh nodes")
omega_ccx = float(np.rad2deg(t._ring_kinematics(coords, u_ccx, tip)["omega"]))
print(f"CalculiX  (S8R,   same vector)  tip section rotation {omega_ccx:+9.4f} deg")
print(f"matched {matched} of {n_nodes} nodes")
rel = abs(omega_ccx - omega_ae) / max(abs(omega_ae), 1e-30)
same = "SAME SIGN" if np.sign(omega_ccx) == np.sign(omega_ae) else "OPPOSITE SIGNS"
print(f"difference {rel:.2%}   {same}")
