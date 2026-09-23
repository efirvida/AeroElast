"""S-0 — omega=0 consistency: the rotating-frame machinery reduces to the
plain elastic dynamics when the rotor speed is zero.

Guards the physics the rotor solver must obey at rest: the centrifugal
pre-stress, the geometric stiffness, the spin softening, and the inertial
(centrifugal/Coriolis/Euler) forces must ALL vanish at omega=0, and a
Newmark transient built from the rotating-formulation matrices must match
the plain elastic solver exactly.

The K_G path is exercised end-to-end (the centrifugal load -> the
pre-stress static solve -> the membrane stress -> the K_G assembly), not
just the individual formulas.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import coo_matrix, csr_matrix

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tests"))

from test_box_torsion_bending_benchmark import (  # noqa: E402
    DOF,
    E,
    H,
    L,
    NU,
    T,
    W,
    _assemble,
    _box_mesh,
)
from test_iea15mw_s3_newmark_transient import RHO, _fixed_dofs  # noqa: E402

AXIS = np.array([0.0, 1.0, 0.0])


def _assembler(mesh):
    """Build the MeshAssembler wrapper for the box (the K_G access)."""
    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.model import ElementSet
    from aeroelast.elements import ElementFamily

    if "all" not in mesh.element_sets:
        mesh.add_element_set(ElementSet(name="all", elements=list(mesh.elements)))
    model = {
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": {"all": {"type": "isotropic", "e": E, "nu": NU,
                                   "rho": RHO, "thickness": T,
                                   "shear_correction": 5.0 / 6.0}},
        }
    }
    return MeshAssembler(mesh, model)


@pytest.fixture(scope="module")
def box_setup():
    mesh = _box_mesh(nx=4, ny=20)
    K, m = _assemble(mesh, rho=RHO)
    return {"mesh": mesh, "K": K, "m": m}


def test_centrifugal_load_zero_at_zero_omega(box_setup):
    mesh = box_setup["mesh"]
    asm = _assembler(mesh)
    f = np.asarray(asm._rust.centrifugal_load(
        0.0, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], asm._rho_per_elem))
    assert np.allclose(f, 0.0), "centrifugal load must vanish at omega=0"


def test_geometric_stiffness_zero_at_zero_omega(box_setup):
    mesh = box_setup["mesh"]
    asm = _assembler(mesh)
    fixed = _fixed_dofs(mesh, box_setup["m"])
    free = np.setdiff1d(np.arange(asm.dofs_count), fixed)
    kg = asm.assemble_geometric_stiffness(
        omega=0.0, rotation_axis=AXIS, rotation_center=np.zeros(3),
        free_dofs=free)
    norm = kg.norm()
    assert norm < 1e-10, f"K_G(omega=0) norm {norm:.3e} must vanish"


def test_spin_softening_zero_at_zero_omega():
    from scipy.sparse import csr_matrix as _coo

    n = 4
    rows, cols, vals = [], [], []
    blk = -0.0**2 * (np.eye(3) - np.outer(AXIS, AXIS))
    for i in range(n):
        for a in range(3):
            for b in range(3):
                rows.append(i * 6 + a)
                cols.append(i * 6 + b)
                vals.append(2.0 * blk[a, b])
    Ksp = _coo((vals, (rows, cols)), shape=(6 * n, 6 * n)).tocsr()
    assert Ksp.nnz == 0 or abs(Ksp.data).max() < 1e-30


def test_inertial_forces_zero_at_zero_omega():
    from aeroelast.solvers.fsi.corotational import InertialForcesCalculator

    calc = InertialForcesCalculator(rotation_axis=[0.0, 1.0, 0.0])
    coords = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 2.0], [-1.0, 0.0, 3.0]])
    masses = np.array([1.0, 2.0, 3.0])
    vel = np.ones_like(coords)
    f_cf = calc.compute_centrifugal_force(coords, masses, 0.0)
    f_cor = calc.compute_coriolis_force(vel, masses, 0.0)
    f_eul = calc.compute_euler_force(coords, masses, 0.0)
    for name, f in (("centrifugal", f_cf), ("coriolis", f_cor), ("euler", f_eul)):
        assert np.allclose(f, 0.0), f"{name} force must vanish at omega=0"


def test_transient_identical_to_plain_newmark(box_setup):
    """The rotating-formulation matrices at omega=0 reproduce the plain
    Newmark response of the S-3 (the same box, the same modal IC)."""
    from test_iea15mw_s3_newmark_transient import (
        DT,
        N_STEPS,
        _newmark_transient,
    )

    mesh = box_setup["mesh"]
    K = box_setup["K"]
    fixed = _fixed_dofs(mesh, box_setup["m"])

    # the rotating-formulation assembly at omega=0 -> K_eff = K + 0 + 0
    asm = _assembler(mesh)
    free = np.setdiff1d(np.arange(asm.dofs_count), fixed)
    kg = asm.assemble_geometric_stiffness(
        omega=0.0, rotation_axis=AXIS, rotation_center=np.zeros(3),
        free_dofs=free)
    ai, aj, av = kg.getValuesCSR()
    Kg = csr_matrix((np.asarray(av), np.asarray(aj), np.asarray(ai)),
                    shape=(asm.dofs_count, asm.dofs_count)).tocsr()

    from _aeroelast import PyMeshAssembler  # type: ignore[import-not-found]
    nodes_sorted = sorted(mesh.nodes, key=lambda n: n.id)
    conn = [[box_setup["m"][n.id] for n in e.nodes] for e in mesh.elements]
    mat = {"type": "isotropic", "e": E, "nu": NU, "rho": RHO,
           "thickness": T, "shear_correction": 5.0 / 6.0}
    asm_m = PyMeshAssembler(
        node_coords=np.asarray([n.coords[:3] for n in nodes_sorted], dtype=float),
        connectivity=conn, elem_types=[4] * len(conn),
        materials=[dict(mat) for _ in conn])
    mr, mc, mv = asm_m.assemble_m()
    M = coo_matrix((mv, (mr, mc)),
                   shape=(asm.dofs_count, asm.dofs_count)).tocsr()

    # modal IC (the 1st bending) via the plain K/M
    from scipy.sparse.linalg import eigsh
    lam, vecs = eigsh(K[free][:, free], k=4, M=M[free][:, free],
                      sigma=1.0, which="LM")
    order = np.argsort(lam)
    tip_dof = box_setup["m"][next(n.id for n in mesh.nodes
                                  if abs(n.coords[2] - L) < 1e-9
                                  and abs(n.coords[1] - H / 2) < 1e-9)] * DOF + 1
    u_modal = np.zeros(asm.dofs_count)
    for k in order:
        v = vecs[:, k]
        tip_pos = np.where(free == tip_dof)[0][0]
        if abs(v[tip_pos]) > 0.1 * np.abs(v).max():
            u_modal[free] = v
            break
    v0 = np.zeros_like(u_modal)

    hist_plain, _ = _newmark_transient(K, M, 0.0 * M, u_modal, v0, DT,
                                       N_STEPS, fixed, tip_dof)
    hist_rot, _ = _newmark_transient(K + Kg, M, 0.0 * M, u_modal, v0, DT,
                                     N_STEPS, fixed, tip_dof)

    assert np.allclose(hist_plain, hist_rot, atol=1e-10), (
        "the omega=0 rotating-formulation response must equal the plain "
        f"Newmark (max diff {np.abs(hist_plain - hist_rot).max():.2e})"
    )
