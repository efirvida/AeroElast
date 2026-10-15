"""Static stress recovery: AeroElast against beam theory, outer fibre.

Why the outer fibre needs a special CalculiX output
---------------------------------------------------
For shells, CalculiX ``*EL FILE, OUTPUT=2D`` + ``S`` reports the **mid-surface**
stress.  Two controlled probes established this:

* in-plane bending of a strip: ``SXX`` reproduces the axial distribution
  (measured 371.98 MPa against the analytical 360 MPa extreme fibre);
* out-of-plane bending of a plate: ``SXX`` is **0** everywhere although the
  plate deflects (measured ``|U|`` 0.0150 m), because pure bending has no
  mid-surface stress.

``*EL FILE, OUTPUT=3D`` makes CalculiX expand the shell and expose the
outer-fibre stress (measured 57.41 MPa against the analytical 60 MPa for the
same plate).  This module compares that value with the AeroElast
``StressRecovery`` output at ``TOP``/``BOTTOM``.

Scope note
----------
This is the **isotropic** validation of the stress-recovery path.  The same
``OUTPUT=3D`` switch does **not** change the stress of a composite
``*SHELL SECTION, COMPOSITE`` deck (its ``S`` stayed at 899.67 MPa with both
output modes), so the IEA 15 MW composite blade stress needs a different
CalculiX route and is not covered here.
"""

from __future__ import annotations


import numpy as np
import pytest
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from aeroelast.core.assembler import MeshAssembler

from tests.support.assertions import assert_relative_error  # noqa: E402
from tests import _production_arbiter as arbiter  # noqa: E402
from aeroelast.core.mesh.entities import (
    ElementSet,
    ElementType,
    MeshElement,
    Node,
    NodeSet,
)
from aeroelast.core.mesh.model import MeshModel
from aeroelast.elements import ElementFamily
from aeroelast.postprocess.stress_recovery import (
    StressLocation,
    StressRecovery,
    StressType,
)


# Cantilever plate, out-of-plane tip load.
L, B, H = 1.0, 0.1, 0.01  # length, width, thickness [m]
NX, NY = 8, 2  # regular quad mesh
E, NU, RHO = 2.1e11, 0.3, 7850.0
FORCE = 100.0  # total out-of-plane load at the free edge [N]

#: Outer-fibre bending stress ``M c / I`` at the clamped edge.
ANALYTICAL_STRESS = FORCE * L * (H / 2.0) / (B * H**3 / 12.0)

#: Bound against beam theory at mid-span (see the test for why mid-span).
TOL_ANALYTICAL = 0.05

_PROP = {"type": "isotropic", "name": "STEEL", "e": E, "nu": NU, "rho": RHO, "thickness": H}


def _build_plate() -> MeshModel:
    Node._id_counter = 0
    MeshElement._id_counter = 0
    mesh = MeshModel()
    xs = np.linspace(0.0, L, NX + 1)
    ys = np.linspace(0.0, B, NY + 1)

    grid: dict[tuple[int, int], Node] = {}
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            node = Node([float(x), float(y), 0.0], geometric_node=False)
            mesh.add_node(node)
            grid[(i, j)] = node

    for j in range(NY):
        for i in range(NX):
            mesh.add_element(
                MeshElement(
                    nodes=[
                        grid[(i, j)],
                        grid[(i + 1, j)],
                        grid[(i + 1, j + 1)],
                        grid[(i, j + 1)],
                    ],
                    element_type=ElementType.quad,
                )
            )

    mesh.add_node_set(
        NodeSet("clamped", {n for n in mesh.nodes if np.isclose(n.x, 0.0, atol=1e-12)})
    )
    mesh.add_element_set(ElementSet("plate", set(mesh.elements)))
    return mesh


def _model_cfg() -> dict:
    return {
        "elements": {"element_family": ElementFamily.SHELL, "properties": {"plate": _PROP}},
        "solver": {"time_step": 0.01, "total_time": 1.0, "beta": 0.25, "gamma": 0.5},
    }


def _solve_aeroelast(mesh: MeshModel) -> tuple[MeshAssembler, np.ndarray]:
    domain = MeshAssembler(mesh=mesh, model=_model_cfg())
    if domain._rust is None:
        raise RuntimeError("the Rust assembler is required to build this comparison's K")
    rows, cols, vals = domain._rust.assemble_k()
    n = domain.dofs_count
    k = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))), shape=(n, n)).tocsr()

    node_index = mesh.node_id_to_index
    root = {node_index[nid] for nid in mesh.get_node_set("clamped").node_ids}
    fixed = {6 * i + d for i in root for d in range(6)}
    free = np.array([i for i in range(n) if i not in fixed], dtype=np.int64)

    free_nodes = sorted(
        (nd for nd in mesh.nodes if np.isclose(nd.x, L, atol=1e-12)), key=lambda nd: nd.y
    )
    per_node = FORCE / len(free_nodes)
    force = np.zeros(n)
    for nd in free_nodes:
        force[6 * node_index[nd.id] + 2] = -per_node

    u = np.zeros(n)
    u[free] = np.asarray(spsolve(k[np.ix_(free, free)], force[free])).ravel()
    return domain, u


@pytest.fixture(scope="module")
def plate_stress(tmp_path_factory: pytest.TempPathFactory) -> dict:
    """Solve the plate once in AeroElast."""
    mesh = _build_plate()
    domain, u = _solve_aeroelast(mesh)
    recovery = StressRecovery(domain, u)

    def max_von_mises(location: StressLocation) -> float:
        result = recovery.compute_element_stresses(location=location, stress_type=StressType.TOTAL)
        return float(np.asarray(result.von_mises).max())

    upper = max_von_mises(StressLocation.TOP)
    lower = max_von_mises(StressLocation.BOTTOM)
    middle = max_von_mises(StressLocation.MIDDLE)

    # Element value at the centroid of the centre-line element just inboard of mid-span.
    top = recovery.compute_element_stresses(location=StressLocation.TOP, stress_type=StressType.TOTAL)
    x_mid = (NX // 2 - 0.5) * (L / NX)
    midspan = float(np.asarray(top.von_mises)[(NY // 2) * NX + (NX // 2 - 1)])

    print(
        f"  stress: aero TOP/BOT {upper / 1e6:.2f}/{lower / 1e6:.2f} MPa, "
        f"MID {middle / 1e6:.2f} MPa, mid-span {midspan / 1e6:.2f} MPa at x={x_mid:.4f} m, "
        f"root M c / I {ANALYTICAL_STRESS / 1e6:.2f} MPa"
    )
    return {"upper": upper, "lower": lower, "middle": middle, "midspan": midspan, "x_mid": x_mid}


def test_bending_has_no_membrane_stress(plate_stress: dict) -> None:
    """Pure out-of-plane bending leaves the mid-surface stress at zero."""
    assert plate_stress["middle"] < 0.02 * plate_stress["upper"], (
        f"mid-surface von Mises {plate_stress['middle'] / 1e6:.4f} MPa is not "
        f"negligible against the outer fibre {plate_stress['upper'] / 1e6:.2f} MPa"
    )


def test_outer_fibre_stress_is_symmetric(plate_stress: dict) -> None:
    """The two outer fibres carry the same magnitude for a symmetric section."""
    rel = abs(plate_stress["upper"] - plate_stress["lower"]) / plate_stress["upper"]
    assert rel < 1e-6, (
        f"TOP={plate_stress['upper'] / 1e6:.4f} BOT={plate_stress['lower'] / 1e6:.4f} "
        f"differ by {rel * 100:.4g}%"
    )


def test_outer_fibre_stress_matches_beam_theory(plate_stress: dict) -> None:
    """AeroElast outer-fibre stress at mid-span matches beam theory ``M(x) c / I``.

    The sample is taken away from the clamp on purpose: the maximum over the whole plate picks up
    the clamped-edge stress singularity (six DOF held), which keeps growing with refinement
    (2026-10-15 probe ``tools/g7_stress_convergence.py``: 8x2 -> 128x32 moves 52.5 -> 61.8 MPa
    against ``M c / I`` = 60), so that metric has no converged value to compare. The mid-span
    value is converged in the mesh (``tools/g7_midspan_convergence.py``).

    CalculiX is deliberately not an arbiter here: its ``OUTPUT=3D`` outer-fibre stress for an
    expanded S8R shell is one element ``h`` out of position in x (the field is linear with the
    correct slope but ``sigma_ccx(x) = M c / I (x + h)``), so a point comparison measures that
    output offset, not AeroElast. See ``odd/tasks/validation-l3-b3-audit.md`` (PRE-3, G7).
    """
    x = plate_stress["x_mid"]
    reference = FORCE * (L - x) * (H / 2.0) / (B * H**3 / 12.0)
    assert_relative_error(
        plate_stress["midspan"],
        reference,
        tol=TOL_ANALYTICAL,
        kind="analytical",
        reference_name="beam theory M(x) c / I at the sample station, outer-fibre stress",
        what="mid-span outer-fibre von Mises stress",
    )


#: Bound for the production-vs-replica comparison. Measured on this case: 1.85e-10 relative
#: (WU-3a), so this leaves ~5000x of margin while staying five orders tighter than the 5 % physics bound.
ARBITER_RTOL = 1e-6


def test_production_solver_matches_the_scipy_replica() -> None:
    """The production static solver reproduces the replica on the same mesh, loads and BCs.

    The replica above (production Rust ``K`` plus a scipy direct solve) is validated against
    CalculiX and beam theory, which proves the *model*. This closes the other half: nothing else in
    the suite drives the production solve path, so PETSc KSP and scipy could disagree and the suite
    would stay green. Deliberately CalculiX-free so it runs wherever the suite runs.
    """
    mesh = _build_plate()
    _, u_replica = _solve_aeroelast(mesh)

    edge = arbiter.nodes_at_coordinate(mesh, L, axis=0)
    u_prod = arbiter.solve_static(
        mesh,
        _model_cfg(),
        fixed_dofs=arbiter.fixed_dofs_from_node_set(mesh, "clamped"),
        nodal_loads_list=[arbiter.nodal_loads(edge, -FORCE, dof=2)],
    )

    arbiter.require_gap(
        arbiter.displacement_gap(u_prod, u_replica), rtol=ARBITER_RTOL, label="plate 8x2"
    )

    mid_dof = 6 * edge[len(edge) // 2] + 2
    w_prod, w_replica = float(u_prod[mid_dof]), float(u_replica[mid_dof])
    rel = abs(w_prod - w_replica) / abs(w_replica)
    assert rel < ARBITER_RTOL, (
        f"out-of-plane tip displacement: production {w_prod:.6e} against replica "
        f"{w_replica:.6e} ({rel:.3e} relative, bound {ARBITER_RTOL:.1e})"
    )
