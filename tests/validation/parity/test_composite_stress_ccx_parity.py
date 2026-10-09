"""Composite outer-fibre ply stress: AeroElast against CalculiX S8R and CLT.

Why this row closes issue #27
-----------------------------
Issue #27 assumed CalculiX cannot judge a composite shell because
``*SHELL SECTION, COMPOSITE`` ignores ``OUTPUT=3D``.  Measured on 2026-10-09
(``odd/tasks/composite-ply-stress-recovery.md``), both halves of that premise
are wrong, and the missing half was ours: the recovery in the tree returned the
thickness-mean equivalent stress (E1), not a ply stress.  The ply-resolved
recovery landed in T2/T2b/T3 (``95bfbab``, ``208f220``, ``bd4cc05``); this
module is the converged-vs-converged row against CalculiX with the closed form
as the second, independent reference.

E2 - what the CalculiX composite-shell FRD value is
---------------------------------------------------
The value CalculiX writes for a composite shell is the **outer-fibre ply
stress**, not a thickness mean.  The card expands the section: the FRD carries
a midsurface node (zero stress) plus one node per ply boundary and ply midpoint
(interface nodes duplicated per ply), and the outer-fibre value is the topmost
surface node.  On this deck
that value follows the closure ply (0.3150 MPa at 8x2 and 0.3148 MPa at 16x4
for ``[90/0]s``, 3.6820 / 3.6853 MPa for ``[0/90]s``), while a thickness mean
would be 2.000 MPa for both stacks.  (The feature document's original probe
quotes 0.2991 / 0.2994 and 3.6984 / 3.7005 for the same case; this deck
reproduces the mechanism -- the value follows the closure ply and is not the
mean -- but not those exact splits, which is a finding reported with this
row.)

E3 - ``OUTPUT=3D`` is inert for composite sections
--------------------------------------------------
Measured on a bending-dominated ``[0/90/90/0]`` cantilever with the same deck:
``OUTPUT=2D`` and ``OUTPUT=3D`` coincide to a relative difference of
0.000e+00 (max ``|SZZ|`` 595.588 MPa, max von Mises 574.521 MPa both ways).
This is why this module leaves the deck at ``OUTPUT=2D``; the isotropic sibling
(``test_shell_stress_ccx_parity.py``) still needs ``OUTPUT=3D`` because an
isotropic shell reports the midsurface stress there.  CalculiX requires S8R/S6
for ``*SHELL SECTION, COMPOSITE``; linear S4/S3 does not accept the card.

E4 - a boundary singularity, not the answer
-------------------------------------------
On this case the global maximum ``SXX`` is a load-introduction singularity at
the loaded free edge and **grows with refinement** (measured 5.8838 MPa at 8x2
and 6.5803 MPa at 16x4 for ``[0/90]s``), while the free-field centre stays at
3.685 MPa.  The row therefore reads the **free-field centre**, never the global
maximum, or it would validate a mesh-dependent singularity.  The mesh-stability
of the centre and the non-stability of the maximum are asserted explicitly
below.

Reference independence
----------------------
The closed-form reference is re-implemented in this file from first principles
-- own ``Q``, own ``Qbar(theta)``, own ``A`` (Jones 1999 / Reddy 2004) -- and
never imports ``aeroelast.core.laminate``.  The CalculiX reference is the value
CalculiX writes on its own S8R mesh.  The two are declared separately through
``assert_relative_error`` (``kind="code"`` and ``kind="analytical"``) so the
validation store reads each comparison.

Arbiter
-------
There is deliberately no production-vs-replica arbiter test here.  The
isotropic sibling already drives the production solve path
(``test_production_solver_matches_the_scipy_replica``), and this module reuses
the same ``_production_arbiter.solve_static`` production solve, so a second
arbiter would only duplicate that claim on a different geometry.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests import _production_arbiter as arbiter  # noqa: E402
from tests.support.assertions import assert_relative_error  # noqa: E402
from tests.support.ccx_io import parse_frd_stress  # noqa: E402
from tests.conftest import ccx_bin_or_skip  # noqa: E402

from aeroelast.core.assembler import MeshAssembler  # noqa: E402
from aeroelast.core.laminate import create_laminate_from_angles  # noqa: E402
from aeroelast.core.material import OrthotropicMaterial  # noqa: E402
from aeroelast.core.mesh.entities import (  # noqa: E402
    ElementSet,
    ElementType,
    MeshElement,
    Node,
    NodeSet,
)
from aeroelast.core.mesh.io.writers import write_ccx_mesh  # noqa: E402
from aeroelast.core.mesh.model import MeshModel  # noqa: E402
from aeroelast.core.properties import CompositeShellProperty  # noqa: E402
from aeroelast.elements import ElementFamily  # noqa: E402
from aeroelast.postprocess.stress_recovery import (  # noqa: E402
    StressLocation,
    StressRecovery,
    StressType,
)

# ============================================================================
# Geometry, material, load
# ============================================================================

L = 1.0  # length along x [m]
B = 0.1  # width along y [m]
N_PLY = 4
PLY_THICKNESS = 0.00125  # [m]
H = N_PLY * PLY_THICKNESS  # 0.005 m total

FORCE = 1000.0  # total axial resultant at the free edge, +x [N]
N_X = FORCE / B  # membrane resultant per unit width [N/m] = 1e4

E1, E2, G12, G23, NU12, RHO = 120e9, 10e9, 5e9, 3e9, 0.3, 1000.0

#: Layups under test, ply angles bottom to top.  Symmetric, so B = 0 and the
#: closed form is pure membrane: sigma(z) is constant within each ply, and the
#: outer (TOP) ply is the last angle in the list.
LAYUPS: dict[str, list[float]] = {
    "sym_0_90s": [0.0, 90.0, 90.0, 0.0],
    "sym_90_0s": [90.0, 0.0, 0.0, 90.0],
}

#: Each code solves its own mesh sequence: our QUAD4 strip and CalculiX's S8R
#: strip (``quadratic=True``) at the same two refinements.
MESHES: dict[str, tuple[int, int]] = {"8x2": (8, 2), "16x4": (16, 4)}

#: Bound in the style of the isotropic sibling (``TOL_CCX = 0.05`` there).
#: It is a bound, not a convergence claim.  The measured ours-vs-CCX gap on this
#: deck is 0.29 % ([0/90]s) and 3.41 % ([90/0]s) at 16x4, and the [90/0]s gap
#: barely moves with refinement (3.48 % -> 3.41 %), so agreement below ~3 % on
#: that stack is not demonstrated.  The whole of the gap sits on the CalculiX
#: side (ours matches the closed form to <= 0.0006 %), between two different
#: discretisations: our QUAD4 strip against S8R quadratic with two elements
#: across the width.  5 % is ~1.5x the measured worst case.
TOL_CCX = 0.05
TOL_ANALYTICAL = 0.05

#: Published CLT targets for the outer-fibre sigma_xx (contract, 2026-10-09).
TARGET_OUTER_0_DEG = 3.696e6  # outer 0 deg ply of [0/90]s [Pa]
TARGET_OUTER_90_DEG = 0.304e6  # outer 90 deg ply of [90/0]s [Pa]

#: Measured CCX binary on this tree: ``ccx -v`` -> "This is Version 2.20".
#: Do not copy the isotropic sibling's "2.23": that string is stale here.
CCX_VERSION = "2.20"

#: The free-field centre must not move with refinement; the measured drift is
#: <= 0.09 % (CalculiX) and <= 0.03 % (ours) on both stacks, so 2 % is a
#: stability guard, not a physics tolerance.
TOL_FIELD_STABILITY = 0.02


# ============================================================================
# Independent CLT reference (Jones 1999 / Reddy 2004), re-implemented in-file
# ============================================================================


def _hand_q(E1: float, E2: float, G12: float, nu12: float) -> np.ndarray:
    """Plane-stress reduced stiffness of one ply in its own material axes."""
    nu21 = nu12 * E2 / E1
    denom = 1.0 - nu12 * nu21
    q11 = E1 / denom
    q22 = E2 / denom
    q12 = nu12 * E2 / denom
    q66 = G12
    return np.array([[q11, q12, 0.0], [q12, q22, 0.0], [0.0, 0.0, q66]])


def _hand_qbar(q: np.ndarray, theta_deg: float) -> np.ndarray:
    """Transformed reduced stiffness ``Qbar(theta)`` (explicit component form).

    Standard reduced-stiffness transformation (Jones 1999 eq. 2.72ff / Reddy
    2004): the 0-deg limit reproduces ``q`` and the 90-deg limit swaps q11/q22.
    """
    c = np.cos(np.deg2rad(theta_deg))
    s = np.sin(np.deg2rad(theta_deg))
    q11, q12, q22, q66 = q[0, 0], q[0, 1], q[1, 1], q[2, 2]
    return np.array(
        [
            [
                q11 * c**4 + 2 * (q12 + 2 * q66) * s * s * c * c + q22 * s**4,
                (q11 + q22 - 4 * q66) * s * s * c * c + q12 * (c**4 + s**4),
                (q11 - q12 - 2 * q66) * c**3 * s - (q22 - q12 - 2 * q66) * c * s**3,
            ],
            [
                (q11 + q22 - 4 * q66) * s * s * c * c + q12 * (c**4 + s**4),
                q11 * s**4 + 2 * (q12 + 2 * q66) * s * s * c * c + q22 * c**4,
                (q11 - q12 - 2 * q66) * c * s**3 - (q22 - q12 - 2 * q66) * c**3 * s,
            ],
            [
                (q11 - q12 - 2 * q66) * c**3 * s - (q22 - q12 - 2 * q66) * c * s**3,
                (q11 - q12 - 2 * q66) * c * s**3 - (q22 - q12 - 2 * q66) * c**3 * s,
                (q11 + q22 - 2 * q12 - 2 * q66) * s * s * c * c + q66 * (c**4 + s**4),
            ],
        ]
    )


def _hand_clt_ply_stresses(angles: list[float], n_x: float) -> tuple[list[float], float]:
    """Per-ply ``sigma_xx`` (global x axes, bottom to top) and the thickness mean.

    Classical lamination theory for a symmetric laminate under a pure membrane
    resultant ``N = [n_x, 0, 0]``: ``eps0 = A^-1 N`` and ``sigma_k =
    Qbar(theta_k) . eps0`` (constant within each ply).  ``A = sum Qbar_k
    (z_top - z_bottom)`` is assembled from the same first-principles ``Qbar``
    -- nothing is imported from ``aeroelast.core.laminate``.
    """
    q = _hand_q(E1, E2, G12, NU12)
    n = len(angles)
    z = np.linspace(-0.5 * n * PLY_THICKNESS, 0.5 * n * PLY_THICKNESS, n + 1)

    A = np.zeros((3, 3))
    qbars = []
    for k, theta in enumerate(angles):
        qb = _hand_qbar(q, theta)
        qbars.append(qb)
        A += qb * (z[k + 1] - z[k])

    eps0 = np.linalg.solve(A, np.array([n_x, 0.0, 0.0]))
    sigma_xx = [float((qb @ eps0)[0]) for qb in qbars]
    mean = float((A @ eps0)[0]) / H  # = N_x / h by equilibrium
    return sigma_xx, mean


# ============================================================================
# Production-path model (mirrors test_composite_ply_stress_parity.py)
# ============================================================================


def _material() -> OrthotropicMaterial:
    # OrthotropicMaterial convention: G = (G12, G23, G31), nu = (nu12, nu23, nu31).
    # Only E1, E2, G12, nu12 enter the CLT reference; G13 = G23 and
    # nu23 = nu12 are stand-ins for the unused out-of-plane terms.
    return OrthotropicMaterial(
        name="cfrp",
        E=(E1, E2, E2),
        G=(G12, G23, G23),
        nu=(NU12, NU12, 0.0),
        rho=RHO,
    )


def _laminate(angles: list[float]):
    return create_laminate_from_angles(
        _material(), PLY_THICKNESS, angles, shear_correction_factor=5.0 / 6.0
    )


def _model_cfg(laminate) -> dict:
    return {
        "elements": {
            "element_family": ElementFamily.SHELL,
            "properties": {"plate": CompositeShellProperty(laminate=laminate)},
        },
        "solver": {"time_step": 0.01, "total_time": 1.0, "beta": 0.25, "gamma": 0.5},
    }


def _build_mesh(nx: int, ny: int) -> MeshModel:
    """Flat cantilever strip: x is the span (clamped at x=0), y the width."""
    Node._id_counter = 0
    MeshElement._id_counter = 0

    mesh = MeshModel()
    xs = np.linspace(0.0, L, nx + 1)
    ys = np.linspace(0.0, B, ny + 1)
    grid: dict[tuple[int, int], Node] = {}
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            node = Node([float(x), float(y), 0.0], geometric_node=False)
            mesh.add_node(node)
            grid[(i, j)] = node

    for j in range(ny):
        for i in range(nx):
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
    mesh.add_node_set(
        NodeSet("free_face", {n for n in mesh.nodes if np.isclose(n.x, L, atol=1e-12)})
    )
    mesh.add_element_set(ElementSet("plate", set(mesh.elements)))
    return mesh


def _solve_production(mesh: MeshModel, model_cfg: dict) -> np.ndarray:
    """Production static solve (Rust assembly + PETSc KSP) on the membrane case."""
    edge = arbiter.nodes_at_coordinate(mesh, L, axis=0)
    return arbiter.solve_static(
        mesh,
        model_cfg,
        fixed_dofs=arbiter.fixed_dofs_from_node_set(mesh, "clamped"),
        nodal_loads_list=[arbiter.nodal_loads(edge, FORCE, dof=0)],
    )


def _centre_element_index(domain: MeshAssembler) -> int:
    """Index (into ``domain.elements``) of the element whose centroid is nearest (L/2, B/2)."""
    node_pos = {n.id: (n.x, n.y, n.z) for n in domain.nodes}
    target = np.array([L / 2.0, B / 2.0, 0.0])
    best_i, best_d = -1, np.inf
    for i, elem in enumerate(domain.elements):
        centroid = np.mean([node_pos[nid] for nid in elem.node_ids], axis=0)
        d = float(np.linalg.norm(centroid - target))
        if d < best_d:
            best_i, best_d = i, d
    assert best_i >= 0, "no elements in the assembled domain"
    return best_i


# ============================================================================
# CalculiX S8R composite deck
# ============================================================================


def _parse_frd_coords(frd_file: Path) -> dict[int, np.ndarray]:
    """Parse the FRD nodal-coordinate block by fixed columns (``2C`` record)."""
    coords: dict[int, np.ndarray] = {}
    in_coords = False
    with open(frd_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("2C"):
                in_coords = True
                continue
            if not in_coords:
                continue
            if stripped.startswith("-3"):
                break
            if not stripped.startswith("-1"):
                continue
            try:
                node_id = int(line[3:13])
                coords[node_id] = np.array(
                    [float(line[13:25]), float(line[25:37]), float(line[37:49])],
                    dtype=float,
                )
            except ValueError:
                continue
    return coords


def _run_ccx(tmp_path: Path, mesh: MeshModel, laminate) -> tuple[float, float]:
    """Write and run the S8R composite deck; return (free-field SXX, max SXX).

    ``write_ccx_mesh(..., quadratic=True)`` emits ``*SHELL SECTION, COMPOSITE``
    with per-ply data and converts the QUAD4 mesh to S8R (CalculiX requires
    S8R/S6 for a composite section).  The resultant axial load is distributed
    over the free-face node set by the writer.  ``OUTPUT=2D`` is left as written
    because for a composite section ``OUTPUT=3D`` is inert (E3).
    """
    ccx_bin = ccx_bin_or_skip()
    inp = tmp_path / "composite.inp"
    write_ccx_mesh(
        mesh,
        str(inp),
        properties={"plate": CompositeShellProperty(laminate=laminate)},
        boundary_nodeset="clamped",
        load_nodeset="free_face",
        load_vector=[FORCE, 0.0, 0.0],
        solver_type="LinearStatic",
        quadratic=True,
    )

    proc = subprocess.run([str(ccx_bin), inp.stem], cwd=tmp_path, capture_output=True, text=True)
    if proc.returncode != 0:
        pytest.fail(
            f"CalculiX failed (rc={proc.returncode}):\n{proc.stdout[-2000:]}\n{proc.stderr[-1000:]}"
        )

    frd = tmp_path / "composite.frd"
    stress = parse_frd_stress(frd)
    coords = _parse_frd_coords(frd)
    common = [nid for nid in coords if nid in stress]
    assert common, "no FRD node carries both coordinates and stress"

    target = np.array([L / 2.0, B / 2.0, H / 2.0])
    # The card expands the section: the midsurface node carries zero and the
    # outer-fibre value lives on the topmost surface node.  Target +H/2 (TOP).
    free_node = min(common, key=lambda nid: float(np.linalg.norm(coords[nid] - target)))
    free_sxx = float(stress[free_node][0])
    max_sxx = max(float(v[0]) for v in stress.values())
    return free_sxx, max_sxx


# ============================================================================
# Module fixture: both codes on both meshes, per laminate
# ============================================================================


@pytest.fixture(scope="module")
def cases(tmp_path_factory: pytest.TempPathFactory) -> dict:
    """Solve each layup through the production path and CalculiX on both meshes."""
    out: dict[str, dict] = {}
    for name, angles in LAYUPS.items():
        laminate = _laminate(angles)
        sigma_ref, mean_ref = _hand_clt_ply_stresses(angles, N_X)
        meshes: dict[str, dict] = {}
        for mname, (nx, ny) in MESHES.items():
            mesh = _build_mesh(nx, ny)
            model_cfg = _model_cfg(laminate)

            u = _solve_production(mesh, model_cfg)
            domain = MeshAssembler(mesh=mesh, model=model_cfg)
            recovery = StressRecovery(domain, u)
            idx = _centre_element_index(domain)
            res = recovery.compute_element_stresses(
                location=StressLocation.TOP, stress_type=StressType.TOTAL
            )
            assert res.sigma_xx.shape[0] == len(list(domain.elements)), (
                "compute_element_stresses must return one row per assembled element"
            )
            ours_top = float(res.sigma_xx[idx])
            ours_max = float(np.max(res.sigma_xx))

            workdir = tmp_path_factory.mktemp(f"ccx_{name}_{mname}")
            ccx_free, ccx_max = _run_ccx(workdir, mesh, laminate)

            meshes[mname] = {
                "ours_top": ours_top,
                "ours_max": ours_max,
                "ccx_free": ccx_free,
                "ccx_max": ccx_max,
            }
            print(
                f"  {name} {mname}: ours TOP {ours_top / 1e6:.4f} MPa, "
                f"ours max {ours_max / 1e6:.4f} MPa, "
                f"ccx free {ccx_free / 1e6:.4f} MPa, ccx max {ccx_max / 1e6:.4f} MPa"
            )
        out[name] = {
            "angles": angles,
            "ref_outer": float(sigma_ref[-1]),
            "ref_mean": mean_ref,
            "meshes": meshes,
        }
        print(
            f"  {name}: CLT outer {sigma_ref[-1] / 1e6:.4f} MPa, "
            f"thickness mean {mean_ref / 1e6:.4f} MPa"
        )
    return out


# ============================================================================
# The capability comparison (this is the row)
# ============================================================================


@pytest.mark.parametrize("name", sorted(LAYUPS))
def test_outer_fibre_stress_matches_ccx_and_analytical(cases: dict, name: str) -> None:
    """Converged free-field outer-fibre ply stress against CalculiX S8R and CLT.

    The measured quantity is the 16x4 free-field centre TOP ``sigma_xx`` from
    the production path.  It is compared against two independent references:
    the value CalculiX writes for the same physical point (``kind="code"``,
    outer-fibre ply stress of ``*SHELL SECTION, COMPOSITE``) and the CLT
    closed form re-implemented in this file (``kind="analytical"``).  The two
    references are declared separately, one comparison each.
    """
    case = cases[name]
    fine = case["meshes"]["16x4"]
    ours = fine["ours_top"]
    ccx = fine["ccx_free"]
    ref = case["ref_outer"]
    coarse = case["meshes"]["8x2"]
    print(
        f"  {name}: 8x2 ours-vs-ccx "
        f"{abs(coarse['ours_top'] - coarse['ccx_free']) / coarse['ccx_free'] * 100:.4f}%, "
        f"16x4 ours-vs-ccx {abs(ours - ccx) / ccx * 100:.4f}%, "
        f"ours-vs-CLT {abs(ours - ref) / ref * 100:.4f}%"
    )
    assert_relative_error(
        ours,
        ccx,
        tol=TOL_CCX,
        kind="code",
        reference_name=(
            f"CalculiX {CCX_VERSION} S8R, *SHELL SECTION, COMPOSITE, outer-fibre ply stress"
        ),
        what=f"{name} TOP outer-fibre ply sigma_xx",
    )
    assert_relative_error(
        ours,
        ref,
        tol=TOL_ANALYTICAL,
        kind="analytical",
        reference_name=(
            "CLT closed form re-implemented in this test (Jones 1999 / Reddy 2004), "
            "outer-fibre ply sigma_xx"
        ),
        what=f"{name} TOP outer-fibre ply sigma_xx",
    )


# ============================================================================
# Guards (machinery, not rows)
# ============================================================================


def test_reference_reproduces_published_clt_targets() -> None:
    """The in-file closed form reproduces the published 3.696 / 0.304 MPa.

    Pins the hand reference to the contract's target values so a transcription
    error in the reference cannot masquerade as a solver result.
    """
    sigma_0, mean_0 = _hand_clt_ply_stresses(LAYUPS["sym_0_90s"], N_X)
    sigma_90, mean_90 = _hand_clt_ply_stresses(LAYUPS["sym_90_0s"], N_X)
    outer_0 = sigma_0[-1]  # top ply of [0/90]s is 0 deg
    outer_90 = sigma_90[-1]  # top ply of [90/0]s is 90 deg
    print(
        f"  hand CLT: [0/90]s outer {outer_0 / 1e6:.4f} MPa, "
        f"[90/0]s outer {outer_90 / 1e6:.4f} MPa, "
        f"mean {mean_0 / 1e6:.4f} / {mean_90 / 1e6:.4f} MPa"
    )
    assert abs(outer_0 - TARGET_OUTER_0_DEG) / TARGET_OUTER_0_DEG < 0.005
    assert abs(outer_90 - TARGET_OUTER_90_DEG) / TARGET_OUTER_90_DEG < 0.005
    assert abs(mean_0 - N_X / H) < 1e-12 * (N_X / H)
    assert abs(mean_90 - N_X / H) < 1e-12 * (N_X / H)


@pytest.mark.parametrize("name", sorted(LAYUPS))
def test_top_is_not_the_thickness_mean(cases: dict, name: str) -> None:
    """TOP must differ from the smeared thickness mean (2.000 MPa would fail both references).

    A thickness-mean value (``N_x / h = 2.000 MPa``) is within neither 5 %
    bound: it is 45.9 % below the 3.696 MPa outer ply and 558 % above the
    0.304 MPa one.  This guard keeps the row from silently regressing to the
    pre-T2 smear.
    """
    case = cases[name]
    top = case["meshes"]["16x4"]["ours_top"]
    mean = case["ref_mean"]
    dev = abs(top - mean) / abs(mean)
    print(
        f"  {name}: TOP {top / 1e6:.4f} MPa vs thickness mean {mean / 1e6:.4f} MPa "
        f"(deviation {dev:.4%})"
    )
    assert dev > 0.10, (
        f"{name}: TOP sigma_xx {top:.6e} equals the thickness mean {mean:.6e} "
        f"(deviation {dev:.4%}) -- the recovery is smeared, not ply-resolved"
    )


@pytest.mark.parametrize("name", sorted(LAYUPS))
def test_free_field_is_mesh_stable_but_the_global_max_is_not(cases: dict, name: str) -> None:
    """E4: the free-field centre converges while the global maximum grows.

    The global maximum is a load-introduction singularity at the loaded free
    edge (x = L; measured on node (1.0, 0.0, -H/2) for ``[0/90]s``), not the
    clamped edge at x = 0, whose value is not singular.  It is therefore
    mesh-dependent and the row must read the centre.  This asserts the
    stability of both codes' free-field values and the growth of the CalculiX
    maximum; the production maximum also grows (3.6968 -> 4.4201 MPa for
    ``[0/90]s``) but is not asserted here.
    """
    case = cases[name]
    coarse = case["meshes"]["8x2"]
    fine = case["meshes"]["16x4"]

    ours_stab = abs(fine["ours_top"] - coarse["ours_top"]) / abs(fine["ours_top"])
    ccx_stab = abs(fine["ccx_free"] - coarse["ccx_free"]) / abs(fine["ccx_free"])
    print(
        f"  {name}: free-field drift ours {ours_stab:.4%}, ccx {ccx_stab:.4%} "
        f"(bound {TOL_FIELD_STABILITY:.2%}); "
        f"ccx max {coarse['ccx_max'] / 1e6:.4f} -> {fine['ccx_max'] / 1e6:.4f} MPa "
        f"vs free field {fine['ccx_free'] / 1e6:.4f} MPa"
    )
    assert ours_stab < TOL_FIELD_STABILITY
    assert ccx_stab < TOL_FIELD_STABILITY
    # The loaded-edge singularity is materially above the free field and grows
    # with refinement: reading the global maximum would validate a mesh effect.
    assert coarse["ccx_max"] > coarse["ccx_free"] * 1.05
    assert fine["ccx_max"] > coarse["ccx_max"]
