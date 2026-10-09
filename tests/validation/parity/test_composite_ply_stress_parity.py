"""Composite outer-fibre ply stress: AeroElast against a CLT closed form.

Why this test exists
--------------------
For a composite shell, the stress at a through-thickness location is the
**outer ply's stress**: within each ply the membrane stress is constant, but it
differs between plies whose fibre angles differ.  The recovery currently in the
tree does not deliver that: it returns the homogenised, thickness-averaged
equivalent stress of the section (``sigma_m = (C/h) . B_m u`` with
``C = Cm()/h`` and ``Cm()`` the laminate's integrated ``A`` — the documented
``h`` trap), i.e. the smeared value ``N_x / h`` at every location.

The discriminating case is a symmetric cross-ply strip in membrane tension,
where CLT gives distinct ply stresses:

===============  ===========  ===========  ===============
layup            outer ply    second ply   thickness mean
===============  ===========  ===========  ===============
``[0/90]s``      3.696 MPa    0.304 MPa    2.000 MPa
``[90/0]s``      0.304 MPa    3.696 MPa    2.000 MPa
===============  ===========  ===========  ===============

Measured on this tree (8x2 and 16x4, stable): ``StressRecovery`` reports
2.0000 MPa at TOP / MIDDLE / BOTTOM for **both** stacks — the smeared mean.
The error against the true outer fibre is **-45.9 %** for ``[0/90]s`` and
**+558 %** for ``[90/0]s``, far outside the 5 % bound below.  This test is the
RED half of the ply-stress capability (``odd/tasks/composite-ply-stress-recovery.md``,
T1); the fix (T2) resolves the ply containing the requested ``z`` and evaluates
``sigma = Qbar(theta_ply) . eps(z)``.

Reference independence
----------------------
The reference is re-implemented in this file from first principles — own ``Q``,
own ``Qbar(theta)``, own ``A`` (Jones 1999 / Reddy 2004 conventions), following
the ``_hand_clt_abd`` pattern of ``test_composite_beam_parity.py``.  Nothing is
imported from ``aeroelast.core.laminate`` for the reference value; the library
is only used to build the model under test.

Restraint singularity
---------------------
The clamped edge is a restraint singularity whose interior maximum grows under
refinement, so the comparison reads the **centre element** (nearest centroid to
``(L/2, B/2)``), never the global maximum.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests import _production_arbiter as arbiter  # noqa: E402
from tests.support.assertions import assert_relative_error  # noqa: E402
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

NX, NY = 8, 2  # 8 along x, 2 across y (the E1 probe mesh; stable vs 16x4)
FORCE = 1000.0  # total axial resultant at the free edge, +x [N]
N_X = FORCE / B  # membrane resultant per unit width [N/m] = 1e4

E1, E2, G12, G23, NU12, RHO = 120e9, 10e9, 5e9, 3e9, 0.3, 1000.0

#: Layups under test, ply angles bottom to top.  Symmetric, so B = 0 and the
#: closed form is pure membrane: sigma(z) is constant within each ply.
LAYUPS: dict[str, list[float]] = {
    "sym_0_90s": [0.0, 90.0, 90.0, 0.0],
    "sym_90_0s": [90.0, 0.0, 0.0, 90.0],
}

#: Bound in the style of the isotropic sibling (``TOL_CCX = 0.05`` in
#: ``test_shell_stress_ccx_parity.py``).  The pre-fix smear misses it by
#: -45.9 % ([0/90]s) and +558 % ([90/0]s) — measured on this tree.
TOL_CLT = 0.05

#: Published CLT targets for the reference self-check (contract, 2026-10-09).
TARGET_OUTER_0_DEG = 3.696e6  # sigma_xx of a 0 deg ply in [0/90]s [Pa]
TARGET_OUTER_90_DEG = 0.304e6  # sigma_xx of a 90 deg ply in [0/90]s [Pa]


# ============================================================================
# Independent CLT reference (Jones 1999 / Reddy 2004), re-implemented
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


def _hand_clt_ply_stresses(
    E1: float, E2: float, G12: float, nu12: float, angles: list[float], n_x: float
) -> tuple[list[float], float]:
    """Per-ply ``sigma_xx`` (global x axes, bottom to top) and the thickness mean.

    Classical lamination theory for a symmetric laminate under a pure membrane
    resultant ``N = [n_x, 0, 0]``: ``eps0 = A^-1 N`` and ``sigma_k =
    Qbar(theta_k) . eps0`` (constant within each ply).  ``A = sum Qbar_k
    (z_top - z_bottom)`` is assembled from the same first-principles ``Qbar``
    — nothing is imported from ``aeroelast.core.laminate``.
    """
    q = _hand_q(E1, E2, G12, nu12)
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


def test_reference_reproduces_published_clt_targets() -> None:
    """Guard: the in-file closed form reproduces the published 3.696 / 0.304 MPa.

    Machinery, not a row: it pins the hand reference to the contract's target
    values so a transcription error in the reference cannot masquerade as a
    solver result.
    """
    sigma_xx, mean = _hand_clt_ply_stresses(E1, E2, G12, NU12, LAYUPS["sym_0_90s"], N_X)
    # [0/90]s bottom-to-top: 0 deg ply is both the bottom and the top (outer) ply.
    outer_0 = sigma_xx[-1]
    outer_90 = sigma_xx[1]  # second ply from the bottom is a 90 deg ply
    rel_0 = abs(outer_0 - TARGET_OUTER_0_DEG) / TARGET_OUTER_0_DEG
    rel_90 = abs(outer_90 - TARGET_OUTER_90_DEG) / TARGET_OUTER_90_DEG
    print(
        f"  hand CLT: 0 deg ply {outer_0 / 1e6:.4f} MPa, 90 deg ply {outer_90 / 1e6:.4f} MPa, "
        f"thickness mean {mean / 1e6:.4f} MPa (targets 3.696 / 0.304 / 2.000 MPa)"
    )
    assert rel_0 < 0.005, f"hand CLT 0 deg ply {outer_0:.6e} vs target 3.696e6 ({rel_0:.3%})"
    assert rel_90 < 0.005, f"hand CLT 90 deg ply {outer_90:.6e} vs target 0.304e6 ({rel_90:.3%})"
    assert abs(mean - N_X / H) < 1e-12 * (N_X / H), (
        f"thickness mean {mean:.6e} must equal N_x/h = {N_X / H:.6e} by equilibrium"
    )


# ============================================================================
# Production-path model
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
    # shear_correction_factor=5/6 matches the value ``build_rust_properties``
    # hardcodes, so the CompositeShellProperty form and the production Rust
    # ``Laminate`` form define the *identical* section (same A, B, D and Cs) —
    # required by the two-forms-agree cross-check below.
    return create_laminate_from_angles(
        _material(), PLY_THICKNESS, angles, shear_correction_factor=5.0 / 6.0
    )


def _model_cfg(laminate) -> dict:
    return _model_cfg_from_property(CompositeShellProperty(laminate=laminate))


def _model_cfg_from_property(prop) -> dict:
    """Model config carrying an arbitrary section-property form under ``'plate'``.

    Used with ``CompositeShellProperty`` (test form), the production Rust
    ``_aeroelast.Laminate`` from ``build_rust_properties`` and the legacy raw
    ABD dict — the assembler accepts all three in ``properties``.
    """
    return {
        "elements": {
            "element_family": ElementFamily.SHELL,
            "properties": {"plate": prop},
        },
        "solver": {"time_step": 0.01, "total_time": 1.0, "beta": 0.25, "gamma": 0.5},
    }


def _build_mesh(nx: int = NX, ny: int = NY) -> MeshModel:
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
    """Production static solve (Rust assembly + PETSc KSP), driven as the arbiter drives it.

    Mirrors ``test_shell_stress_ccx_parity.py::test_production_solver_matches_the_scipy_replica``:
    ``arbiter.solve_static`` wraps the production ``StaticLinearSolver`` on the
    same mesh, BCs and nodal loads — no test-only scipy solve anywhere.
    """
    edge = arbiter.nodes_at_coordinate(mesh, L, axis=0)
    return arbiter.solve_static(
        mesh,
        model_cfg,
        fixed_dofs=arbiter.fixed_dofs_from_node_set(mesh, "clamped"),
        nodal_loads_list=[arbiter.nodal_loads(edge, FORCE, dof=0)],
    )


def _centre_element_index(domain: MeshAssembler) -> int:
    """Index (into ``domain.elements``) of the element whose centroid is nearest (L/2, B/2).

    The stress rows returned by ``compute_element_stresses`` follow the Rust
    assembler's element order, which is the same order ``domain.elements``
    iterates — the assumption ``compute_nodal_stresses`` already relies on
    (``stress_recovery.py``: it zips ``enumerate(self.domain.elements)`` with
    the element-major Rust array).
    """
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


@pytest.fixture(scope="module")
def cases() -> dict:
    """Solve each layup once through the production path; recover stresses at the centre element."""
    out: dict[str, dict] = {}
    for name, angles in LAYUPS.items():
        mesh = _build_mesh()
        model_cfg = _model_cfg(_laminate(angles))

        u = _solve_production(mesh, model_cfg)
        # The production assembler (same Rust ``from_model`` path the solver
        # builds internally) provides the element map for StressRecovery.
        domain = MeshAssembler(mesh=mesh, model=model_cfg)
        recovery = StressRecovery(domain, u)

        idx = _centre_element_index(domain)
        sxx = {}
        sxx_arr = {}
        for key, loc in (
            ("top", StressLocation.TOP),
            ("middle", StressLocation.MIDDLE),
            ("bottom", StressLocation.BOTTOM),
        ):
            res = recovery.compute_element_stresses(location=loc, stress_type=StressType.TOTAL)
            assert res.sigma_xx.shape[0] == len(list(domain.elements)), (
                "compute_element_stresses must return one row per assembled element"
            )
            sxx[key] = float(res.sigma_xx[idx])
            sxx_arr[key] = res.sigma_xx.copy()

        edge = arbiter.nodes_at_coordinate(mesh, L, axis=0)
        tip_ux = float(np.mean([u[6 * i + 0] for i in edge]))

        sigma_xx_ref, mean_ref = _hand_clt_ply_stresses(E1, E2, G12, NU12, angles, N_X)
        out[name] = {
            "angles": angles,
            "tip_ux": tip_ux,
            "sxx": sxx,
            "sxx_arr": sxx_arr,
            "sigma_xx_ref": sigma_xx_ref,  # per ply, bottom to top
            "mean_ref": mean_ref,
            "centre_index": idx,
        }
        print(
            f"{name}: tip u_x = {tip_ux:.6e} m; centre element #{idx} sigma_xx "
            f"TOP/MID/BOT = {sxx['top'] / 1e6:.4f} / {sxx['middle'] / 1e6:.4f} / "
            f"{sxx['bottom'] / 1e6:.4f} MPa; hand CLT per ply (bottom→top) = "
            f"{[f'{v / 1e6:.4f}' for v in sigma_xx_ref]} MPa, mean {mean_ref / 1e6:.4f} MPa"
        )
    return out


# ============================================================================
# The capability comparisons
# ============================================================================


@pytest.mark.parametrize("name", sorted(LAYUPS))
def test_top_fibre_ply_stress_matches_clt(cases: dict, name: str) -> None:
    """TOP stress of a composite shell is the topmost ply's stress, not the smeared mean."""
    case = cases[name]
    # TOP is z = +h/2: the topmost ply of the stack carries it.
    ref_outer = case["sigma_xx_ref"][-1]
    assert_relative_error(
        case["sxx"]["top"],
        ref_outer,
        tol=TOL_CLT,
        kind="analytical",
        reference_name=(
            "CLT closed form re-implemented in this test (Jones 1999 / Reddy 2004), "
            "outer-fibre ply sigma_xx"
        ),
        what=f"{name} TOP outer-fibre ply sigma_xx",
    )


@pytest.mark.parametrize("name", sorted(LAYUPS))
def test_bottom_fibre_ply_stress_matches_clt(cases: dict, name: str) -> None:
    """BOTTOM stress is the bottommost ply's stress — the symmetric-twin claim of the TOP one."""
    case = cases[name]
    # BOTTOM is z = -h/2: the bottom ply of the stack carries it.
    ref_bottom = case["sigma_xx_ref"][0]
    assert_relative_error(
        case["sxx"]["bottom"],
        ref_bottom,
        tol=TOL_CLT,
        kind="analytical",
        reference_name=(
            "CLT closed form re-implemented in this test (Jones 1999 / Reddy 2004), "
            "outer-fibre ply sigma_xx"
        ),
        what=f"{name} BOTTOM outer-fibre ply sigma_xx",
    )


@pytest.mark.parametrize("name", sorted(LAYUPS))
def test_top_is_not_the_thickness_mean(cases: dict, name: str) -> None:
    """Through-thickness shape the closed form requires: TOP differs from the smeared mean.

    Under membrane load the stress is constant within each ply but differs
    between the 0 and 90 deg plies, so TOP must NOT equal the thickness mean
    ``N_x / h``.  In the closed form the deviation is 84.8 % of the mean; the
    pre-fix smear delivers exactly the mean (deviation 0 %).
    """
    case = cases[name]
    dev = abs(case["sxx"]["top"] - case["mean_ref"]) / abs(case["mean_ref"])
    print(
        f"{name}: TOP {case['sxx']['top']:.6e} vs thickness mean {case['mean_ref']:.6e} "
        f"(deviation {dev:.4%}; closed form requires {abs(case['sigma_xx_ref'][-1] - case['mean_ref']) / case['mean_ref']:.4%})"
    )
    assert dev > 0.10, (
        f"{name}: TOP sigma_xx {case['sxx']['top']:.6e} equals the thickness-mean equivalent "
        f"{case['mean_ref']:.6e} (deviation {dev:.4%}) — the recovery is smeared, not ply-resolved"
    )


# ============================================================================
# T2b — the PRODUCTION section form (Rust Laminate from build_rust_properties)
# ============================================================================


def _numad_data(angles: list[float]) -> dict:
    """Minimal NuMAD payload built the way production builds sections.

    The layup ``[[mat_name, thickness, angle], ...]`` reproduces the same
    ``[0/90]s`` / ``[90/0]s`` stacks (bottom to top) as the ``CompositeShellProperty``
    cases, with the same material constants — so the two section forms must give
    the same numbers everywhere.
    """
    return {
        "materials": [
            {
                "name": "cfrp",
                "density": RHO,
                "elastic": {
                    "E": [E1, E2, E2],
                    "G": [G12, G23, G23],
                    "nu": [NU12, NU12, 0.0],
                },
            }
        ],
        "sections": [
            {
                "elementSet": "plate",
                "layup": [["cfrp", PLY_THICKNESS, a] for a in angles],
            }
        ],
    }


@pytest.fixture(scope="module")
def rust_cases(cases: dict) -> dict:
    """Solve each layup through the PRODUCTION section form (Rust ``Laminate``).

    ``build_rust_properties`` is the map builder ``runner._extract_blade_properties``
    delegates to — its values are ``_aeroelast.Laminate`` (or isotropic dicts), never
    a Python ``Laminate``/``CompositeShellProperty``.  This fixture is the guard that
    the ply resolution reaches that production form.
    """
    import _aeroelast

    from aeroelast.models.blade.model import build_rust_properties

    out: dict[str, dict] = {}
    for name, angles in LAYUPS.items():
        props = build_rust_properties(_numad_data(angles))
        prop = props["plate"]
        # Proof the production-shaped map really contains the Rust type:
        assert type(prop) is _aeroelast.Laminate, (
            f"build_rust_properties returned {type(prop)!r}, expected _aeroelast.Laminate"
        )
        print(f"{name}: properties['plate'] type = {type(prop)}")

        mesh = _build_mesh()
        model_cfg = _model_cfg_from_property(prop)
        u = _solve_production(mesh, model_cfg)
        domain = MeshAssembler(mesh=mesh, model=model_cfg)
        recovery = StressRecovery(domain, u)

        idx = _centre_element_index(domain)
        sxx = {}
        for key, loc in (
            ("top", StressLocation.TOP),
            ("middle", StressLocation.MIDDLE),
            ("bottom", StressLocation.BOTTOM),
        ):
            res = recovery.compute_element_stresses(location=loc, stress_type=StressType.TOTAL)
            sxx[key] = res.sigma_xx.copy()

        sigma_xx_ref, mean_ref = _hand_clt_ply_stresses(E1, E2, G12, NU12, angles, N_X)
        out[name] = {
            "sxx": sxx,
            "sigma_xx_ref": sigma_xx_ref,
            "mean_ref": mean_ref,
            "centre_index": idx,
            "tip_ux": float(
                np.mean([u[6 * i + 0] for i in arbiter.nodes_at_coordinate(mesh, L, axis=0)])
            ),
            "prop": prop,
        }
        print(
            f"{name} [Rust Laminate]: tip u_x = {out[name]['tip_ux']:.6e} m; centre element "
            f"#{idx} sigma_xx TOP/MID/BOT = {sxx['top'][idx] / 1e6:.4f} / "
            f"{sxx['middle'][idx] / 1e6:.4f} / {sxx['bottom'][idx] / 1e6:.4f} MPa"
        )
    return out


@pytest.mark.parametrize("name", sorted(LAYUPS))
def test_production_rust_laminate_top_fibre_matches_clt(rust_cases: dict, name: str) -> None:
    """The production section form (Rust ``Laminate``) resolves the same outer ply.

    Same closed-form targets as the ``CompositeShellProperty`` cases: the ply
    resolution must reach the production path, not only the test construction.
    """
    case = rust_cases[name]
    assert_relative_error(
        float(case["sxx"]["top"][case["centre_index"]]),
        case["sigma_xx_ref"][-1],
        tol=TOL_CLT,
        kind="analytical",
        reference_name=(
            "CLT closed form re-implemented in this test (Jones 1999 / Reddy 2004), "
            "outer-fibre ply sigma_xx, production Rust-Laminate section"
        ),
        what=f"{name} TOP outer-fibre ply sigma_xx (production Rust Laminate)",
    )


@pytest.mark.parametrize("name", sorted(LAYUPS))
def test_production_rust_laminate_bottom_fibre_matches_clt(rust_cases: dict, name: str) -> None:
    """BOTTOM of the production form is the bottommost ply's stress (same claim, other edge)."""
    case = rust_cases[name]
    assert_relative_error(
        float(case["sxx"]["bottom"][case["centre_index"]]),
        case["sigma_xx_ref"][0],
        tol=TOL_CLT,
        kind="analytical",
        reference_name=(
            "CLT closed form re-implemented in this test (Jones 1999 / Reddy 2004), "
            "outer-fibre ply sigma_xx, production Rust-Laminate section"
        ),
        what=f"{name} BOTTOM outer-fibre ply sigma_xx (production Rust Laminate)",
    )


@pytest.mark.parametrize("name", sorted(LAYUPS))
def test_production_and_test_section_forms_agree(cases: dict, rust_cases: dict, name: str) -> None:
    """Cross-check: ``CompositeShellProperty`` and Rust ``Laminate`` give equal numbers.

    The two forms are normalised into ONE ply-resolution path; this is the guard
    that the production route is not a second implementation.  Both solves run the
    identical system (same material constants, same stack, same mesh), so the full
    per-element ``sigma_xx`` vectors must agree to solver round-off.
    """
    case, rcase = cases[name], rust_cases[name]
    assert rcase["tip_ux"] == pytest.approx(case["tip_ux"], rel=1e-12), (
        "the two section forms must assemble the identical stiffness (same tip displacement)"
    )
    for key in ("top", "middle", "bottom"):
        np.testing.assert_allclose(
            rcase["sxx"][key],
            case["sxx_arr"][key],
            rtol=1e-10,
            atol=0.0,
            err_msg=(
                f"{name} {key}: production Rust-Laminate form diverges from the "
                "CompositeShellProperty form — the ply resolution split into a second path"
            ),
        )
        print(
            f"{name} {key}: centre element {case['sxx_arr'][key][case['centre_index']]:.6e} Pa "
            f"(both forms, max |diff| over elements "
            f"{np.max(np.abs(rcase['sxx'][key] - case['sxx_arr'][key])):.3e} Pa)"
        )


# ============================================================================
# T2b — raw composite ABD dicts: pinned smeared limitation + one-time warning
# ============================================================================


def _raw_abd_property(angles: list[float]) -> dict:
    """Legacy raw composite ABD dict form (no ply stack), built from the same laminate.

    Mirrors the dict shape ``assembler._orthotropic_shell_to_composite_dict`` emits
    and the Rust ``parse_material`` composite branch consumes.
    """
    from aeroelast.core.laminate import create_laminate_from_angles as _clf

    lam = _clf(_material(), PLY_THICKNESS, angles)
    t = float(lam.total_thickness)
    return {
        "type": "composite",
        "cm": lam.A.ravel().tolist(),
        "b_coupling": lam.B.ravel().tolist(),
        "cb": lam.D.ravel().tolist(),
        "cs": lam.Cs.ravel().tolist(),
        "thickness": t,
        "e_equiv": float(np.trace(lam.A) / (3.0 * t)),
        "mass_per_area": float(sum(p.material.rho * p.thickness for p in lam.plies)),
        "rotational_inertia": float(
            sum(p.material.rho * (p.z_top**3 - p.z_bottom**3) / 3.0 for p in lam.plies)
        ),
    }


def test_raw_abd_section_is_smeared_and_warns(caplog) -> None:
    """Raw ABD dicts carry no ply stack: the smeared thickness mean is pinned, not silent.

    Characterization of the documented limitation: for a raw composite ABD dict the
    recovery cannot resolve plies, keeps the homogenised value (``N_x / h`` at every
    through-thickness location) and emits a one-time warning naming the element set —
    the caller is never left to mistake the thickness mean for a ply stress.
    """
    mesh = _build_mesh()
    model_cfg = _model_cfg_from_property(_raw_abd_property(LAYUPS["sym_0_90s"]))
    u = _solve_production(mesh, model_cfg)
    domain = MeshAssembler(mesh=mesh, model=model_cfg)
    recovery = StressRecovery(domain, u)

    with caplog.at_level(logging.WARNING, logger="aeroelast.postprocess.stress_recovery"):
        res = recovery.compute_element_stresses(
            location=StressLocation.TOP, stress_type=StressType.TOTAL
        )

    idx = _centre_element_index(domain)
    top = float(res.sigma_xx[idx])
    mean = N_X / H
    print(f"raw ABD [0/90]s: TOP sigma_xx = {top:.6e} Pa, thickness mean = {mean:.6e} Pa")
    # The smeared value IS the thickness mean (the pre-T2 behaviour for every form).
    assert_relative_error(
        top,
        mean,
        tol=0.01,
        kind="analytical",
        reference_name="N_x / h thickness-mean equivalent (raw ABD dict, no ply stack)",
        what="[0/90]s TOP sigma_xx (raw composite ABD dict)",
    )
    # And the warning names the element set and the limitation.
    messages = [r.getMessage() for r in caplog.records]
    assert any("plate" in m and "thickness-mean" in m and "no ply stack" in m for m in messages), (
        f"expected a one-time raw-ABD warning naming 'plate'; got {messages}"
    )
