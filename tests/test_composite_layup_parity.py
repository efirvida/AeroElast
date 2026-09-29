"""Composite layup parity against CalculiX S8R.

``test_composite_beam_parity.py`` validates one layup (``[0/90/45/-45]s``) under
axial, isotropic-equivalent and bending loads.  This module adds the layups that
had no CalculiX reference, so the composite **material** behaviour is compared
against an independent implementation and not only against the CLT formula the
code itself uses:

- ``[0]`` and ``[90]``: unidirectional, the anisotropy limit;
- ``[0/90]s``: symmetric, B = 0 (the control that coupling does not fire);
- ``[0/90]``: asymmetric, B != 0, the membrane-bending coupling signature;
- ``[0/45/-45/90]s``: quasi-isotropic.

Loading is the **total resultant distributed equally over the free edge**, and
the measured quantity is the **mean over that edge**, so the comparison is a
structural response instead of a single-node local effect.  Measured on this
tree at the 4x10 mesh: the mean axial extension agrees with CCX within 1.83%
(worst case ``[90]``), the mean transverse bending within 0.51%, and the
asymmetric B-coupling within 0.37%.  The tolerances below carry those numbers,
and every test prints the measured value.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from conftest import ccx_bin_or_skip

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler  # noqa: E402

from _ccx_io import fail_ccx, parse_frd_disp, run_ccx  # noqa: E402
from aeroelast.core.laminate import create_laminate_from_angles  # noqa: E402
from aeroelast.core.material import Material  # noqa: E402
from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node, NodeSet  # noqa: E402
from aeroelast.core.mesh.io.writers import write_ccx_mesh  # noqa: E402
from aeroelast.core.mesh.model import MeshModel  # noqa: E402
from aeroelast.core.properties import CompositeShellProperty  # noqa: E402

# ============================================================================
# Geometry and material (same cantilever as test_composite_beam_parity)
# ============================================================================

L = 1.0  # length along Y (cantilever axis)
B = 0.1  # width along X
THICKNESS = 0.005  # total laminate thickness

E1, E2, G12, NU12 = 120e9, 10e9, 5e9, 0.3

#: Layups under test, as fibre angles bottom to top.
LAYUPS: dict[str, list[float]] = {
    "uni_0": [0.0],
    "uni_90": [90.0],
    "sym_0_90s": [0.0, 90.0, 90.0, 0.0],
    "asym_0_90": [0.0, 90.0],
    "quasi_iso": [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0],
}

#: (load DOF, total resultant, measurement DOF).  DOF 0=X, 1=Y, 2=Z.
CASES: dict[str, tuple[int, float, int]] = {
    "axial": (1, 1000.0, 1),
    "bending": (2, 100.0, 2),
    "b_coupling": (1, 1000.0, 2),
}

# Tolerances from the measured margins printed by this module (see the docstring).
AXIAL_TOL = 0.025  # measured max 1.83% (uni_90, 4x10)
BENDING_TOL = 0.01  # measured max 0.51% (uni_0, 4x10)
B_COUPLING_TOL = 0.02  # measured 0.37% (asym_0_90, 4x10)
SYMMETRIC_B_ABS = 1e-11  # symmetric layups must give |w| below this, aero and CCX


# ============================================================================
# Helpers
# ============================================================================


def _material() -> Material:
    return Material(name="cfrp", E=(E1, E2, E2), G=(G12, G12, G12), nu=(NU12, NU12, 0.0), rho=0.0)


def _laminate(angles: list[float]):
    return create_laminate_from_angles(_material(), THICKNESS / len(angles), angles)


def _aero_material_dict(laminate) -> dict:
    """Flat ABD material dict expected by the Rust ``PyMeshAssembler``."""
    abd = laminate.get_ABD_matrix()
    A, B_mat, D = abd[:3, :3], abd[:3, 3:], abd[3:, 3:]
    h = laminate.total_thickness
    return {
        "type": "composite",
        "cm": A.ravel().tolist(),
        "b_coupling": B_mat.ravel().tolist(),
        "cb": D.ravel().tolist(),
        "cs": laminate.Cs.ravel().tolist(),
        # ADR-1 (amended): the element consumes the UNCORRECTED section shear.
        "cs_uncorrected": laminate.shear_stiffness_uncorrected().ravel().tolist(),
        "thickness": h,
        "e_equiv": A[0, 0] / h,
        "mass_per_area": 0.0,
        "rotational_inertia": 0.0,
    }


def _build_mesh(nx: int = 4, ny: int = 10) -> MeshModel:
    """Cantilever shell strip in the XY plane: X width, Y length, Z out of plane."""
    Node._id_counter = 0
    MeshElement._id_counter = 0

    mesh = MeshModel()
    xs = np.linspace(0, B, nx + 1)
    ys = np.linspace(0, L, ny + 1)
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
                    nodes=[grid[(i, j)], grid[(i + 1, j)], grid[(i + 1, j + 1)], grid[(i, j + 1)]],
                    element_type=ElementType.quad,
                )
            )

    clamped = {n for n in mesh.nodes if np.isclose(n.y, 0.0, atol=1e-12)}
    free_face = {n for n in mesh.nodes if np.isclose(n.y, L, atol=1e-12)}
    mesh.add_node_set(NodeSet("clamped", clamped))
    mesh.add_node_set(NodeSet("free_face", free_face))
    mesh.add_element_set(ElementSet("plate", set(mesh.elements)))
    return mesh


def _aero_solution(mesh: MeshModel, laminate, load_dof: int, load: float) -> np.ndarray:
    """Solve the strip with a total resultant distributed over ``free_face``."""
    coords = np.asarray([[n.x, n.y, n.z] for n in mesh.nodes], dtype=float)
    conn = [[mesh.node_id_to_index[nid] for nid in el.node_ids] for el in mesh.elements]
    mats = [_aero_material_dict(laminate)] * len(mesh.elements)

    asm = PyMeshAssembler(
        node_coords=coords, connectivity=conn, elem_types=[4] * len(conn), materials=mats
    )
    rows, cols, vals = asm.assemble_k()
    n = asm.dofs_count
    K = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()

    face = [mesh.node_id_to_index[n_.id] for n_ in mesh.get_node_set("free_face").nodes.values()]
    f = np.zeros(n, dtype=float)
    for nd in face:
        f[6 * nd + load_dof] = load / len(face)

    clamped = {
        mesh.node_id_to_index[n_.id] * 6 + i
        for n_ in mesh.get_node_set("clamped").nodes.values()
        for i in range(6)
    }
    free_mask = np.ones(n, dtype=bool)
    for dof in clamped:
        free_mask[dof] = False
    free = np.where(free_mask)[0]

    u = np.zeros(n, dtype=float)
    u[free] = spsolve(K[np.ix_(free, free)], f[free])
    return u


def _ccx_solution(
    mesh: MeshModel, laminate, load_dof: int, load: float, ccx_bin: str, workdir: Path
):
    """Run CCX S8R with the same resultant on the same node set; return mean per DOF."""
    inp_path = workdir / "layup.inp"
    load_vector = [0.0, 0.0, 0.0]
    load_vector[load_dof] = load

    write_ccx_mesh(
        mesh,
        str(inp_path),
        properties={"plate": CompositeShellProperty(laminate=laminate)},
        boundary_nodeset="clamped",
        load_nodeset="free_face",
        load_vector=load_vector,
        solver_type="LinearStatic",
        quadratic=True,
    )

    node_ids = [
        mesh.node_id_to_index[n_.id] + 1 for n_ in mesh.get_node_set("free_face").nodes.values()
    ]
    result = run_ccx(inp_path, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp_path)

    frd_path = inp_path.with_suffix(".frd")
    if not frd_path.exists():
        pytest.fail(f"No FRD output from CCX for {inp_path.name}")
    disp = parse_frd_disp(frd_path, node_ids)
    if not disp:
        pytest.fail(f"Could not parse CCX displacements from {frd_path.name}")
    return disp, node_ids


@pytest.fixture(scope="module")
def parity(tmp_path_factory: pytest.TempPathFactory) -> dict:
    """Run every (layup, case) once: {key: (mean_aero, mean_ccx)}."""
    ccx_bin = ccx_bin_or_skip()
    workdir = tmp_path_factory.mktemp("composite_layup_parity")
    mesh = _build_mesh()
    face = [mesh.node_id_to_index[n_.id] for n_ in mesh.get_node_set("free_face").nodes.values()]

    results: dict[tuple[str, str], tuple[float, float]] = {}
    for name, angles in LAYUPS.items():
        laminate = _laminate(angles)
        for case, (load_dof, load, measure_dof) in CASES.items():
            u = _aero_solution(mesh, laminate, load_dof, load)
            aero_mean = float(np.mean([u[6 * nd + measure_dof] for nd in face]))

            with tempfile.TemporaryDirectory(dir=workdir) as td:
                disp, node_ids = _ccx_solution(mesh, laminate, load_dof, load, ccx_bin, Path(td))
                ccx_mean = float(np.mean([disp[nid][measure_dof] for nid in node_ids]))

            results[(name, case)] = (aero_mean, ccx_mean)
            rel = abs(aero_mean - ccx_mean) / max(abs(ccx_mean), 1e-30)
            print(
                f"{name:10s} {case:11s} aero={aero_mean: .6e} ccx={ccx_mean: .6e} "
                f"rel_err={rel * 100:5.2f}%"
            )
    return results


# ============================================================================
# Tests
# ============================================================================


@pytest.mark.parametrize("layup", list(LAYUPS))
def test_axial_extension_matches_ccx(parity: dict, layup: str) -> None:
    """Mean axial extension under an edge resultant matches CCX S8R."""
    aero, ccx = parity[(layup, "axial")]
    rel = abs(aero - ccx) / max(abs(ccx), 1e-30)
    assert rel < AXIAL_TOL, (
        f"{layup}: axial extension aero={aero:.6e} ccx={ccx:.6e} rel_err={rel * 100:.2f}% "
        f"(tol {AXIAL_TOL * 100:.1f}%)"
    )


@pytest.mark.parametrize("layup", list(LAYUPS))
def test_transverse_bending_matches_ccx(parity: dict, layup: str) -> None:
    """Mean out-of-plane bending deflection matches CCX S8R."""
    aero, ccx = parity[(layup, "bending")]
    rel = abs(aero - ccx) / max(abs(ccx), 1e-30)
    assert rel < BENDING_TOL, (
        f"{layup}: bending aero={aero:.6e} ccx={ccx:.6e} rel_err={rel * 100:.2f}% "
        f"(tol {BENDING_TOL * 100:.1f}%)"
    )


def test_asymmetric_b_coupling_matches_ccx(parity: dict) -> None:
    """The asymmetric [0/90] strip bends out of plane under axial load, as in CCX.

    This is the membrane-bending coupling signature. The test is non-vacuous:
    the CCX reference itself has to be macroscopic, not round-off.
    """
    aero, ccx = parity[("asym_0_90", "b_coupling")]
    assert abs(ccx) > 1e-3, f"CCX reference is not macroscopic: {ccx:.3e}"
    rel = abs(aero - ccx) / abs(ccx)
    assert rel < B_COUPLING_TOL, (
        f"asym_0_90: B-coupling aero={aero:.6e} ccx={ccx:.6e} rel_err={rel * 100:.2f}% "
        f"(tol {B_COUPLING_TOL * 100:.1f}%)"
    )


@pytest.mark.parametrize("layup", ["sym_0_90s", "quasi_iso"])
def test_symmetric_laminates_have_no_b_coupling(parity: dict, layup: str) -> None:
    """Symmetric laminates stay flat under axial load, in both codes (B = 0)."""
    aero, ccx = parity[(layup, "b_coupling")]
    assert abs(aero) < SYMMETRIC_B_ABS, f"{layup}: aero w={aero:.3e} (expected ~0)"
    assert abs(ccx) < SYMMETRIC_B_ABS, f"{layup}: ccx w={ccx:.3e} (expected ~0)"
