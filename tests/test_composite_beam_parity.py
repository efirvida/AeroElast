"""Composite Shell Cantilever Parity: AeroElast vs CalculiX (CCX).

Geometry:
- Cantilever plate/shell clamped at Z=0, free at Z=L
- Uses composite shell elements with multiple layers
- Different fiber orientations (typical layup: [0/90/45/-45])

Load Cases:
1. Axial tension: +Fz (extension along Z)
2. Bending: Fx applied at free end (bending about Y-axis)

This test validates composite shell formulation against CalculiX.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from conftest import ccx_bin_or_skip

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler

from aeroelast.core.laminate import create_laminate_from_angles
from aeroelast.core.material import Material
from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node, NodeSet
from aeroelast.core.mesh.io.writers import write_ccx_mesh
from aeroelast.core.mesh.model import MeshModel
from aeroelast.core.properties import CompositeShellProperty

pytestmark = [pytest.mark.slow]


# ============================================================================
# Geometry and Material Constants
# ============================================================================

L = 1.0  # Beam/plate length along Z (m)
B = 0.1  # Width in X (m)
thickness = 0.005  # Total laminate thickness (5mm)

# Composite layup: typical [0/90/45/-45]s configuration
E1 = 120e9  # Longitudinal modulus (Pa)
E2 = 10e9  # Transverse modulus
G12 = 5e9  # Shear modulus
nu12 = 0.3  # Poisson ratio

# ============================================================================
# CCX Helpers
# ============================================================================


def _parse_frd_disp(frd_file: Path, node_ids: list[int]) -> dict[int, np.ndarray]:
    """Parse the last FRD displacement block, returning only requested nodes.

    CalculiX writes nodal records in fixed-width columns
    (``(1X,'-1',I10,3E12.5)``): a negative component abuts the previous one
    with no separator, so a greedy character-class regex silently drops every
    node that has a negative component.  The node number is read from columns
    4-13 and each component from the following 12-column fields.
    """
    wanted = set(node_ids)
    disps: dict[int, np.ndarray] = {}

    with open(frd_file, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    last_disp_start = -1
    for i, line in enumerate(lines):
        if "-4" in line and "DISP" in line.upper():
            last_disp_start = i

    if last_disp_start == -1:
        return disps

    for line in lines[last_disp_start + 1 :]:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("-3") or stripped.startswith("*") or "STEP" in stripped.upper():
            break
        if not stripped.startswith("-1"):
            continue
        try:
            node_id = int(line[3:13])
            if node_id not in wanted:
                continue
            disps[node_id] = np.array(
                [float(line[13:25]), float(line[25:37]), float(line[37:49])],
                dtype=float,
            )
        except ValueError:
            continue

    return disps


def _run_ccx(inp_path: Path, ccx_bin: str) -> subprocess.CompletedProcess:
    """Run CalculiX on input file, inside its own directory.

    ``cwd`` is passed to the subprocess instead of mutating the process-wide
    working directory with ``os.chdir`` (which leaked into later tests).
    """
    return subprocess.run(
        [ccx_bin, inp_path.stem],
        cwd=inp_path.parent,
        capture_output=True,
        text=True,
    )


def _fail_ccx(result: subprocess.CompletedProcess, inp_path: Path) -> None:
    """Hard-fail on a non-zero CCX exit; ccx writes its diagnostics to stdout."""
    pytest.fail(
        f"CCX failed for {inp_path.name} (rc={result.returncode}):\n"
        f"STDOUT tail:\n{result.stdout[-2000:]}\n"
        f"STDERR tail:\n{result.stderr[-2000:]}"
    )


# ============================================================================
# Build Mesh - Cantilever Plate with Composite Shells
# ============================================================================


def _build_composite_plate_mesh(nx: int = 4, ny: int = 10) -> MeshModel:
    """Build cantilever plate mesh for composite shell test.

    ``nx`` and ``ny`` are the element counts along X (width) and Y (length);
    the defaults reproduce the original 4x10 mesh.  Other densities are used
    by the laminate mesh-dependence study.

    The mesh is a proper 2D shell cantilever in the XY plane:
    - X: width direction (-B/2 to +B/2)
    - Y: length direction (0 to L, cantilever axis)
    - Z = 0: mid-surface of the shell
    - Elements span X and Y, Z is fixed

    This matches the CCX *SHELL SECTION setup where the shell lies in
    the XY plane and deforms out-of-plane (Z displacement).
    """
    # Reset global counters so node ids are 0-based and contiguous.  The CCX
    # quadratic upgrade indexes coordinates by node id, so a stale counter
    # (from a previous mesh built in the same process) breaks it.
    Node._id_counter = 0
    MeshElement._id_counter = 0

    mesh = MeshModel()

    # Create node grid in XY plane (Z=0 everywhere)
    xs = np.linspace(0, B, nx + 1)
    ys = np.linspace(0, L, ny + 1)
    # Z is fixed at 0 (mid-surface)
    # We use a 2D grid (X, Y) — one node per (i, j) pair
    grid: dict[tuple[int, int], Node] = {}
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            node = Node([float(x), float(y), 0.0], geometric_node=False)
            mesh.add_node(node)
            grid[(i, j)] = node

    # Create quadrilateral shell elements in the XY plane.
    # Nodes are numbered consecutively: node_id = j * (nx + 1) + i
    for j in range(ny):
        for i in range(nx):
            n00 = grid[(i, j)]
            n10 = grid[(i + 1, j)]
            n11 = grid[(i + 1, j + 1)]
            n01 = grid[(i, j + 1)]
            mesh.add_element(
                MeshElement(
                    nodes=[n00, n10, n11, n01],
                    element_type=ElementType.quad,
                )
            )

    # Node sets
    clamped = {n for n in mesh.nodes if np.isclose(n.y, 0.0, atol=1e-12)}
    free_face = {n for n in mesh.nodes if np.isclose(n.y, L, atol=1e-12)}
    free_center = min(free_face, key=lambda n: abs(float(n.x) - B / 2))

    mesh.add_node_set(NodeSet("clamped", clamped))
    mesh.add_node_set(NodeSet("free_face", free_face))
    mesh.add_node_set(NodeSet("free_center", {free_center}))
    mesh.add_element_set(ElementSet("plate", set(mesh.elements)))

    return mesh


# ============================================================================
# CCX Input Writer - Composite Shells
# ============================================================================


def _write_composite_ccx_inp(
    inp_path: Path,
    mesh: MeshModel,
    load_vector: tuple[float, float, float],
    laminate=None,
) -> None:
    """Write CCX input for a laminate (S8R + COMPOSITE section).

    ``laminate`` defaults to the real 8-ply symmetric [0/90/45/-45]s laminate,
    preserving the original behaviour.  Pass an explicit laminate to override
    it -- e.g. a single-ply isotropic laminate used to isolate the
    element-formulation difference from the material difference.
    """

    if laminate is None:
        # Build the real 8-ply symmetric [0/90/45/-45]s laminate.
        from aeroelast.core.laminate import create_laminate_from_angles
        from aeroelast.core.material import Material

        mat = Material(
            name="comp",
            E=(E1, E2, E2),
            G=(G12, G12, G12),
            nu=(nu12, nu12, 0.0),
            rho=0.0,
        )
        ply_t = thickness / 8
        angles = [0.0, 90.0, 45.0, -45.0, -45.0, 45.0, 90.0, 0.0]
        laminate = create_laminate_from_angles(mat, ply_t, angles)

    # CCX needs S8R + *SHELL SECTION, COMPOSITE for a real laminate; linear S4
    # does not support composite sections, so quadratic=True is required.
    props = {"plate": CompositeShellProperty(laminate=laminate)}

    write_ccx_mesh(
        mesh,
        str(inp_path),
        properties=props,
        load_nodeset="free_center",
        load_vector=list(load_vector),
        solver_type="LinearStatic",
        quadratic=True,
    )


# ============================================================================
# Laminate Material Helpers
# ============================================================================


def _make_laminate_mat(E1: float, E2: float, G12: float, nu12: float, thickness: float) -> dict:
    """Compute ABD matrices and return a material dict for PyMeshAssembler.

    Uses the Python CLT to build the 8-ply symmetric [0/90/45/-45]s layup,
    then exposes the result in the flat format expected by the Rust binding.
    """
    # Engineering constants -> orthotropic material
    mat = Material(
        name="comp",
        E=(E1, E2, E2),
        G=(G12, G12, G12),
        nu=(nu12, nu12, 0.0),
        rho=0.0,
    )

    # 8-ply symmetric: [0/90/45/-45/-45/45/90/0]
    ply_t = thickness / 8
    angles = [0.0, 90.0, 45.0, -45.0, -45.0, 45.0, 90.0, 0.0]
    laminate = create_laminate_from_angles(mat, ply_t, angles)

    ABD = laminate.get_ABD_matrix()  # 6x6 [[A,B],[B,D]]
    A = ABD[:3, :3]
    B = ABD[:3, 3:]
    D = ABD[3:, 3:]
    Cs = laminate.Cs

    # Flatten for Rust binding: row-major
    cm = A.ravel()  # 9 elements
    cb_coupling = B.ravel()  # 9 elements
    cb = D.ravel()  # 9 elements
    cs = Cs.ravel()  # 4 elements

    h = laminate.total_thickness
    rho_eq = 0.0  # no density needed for stiffness test
    mass_per_area = rho_eq * h
    rotational_inertia = rho_eq * h**3 / 12.0

    return {
        "type": "composite",
        "cm": cm.tolist(),
        "b_coupling": cb_coupling.tolist(),
        "cb": cb.tolist(),
        "cs": cs.tolist(),
        "thickness": h,
        "e_equiv": A[0, 0] / h,
        "mass_per_area": mass_per_area,
        "rotational_inertia": rotational_inertia,
    }


# ============================================================================
# Tests
# ============================================================================


def _ccx_bin_or_skip() -> str:
    """Find CCX binary or skip (shared resolver in ``conftest``)."""
    return ccx_bin_or_skip()


# =============================================================================
# Unit tests (no CCX needed)
# =============================================================================


class TestCompositeMaterial:
    """Test composite material and ABD matrix computation."""

    def test_laminate_abd_matrices(self):
        """Verify ABD matrices are computed correctly."""
        E1, E2, G12, nu12, t = 120e9, 10e9, 5e9, 0.3, 0.005
        mat_dict = _make_laminate_mat(E1, E2, G12, nu12, t)

        # Check keys exist
        assert "cm" in mat_dict
        assert "b_coupling" in mat_dict
        assert "cb" in mat_dict
        assert "cs" in mat_dict
        assert mat_dict["thickness"] == t

        # Check shapes (3x3 matrices flattened to 9 elements)
        assert len(mat_dict["cm"]) == 9  # A matrix
        assert len(mat_dict["b_coupling"]) == 9  # B matrix
        assert len(mat_dict["cb"]) == 9  # D matrix

    def test_mesh_connectivity(self):
        """Verify mesh has correct node sets and connectivity."""
        mesh = _build_composite_plate_mesh()

        # Check node sets exist
        assert mesh.get_node_set("clamped") is not None
        assert mesh.get_node_set("free_face") is not None
        assert mesh.get_node_set("free_center") is not None

        # Check clamped at y=0
        clamped = list(mesh.get_node_set("clamped").nodes.values())
        for n in clamped:
            assert np.isclose(n.y, 0.0, atol=1e-12)

        # Check free face at y=L
        free = list(mesh.get_node_set("free_face").nodes.values())
        for n in free:
            assert np.isclose(n.y, L, atol=1e-12)


# =============================================================================
# Integration tests (require CCX)
# =============================================================================


# Test de composite shell - tensión axial
def test_composite_axial_tension(tmp_path: Path):
    """Test composite shell under axial tension using full ABD matrices."""
    ccx_bin = _ccx_bin_or_skip()
    mesh = _build_composite_plate_mesh()

    # Load: 1000N in Y direction (along cantilever axis, not Z which is out-of-plane)
    # Mesh is in XY plane: X=width, Y=length (cantilever axis), Z=0
    # For axial tension, load should be along Y (length direction)
    load = (0.0, 1000.0, 0.0)

    # AeroElast solution (MITC4: 6 DOFs/node)
    node_coords = np.asarray([[n.x, n.y, n.z] for n in mesh.nodes], dtype=float)
    conn = [[mesh.node_id_to_index[nid] for nid in el.node_ids] for el in mesh.elements]
    elem_types = [4] * len(mesh.elements)  # MITC4

    # Composite material - pre-computed ABD matrices from CLT
    mat_dict = _make_laminate_mat(E1, E2, G12, nu12, thickness)
    mats = [mat_dict] * len(mesh.elements)

    asm = PyMeshAssembler(
        node_coords=node_coords, connectivity=conn, elem_types=elem_types, materials=mats
    )
    # ... rest of test stays the same

    # AeroElast solution (MITC4: 6 DOFs/node)
    node_coords = np.asarray([[n.x, n.y, n.z] for n in mesh.nodes], dtype=float)
    conn = [[mesh.node_id_to_index[nid] for nid in el.node_ids] for el in mesh.elements]
    elem_types = [4] * len(mesh.elements)  # MITC4

    # Composite material - pre-computed ABD matrices from CLT
    mat_dict = _make_laminate_mat(E1, E2, G12, nu12, thickness)
    mats = [mat_dict] * len(mesh.elements)

    asm = PyMeshAssembler(
        node_coords=node_coords, connectivity=conn, elem_types=elem_types, materials=mats
    )
    rows, cols, vals = asm.assemble_k()
    n = asm.dofs_count
    K = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()

    # Load vector (MITC4: 6 DOFs/node)
    f = np.zeros(n, dtype=float)
    center = next(iter(mesh.get_node_set("free_center").nodes.values()))
    i0 = mesh.node_id_to_index[center.id] * 6  # 6 DOFs per node for MITC4
    f[i0 + 1] = load[1]  # Y direction (along cantilever axis)

    # Boundary conditions (MITC4: 6 DOFs/node)
    clamped = {
        mesh.node_id_to_index[n.id] * 6 + i
        for n in mesh.get_node_set("clamped").nodes.values()
        for i in range(6)
    }
    free_mask = np.ones(n, dtype=bool)
    for dof in clamped:
        free_mask[dof] = False
    free = np.where(free_mask)[0]

    K_ff = K[np.ix_(free, free)]
    f_free = f[free]

    u_free = spsolve(K_ff, f_free)
    u = np.zeros(n, dtype=float)
    u[free] = u_free

    aero_disp = u[i0 : i0 + 6]

    # CCX solution
    inp_path = tmp_path / "composite_axial.inp"
    _write_composite_ccx_inp(inp_path, mesh, load)

    # Get node IDs for CCX output
    node_ids = [
        mesh.node_id_to_index[n.id] + 1 for n in mesh.get_node_set("free_center").nodes.values()
    ]

    result = _run_ccx(inp_path, ccx_bin)
    if result.returncode != 0:
        _fail_ccx(result, inp_path)

    frd_path = inp_path.with_suffix(".frd")
    if not frd_path.exists():
        pytest.fail(f"No FRD output from CCX for {inp_path.name}")

    ccx_disp = _parse_frd_disp(frd_path, node_ids)
    if not ccx_disp:
        pytest.fail(f"Could not parse CCX displacements from {frd_path.name}")

    # Compare the same free-centre node on both sides
    ccx_uy = abs(ccx_disp[node_ids[0]][1])

    aero_uy = abs(aero_disp[1])
    rel_error = abs(aero_uy - ccx_uy) / max(ccx_uy, 1e-10)
    print(f"AeroElast UY: {aero_uy * 1e6:.2f} um, CCX: {ccx_uy * 1e6:.2f} um")
    print(f"Relative error: {rel_error * 100:.2f}%")

    assert rel_error < 0.1, f"Composite axial: {rel_error * 100:.1f}% error (max 10%)"


# Test de composite shell - isotrópico equivalente
def test_composite_isotropic_equiv(tmp_path: Path):
    """Test composite shell using isotropic equivalent (E_avg) for comparison."""
    ccx_bin = _ccx_bin_or_skip()
    mesh = _build_composite_plate_mesh()

    # Load: 1000N in Y direction
    load = (0.0, 1000.0, 0.0)

    # AeroElast solution with ISOTROPIC equivalent (same as CCX uses)
    node_coords = np.asarray([[n.x, n.y, n.z] for n in mesh.nodes], dtype=float)
    conn = [[mesh.node_id_to_index[nid] for nid in el.node_ids] for el in mesh.elements]
    elem_types = [4] * len(mesh.elements)

    # Isotropic equivalent material (average E, same as CCX)
    e_equiv = (E1 + E2) / 2
    mats_iso = [
        {
            "type": "isotropic",
            "e": e_equiv,
            "nu": nu12,
            "rho": 0.0,
            "thickness": thickness,
        }
    ] * len(mesh.elements)

    asm_iso = PyMeshAssembler(
        node_coords=node_coords, connectivity=conn, elem_types=elem_types, materials=mats_iso
    )
    rows, cols, vals = asm_iso.assemble_k()
    n = asm_iso.dofs_count
    K_iso = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()

    # Load vector
    f = np.zeros(n, dtype=float)
    center = next(iter(mesh.get_node_set("free_center").nodes.values()))
    i0 = mesh.node_id_to_index[center.id] * 6
    f[i0 + 1] = load[1]

    # Boundary conditions
    clamped = {
        mesh.node_id_to_index[n.id] * 6 + i
        for n in mesh.get_node_set("clamped").nodes.values()
        for i in range(6)
    }
    free_mask = np.ones(n, dtype=bool)
    for dof in clamped:
        free_mask[dof] = False
    free = np.where(free_mask)[0]

    K_ff = K_iso[np.ix_(free, free)]
    f_free = f[free]
    u_free = spsolve(K_ff, f_free)
    u = np.zeros(n, dtype=float)
    u[free] = u_free

    aero_iso_disp = u[i0 : i0 + 6]

    # CCX solution -- the same isotropic material expressed as a single-ply
    # laminate.  A single isotropic ply gives A = Q*t and D = Q*t^3/12 with
    # Q11 = E/(1 - nu^2), i.e. exactly an isotropic shell, so the only
    # remaining difference is the element formulation (MITC4 vs S8R).
    g_iso = e_equiv / (2.0 * (1.0 + nu12))
    mat_iso = Material(
        name="iso",
        E=(e_equiv, e_equiv, e_equiv),
        G=(g_iso, g_iso, g_iso),
        nu=(nu12, nu12, nu12),
        rho=0.0,
    )
    lam_iso = create_laminate_from_angles(mat_iso, thickness, [0.0])
    inp_path = tmp_path / "composite_iso.inp"
    _write_composite_ccx_inp(inp_path, mesh, load, laminate=lam_iso)

    result = _run_ccx(inp_path, ccx_bin)
    if result.returncode != 0:
        _fail_ccx(result, inp_path)

    frd_path = inp_path.with_suffix(".frd")
    if not frd_path.exists():
        pytest.fail(f"No FRD output from CCX for {inp_path.name}")

    node_ids = [
        mesh.node_id_to_index[n.id] + 1 for n in mesh.get_node_set("free_center").nodes.values()
    ]
    ccx_disp = _parse_frd_disp(frd_path, node_ids)
    if not ccx_disp:
        pytest.fail(f"Could not parse CCX displacements from {frd_path.name}")

    ccx_uy = abs(ccx_disp[node_ids[0]][1])

    # Compare
    aero_uy = abs(aero_iso_disp[1])
    rel_error = abs(aero_uy - ccx_uy) / max(ccx_uy, 1e-10)
    print(f"Isotropic equiv UY: {aero_uy * 1e6:.2f} um, CCX: {ccx_uy * 1e6:.2f} um")
    print(f"Relative error: {rel_error * 100:.2f}%")

    assert rel_error < 0.1, f"Isotropic equiv: {rel_error * 100:.1f}% error (max 10%)"


# Test de composite shell - bending
def test_composite_bending(tmp_path: Path):
    """Test composite shell under transverse bending."""
    ccx_bin = _ccx_bin_or_skip()
    mesh = _build_composite_plate_mesh()

    # Load: 100N in X direction at free end
    load = (100.0, 0.0, 0.0)

    # AeroElast solution (MITC4: 6 DOFs/node)
    node_coords = np.asarray([[n.x, n.y, n.z] for n in mesh.nodes], dtype=float)
    conn = [[mesh.node_id_to_index[nid] for nid in el.node_ids] for el in mesh.elements]
    elem_types = [4] * len(mesh.elements)

    # Composite material - pre-computed ABD matrices from CLT
    mat_dict = _make_laminate_mat(E1, E2, G12, nu12, thickness)
    mats = [mat_dict] * len(mesh.elements)

    asm = PyMeshAssembler(
        node_coords=node_coords, connectivity=conn, elem_types=elem_types, materials=mats
    )
    rows, cols, vals = asm.assemble_k()
    n = asm.dofs_count
    K = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()

    # Load vector (MITC4: 6 DOFs/node)
    f = np.zeros(n, dtype=float)
    center = next(iter(mesh.get_node_set("free_center").nodes.values()))
    i0 = mesh.node_id_to_index[center.id] * 6  # 6 DOFs per node for MITC4
    f[i0] = load[0]  # X direction

    # Boundary conditions (MITC4: 6 DOFs/node)
    clamped = {
        mesh.node_id_to_index[n.id] * 6 + i
        for n in mesh.get_node_set("clamped").nodes.values()
        for i in range(6)
    }
    free_mask = np.ones(n, dtype=bool)
    for dof in clamped:
        free_mask[dof] = False
    free = np.where(free_mask)[0]

    K_ff = K[np.ix_(free, free)]
    f_free = f[free]

    u_free = spsolve(K_ff, f_free)
    u = np.zeros(n, dtype=float)
    u[free] = u_free

    aero_disp = u[i0 : i0 + 6]

    # CCX solution
    inp_path = tmp_path / "composite_bending.inp"
    _write_composite_ccx_inp(inp_path, mesh, load)

    node_ids = [
        mesh.node_id_to_index[n.id] + 1 for n in mesh.get_node_set("free_center").nodes.values()
    ]

    result = _run_ccx(inp_path, ccx_bin)
    if result.returncode != 0:
        _fail_ccx(result, inp_path)

    frd_path = inp_path.with_suffix(".frd")
    if not frd_path.exists():
        pytest.fail(f"No FRD output from CCX for {inp_path.name}")

    ccx_disp = _parse_frd_disp(frd_path, node_ids)
    if not ccx_disp:
        pytest.fail(f"Could not parse CCX displacements from {frd_path.name}")

    ccx_vals = ccx_disp[node_ids[0]]

    # Compare (CCX has 3 DOFs, MITC4 has 6)
    rel_error = abs(aero_disp[0] - ccx_vals[0]) / max(abs(ccx_vals[0]), 1e-10)
    print(f"AeroElast X: {aero_disp[0] * 1e6:.2f} um, CCX: {ccx_vals[0] * 1e6:.2f} um")
    print(f"Relative error: {rel_error * 100:.2f}%")

    assert rel_error < 0.1, f"Composite bending: {rel_error * 100:.1f}% error (max 10%)"
